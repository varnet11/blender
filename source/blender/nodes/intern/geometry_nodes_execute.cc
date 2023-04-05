/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup nodes
 *
 * This file mainly converts a #bNodeTree into a lazy-function graph. This generally works by
 * creating a lazy-function for every node, which is then put into the lazy-function graph. Then
 * the nodes in the new graph are linked based on links in the original #bNodeTree. Some additional
 * nodes are inserted for things like type conversions and multi-input sockets.
 *
 * Currently, lazy-functions are even created for nodes that don't strictly require it, like
 * reroutes or muted nodes. In the future we could avoid that at the cost of additional code
 * complexity. So far, this does not seem to be a performance issue.
 */

#include "NOD_geometry_exec.hh"
#include "NOD_geometry_nodes_lazy_function.hh"
#include "NOD_multi_function.hh"
#include "NOD_node_declaration.hh"

#include "BLI_cpp_types.hh"
#include "BLI_dot_export.hh"
#include "BLI_hash.h"
#include "BLI_lazy_threading.hh"
#include "BLI_map.hh"

#include "DNA_ID.h"

#include "BKE_compute_contexts.hh"
#include "BKE_geometry_set.hh"
#include "BKE_idprop.hh"
#include "BKE_type_conversions.hh"

#include "FN_field_cpp_type.hh"
#include "FN_lazy_function_execute.hh"
#include "FN_lazy_function_graph_executor.hh"

#include "DEG_depsgraph_query.h"

namespace blender::nodes {

/**
 * \return Whether using an attribute to input values of this type is supported, and the node
 * group's input for this socket accepts a field rather than just single values.
 */
static bool input_has_attribute_toggle(const bNodeTree &node_tree, const int socket_index)
{
  BLI_assert(node_tree.runtime->field_inferencing_interface);
  const nodes::FieldInferencingInterface &field_interface =
      *node_tree.runtime->field_inferencing_interface;
  return field_interface.inputs[socket_index] != nodes::InputSocketFieldType::None;
}

static std::unique_ptr<IDProperty, bke::idprop::IDPropertyDeleter> id_property_create_from_socket(
    const bNodeSocket &socket)
{
  switch (socket.type) {
    case SOCK_FLOAT: {
      const bNodeSocketValueFloat *value = static_cast<const bNodeSocketValueFloat *>(
          socket.default_value);
      auto property = bke::idprop::create(socket.identifier, value->value);
      IDPropertyUIDataFloat *ui_data = (IDPropertyUIDataFloat *)IDP_ui_data_ensure(property.get());
      ui_data->base.rna_subtype = value->subtype;
      ui_data->soft_min = double(value->min);
      ui_data->soft_max = double(value->max);
      ui_data->default_value = value->value;
      return property;
    }
    case SOCK_INT: {
      const bNodeSocketValueInt *value = static_cast<const bNodeSocketValueInt *>(
          socket.default_value);
      auto property = bke::idprop::create(socket.identifier, value->value);
      IDPropertyUIDataInt *ui_data = (IDPropertyUIDataInt *)IDP_ui_data_ensure(property.get());
      ui_data->base.rna_subtype = value->subtype;
      ui_data->soft_min = value->min;
      ui_data->soft_max = value->max;
      ui_data->default_value = value->value;
      return property;
    }
    case SOCK_VECTOR: {
      const bNodeSocketValueVector *value = static_cast<const bNodeSocketValueVector *>(
          socket.default_value);
      auto property = bke::idprop::create(
          socket.identifier, Span<float>{value->value[0], value->value[1], value->value[2]});
      IDPropertyUIDataFloat *ui_data = (IDPropertyUIDataFloat *)IDP_ui_data_ensure(property.get());
      ui_data->base.rna_subtype = value->subtype;
      ui_data->soft_min = double(value->min);
      ui_data->soft_max = double(value->max);
      ui_data->default_array = (double *)MEM_mallocN(sizeof(double[3]), "mod_prop_default");
      ui_data->default_array_len = 3;
      for (const int i : IndexRange(3)) {
        ui_data->default_array[i] = double(value->value[i]);
      }
      return property;
    }
    case SOCK_RGBA: {
      const bNodeSocketValueRGBA *value = static_cast<const bNodeSocketValueRGBA *>(
          socket.default_value);
      auto property = bke::idprop::create(
          socket.identifier,
          Span<float>{value->value[0], value->value[1], value->value[2], value->value[3]});
      IDPropertyUIDataFloat *ui_data = (IDPropertyUIDataFloat *)IDP_ui_data_ensure(property.get());
      ui_data->base.rna_subtype = PROP_COLOR;
      ui_data->default_array = (double *)MEM_mallocN(sizeof(double[4]), __func__);
      ui_data->default_array_len = 4;
      ui_data->min = 0.0;
      ui_data->max = FLT_MAX;
      ui_data->soft_min = 0.0;
      ui_data->soft_max = 1.0;
      for (const int i : IndexRange(4)) {
        ui_data->default_array[i] = double(value->value[i]);
      }
      return property;
    }
    case SOCK_BOOLEAN: {
      const bNodeSocketValueBoolean *value = static_cast<const bNodeSocketValueBoolean *>(
          socket.default_value);
      auto property = bke::idprop::create_bool(socket.identifier, value->value);
      IDPropertyUIDataBool *ui_data = (IDPropertyUIDataBool *)IDP_ui_data_ensure(property.get());
      ui_data->default_value = value->value != 0;
      return property;
    }
    case SOCK_STRING: {
      const bNodeSocketValueString *value = static_cast<const bNodeSocketValueString *>(
          socket.default_value);
      auto property = bke::idprop::create(socket.identifier, value->value);
      IDPropertyUIDataString *ui_data = (IDPropertyUIDataString *)IDP_ui_data_ensure(
          property.get());
      ui_data->default_value = BLI_strdup(value->value);
      return property;
    }
    case SOCK_OBJECT: {
      const bNodeSocketValueObject *value = static_cast<const bNodeSocketValueObject *>(
          socket.default_value);
      auto property = bke::idprop::create(socket.identifier, reinterpret_cast<ID *>(value->value));
      IDPropertyUIDataID *ui_data = (IDPropertyUIDataID *)IDP_ui_data_ensure(property.get());
      ui_data->id_type = ID_OB;
      return property;
    }
    case SOCK_COLLECTION: {
      const bNodeSocketValueCollection *value = static_cast<const bNodeSocketValueCollection *>(
          socket.default_value);
      return bke::idprop::create(socket.identifier, reinterpret_cast<ID *>(value->value));
    }
    case SOCK_TEXTURE: {
      const bNodeSocketValueTexture *value = static_cast<const bNodeSocketValueTexture *>(
          socket.default_value);
      return bke::idprop::create(socket.identifier, reinterpret_cast<ID *>(value->value));
    }
    case SOCK_IMAGE: {
      const bNodeSocketValueImage *value = static_cast<const bNodeSocketValueImage *>(
          socket.default_value);
      return bke::idprop::create(socket.identifier, reinterpret_cast<ID *>(value->value));
    }
    case SOCK_MATERIAL: {
      const bNodeSocketValueMaterial *value = static_cast<const bNodeSocketValueMaterial *>(
          socket.default_value);
      return bke::idprop::create(socket.identifier, reinterpret_cast<ID *>(value->value));
    }
  }
  return nullptr;
}

static bool id_property_type_matches_socket(const bNodeSocket &socket, const IDProperty &property)
{
  switch (socket.type) {
    case SOCK_FLOAT:
      return ELEM(property.type, IDP_FLOAT, IDP_DOUBLE);
    case SOCK_INT:
      return property.type == IDP_INT;
    case SOCK_VECTOR:
      return property.type == IDP_ARRAY && property.subtype == IDP_FLOAT && property.len == 3;
    case SOCK_RGBA:
      return property.type == IDP_ARRAY && property.subtype == IDP_FLOAT && property.len == 4;
    case SOCK_BOOLEAN:
      return property.type == IDP_BOOLEAN;
    case SOCK_STRING:
      return property.type == IDP_STRING;
    case SOCK_OBJECT:
    case SOCK_COLLECTION:
    case SOCK_TEXTURE:
    case SOCK_IMAGE:
    case SOCK_MATERIAL:
      return property.type == IDP_ID;
  }
  BLI_assert_unreachable();
  return false;
}

void initialize_group_input(const bNodeTree &tree,
                            const IDProperty *properties,
                            const int input_index,
                            void *r_value)
{
  const bNodeSocket &io_input = *tree.interface_inputs()[input_index];
  const bNodeSocketType &socket_type = *io_input.typeinfo;
  const eNodeSocketDatatype socket_data_type = static_cast<eNodeSocketDatatype>(io_input.type);
  if (properties == nullptr) {
    socket_type.get_geometry_nodes_cpp_value(io_input, r_value);
    return;
  }
  const IDProperty *property = IDP_GetPropertyFromGroup(properties, io_input.identifier);
  if (property == nullptr) {
    socket_type.get_geometry_nodes_cpp_value(io_input, r_value);
    return;
  }
  if (!id_property_type_matches_socket(io_input, *property)) {
    socket_type.get_geometry_nodes_cpp_value(io_input, r_value);
    return;
  }

  if (!input_has_attribute_toggle(tree, input_index)) {
    init_socket_cpp_value_from_property(*property, socket_data_type, r_value);
    return;
  }

  const IDProperty *property_use_attribute = IDP_GetPropertyFromGroup(
      properties, (io_input.identifier + use_attribute_suffix).c_str());
  const IDProperty *property_attribute_name = IDP_GetPropertyFromGroup(
      properties, (io_input.identifier + attribute_name_suffix).c_str());
  if (property_use_attribute == nullptr || property_attribute_name == nullptr) {
    init_socket_cpp_value_from_property(*property, socket_data_type, r_value);
    return;
  }

  const bool use_attribute = IDP_Int(property_use_attribute) != 0;
  if (use_attribute) {
    const StringRef attribute_name{IDP_String(property_attribute_name)};
    if (!bke::allow_procedural_attribute_access(attribute_name)) {
      init_socket_cpp_value_from_property(*property, socket_data_type, r_value);
      return;
    }
    fn::GField attribute_field = bke::AttributeFieldInput::Create(attribute_name,
                                                                  *socket_type.base_cpp_type);
    const auto *value_or_field_cpp_type = fn::ValueOrFieldCPPType::get_from_self(
        *socket_type.geometry_nodes_cpp_type);
    BLI_assert(value_or_field_cpp_type != nullptr);
    value_or_field_cpp_type->construct_from_field(r_value, std::move(attribute_field));
  }
  else {
    init_socket_cpp_value_from_property(*property, socket_data_type, r_value);
  }
}

static void init_socket_cpp_value_from_property(const IDProperty &property,
                                                const eNodeSocketDatatype socket_value_type,
                                                void *r_value)
{
  switch (socket_value_type) {
    case SOCK_FLOAT: {
      float value = 0.0f;
      if (property.type == IDP_FLOAT) {
        value = IDP_Float(&property);
      }
      else if (property.type == IDP_DOUBLE) {
        value = float(IDP_Double(&property));
      }
      new (r_value) fn::ValueOrField<float>(value);
      break;
    }
    case SOCK_INT: {
      int value = IDP_Int(&property);
      new (r_value) fn::ValueOrField<int>(value);
      break;
    }
    case SOCK_VECTOR: {
      float3 value = (const float *)IDP_Array(&property);
      new (r_value) fn::ValueOrField<float3>(value);
      break;
    }
    case SOCK_RGBA: {
      ColorGeometry4f value = (const float *)IDP_Array(&property);
      new (r_value) fn::ValueOrField<ColorGeometry4f>(value);
      break;
    }
    case SOCK_BOOLEAN: {
      const bool value = IDP_Bool(&property);
      new (r_value) fn::ValueOrField<bool>(value);
      break;
    }
    case SOCK_STRING: {
      std::string value = IDP_String(&property);
      new (r_value) fn::ValueOrField<std::string>(std::move(value));
      break;
    }
    case SOCK_OBJECT: {
      ID *id = IDP_Id(&property);
      Object *object = (id && GS(id->name) == ID_OB) ? (Object *)id : nullptr;
      *(Object **)r_value = object;
      break;
    }
    case SOCK_COLLECTION: {
      ID *id = IDP_Id(&property);
      Collection *collection = (id && GS(id->name) == ID_GR) ? (Collection *)id : nullptr;
      *(Collection **)r_value = collection;
      break;
    }
    case SOCK_TEXTURE: {
      ID *id = IDP_Id(&property);
      Tex *texture = (id && GS(id->name) == ID_TE) ? (Tex *)id : nullptr;
      *(Tex **)r_value = texture;
      break;
    }
    case SOCK_IMAGE: {
      ID *id = IDP_Id(&property);
      Image *image = (id && GS(id->name) == ID_IM) ? (Image *)id : nullptr;
      *(Image **)r_value = image;
      break;
    }
    case SOCK_MATERIAL: {
      ID *id = IDP_Id(&property);
      Material *material = (id && GS(id->name) == ID_MA) ? (Material *)id : nullptr;
      *(Material **)r_value = material;
      break;
    }
    default: {
      BLI_assert_unreachable();
      break;
    }
  }
}

static GeometrySet compute_geometry(const bNodeTree &btree,
                                    const nodes::GeometryNodesLazyFunctionGraphInfo &lf_graph_info,
                                    GeometrySet input_geometry_set,
                                    NodesModifierData *nmd,
                                    const ModifierEvalContext *ctx)
{

  const nodes::GeometryNodeLazyFunctionGraphMapping &mapping = lf_graph_info.mapping;

  Vector<const lf::OutputSocket *> graph_inputs = mapping.group_input_sockets;
  graph_inputs.extend(mapping.group_output_used_sockets);
  graph_inputs.extend(mapping.attribute_set_by_geometry_output.values().begin(),
                      mapping.attribute_set_by_geometry_output.values().end());
  Vector<const lf::InputSocket *> graph_outputs = mapping.standard_group_output_sockets;

  Array<GMutablePointer> param_inputs(graph_inputs.size());
  Array<GMutablePointer> param_outputs(graph_outputs.size());
  Array<std::optional<lf::ValueUsage>> param_input_usages(graph_inputs.size());
  Array<lf::ValueUsage> param_output_usages(graph_outputs.size(), lf::ValueUsage::Used);
  Array<bool> param_set_outputs(graph_outputs.size(), false);

  nodes::GeometryNodesLazyFunctionLogger lf_logger(lf_graph_info);
  nodes::GeometryNodesLazyFunctionSideEffectProvider lf_side_effect_provider;

  lf::GraphExecutor graph_executor{
      lf_graph_info.graph, graph_inputs, graph_outputs, &lf_logger, &lf_side_effect_provider};

  nodes::GeoNodesModifierData geo_nodes_modifier_data;
  geo_nodes_modifier_data.depsgraph = ctx->depsgraph;
  geo_nodes_modifier_data.self_object = ctx->object;
  auto eval_log = std::make_unique<geo_log::GeoModifierLog>();

  Set<ComputeContextHash> socket_log_contexts;
  if (logging_enabled(ctx)) {
    geo_nodes_modifier_data.eval_log = eval_log.get();

    find_socket_log_contexts(*nmd, *ctx, socket_log_contexts);
    geo_nodes_modifier_data.socket_log_contexts = &socket_log_contexts;
  }
  MultiValueMap<ComputeContextHash, const lf::FunctionNode *> r_side_effect_nodes;
  find_side_effect_nodes(*nmd, *ctx, r_side_effect_nodes);
  geo_nodes_modifier_data.side_effect_nodes = &r_side_effect_nodes;
  nodes::GeoNodesLFUserData user_data;
  user_data.modifier_data = &geo_nodes_modifier_data;
  bke::ModifierComputeContext modifier_compute_context{nullptr, nmd->modifier.name};
  user_data.compute_context = &modifier_compute_context;

  LinearAllocator<> allocator;
  Vector<GMutablePointer> inputs_to_destruct;

  const IDProperty *properties = nmd->settings.properties;
  int input_index = -1;
  for (const int i : btree.interface_inputs().index_range()) {
    input_index++;
    const bNodeSocket &interface_socket = *btree.interface_inputs()[i];
    if (interface_socket.type == SOCK_GEOMETRY && input_index == 0) {
      param_inputs[input_index] = &input_geometry_set;
      continue;
    }

    const CPPType *type = interface_socket.typeinfo->geometry_nodes_cpp_type;
    BLI_assert(type != nullptr);
    void *value = allocator.allocate(type->size(), type->alignment());
    initialize_group_input(btree, properties, i, value);
    param_inputs[input_index] = {type, value};
    inputs_to_destruct.append({type, value});
  }

  Array<bool> output_used_inputs(btree.interface_outputs().size(), true);
  for (const int i : btree.interface_outputs().index_range()) {
    input_index++;
    param_inputs[input_index] = &output_used_inputs[i];
  }

  Array<bke::AnonymousAttributeSet> attributes_to_propagate(
      mapping.attribute_set_by_geometry_output.size());
  for (const int i : attributes_to_propagate.index_range()) {
    input_index++;
    param_inputs[input_index] = &attributes_to_propagate[i];
  }

  for (const int i : graph_outputs.index_range()) {
    const lf::InputSocket &socket = *graph_outputs[i];
    const CPPType &type = socket.type();
    void *buffer = allocator.allocate(type.size(), type.alignment());
    param_outputs[i] = {type, buffer};
  }

  lf::Context lf_context;
  lf_context.storage = graph_executor.init_storage(allocator);
  lf_context.user_data = &user_data;
  lf::BasicParams lf_params{graph_executor,
                            param_inputs,
                            param_outputs,
                            param_input_usages,
                            param_output_usages,
                            param_set_outputs};
  graph_executor.execute(lf_params, lf_context);
  graph_executor.destruct_storage(lf_context.storage);

  for (GMutablePointer &ptr : inputs_to_destruct) {
    ptr.destruct();
  }

  GeometrySet output_geometry_set = std::move(*static_cast<GeometrySet *>(param_outputs[0].get()));
  store_output_attributes(output_geometry_set, btree, properties, output_node, param_outputs);

  for (GMutablePointer &ptr : param_outputs) {
    ptr.destruct();
  }

  if (logging_enabled(ctx)) {
    NodesModifierData *nmd_orig = reinterpret_cast<NodesModifierData *>(
        BKE_modifier_get_original(ctx->object, &nmd->modifier));
    delete static_cast<geo_log::GeoModifierLog *>(nmd_orig->runtime_eval_log);
    nmd_orig->runtime_eval_log = eval_log.release();
  }

  return output_geometry_set;
}

}  // namespace blender::nodes
