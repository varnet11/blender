/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup edcurves
 */

#include "ED_curves.h"
#include "ED_object.h"
#include "ED_screen.h"
#include "ED_select_utils.h"
#include "ED_view3d.h"

#include "WM_api.h"

#include "BKE_attribute_math.hh"
#include "BKE_compute_contexts.hh"
#include "BKE_context.h"
#include "BKE_curves.hh"
#include "BKE_geometry_set.hh"
#include "BKE_layer.h"
#include "BKE_lib_id.h"
#include "BKE_node_runtime.hh"
#include "BKE_object.h"
#include "BKE_report.h"

#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "DEG_depsgraph.h"
#include "DEG_depsgraph_query.h"

#include "RNA_access.h"
#include "RNA_define.h"
#include "RNA_enum_types.h"
#include "RNA_prototypes.h"

#include "UI_interface.h"
#include "UI_resources.h"

#include "FN_lazy_function_execute.hh"

#include "NOD_geometry_nodes_execute.hh"
#include "NOD_geometry_nodes_lazy_function.hh"

namespace blender::ed::curves {

static GeometrySet compute_geometry(const bNodeTree &btree,
                                    const nodes::GeometryNodesLazyFunctionGraphInfo &lf_graph_info,
                                    Object &object,
                                    Depsgraph &depsgraph,
                                    const IDProperty *properties,
                                    const bNode &output_node,
                                    GeometrySet input_geometry_set)
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
  geo_nodes_modifier_data.depsgraph = &depsgraph;
  geo_nodes_modifier_data.self_object = &object;

  Set<ComputeContextHash> socket_log_contexts;

  MultiValueMap<ComputeContextHash, const lf::FunctionNode *> r_side_effect_nodes;
  nodes::GeoNodesLFUserData user_data;
  user_data.modifier_data = &geo_nodes_modifier_data;
  bke::ModifierComputeContext modifier_compute_context{nullptr, "actually a node group"};
  user_data.compute_context = &modifier_compute_context;

  LinearAllocator<> allocator;
  Vector<GMutablePointer> inputs_to_destruct;

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
    nodes::initialize_group_input(btree, properties, i, value);
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
  // store_output_attributes(output_geometry_set, btree, properties, output_node, param_outputs);

  for (GMutablePointer &ptr : param_outputs) {
    ptr.destruct();
  }

  return output_geometry_set;
}

static int run_node_group_exec(bContext *C, wmOperator *op)
{
  Scene *scene = CTX_data_scene(C);
  ViewLayer *view_layer = CTX_data_view_layer(C);
  Object *object = CTX_data_active_object(C);
  if (!object) {
    return OPERATOR_CANCELLED;
  }
  eObjectMode mode = eObjectMode(object->mode);

  const uint32_t session_uuid = RNA_int_get(op->ptr, "session_uuid");
  const ID *id = BKE_libblock_find_session_uuid(CTX_data_main(C), ID_NT, session_uuid);
  if (!id) {
    return OPERATOR_CANCELLED;
  }
  const bNodeTree &node_tree = reinterpret_cast<const bNodeTree &>(*id);
  const nodes::GeometryNodesLazyFunctionGraphInfo *lf_graph_info =
      nodes::ensure_geometry_nodes_lazy_function_graph(node_tree);
  if (lf_graph_info == nullptr) {
    BKE_report(op->reports, RPT_ERROR, "Cannot evaluate node group");
    return;
  }

  uint objects_len = 0;
  Object **objects = BKE_view_layer_array_from_objects_in_mode_unique_data(
      scene, view_layer, CTX_wm_view3d(C), &objects_len, mode);

  for (Object *object : Span(objects, objects_len)) {
    Curves &curves = *static_cast<Curves *>(object->data);
    GeometrySet geometry_set = GeometrySet::create_with_curves(&curves,
                                                               GeometryOwnershipType::Editable);
  }

  return OPERATOR_FINISHED;
}

void CURVES_OT_node_group(wmOperatorType *ot)
{
  ot->name = "Run Node Group";
  ot->idname = __func__;
  ot->description = "Execute a node group on curves";  // TODO: Retrieve from node group.

  ot->exec = run_node_group_exec;
  ot->poll = editable_curves_poll;

  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;

  PropertyRNA *prop = RNA_def_int(ot->srna,
                                  "session_uuid",
                                  0,
                                  INT32_MIN,
                                  INT32_MAX,
                                  "Session UUID",
                                  "Session UUID of the node group",
                                  INT32_MIN,
                                  INT32_MAX);
  RNA_def_property_flag(prop, PropertyFlag(PROP_SKIP_SAVE | PROP_HIDDEN));
}

}  // namespace blender::ed::curves
