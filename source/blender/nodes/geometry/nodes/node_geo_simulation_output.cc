/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "BLI_string_utils.h"

#include "BKE_compute_contexts.hh"
#include "BKE_scene.h"

#include "DEG_depsgraph_query.h"

#include "UI_interface.h"
#include "UI_resources.h"

#include "NOD_common.h"
#include "NOD_socket.h"

#include "node_geometry_util.hh"

namespace blender::nodes {

void socket_declarations_for_simulation_items(const Span<NodeSimulationItem> items,
                                              NodeDeclaration &r_declaration)
{
  for (const NodeSimulationItem &item : items) {
    switch (eNodeSocketDatatype(item.socket_type)) {
      case SOCK_GEOMETRY: {
        {
          std::unique_ptr<decl::Geometry> decl = std::make_unique<decl::Geometry>();
          decl->name = item.name;
          decl->identifier = item.name;
          decl->in_out = SOCK_IN;
          r_declaration.inputs.append(std::move(decl));
        }
        {
          std::unique_ptr<decl::Geometry> decl = std::make_unique<decl::Geometry>();
          decl->name = item.name;
          decl->identifier = item.name;
          decl->in_out = SOCK_OUT;
          r_declaration.outputs.append(std::move(decl));
        }
        break;
      }
      default:
        BLI_assert_unreachable();
    }
  }
  r_declaration.inputs.append(decl::create_extend_declaration(SOCK_IN));
  r_declaration.outputs.append(decl::create_extend_declaration(SOCK_OUT));
}

static bool simulation_items_unique_name_check(void *arg, const char *name)
{
  const NodeGeometrySimulationOutput &storage = *static_cast<const NodeGeometrySimulationOutput *>(
      arg);
  for (const NodeSimulationItem &item : Span(storage.items, storage.items_num)) {
    if (STREQ(item.name, name)) {
      return true;
    }
  }
  if (STREQ(name, "Delta Time")) {
    return true;
  }
  return false;
}

NodeSimulationItem *simulation_item_add_from_socket(NodeGeometrySimulationOutput &storage,
                                                    const bNodeSocket &socket)
{
  if (socket.type != SOCK_GEOMETRY) {
    return nullptr;
  }
  char unique_name[MAX_NAME + 4] = "";
  BLI_uniquename_cb(simulation_items_unique_name_check,
                    &storage,
                    socket.name,
                    '.',
                    unique_name,
                    ARRAY_SIZE(unique_name));

  NodeSimulationItem *old_items = storage.items;
  storage.items = MEM_cnew_array<NodeSimulationItem>(storage.items_num + 1, __func__);
  for (const int i : IndexRange(storage.items_num)) {
    storage.items[i].name = old_items[i].name;
    storage.items[i].socket_type = old_items[i].socket_type;
  }

  NodeSimulationItem &added_item = storage.items[storage.items_num];
  added_item.name = BLI_strdup(unique_name);
  added_item.socket_type = socket.type;

  storage.items_num++;
  MEM_SAFE_FREE(old_items);

  return &added_item;
}

const CPPType &get_simulation_item_cpp_type(const NodeSimulationItem &item)
{
  switch (item.socket_type) {
    case SOCK_GEOMETRY:
      return CPPType::get<GeometrySet>();
    default:
      BLI_assert_unreachable();
      return CPPType::get<GeometrySet>();
  }
}

}  // namespace blender::nodes

namespace blender::nodes::node_geo_simulation_output_cc {

NODE_STORAGE_FUNCS(NodeGeometrySimulationOutput);

struct EvalData {
  bool is_first_evaluation = true;
};

class LazyFunctionForSimulationOutputNode final : public LazyFunction {
  int32_t node_id_;
  Span<NodeSimulationItem> simulation_items_;

 public:
  LazyFunctionForSimulationOutputNode(const bNode &node) : node_id_(node.identifier)
  {
    const NodeGeometrySimulationOutput &storage = node_storage(node);
    simulation_items_ = {storage.items, storage.items_num};
    for (const NodeSimulationItem &item : Span(storage.items, storage.items_num)) {
      const CPPType &type = get_simulation_item_cpp_type(item);
      inputs_.append_as(item.name, type, lf::ValueUsage::Maybe);
      outputs_.append_as(item.name, type);
    }
  }

  void *init_storage(LinearAllocator<> &allocator) const
  {
    return allocator.construct<EvalData>().get();
  }

  void destruct_storage(void *storage) const
  {
    std::destroy_at(static_cast<EvalData *>(storage));
  }

  void execute_impl(lf::Params &params, const lf::Context &context) const final
  {
    GeoNodesLFUserData &user_data = *static_cast<GeoNodesLFUserData *>(context.user_data);
    GeoNodesModifierData &modifier_data = *user_data.modifier_data;
    EvalData &eval_data = *static_cast<EvalData *>(context.storage);
    BLI_SCOPED_DEFER([&]() { eval_data.is_first_evaluation = false; });

    const bke::sim::SimulationZoneID zone_id = get_simulation_zone_id(*user_data.compute_context,
                                                                      node_id_);

    const bke::sim::SimulationZoneState *current_zone_state =
        modifier_data.current_simulation_state ?
            modifier_data.current_simulation_state->get_zone_state(zone_id) :
            nullptr;
    if (eval_data.is_first_evaluation && current_zone_state != nullptr) {
      this->output_cached_state(params, *current_zone_state);
      return;
    }

    if (modifier_data.current_simulation_state_for_write == nullptr) {
      const bke::sim::SimulationZoneState *prev_zone_state =
          modifier_data.prev_simulation_state ?
              modifier_data.prev_simulation_state->get_zone_state(zone_id) :
              nullptr;
      if (prev_zone_state == nullptr) {
        params.set_default_remaining_outputs();
        return;
      }
      const bke::sim::SimulationZoneState *next_zone_state =
          modifier_data.next_simulation_state ?
              modifier_data.next_simulation_state->get_zone_state(zone_id) :
              nullptr;
      if (next_zone_state == nullptr) {
        this->output_cached_state(params, *prev_zone_state);
        return;
      }
      this->output_mixed_cached_state(
          params, *prev_zone_state, *next_zone_state, modifier_data.simulation_state_mix_factor);
      return;
    }

    bke::sim::SimulationZoneState &new_zone_state =
        modifier_data.current_simulation_state_for_write->get_zone_state_for_write(zone_id);
    new_zone_state.items.reinitialize(simulation_items_.size());

    bool all_available = true;
    for (const int i : simulation_items_.index_range()) {
      GeometrySet *input_geometry = params.try_get_input_data_ptr_or_request<GeometrySet>(i);
      if (input_geometry == nullptr) {
        all_available = false;
        continue;
      }
      input_geometry->ensure_owns_direct_data();
      new_zone_state.items[i] = std::make_unique<bke::sim::GeometrySimulationStateItem>(
          std::move(*input_geometry));
    }

    if (all_available) {
      this->output_cached_state(params, new_zone_state);
    }
  }

  void output_cached_state(lf::Params &params, const bke::sim::SimulationZoneState &state) const
  {
    for (const int i : simulation_items_.index_range()) {
      if (i >= state.items.size()) {
        continue;
      }
      const bke::sim::SimulationStateItem *item = state.items[i].get();
      if (item == nullptr) {
        continue;
      }
      if (auto *geometry_item = dynamic_cast<const bke::sim::GeometrySimulationStateItem *>(
              item)) {
        params.set_output(i, geometry_item->geometry());
      }
    }
    params.set_default_remaining_outputs();
  }

  void output_mixed_cached_state(lf::Params &params,
                                 const bke::sim::SimulationZoneState &prev_state,
                                 const bke::sim::SimulationZoneState &next_state,
                                 const float mix_factor) const
  {
    /* TODO: Implement subframe mixing. */
    this->output_cached_state(params, prev_state);
    UNUSED_VARS(next_state, mix_factor);
  }
};

}  // namespace blender::nodes::node_geo_simulation_output_cc

namespace blender::nodes {

std::unique_ptr<LazyFunction> get_simulation_output_lazy_function(const bNode &node)
{
  namespace file_ns = blender::nodes::node_geo_simulation_output_cc;
  BLI_assert(node.type == GEO_NODE_SIMULATION_OUTPUT);
  return std::make_unique<file_ns::LazyFunctionForSimulationOutputNode>(node);
}

bke::sim::SimulationZoneID get_simulation_zone_id(const ComputeContext &compute_context,
                                                  const int output_node_id)
{
  bke::sim::SimulationZoneID zone_id;
  for (const ComputeContext *context = &compute_context; context != nullptr;
       context = context->parent()) {
    if (const auto *node_context = dynamic_cast<const bke::NodeGroupComputeContext *>(context)) {
      zone_id.node_ids.append(node_context->node_id());
    }
  }
  std::reverse(zone_id.node_ids.begin(), zone_id.node_ids.end());
  zone_id.node_ids.append(output_node_id);
  return zone_id;
}

}  // namespace blender::nodes

namespace blender::nodes::node_geo_simulation_output_cc {

static void node_declare_dynamic(const bNodeTree & /*node_tree*/,
                                 const bNode &node,
                                 NodeDeclaration &r_declaration)
{
  const NodeGeometrySimulationOutput &storage = node_storage(node);
  socket_declarations_for_simulation_items({storage.items, storage.items_num}, r_declaration);
}

static void node_init(bNodeTree * /*tree*/, bNode *node)
{
  NodeGeometrySimulationOutput *data = MEM_cnew<NodeGeometrySimulationOutput>(__func__);
  data->items = MEM_cnew_array<NodeSimulationItem>(1, __func__);
  data->items[0].name = BLI_strdup(DATA_("Geometry"));
  data->items[0].socket_type = SOCK_GEOMETRY;
  data->items_num = 1;
  node->storage = data;
}

static void node_free_storage(bNode *node)
{
  if (!node->storage) {
    return;
  }
  NodeGeometrySimulationOutput &storage = node_storage(*node);
  for (NodeSimulationItem &item : MutableSpan(storage.items, storage.items_num)) {
    MEM_SAFE_FREE(item.name);
  }
  MEM_SAFE_FREE(storage.items);
  MEM_freeN(node->storage);
}

static void node_copy_storage(bNodeTree * /*dst_tree*/, bNode *dst_node, const bNode *src_node)
{
  const NodeGeometrySimulationOutput &src_storage = node_storage(*src_node);
  NodeGeometrySimulationOutput *dst_storage = MEM_cnew<NodeGeometrySimulationOutput>(__func__);

  dst_storage->items = MEM_cnew_array<NodeSimulationItem>(src_storage.items_num, __func__);
  dst_storage->items_num = src_storage.items_num;
  for (const int i : IndexRange(src_storage.items_num)) {
    if (char *name = src_storage.items[i].name) {
      dst_storage->items[i].name = BLI_strdup(name);
      dst_storage->items[i].socket_type = src_storage.items[i].socket_type;
    }
  }

  dst_node->storage = dst_storage;
}

static bool node_insert_link(bNodeTree *ntree, bNode *node, bNodeLink *link)
{
  NodeGeometrySimulationOutput &storage = node_storage(*node);
  if (link->tonode == node) {
    if (link->tosock->identifier == StringRef("__extend__")) {
      if (const NodeSimulationItem *item = simulation_item_add_from_socket(storage,
                                                                           *link->fromsock)) {
        update_node_declaration_and_sockets(*ntree, *node);
        link->tosock = nodeFindSocket(node, SOCK_IN, item->name);
      }
      else {
        return false;
      }
    }
  }
  else {
    BLI_assert(link->fromnode == node);
    if (link->fromsock->identifier == StringRef("__extend__")) {
      if (const NodeSimulationItem *item = simulation_item_add_from_socket(storage,
                                                                           *link->tosock)) {
        update_node_declaration_and_sockets(*ntree, *node);
        link->fromsock = nodeFindSocket(node, SOCK_OUT, item->name);
      }
      else {
        return false;
      }
    }
  }
  return true;
}

}  // namespace blender::nodes::node_geo_simulation_output_cc

void register_node_type_geo_simulation_output()
{
  namespace file_ns = blender::nodes::node_geo_simulation_output_cc;

  static bNodeType ntype;

  geo_node_type_base(
      &ntype, GEO_NODE_SIMULATION_OUTPUT, "Simulation Output", NODE_CLASS_INTERFACE);
  ntype.initfunc = file_ns::node_init;
  ntype.declare_dynamic = file_ns::node_declare_dynamic;
  ntype.insert_link = file_ns::node_insert_link;
  node_type_storage(&ntype,
                    "NodeGeometrySimulationOutput",
                    file_ns::node_free_storage,
                    file_ns::node_copy_storage);
  nodeRegisterType(&ntype);
}
