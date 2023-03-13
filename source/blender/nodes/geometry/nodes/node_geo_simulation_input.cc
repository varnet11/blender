/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "BKE_compute_cache.hh"
#include "BKE_compute_contexts.hh"
#include "BKE_scene.h"

#include "DEG_depsgraph_query.h"

#include "UI_interface.h"
#include "UI_resources.h"

#include "NOD_socket.h"

#include "node_geometry_util.hh"

namespace blender::nodes::node_geo_simulation_input_cc {

NODE_STORAGE_FUNCS(NodeGeometrySimulationInput);

class LazyFunctionForSimulationInputNode final : public LazyFunction {
  int32_t output_node_id_;
  Span<NodeSimulationItem> simulation_items_;

 public:
  LazyFunctionForSimulationInputNode(const bNodeTree &node_tree, const bNode &node)
  {
    output_node_id_ = node_storage(node).output_node_id;
    const bNode &output_node = *node_tree.node_by_id(output_node_id_);
    const NodeGeometrySimulationOutput &storage = *static_cast<NodeGeometrySimulationOutput *>(
        output_node.storage);
    simulation_items_ = {storage.items, storage.items_num};
    outputs_.append_as("Delta Time", CPPType::get<float>());
    for (const NodeSimulationItem &item : Span(storage.items, storage.items_num)) {
      const CPPType &type = get_simulation_item_cpp_type(item);
      inputs_.append_as(item.name, type, lf::ValueUsage::Maybe);
      outputs_.append_as(item.name, type);
    }
  }

  void execute_impl(lf::Params &params, const lf::Context &context) const final
  {
    const GeoNodesLFUserData &user_data = *dynamic_cast<const GeoNodesLFUserData *>(
        context.user_data);
    const Scene *scene = DEG_get_input_scene(user_data.modifier_data->depsgraph);
    const float scene_ctime = BKE_scene_ctime_get(scene);
    const bke::sim::TimePoint time{int(scene_ctime), scene_ctime};

    bke::sim::ComputeCaches &all_caches = *user_data.modifier_data->cache_per_frame;
    const bke::NodeGroupComputeContext cache_context(user_data.compute_context, output_node_id_);
    bke::sim::SimulationCache *cache = all_caches.lookup_context(cache_context.hash());

    if (!params.output_was_set(0) && params.get_output_usage(0) == lf::ValueUsage::Used) {
      if (cache && !cache->is_empty()) {
        const float time_diff = scene_ctime - cache->last_run_time()->time;
        const double frame_rate = double(scene->r.frs_sec) / double(scene->r.frs_sec_base);
        params.set_output<float>(0, float(std::max(0.0f, time_diff) / frame_rate));
      }
      else {
        params.set_output<float>(0, 0.0f);
      }
    }

    for (const int i : inputs_.index_range()) {
      const int output_i = i + 1;
      if (params.get_output_usage(output_i) != lf::ValueUsage::Used) {
        continue;
      }
      if (params.output_was_set(output_i)) {
        continue;
      }
      const StringRef name = simulation_items_[i].name;
      const CPPType &type = *inputs_[i].type;
      void *output = params.get_output_data_ptr(output_i);
      if (cache) {
        /* Try to retrieve a cached value for the previous evaluation. */
        if (std::optional<GeometrySet> cached_value = cache->value_before_time(name, time)) {
          type.move_construct(&*cached_value, output);
          params.output_set(output_i);
          continue;
        }
      }

      /* There was no cached value, so try to retrieve the input value. */
      void *value = params.try_get_input_data_ptr_or_request(i);
      if (!value) {
        continue;
      }
      type.move_construct(value, output);
      params.output_set(output_i);
    }
  }
};

}  // namespace blender::nodes::node_geo_simulation_input_cc

namespace blender::nodes {

std::unique_ptr<LazyFunction> get_simulation_input_lazy_function(const bNodeTree &node_tree,
                                                                 const bNode &node)
{
  namespace file_ns = blender::nodes::node_geo_simulation_input_cc;
  BLI_assert(node.type == GEO_NODE_SIMULATION_INPUT);
  return std::make_unique<file_ns::LazyFunctionForSimulationInputNode>(node_tree, node);
}

}  // namespace blender::nodes

namespace blender::nodes::node_geo_simulation_input_cc {

static void node_declare_dynamic(const bNodeTree &node_tree,
                                 const bNode &node,
                                 NodeDeclaration &r_declaration)
{
  const bNode *output_node = node_tree.node_by_id(node_storage(node).output_node_id);
  if (!output_node) {
    return;
  }

  std::unique_ptr<decl::Float> delta_time = std::make_unique<decl::Float>();
  delta_time->identifier = N_("Delta Time");
  delta_time->name = delta_time->identifier;
  delta_time->in_out = SOCK_OUT;
  r_declaration.outputs.append(std::move(delta_time));

  const NodeGeometrySimulationOutput &storage = *static_cast<const NodeGeometrySimulationOutput *>(
      output_node->storage);
  socket_declarations_for_simulation_items({storage.items, storage.items_num}, r_declaration);
}

static void node_init(bNodeTree *tree, bNode *node)
{
  NodeGeometrySimulationInput *data = MEM_cnew<NodeGeometrySimulationInput>(__func__);

  VectorSet<int32_t> sim_output_ids;
  Set<int32_t> sim_input_output_ids;
  for (bNode *other_node : tree->all_nodes()) {
    if (other_node->type == GEO_NODE_SIMULATION_INPUT && other_node != node &&
        other_node->storage) {
      const NodeGeometrySimulationInput &storage = node_storage(*other_node);
      sim_input_output_ids.add_new(storage.output_node_id);
    }
    else if (other_node->type == GEO_NODE_SIMULATION_OUTPUT) {
      sim_output_ids.add_new(other_node->identifier);
    }
  }

  sim_output_ids.remove_if(
      [&](const int32_t identifier) { return sim_input_output_ids.contains(identifier); });

  if (sim_output_ids.size() == 1) {
    data->output_node_id = sim_output_ids[0];
  }
  else {
    data->output_node_id = 0;
  }

  node->storage = data;
}

static bool node_insert_link(bNodeTree *ntree, bNode *node, bNodeLink *link)
{
  const bNode *output_node = ntree->node_by_id(node_storage(*node).output_node_id);
  if (!output_node) {
    return true;
  }

  NodeGeometrySimulationOutput &storage = *static_cast<NodeGeometrySimulationOutput *>(
      output_node->storage);

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

}  // namespace blender::nodes::node_geo_simulation_input_cc

void register_node_type_geo_simulation_input()
{
  namespace file_ns = blender::nodes::node_geo_simulation_input_cc;

  static bNodeType ntype;
  geo_node_type_base(&ntype, GEO_NODE_SIMULATION_INPUT, "Simulation Input", NODE_CLASS_INTERFACE);
  ntype.initfunc = file_ns::node_init;
  ntype.declare_dynamic = file_ns::node_declare_dynamic;
  ntype.insert_link = file_ns::node_insert_link;
  node_type_storage(&ntype,
                    "NodeGeometrySimulationInput",
                    node_free_standard_storage,
                    node_copy_standard_storage);
  nodeRegisterType(&ntype);
}
