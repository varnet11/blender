/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "BKE_compute_cache.hh"
#include "BKE_compute_contexts.hh"
#include "BKE_scene.h"

#include "DEG_depsgraph_query.h"

#include "UI_interface.h"
#include "UI_resources.h"

#include "node_geometry_util.hh"

namespace blender::nodes::node_geo_simulation_input_cc {

NODE_STORAGE_FUNCS(NodeGeometrySimulationInput);

static void node_declare(NodeDeclarationBuilder &b)
{
  b.add_output<decl::Float>(N_("Delta Time"));

  const bNode &node = b.node();
  const NodeGeometrySimulationInput &storage = node_storage(node);
  const int32_t sim_output_node_id = storage.output_node_id;
  /* TODO: Add node tree and ndoe as arguments to new dynamic declaration function. */
  node.owner_tree().ensure_topology_cache();
  const bNode *sim_output_node = node.owner_tree().node_by_id(sim_output_node_id);
  if (!sim_output_node) {
    return;
  }
  const NodeGeometrySimulationOutput &output_storage =
      *static_cast<const NodeGeometrySimulationOutput *>(sim_output_node->storage);
  const Span<SimulationStateItem> items(output_storage.state_items,
                                        output_storage.state_items_num);

  for (const int i : items.index_range()) {
    const SimulationStateItem &item = items[i];
    switch (item.data_type) {
      case SOCK_FLOAT:
        b.add_input<decl::Float>(item.name).supports_field();
        b.add_output<decl::Float>(item.name).dependent_field({i});
        break;
      case SOCK_VECTOR:
        b.add_input<decl::Vector>(item.name).supports_field();
        b.add_output<decl::Vector>(item.name).dependent_field({i});
        break;
      case SOCK_RGBA:
        b.add_input<decl::Color>(item.name).supports_field();
        b.add_output<decl::Color>(item.name).dependent_field({i});
        break;
      case SOCK_BOOLEAN:
        b.add_input<decl::Bool>(item.name).supports_field();
        b.add_output<decl::Bool>(item.name).dependent_field({i});
        break;
      case SOCK_INT:
        b.add_input<decl::Int>(item.name).supports_field();
        b.add_output<decl::Int>(item.name).dependent_field({i});
        break;
      case SOCK_STRING:
        b.add_input<decl::String>(item.name);
        b.add_output<decl::String>(item.name);
        break;
      case SOCK_GEOMETRY:
        b.add_input<decl::Geometry>(item.name);
        b.add_output<decl::Geometry>(item.name);
        break;
    }
  }

  b.add_input<decl::Extend>("", N_("__extend__"));
  b.add_output<decl::Extend>("", N_("__extend__"));
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

static void node_geo_exec(GeoNodeExecParams params)
{
  const bNode &node = params.node();
  const NodeGeometrySimulationInput &storage = node_storage(node);
  const int32_t sim_output_node_id = storage.output_node_id;
  const bNode *sim_output_node = node.owner_tree().node_by_id(sim_output_node_id);
  if (!sim_output_node) {
    params.error_message_add(NodeWarningType::Error, TIP_("Missing simulation output node"));
    params.set_default_remaining_outputs();
    return;
  }
  const NodeGeometrySimulationOutput &output_storage =
      *static_cast<const NodeGeometrySimulationOutput *>(sim_output_node->storage);
  const Span<SimulationStateItem> items(output_storage.state_items,
                                        output_storage.state_items_num);

  const Scene *scene = DEG_get_input_scene(params.depsgraph());
  const float scene_ctime = BKE_scene_ctime_get(scene);
  const bke::sim::TimePoint time{int(scene_ctime), scene_ctime};

  const GeoNodesLFUserData &lf_data = *params.user_data();
  bke::sim::ComputeCaches &all_caches = *lf_data.modifier_data->cache_per_frame;

  const bke::NodeGroupComputeContext cache_context(lf_data.compute_context, sim_output_node_id);
  bke::sim::SimulationCache *cache = all_caches.lookup_context(cache_context.hash());

  if (params.lazy_output_is_required("Delta Time")) {
    if (cache) {
      const float time_diff = cache->is_empty() ? 0.0f :
                                                  scene_ctime - cache->last_run_time()->time;
      const double frame_rate = (double(scene->r.frs_sec) / double(scene->r.frs_sec_base));
      const float delta_time = float(std::max(0.0f, time_diff) / frame_rate);
      params.set_output("Delta Time", delta_time);
    }
    else {
      params.set_output("Delta Time", 0.0f);
    }
  }

  for (const SimulationStateItem &item : items) {
    /* TODO: Generic data type. */
    if (!cache) {
      params.set_output(item.name, params.extract_input<GeometrySet>(item.name));
      continue;
    }

    if (std::optional<GeometrySet> cached_value = cache->value_before_time(item.name, time)) {
      if (params.lazy_output_is_required(item.name)) {
        params.set_output(item.name, std::move(*cached_value));
      }
      continue;
    }

    if (params.lazy_require_input(item.name)) {
      continue;
    }

    GeometrySet geometry_set = params.extract_input<GeometrySet>(item.name);
    params.set_output(item.name, std::move(geometry_set));
  }
}

}  // namespace blender::nodes::node_geo_simulation_input_cc

void register_node_type_geo_simulation_input()
{
  namespace file_ns = blender::nodes::node_geo_simulation_input_cc;

  static bNodeType ntype;
  geo_node_type_base(&ntype, GEO_NODE_SIMULATION_INPUT, "Simulation Input", NODE_CLASS_INTERFACE);
  ntype.initfunc = file_ns::node_init;
  ntype.geometry_node_execute = file_ns::node_geo_exec;
  ntype.declare = file_ns::node_declare;
  node_type_storage(&ntype,
                    "NodeGeometrySimulationInput",
                    node_free_standard_storage,
                    node_copy_standard_storage);

  ntype.geometry_node_execute_supports_laziness = true;
  ntype.declaration_is_dynamic = true;
  nodeRegisterType(&ntype);
}
