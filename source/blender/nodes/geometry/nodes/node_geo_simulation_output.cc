/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "BKE_compute_cache.hh"
#include "BKE_compute_contexts.hh"
#include "BKE_scene.h"

#include "DEG_depsgraph_query.h"

#include "UI_interface.h"
#include "UI_resources.h"

#include "node_geometry_util.hh"

namespace blender::nodes::node_geo_simulation_output_cc {

NODE_STORAGE_FUNCS(NodeGeometrySimulationOutput);

static void node_declare(NodeDeclarationBuilder &b)
{
  const NodeGeometrySimulationOutput &storage = node_storage(b.node());
  const Span<SimulationStateItem> items(storage.state_items, storage.state_items_num);
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

static void node_init(bNodeTree * /*tree*/, bNode *node)
{
  NodeGeometrySimulationOutput *data = MEM_cnew<NodeGeometrySimulationOutput>(__func__);
  data->state_items = MEM_cnew_array<SimulationStateItem>(1, __func__);
  data->state_items[0].name = BLI_strdup(DATA_("Geometry"));
  data->state_items[0].data_type = SOCK_GEOMETRY;
  data->state_items_num = 1;
  data->use_persistent_cache = false;
  node->storage = data;
}

static void node_free_storage(bNode *node)
{
  NodeGeometrySimulationOutput &storage = node_storage(*node);
  MEM_SAFE_FREE(storage.state_items);
}

void node_copy_storage(bNodeTree * /*dest_ntree*/, bNode *dst_node, const bNode *src_node)
{
  const NodeGeometrySimulationOutput &src = node_storage(*src_node);
  NodeGeometrySimulationOutput &dst = node_storage(*dst_node);
  MEM_SAFE_FREE(dst.state_items);
  dst.state_items = MEM_cnew_array<SimulationStateItem>(src.state_items_num, __func__);
  dst.state_items_num = src.state_items_num;
}

static void node_geo_exec(GeoNodeExecParams params)
{
  const bNode &node = params.node();
  const NodeGeometrySimulationOutput &storage = node_storage(node);
  const Scene *scene = DEG_get_input_scene(params.depsgraph());
  const float scene_ctime = BKE_scene_ctime_get(scene);
  const bke::sim::TimePoint time{int(scene_ctime), scene_ctime};

  const GeoNodesLFUserData &lf_data = *params.user_data();
  bke::sim::ComputeCaches &all_caches = *lf_data.modifier_data->cache_per_frame;
  const bke::NodeGroupComputeContext cache_context(lf_data.compute_context, node.identifier);
  bke::sim::SimulationCache &cache = all_caches.ensure_for_context(cache_context.hash());

  for (const SimulationStateItem &item : Span(storage.state_items, storage.state_items_num)) {
    /* TODO: Generic data types. */
    if (std::optional<GeometrySet> value = cache.value_at_time(item.name, time)) {
      params.set_output(item.name, std::move(*value));
      params.set_input_unused(item.name);
      continue;
    }

    if (params.lazy_require_input(item.name)) {
      continue;
    }

    GeometrySet geometry_set = params.extract_input<GeometrySet>(item.name);
    geometry_set.ensure_owns_direct_data();
    if (storage.use_persistent_cache) {
      if (!cache.is_empty()) {
        if (time.time < cache.start_time()->time) {
          cache.clear();
        }
      }
      /* If using the cache or there is no cached data yet, write the input in a new cache
       * value.
       */
      cache.store_persistent(item.name, time, geometry_set);
    }
    else {
      /* TODO: Maybe don't clear the whole cache here. */
      /* TODO: Move the geometry set here if the output isn't needed. */
      cache.clear();
      cache.store_temporary(item.name, time, geometry_set);
    }
    params.set_output(item.name, std::move(geometry_set));
  }
}

}  // namespace blender::nodes::node_geo_simulation_output_cc

void register_node_type_geo_simulation_output()
{
  namespace file_ns = blender::nodes::node_geo_simulation_output_cc;

  static bNodeType ntype;

  geo_node_type_base(
      &ntype, GEO_NODE_SIMULATION_OUTPUT, "Simulation Output", NODE_CLASS_INTERFACE);
  ntype.initfunc = file_ns::node_init;
  ntype.geometry_node_execute = file_ns::node_geo_exec;
  ntype.declare = file_ns::node_declare;
  node_type_storage(&ntype,
                    "NodeGeometrySimulationOutput",
                    file_ns::node_free_storage,
                    file_ns::node_copy_storage);
  ntype.geometry_node_execute_supports_laziness = true;
  ntype.declaration_is_dynamic = true;
  nodeRegisterType(&ntype);
}
