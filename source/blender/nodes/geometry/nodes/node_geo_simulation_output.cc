/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "BLI_string_utils.h"

#include "BKE_compute_cache.hh"
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

}  // namespace blender::nodes

namespace blender::nodes::node_geo_simulation_output_cc {

NODE_STORAGE_FUNCS(NodeGeometrySimulationOutput);

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
  data->use_persistent_cache = false;
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

  /* TODO: All items. */
  if (std::optional<GeometrySet> value = cache.value_at_time("Geometry", time)) {
    params.set_output("Geometry", std::move(*value));
    params.set_input_unused("Geometry");
    return;
  }

  for (const NodeSimulationItem &item : Span(storage.items, storage.items_num)) {
    if (params.lazy_require_input(item.name)) {
      return;
    }
  }

  /* TODO: All items. */
  GeometrySet geometry_set = params.extract_input<GeometrySet>("Geometry");
  geometry_set.ensure_owns_direct_data();
  if (storage.use_persistent_cache) {
    if (!cache.is_empty()) {
      if (time.time < cache.start_time()->time) {
        cache.clear();
      }
    }
    /* If using the cache or there is no cached data yet, write the input in a new cache value. */
    cache.store_persistent("Geometry", time, geometry_set);
  }
  else {
    /* TODO: Maybe don't clear the whole cache here. */
    /* TODO: Move the geometry set here if the output isn't needed. */
    cache.clear();
    cache.store_temporary("Geometry", time, geometry_set);
  }

  params.set_output("Geometry", std::move(geometry_set));
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
  ntype.declare_dynamic = file_ns::node_declare_dynamic;
  ntype.insert_link = file_ns::node_insert_link;
  node_type_storage(&ntype,
                    "NodeGeometrySimulationOutput",
                    file_ns::node_free_storage,
                    file_ns::node_copy_storage);
  ntype.geometry_node_execute_supports_laziness = true;
  nodeRegisterType(&ntype);
}
