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

  void execute_impl(lf::Params &params, const lf::Context &context) const final
  {
    const GeoNodesLFUserData &user_data = *dynamic_cast<const GeoNodesLFUserData *>(
        context.user_data);
    const Scene *scene = DEG_get_input_scene(user_data.modifier_data->depsgraph);
    const float scene_ctime = BKE_scene_ctime_get(scene);
    const bke::sim::TimePoint time{int(scene_ctime), scene_ctime};

    bke::sim::ComputeCaches &all_caches = *user_data.modifier_data->cache_per_frame;
    const bke::NodeGroupComputeContext cache_context(user_data.compute_context, node_id_);
    bke::sim::SimulationCache &cache = all_caches.ensure_for_context(cache_context.hash());

    /* Retrieve cached items for the current exact time. If all items are cached, the cache is
     * considered valid and no more evaluation has to be done. Otherwise it must be recreated. */
    bool all_items_cached = true;
    for (const int i : inputs_.index_range()) {
      const StringRef name = simulation_items_[i].name;
      const CPPType &type = *inputs_[i].type;
      std::optional<GeometrySet> value = cache.value_at_time(name, time);
      if (!value) {
        all_items_cached = false;
        continue;
      }
      void *output = params.get_output_data_ptr(i);
      type.move_construct(&*value, output);
      params.set_input_unused(i);
      params.output_set(i);
    }
    if (all_items_cached) {
      return;
    }

    for (const int i : inputs_.index_range()) {
      if (params.get_output_usage(i) != lf::ValueUsage::Used) {
        continue;
      }
      if (params.output_was_set(i)) {
        continue;
      }
      void *value = params.try_get_input_data_ptr_or_request(i);
      if (!value) {
        continue;
      }

      const StringRef name = simulation_items_[i].name;
      const CPPType &type = *inputs_[i].type;

      GeometrySet &geometry_set = *static_cast<GeometrySet *>(value);
      geometry_set.ensure_owns_direct_data();

      cache.remove(name);
      cache.store_temporary(name, time, geometry_set);

      void *output = params.get_output_data_ptr(i);
      type.move_construct(&geometry_set, output);
      params.output_set(i);
    }
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
