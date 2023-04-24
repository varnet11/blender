/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2013 Blender Foundation */

/** \file
 * \ingroup depsgraph
 */

#include "intern/node/deg_node_time.h"

#include "DNA_scene_types.h"

#include "intern/depsgraph.h"
#include "intern/depsgraph_relation.h"
#include "intern/node/deg_node_factory.h"

namespace blender::deg {

const char *timeSourceTypeAsString(eTimeSourceType source_type)
{
  switch (source_type) {
    case eTimeSourceType::DEG_TIME_SOURCE_SCENE:
      return "SCENE";
    case eTimeSourceType::DEG_TIME_SOURCE_REALTIME:
      return "REALTIME";
  }
  BLI_assert_msg(0, "Unhandled time source type, should never happen.");
  return "UNKNOWN";
}

TimeSourceNode::TimeSourceNode() : source_type(eTimeSourceType::DEG_TIME_SOURCE_SCENE)
{
}

string TimeSourceNode::identifier() const
{
  const string type_name = type_get_factory(type)->type_name();
  const string name_part = name[0] ? (string(" '") + name + "'") : "";

  return "[" + type_name + "]" + name_part + " : " +
         "(source_type: " + timeSourceTypeAsString(source_type) + ")";
}

void TimeSourceNode::tag_update(Depsgraph * /*graph*/, eUpdateSource /*source*/)
{
  tagged_for_update = true;
}

void TimeSourceNode::flush_update_tag(Depsgraph *graph)
{
  if (!tagged_for_update) {
    return;
  }
  for (Relation *rel : outlinks) {
    Node *node = rel->to;
    node->tag_update(graph, DEG_UPDATE_SOURCE_TIME);
  }
}

}  // namespace blender::deg
