/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2013 Blender Foundation. All rights reserved. */

/** \file
 * \ingroup depsgraph
 */

#include "intern/node/deg_node_time.h"

#include "DNA_scene_types.h"

#include "intern/depsgraph.h"
#include "intern/depsgraph_relation.h"

namespace blender::deg {

const char *timeSourceTypeAsString(eTimeSourceType source_type)
{
  switch (source_type) {
    case eTimeSourceType::DEG_TIME_SOURCE_SCENE:
      return "SCENE_TIMELINE";
    case eTimeSourceType::DEG_TIME_SOURCE_REALTIME:
      return "REALTIME_CLOCK";
  }
  BLI_assert_msg(0, "Unhandled time source type, should never happen.");
  return "UNKNOWN";
}

TimeSourceNode::TimeSourceNode() : source_type(eTimeSourceType::DEG_TIME_SOURCE_SCENE)
{
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
