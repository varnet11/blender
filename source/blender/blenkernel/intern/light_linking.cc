/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2001-2023 Blender Foundation. */

#include "BKE_light_linking.h"

#include "DNA_ID.h"
#include "DNA_collection_types.h"
#include "DNA_object_types.h"

#include "BKE_collection.h"
#include "BKE_lib_id.h"
#include "BKE_report.h"

#include "BLT_translation.h"

#include "DEG_depsgraph.h"
#include "DEG_depsgraph_build.h"

static LightLinking &light_linking_get(Object *object)
{
  return object->light_linking;
}

/* Assign receiver collection from the given light linking.
 * Decreases the counter of the receiver collection.
 * If the light linking has no receiver collection nothing happens. */
static void receiver_collection_unassign(LightLinking &light_linking)
{
  if (!light_linking.receiver_collection) {
    return;
  }

  id_us_min(&light_linking.receiver_collection->id);
  light_linking.receiver_collection = nullptr;
}

/* Unassign current receiver collection (if any) and assign the new one.
 * The user counter is decreased for the old receiver collection, and is increased for the new
 * receiver collection. */
static void receiver_collection_assign(LightLinking &light_linking,
                                       Collection *new_receiver_collection)
{
  receiver_collection_unassign(light_linking);

  id_us_plus(&new_receiver_collection->id);
  light_linking.receiver_collection = new_receiver_collection;
}

Collection *BKE_light_linking_receiver_collection_new(struct Main *bmain, Object *object)
{
  LightLinking &light_linking = light_linking_get(object);

  Collection *new_receiver_collection = BKE_collection_add(
      bmain, nullptr, DATA_("Receiver Collection"));

  receiver_collection_assign(light_linking, new_receiver_collection);

  return new_receiver_collection;
}

bool BKE_light_linking_unlink_id_from_receiver_collection(Main *bmain,
                                                          Object *object,
                                                          ID *id,
                                                          ReportList *reports)
{
  LightLinking &light_linking = light_linking_get(object);

  Collection *receiver_collection = light_linking.receiver_collection;
  if (!receiver_collection) {
    if (reports) {
      BKE_reportf(
          reports, RPT_ERROR, "No light linking collection for object '%s'", object->id.name);
    }
    return false;
  }

  if (ID_IS_LINKED(&receiver_collection->id) || ID_IS_OVERRIDE_LIBRARY(&receiver_collection->id)) {
    BKE_reportf(reports,
                RPT_ERROR,
                "Cannot unlink '%s' from linked receiver collection '%s'",
                id->name + 2,
                receiver_collection->id.name + 2);

    return false;
  }

  const ID_Type id_type = GS(id->name);

  if (id_type == ID_OB) {
    BKE_collection_object_remove(
        bmain, receiver_collection, reinterpret_cast<Object *>(id), false);
  }
  else if (id_type == ID_GR) {
    BKE_collection_child_remove(bmain, receiver_collection, reinterpret_cast<Collection *>(id));
  }

  DEG_id_tag_update(&receiver_collection->id, ID_RECALC_HIERARCHY);

  DEG_relations_tag_update(bmain);

  return true;
}
