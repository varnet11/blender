/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2001-2023 Blender Foundation. */

#include "BKE_light_linking.h"

#include "DNA_ID.h"
#include "DNA_collection_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "BKE_collection.h"
#include "BKE_layer.h"
#include "BKE_lib_id.h"
#include "BKE_report.h"

#include "BLT_translation.h"

#include "DEG_depsgraph.h"
#include "DEG_depsgraph_build.h"

static LightLinking &light_linking_get(Object *object)
{
  return object->light_linking;
}

/* Un-assign light linking collection.
 *
 * Decreases the counter of the receiver collection.
 * If the light linking has no collection nothing happens.
 *
 * NOTE: It is up to the caller to tag the dependency graph for relation update and object for the
 * shading update. */
static void collection_unassign(LightLinking &light_linking)
{
  if (!light_linking.collection) {
    return;
  }

  id_us_min(&light_linking.collection->id);
  light_linking.collection = nullptr;
}

/* Unassign current light linking collection (if any) and assign the new one.
 *
 * The user counter is decreased for the old collection, and is increased for the.
 *
 * NOTE: It is up to the caller to tag the dependency graph for relation update and object for the
 * shading update. */
static void collection_assign(LightLinking &light_linking, Collection *new_collection)
{
  collection_unassign(light_linking);

  id_us_plus(&new_collection->id);
  light_linking.collection = new_collection;
}

Collection *BKE_light_linking_collection_new(struct Main *bmain, Object *object)
{
  LightLinking &light_linking = light_linking_get(object);

  Collection *new_collection = BKE_collection_add(
      bmain, nullptr, DATA_("Light Linking Collection"));

  collection_assign(light_linking, new_collection);

  DEG_id_tag_update(&object->id, ID_RECALC_COPY_ON_WRITE | ID_RECALC_SHADING);
  DEG_relations_tag_update(bmain);

  return new_collection;
}

bool BKE_light_linking_unlink_id_from_collection(Main *bmain,
                                                 Object *object,
                                                 ID *id,
                                                 ReportList *reports)
{
  LightLinking &light_linking = light_linking_get(object);

  Collection *collection = light_linking.collection;
  if (!collection) {
    BKE_reportf(
        reports, RPT_ERROR, "No light linking collection for object '%s'", object->id.name + 2);
    return false;
  }

  if (ID_IS_LINKED(&collection->id) || ID_IS_OVERRIDE_LIBRARY(&collection->id)) {
    BKE_reportf(reports,
                RPT_ERROR,
                "Cannot unlink '%s' from linked receiver collection '%s'",
                id->name + 2,
                collection->id.name + 2);
    return false;
  }

  const ID_Type id_type = GS(id->name);

  if (id_type == ID_OB) {
    BKE_collection_object_remove(bmain, collection, reinterpret_cast<Object *>(id), false);
  }
  else if (id_type == ID_GR) {
    BKE_collection_child_remove(bmain, collection, reinterpret_cast<Collection *>(id));
  }
  else {
    BKE_reportf(reports,
                RPT_ERROR,
                "Cannot unlink unsupported '%s' from light linking collection '%s'",
                id->name + 2,
                collection->id.name + 2);
    return false;
  }

  DEG_id_tag_update(&collection->id, ID_RECALC_HIERARCHY);

  DEG_relations_tag_update(bmain);

  return true;
}

void BKE_light_linking_select_receivers_of_emitter(Scene *scene,
                                                   ViewLayer *view_layer,
                                                   Object *emitter)
{
  LightLinking &light_linking = light_linking_get(emitter);

  Collection *collection = light_linking.collection;
  if (!collection) {
    return;
  }

  BKE_view_layer_synced_ensure(scene, view_layer);

  /* Deselect all currently selected objects in the view layer, but keep the emitter selected.
   * This is because the operation is called from the emitter being active, and it will be
   * confusing to deselect it but keep active. */
  LISTBASE_FOREACH (Base *, base, BKE_view_layer_object_bases_get(view_layer)) {
    if (base->object == emitter) {
      continue;
    }
    base->flag &= ~BASE_SELECTED;
  }

  /* Select objects which are reachable via the receiver collection hierarchy. */
  LISTBASE_FOREACH (CollectionObject *, cob, &collection->gobject) {
    Base *base = BKE_view_layer_base_find(view_layer, cob->ob);
    if (!base) {
      continue;
    }

    /* TODO(sergey): Check whether the object is configured to receive light. */

    base->flag |= BASE_SELECTED;
  }

  DEG_id_tag_update(&scene->id, ID_RECALC_SELECT);
}

void BKE_light_linking_receiver_to_emitter(Main *bmain, Object *emitter, Object *receiver)
{
  LightLinking &light_linking = light_linking_get(emitter);
  Collection *collection = light_linking.collection;

  if (!collection) {
    collection = BKE_light_linking_collection_new(bmain, emitter);
  }

  BKE_light_linking_receiver_to_collection(bmain, collection, &receiver->id);
}

/* Add object to the light linking collection and return corresponding CollectionLightLinking
 * settings.
 *
 * If the object is already in the collection then the content of the collection is not modified,
 * and the existing light linking settings are returned. */
static CollectionLightLinking *light_linking_collection_add_object(Main *bmain,
                                                                   Collection *collection,
                                                                   Object *object)
{
  BKE_collection_object_add(bmain, collection, object);

  LISTBASE_FOREACH (CollectionObject *, collection_object, &collection->gobject) {
    if (collection_object->ob == object) {
      return &collection_object->light_linking;
    }
  }

  BLI_assert_msg(0, "Object was not found after added to the light linking collection");

  return nullptr;
}

/* Add child collection to the light linking collection and return corresponding
 * CollectionLightLinking settings.
 *
 * If the child collection is already in the collection then the content of the collection is not
 * modified, and the existing light linking settings are returned. */
static CollectionLightLinking *light_linking_collection_add_collection(Main *bmain,
                                                                       Collection *collection,
                                                                       Collection *child)
{
  BKE_collection_child_add(bmain, collection, child);

  LISTBASE_FOREACH (CollectionChild *, collection_child, &collection->children) {
    if (collection_child->collection == child) {
      return &collection_child->light_linking;
    }
  }

  BLI_assert_msg(0, "Collection was not found after added to the light linking collection");

  return nullptr;
}

void BKE_light_linking_receiver_to_collection(Main *bmain, Collection *collection, ID *receiver)
{
  const ID_Type id_type = GS(receiver->name);

  CollectionLightLinking *collection_light_linking = nullptr;

  if (id_type == ID_OB) {
    collection_light_linking = light_linking_collection_add_object(
        bmain, collection, reinterpret_cast<Object *>(receiver));
  }
  else if (id_type == ID_GR) {
    collection_light_linking = light_linking_collection_add_collection(
        bmain, collection, reinterpret_cast<Collection *>(receiver));
  }
  else {
    return;
  }

  if (!collection_light_linking) {
    return;
  }

  collection_light_linking->light_state = COLLECTION_LIGHT_LINKING_STATE_INCLUDE;

  DEG_id_tag_update(&collection->id, ID_RECALC_HIERARCHY);
  DEG_id_tag_update(receiver, ID_RECALC_SHADING);

  DEG_relations_tag_update(bmain);
}
