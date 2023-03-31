/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2023 Blender Foundation */

/** \file
 * \ingroup depsgraph
 *
 * Light linking utilities. */

#include "intern/depsgraph_light_linking.h"

#include "MEM_guardedalloc.h"

#include "BLI_hash.hh"
#include "BLI_listbase.h"
#include "BLI_map.hh"
#include "BLI_utildefines.h"

#include "BKE_collection.h"

#include "DNA_collection_types.h"
#include "DNA_layer_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "DEG_depsgraph_light_linking.h"
#include "DEG_depsgraph_light_linking.hh"
#include "DEG_depsgraph_query.h"

#include "intern/depsgraph.h"

namespace deg = blender::deg;

/* -------------------------------------------------------------------- */
/** \name Public C++ API
 * \{ */

namespace blender::deg::light_linking {

void eval_runtime_data(const ::Depsgraph *depsgraph, Object &object_eval)
{
  const deg::Depsgraph *deg_graph = reinterpret_cast<const deg::Depsgraph *>(depsgraph);
  deg_graph->light_linking_cache.eval_runtime_data(object_eval);
}

}  // namespace blender::deg::light_linking

/* \} */

/* -------------------------------------------------------------------- */
/** \name Internal builder API
 * \{ */

namespace blender::deg::light_linking {

namespace internal {

bool LightSet::operator==(const LightSet &other) const
{
  return include_collection_mask == other.include_collection_mask &&
         exclude_collection_mask == other.exclude_collection_mask;
}

uint64_t LightSet::hash() const
{
  return get_default_hash_2(get_default_hash(include_collection_mask),
                            get_default_hash(exclude_collection_mask));
}

uint64_t EmitterData::get_set_membership() const
{
  const uint64_t effective_included_mask = included_sets_mask ? included_sets_mask :
                                                                SET_MEMBERSHIP_ALL;
  return effective_included_mask & ~excluded_sets_mask;
}

}  // namespace internal

namespace {

/* TODO(sergey): Move to a public API, solving the const-correctness. */
template<class T> static inline const T *get_original(const T *id)
{
  if (!id) {
    return nullptr;
  }
  return reinterpret_cast<T *>(DEG_get_original_id(const_cast<ID *>(&id->id)));
}

/* Iterate over all objects of the collection and invoke the given callback with two arguments:
 * the given collection light linking settings, and the object (passed as reference).
 *
 * Note that if an object is reachable from multiple children collection the callback is invoked
 * for all of them. */
template<class Proc>
static void foreach_light_collection_object_inner(
    const CollectionLightLinking &collection_light_linking,
    const Collection &collection,
    Proc &&callback)
{
  LISTBASE_FOREACH (const CollectionChild *, collection_child, &collection.children) {
    foreach_light_collection_object_inner(
        collection_light_linking, *collection_child->collection, callback);
  }

  LISTBASE_FOREACH (const CollectionObject *, collection_object, &collection.gobject) {
    callback(collection_light_linking, *collection_object->ob);
  }
}

/* Iterate over all objects of the collection and invoke the given callback with two arguments:
 * CollectionLightLinking and the actual Object (passed as reference).
 *
 * The CollectionLightLinking denotes the effective light linking settings of the object. It comes
 * from the first level of hierarchy from the given collection.
 *
 * Note that if an object is reachable from multiple children collection the callback is invoked
 * for all of them. */
template<class Proc>
static void foreach_light_collection_object(const Collection &collection, Proc &&callback)
{
  LISTBASE_FOREACH (const CollectionChild *, collection_child, &collection.children) {
    foreach_light_collection_object_inner(
        collection_child->light_linking, *collection_child->collection, callback);
  }

  LISTBASE_FOREACH (const CollectionObject *, collection_object, &collection.gobject) {
    callback(collection_object->light_linking, *collection_object->ob);
  }
}

/* Helper class which takes care of allocating an unique light set IDs, performing checks for
 * overflows. */
class LightSetIDManager {
  using LightSet = internal::LightSet;

 public:
  explicit LightSetIDManager(const Scene &scene) : scene_(scene) {}

  bool get(const LightSet &light_set, uint64_t &id)
  {
    /* Performance note.
     *
     * Always ensure the light set data exists in the map, even when an overflow happens. This has
     * a downside of potentially higher memory usage and when there are many emitters with light
     * linking, but it avoids distinct lookup + add fore the normal cases. */

    const uint64_t light_set_id = light_set_id_map_.lookup_or_add_cb(light_set, [&]() {
      const uint64_t new_light_set_id = next_light_set_id_++;

      if (new_light_set_id == LightSet::MAX_ID + 1) {
        printf("Maximum number of light linking sets (%d) exceeded scene \"%s\".\n",
               LightSet::MAX_ID + 1,
               scene_.id.name + 2);
      }

      return new_light_set_id;
    });

    id = light_set_id;

    return id <= LightSet::MAX_ID;
  }

 private:
  const Scene &scene_;

  /* Next unique ID of a light set. */
  uint64_t next_light_set_id_ = LightSet::DEFAULT_ID + 1;

  /* Map from a link set to its original id. */
  Map<internal::LightSet, uint64_t> light_set_id_map_;
};

}  // namespace

void Cache::clear()
{
  light_linked_sets_.clear();
  emitter_data_map_.clear();
  receiver_light_sets_.clear();
  next_collection_id_ = 0;
}

void Cache::add_emitter(const Scene &scene, const Object &emitter)
{
  BLI_assert(DEG_is_original_id(&emitter.id));

  if (can_skip_emitter(emitter)) {
    return;
  }

  const LightLinking &light_linking = emitter.light_linking;
  Collection *collection = light_linking.collection;

  /* Add collection bit to all receivers affected by this emitter or any emitter with the same
   * receiver collection. */
  /* TODO: optimize by doing this non-recursively, and only going recursive in end_build()? */
  foreach_light_collection_object(
      *collection,
      [&](const CollectionLightLinking &collection_light_linking, const Object &receiver) {
        add_receiver_object(scene, emitter, collection_light_linking, receiver);
      });
}

void Cache::add_receiver_object(const Scene &scene,
                                const Object &emitter,
                                const CollectionLightLinking &collection_light_linking,
                                const Object &receiver)
{
  BLI_assert(DEG_is_original_id(&emitter->id));
  BLI_assert(DEG_is_original_id(&receiver->id));

  if (receiver.light_linking.collection != nullptr) {
    /* If the object has receiver collection configure do not consider it as a receiver, avoiding
     * dependency cycles. */
    return;
  }

  /* Light linking membership. */
  if (collection_light_linking.light_state != COLLECTION_LIGHT_LINKING_STATE_NONE) {
    const EmitterData *emitter_data = ensure_emitter_data_if_possible(scene, emitter);

    if (emitter_data) {
      LightSet &light_set = ensure_light_linked_set(receiver);

      switch (collection_light_linking.light_state) {
        case COLLECTION_LIGHT_LINKING_STATE_NONE:
          break;

        case COLLECTION_LIGHT_LINKING_STATE_INCLUDE:
          light_set.include_collection_mask |= emitter_data->collection_mask;
          light_set.exclude_collection_mask &= ~emitter_data->collection_mask;
          break;

        case COLLECTION_LIGHT_LINKING_STATE_EXCLUDE:
          light_set.exclude_collection_mask |= emitter_data->collection_mask;
          light_set.include_collection_mask &= ~emitter_data->collection_mask;
          break;
      }
    }
  }

  /* TODO(sergey): Handle shadow linking. */
}

void Cache::end_build(const Scene &scene)
{
  if (!has_light_linking()) {
    return;
  }

  LightSetIDManager light_set_id_manager(scene);

  for (const auto it : light_linked_sets_.items()) {
    const Object *receiver = it.key;
    LightSet &light_set = it.value;

    uint64_t light_set_id;
    if (!light_set_id_manager.get(light_set, light_set_id)) {
      continue;
    }

    const uint64_t light_set_mask = uint64_t(1) << light_set_id;

    receiver_light_sets_.add(receiver, light_set_id);

    for (EmitterData &emitter_data : emitter_data_map_.values()) {
      if (emitter_data.collection_mask & light_set.include_collection_mask) {
        emitter_data.included_sets_mask |= light_set_mask;
      }
      if (emitter_data.collection_mask & light_set.exclude_collection_mask) {
        emitter_data.excluded_sets_mask |= light_set_mask;
      }
    }
  }

  light_linked_sets_.clear_and_shrink();
}

/* Set runtime data in light linking. */
void Cache::eval_runtime_data(Object &object_eval) const
{
  LightLinking &light_linking = object_eval.light_linking;

  if (!has_light_linking()) {
    /* No light linking used in the scene. */
    light_linking.runtime.receiver_set = 0;
    light_linking.runtime.set_membership = 0;
    return;
  }

  /* Receiver configuration. */
  const Object *object_orig = get_original(&object_eval);
  light_linking.runtime.receiver_set = receiver_light_sets_.lookup_default(object_orig,
                                                                           LightSet::DEFAULT_ID);

  /* Emitter configuration. */
  const EmitterData *emitter_data = get_emitter_data(object_eval);
  if (emitter_data) {
    light_linking.runtime.set_membership = emitter_data->get_set_membership();
  }
  else {
    light_linking.runtime.set_membership = EmitterData::SET_MEMBERSHIP_ALL;
  }
}

internal::LightSet &Cache::ensure_light_linked_set(const Object &receiver)
{
  BLI_assert(DEG_is_original_id(&receiver.id));

  return light_linked_sets_.lookup_or_add_as(&receiver);
}

bool Cache::can_skip_emitter(const Object &emitter) const
{
  BLI_assert(DEG_is_original_id(&emitter.id));

  const LightLinking &light_linking = emitter.light_linking;
  const Collection *collection = light_linking.collection;

  if (!collection) {
    return true;
  }

  return emitter_data_map_.contains(collection);
}

internal::EmitterData *Cache::ensure_emitter_data_if_possible(const Scene &scene,
                                                              const Object &emitter)
{
  BLI_assert(DEG_is_original_id(&emitter.id));
  BLI_assert(emitter.light_linking.collection);

  const LightLinking &light_linking = emitter.light_linking;
  const Collection *collection = light_linking.collection;

  /* Performance note.
   *
   * Always ensure the emitter data exists in the map, even when an overflow happens. This has a
   * downside of potentially higher memory usage when there are many emitters with light linking,
   * but it avoids distinct lookup + add fore the normal cases.
   *
   * On the API level the function always returns nullptr on overflow, so it is more of an internal
   * behavior. */

  EmitterData &emitter_data = emitter_data_map_.lookup_or_add_cb(collection, [&]() {
    const uint64_t collection_id = next_collection_id_++;

    EmitterData new_emitter_data;

    if (collection_id > EmitterData::MAX_COLLECTION_ID) {
      if (collection_id == EmitterData::MAX_COLLECTION_ID + 1) {
        printf("Maximum number of light linking collections (%d) exceeded in scene \"%s\".\n",
               EmitterData::MAX_COLLECTION_ID + 1,
               scene.id.name + 2);
      }
      new_emitter_data.collection_mask = 0;
    }
    else {
      new_emitter_data.collection_mask = uint64_t(1) << collection_id;
    }

    return new_emitter_data;
  });

  if (!emitter_data.collection_mask) {
    return nullptr;
  }

  return &emitter_data;
}

const internal::EmitterData *Cache::get_emitter_data(const Object &emitter) const
{
  const LightLinking &light_linking = emitter.light_linking;
  const Collection *collection_eval = light_linking.collection;

  if (!collection_eval) {
    return nullptr;
  }

  const Collection *collection_orig = get_original(collection_eval);

  return emitter_data_map_.lookup_ptr(collection_orig);
}

}  // namespace blender::deg::light_linking

/* \} */
