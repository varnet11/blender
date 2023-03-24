/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2023 Blender Foundation. */

/** \file
 * \ingroup depsgraph
 *
 * Light linking utilities. */

#include "intern/depsgraph_light_linking.h"

#include "MEM_guardedalloc.h"

#include "BLI_listbase.h"
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

void eval_runtime_data(const ::Depsgraph *depsgraph, Object *object_eval)
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

void Cache::clear()
{
  receiver_light_sets_.clear();
  emitter_data_map_.clear();
  next_receiver_collection_id_ = 0;
}

void Cache::add_emitter(const Scene *scene, const Object *emitter)
{
  BLI_assert(DEG_is_original_id(&emitter->id));

  const LightLinking &light_linking = emitter->light_linking;
  Collection *receiver_collection = light_linking.receiver_collection;

  if (!receiver_collection) {
    /* (Potential) emitter does not use light linking. */
    return;
  }

  if (emitter_data_map_.contains(receiver_collection)) {
    /* Receiver collection already handled. */
    return;
  }

  /* TODO: auto dedup collections with same content? */
  const uint64_t receiver_collection_id = next_receiver_collection_id_++;
  if (receiver_collection_id >= MAX_RELATION_BITS) {
    if (receiver_collection_id == MAX_RELATION_BITS) {
      printf("Maximum number of light linking collections (%d) exceeded scene \"%s\".",
             MAX_RELATION_BITS,
             scene->id.name + 2);
    }
    return;
  }

  /* Create emitter data for receiver collection, with unique bit. */
  const uint64_t receiver_collection_bit = 1 << receiver_collection_id;
  emitter_data_map_.add_new(receiver_collection, EmitterData(receiver_collection_bit));

  /* Add collection bit to all receivers affected by this emitter or any emitter with the same
   * receiver collection. */
  /* TODO: optimize by doing this non-recursively, and only going recursive in end_build()? */
  FOREACH_COLLECTION_OBJECT_RECURSIVE_BEGIN (receiver_collection, receiver) {
    /* If the object has receiver collection configure do not consider it as a receiver, avoiding
     * dependency cycles. */
    if (receiver->light_linking.receiver_collection == nullptr) {
      receiver_light_sets_.add_or_modify(
          receiver,
          [&](uint64_t *value) { *value = receiver_collection_bit; },
          [&](uint64_t *value) { *value |= receiver_collection_bit; });
    }
  }
  FOREACH_COLLECTION_OBJECT_RECURSIVE_END;
}

void Cache::end_build(const Scene *scene)
{
  /* Compute unique combinations of emitters used by receiver objects. */
  Map<uint64_t, uint> unique_light_sets;
  uint unique_light_set_id = 0;

  /* Default light set for any emitter without light linking. */
  unique_light_sets.add_new(0, unique_light_set_id++);

  /* TODO: parallelize for performance? */
  uint64_t last_mask = 0;
  uint64_t last_set = 0;
  for (uint64_t &receiver_collection_mask : receiver_light_sets_.values()) {
    /* Skip hash map lookup for consecutive objects with same mask for performance. */
    if (receiver_collection_mask != last_mask) {
      last_mask = receiver_collection_mask;
      last_set = unique_light_sets.lookup_or_add_cb(receiver_collection_mask, [&]() {
        /* Assign unique ID for this light set. */
        uint64_t new_id = unique_light_set_id++;
        if (new_id >= MAX_RELATION_BITS) {
          if (new_id == MAX_RELATION_BITS) {
            printf("Maximum number of light linking combinations (%d) exceeded in scene \"%s\".",
                   MAX_RELATION_BITS,
                   scene->id.name + 2);
          }
          return uint64_t(0);
        }

        /* Mark emitters as a member of this light set. */
        for (EmitterData &info : emitter_data_map_.values()) {
          if (receiver_collection_mask & info.receiver_collection_bit) {
            info.light_set_membership |= (1 << new_id);
          }
        }
        return new_id;
      });
    }

    /* Assign set index. */
    receiver_collection_mask = last_set;
  }
}

/* Set runtime data in light linking. */
void Cache::eval_runtime_data(Object *object_eval) const
{
  LightLinking &light_linking = object_eval->light_linking;

  if (receiver_light_sets_.size() == 0) {
    /* No light linking used in the scene. */
    light_linking.runtime.receiver_set = 0;
    light_linking.runtime.set_membership = 0;
    return;
  }

  const Object *object_orig = DEG_get_original_object(object_eval);

  const uint64_t set_membership_all = ~uint64_t(0);
  light_linking.runtime.receiver_set = receiver_light_sets_.lookup_default(object_orig,
                                                                           set_membership_all);

  const Collection *receiver_collection_eval = light_linking.receiver_collection;
  if (receiver_collection_eval) {
    /* TODO(sergey): Fix the const correctness of the depsgraph query API and avoid the const cast.
     */
    const Collection *receiver_collection_orig = reinterpret_cast<const Collection *>(
        DEG_get_original_id(const_cast<ID *>(&receiver_collection_eval->id)));
    light_linking.runtime.set_membership =
        emitter_data_map_.lookup(receiver_collection_orig).light_set_membership;
  }
  else {
    /* Member of all light sets. */
    light_linking.runtime.set_membership = set_membership_all;
  }
}

}  // namespace blender::deg::light_linking

/* \} */
