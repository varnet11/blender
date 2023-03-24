/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2023 Blender Foundation. */

/** \file
 * \ingroup depsgraph
 */

#pragma once

#include <cstdint>

#include "BLI_map.hh"
#include "BLI_span.hh"
#include "BLI_vector.hh"

struct Collection;
struct Object;
struct Scene;

namespace blender::deg::light_linking {

/* Cached light linking evaluation data.
 *
 * This cache is only valid within a specific dependency graph, hence the dependency graph is
 * expected to own this cache.
 *
 * This cache takes care of making it efficient to lookup emitter masks, emitters which affect
 * given receiver and so on. */
class Cache {
 public:
  /* Entirely clear the cache.
   * Should be called whenever the dependency graph is being re-built, in the beginning of the
   * build process. */
  void clear();

  /* Add emitter to the cache.
   *
   * This call does nothing if the emitter does not have light configured linking (as in, if it
   * has receiver_collection set to nullptr).
   *
   * The emitter must be original. This is asserted, but in release builds passing evaluated
   * object leads to an undefined behavior. */
  void add_emitter(const Scene *scene, const Object *emitter);

  /* Compute unique sets of emitters used by receivers.
   *
   * This must be called at the end of depsgraph relations build after all emitters have been
   * added, and before runtime data can be set as part of evaluation. */
  void end_build(const Scene *scene);

  /* Set runtime light linking data on evaluated object. */
  void eval_runtime_data(Object *object_eval) const;

 private:
  /* Maximum number of bits available for light linking relations. */
  const int MAX_RELATION_BITS = 64;

  struct EmitterData {
    EmitterData(const uint64_t receiver_collection_bit)
        : receiver_collection_bit(receiver_collection_bit)
    {
    }

    /* Unique bit for receiver collection, added to receiver objects during build
     * to find which receiver collections affect them. */
    uint64_t receiver_collection_bit;
    /* Bitmask indicating which unique light sets an emitter object is part of. */
    uint64_t light_set_membership = 0;
  };

  /* Map from a receiver collection to data for any emitters using it. */
  Map<const Collection *, EmitterData> emitter_data_map_;

  /* Next unique receiver collection ID. */
  uint64_t next_receiver_collection_id_{0};

  /* Map from an original receiver object.
   *
   * During build: map to bitmask of receiver collections affecting this receiver.
   * After build: map to index of light set for this receiver. */
  Map<const Object *, uint64_t> receiver_light_sets_;
};

}  // namespace blender::deg::light_linking
