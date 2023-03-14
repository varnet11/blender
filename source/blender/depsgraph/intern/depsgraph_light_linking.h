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

struct Object;
struct Collection;

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
  void add_emitter(const Object *emitter);

  /* Get emitters which has light linking configured to affect the given receiver.
   * The receiver can either be evaluated or original. The function returns original objects. */
  Span<const Object *> get_original_emitters_for_receiver(const Object *receiver) const;

  /* Get pre-calculated emission mask of the given emitter.
   *
   * The mask is a bit-field with a single bit set, which corresponds to an identifier of the
   * emitter within the evaluation context.
   *
   * If the given emitter does not have light linking configured (the receiver collection is
   * nullptr) the function returns 0.
   *
   * The emitter can either be evaluated or original. The function returns mask which is only valid
   * within the given dependency graph. */
  uint64_t get_emitter_mask(const Object *emitter) const;

 private:
  /* Collection of emitters, stored in an efficient for traversal manner. */
  using Emitters = Vector<const Object *>;

  /* Add emitter to the cached lookup for the emitters-of-receiver.
   *
   * The receiver_collection can be nullptr, in which case the function does nothing.
   * The collection is non-const because of the BKE API which receivers it as such for the cached
   * object iteration. This call itself does not modify the collection. */
  void add_emitter_to_receivers(const Object *emitter, /*const*/ Collection *receiver_collection);

  /* Cached lookup of emitters affecting a receiver.
   *
   * The map is indexed by an original receiver object, contains a collection of original emitters
   * which has light linking configured to affect the receiver. */
  Map<const Object *, Emitters> emitters_of_receivers_;

  /* A map from an emitter to its emitter_mask. Indexed by an original emitter object. */
  Map<const Collection *, uint64_t> receiver_collection_emitter_mask_map_;

  /* An identifier of an emitter which will be used for the next emitter added to the cache. */
  uint64_t next_unused_id_{0};
};

}  // namespace blender::deg::light_linking
