/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2023 Blender Foundation. */

/** \file
 * \ingroup depsgraph
 */

#pragma once

#include <cstdint>

#include "BLI_map.hh"
#include "BLI_set.hh"

struct Collection;
struct CollectionLightLinking;
struct Object;
struct Scene;

namespace blender::deg::light_linking {

namespace internal {

/* Set of light as seen from a receiver perspective. */
class LightSet {
 public:
  /* Identifier of a light set which is not explicitly linked to anything. */
  static constexpr int DEFAULT_ID = 0;

  bool operator==(const LightSet &other) const;
  bool operator!=(const LightSet &other) const
  {
    return !(*this == other);
  }

  uint64_t hash() const;

  /* Lights which are explicitly included/excluded into the light set. */
  Set<const Object *> include;
  Set<const Object *> exclude;
};

/* Packed information about emitter. */
class EmitterData {
 public:
  static constexpr uint64_t SET_MEMBERSHIP_ALL = ~uint64_t(0);

  EmitterData() = default;

  EmitterData(EmitterData &&other) noexcept = default;
  EmitterData &operator=(EmitterData &&other) = default;

  EmitterData(const EmitterData &other) = delete;
  EmitterData &operator=(const EmitterData &other) = delete;

  /* Get final emitter membership in the light sets, considering its inclusion and exclusion. */
  uint64_t get_set_membership() const;

  /* Bit masks of the emitter membership in the light sets. */
  uint64_t included_sets_mask = 0;
  uint64_t excluded_sets_mask = 0;
};

}  // namespace internal

/* Cached light linking evaluation data.
 *
 * This cache is only valid within a specific dependency graph, hence the dependency graph is
 * expected to own this cache.
 *
 * This cache takes care of making it efficient to lookup emitter masks, emitters which affect
 * given receiver and so on. */
class Cache {
  using LightSet = internal::LightSet;
  using EmitterData = internal::EmitterData;

 public:
  /* Entirely clear the cache.
   * Should be called whenever the dependency graph is being re-built, in the beginning of the
   * build process. */
  void clear();

  /* Add emitter to the cache.
   *
   * This call does nothing if the emitter does not have light configured linking (as in, if it
   * has light linking collection set to nullptr).
   *
   * The emitter must be original. This is asserted, but in release builds passing evaluated
   * object leads to an undefined behavior. */
  void add_emitter(const Scene &scene, const Object &emitter);

  /* Compute unique sets of emitters used by receivers.
   *
   * This must be called at the end of depsgraph relations build after all emitters have been
   * added, and before runtime data can be set as part of evaluation. */
  void end_build(const Scene &scene);

  /* Set runtime light linking data on evaluated object. */
  void eval_runtime_data(Object &object_eval) const;

 private:
  /* Maximum number of bits available for light linking relations. */
  const int MAX_RELATION_BITS = 64;

  /* Add receiver object with the given light linking configuration. */
  void add_receiver_object(const Object &emitter,
                           const CollectionLightLinking &collection_light_linking,
                           const Object &receiver);

  /* Ensure that the light set exists for the given receiver.
   * The receiver must be original. */
  LightSet &ensure_light_linked_set(const Object &receiver);

  /* Returns true if the underlying data of the light linking emitter has been handled, and there
   * is no need to handle the emitter.
   * The emitter must be original object. */
  bool can_skip_emitter(const Object &emitter) const;

  /* Ensure that the data exists for the given emitter.
   * The emitter must be original and have light linking collection. */
  EmitterData &ensure_emitter_data(const Object &emitter);

  /* Get emitter data for the given original or evaluated object.
   * If the light linking is not configured for this emitted nullptr is returned. */
  const EmitterData *get_emitter_data(const Object &emitter) const;

  /* Returns true if there is light linking configuration in the scene. */
  bool has_light_linking() const
  {
    /* Check both collections, as the former one is only non-empty during the build, and the latter
     * one after the build. */
    return !light_linked_sets_.is_empty() || !emitter_data_map_.is_empty();
  }

  /* Receiver-centric view of light sets: indexed by an original receiver object, contains light
   * set which defines from which emitters it receives light.
   *
   * NOTE: Only available during build. */
  Map<const Object *, LightSet> light_linked_sets_;

  /* Emitter-centric information: indexed by an original emitter object, contains accumulated
   * information about emitter. */
  Map<const Collection *, EmitterData> emitter_data_map_;

  /* Map from an original receiver object: map to index of light set for this receiver. */
  Map<const Object *, uint64_t> receiver_light_sets_;
};

}  // namespace blender::deg::light_linking
