/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2023 Blender Foundation. */

/** \file
 * \ingroup depsgraph
 */

#pragma once

#include <cstdint>

#include "BLI_map.hh"

struct Collection;
struct CollectionLightLinking;
struct Object;
struct Scene;

namespace blender::deg::light_linking {

namespace internal {

/* Set of light as seen from a receiver perspective. */
class LightSet {
 public:
  /* Maximum possible identifier of a light set. The identifier is 0-based.
   * The limitation is imposed by the fact that its identifier is converted to a bitmask. */
  static constexpr int MAX_ID = 63;

  /* Identifier of a light set which is not explicitly linked to anything. */
  static constexpr int DEFAULT_ID = 0;

  bool operator==(const LightSet &other) const;
  bool operator!=(const LightSet &other) const
  {
    return !(*this == other);
  }

  uint64_t hash() const;

  /* Lights which are explicitly included/excluded into the light set.
   *
   * The light is denoted as a bit mask of a light linking collection. This mask is allocated for
   * every unique light linking collection on an emitter. */
  uint64_t include_collection_mask;
  uint64_t exclude_collection_mask;
};

/* Packed information about emitter.
 * Emitter is actually corresponding to a light linking collection on an object. */
class EmitterData {
 public:
  /* Maximum possible identifier of a light linking collection. The identifier is 0-based.
   * The limitation is imposed by the fact that its identifier is converted to a bitmask. */
  static constexpr int MAX_COLLECTION_ID = 63;

  /* Bitmask which indicates the emitter belongs to all light sets. */
  static constexpr uint64_t SET_MEMBERSHIP_ALL = ~uint64_t(0);

  /* Get final emitter membership in the light sets, considering its inclusion and exclusion. */
  uint64_t get_set_membership() const;

  /* Mask of a light linking collection this emitter uses in its configuration.
   * A single bit is set in this bitfield which corresponds to an identifier of a light linking
   * collection in the scene. */
  uint64_t collection_mask = 0;

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
  /* Add receiver object with the given light linking configuration. */
  void add_receiver_object(const Scene &scene,
                           const Object &emitter,
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
   * The emitter must be original and have light linking collection.
   *
   * Note that there is limited number of emitters possible within a scene, When this number is
   * exceeded an error is printed and a nullptr is returned. */
  EmitterData *ensure_emitter_data_if_possible(const Scene &scene, const Object &emitter);

  /* Get emitter data for the given original or evaluated object.
   * If the light linking is not configured for this emitted nullptr is returned. */
  const EmitterData *get_emitter_data(const Object &emitter) const;

  /* Returns true if there is light linking configuration in the scene. */
  bool has_light_linking() const
  {
    return !emitter_data_map_.is_empty();
  }

  /* Receiver-centric view of light sets: indexed by an original receiver object, contains light
   * set which defines from which emitters it receives light from.
   *
   * NOTE: Only available during build. */
  Map<const Object *, LightSet> light_linked_sets_;

  /* Emitter-centric information: indexed by an original emitter object, contains accumulated
   * information about emitter. */
  Map<const Collection *, EmitterData> emitter_data_map_;

  /* Map from an original receiver object: map to index of light set for this receiver. */
  Map<const Object *, uint64_t> receiver_light_sets_;

  /* Next unique light linking collection ID. */
  uint64_t next_collection_id_ = 0;
};

}  // namespace blender::deg::light_linking
