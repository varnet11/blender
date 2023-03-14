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

#include "DEG_depsgraph_light_linking.h"
#include "DEG_depsgraph_query.h"

#include "intern/depsgraph.h"

namespace deg = blender::deg;

/* -------------------------------------------------------------------- */
/** \name Public C++ API
 * \{ */

namespace blender::deg::light_linking {

Span<const Object *> get_original_emitters_for_receiver(const ::Depsgraph *depsgraph,
                                                        const Object *receiver)
{
  const deg::Depsgraph *deg_graph = reinterpret_cast<const deg::Depsgraph *>(depsgraph);
  return deg_graph->light_linking_cache.get_original_emitters_for_receiver(receiver);
}

uint64_t get_emitter_mask(const ::Depsgraph *depsgraph, const Object *emitter)
{
  const deg::Depsgraph *deg_graph = reinterpret_cast<const deg::Depsgraph *>(depsgraph);
  return deg_graph->light_linking_cache.get_emitter_mask(emitter);
}

}  // namespace blender::deg::light_linking

/* \} */

/* -------------------------------------------------------------------- */
/** \name Internal builder API
 * \{ */

namespace blender::deg::light_linking {

void Cache::clear()
{
  emitters_of_receivers_.clear();
  receiver_collection_emitter_mask_map_.clear();
  next_unused_id_ = 0;
}

void Cache::add_emitter(const Object *emitter)
{
  BLI_assert(DEG_is_original_id(&emitter->id));

  const LightLinking &light_linking = emitter->light_linking;
  Collection *receiver_collection = light_linking.receiver_collection;

  if (!receiver_collection) {
    /* The given object is not really an emitter.
     * At least, not an emitter which requires special caching. */
    return;
  }

  add_emitter_to_receivers(emitter, receiver_collection);

  receiver_collection_emitter_mask_map_.lookup_or_add_cb(receiver_collection, [&]() {
    const uint64_t id = next_unused_id_++;
    return uint64_t(1) << id;
  });
}

void Cache::add_emitter_to_receivers(const Object *emitter,
                                     /*const*/ Collection *receiver_collection)
{
  if (!receiver_collection) {
    return;
  }

  FOREACH_COLLECTION_OBJECT_RECURSIVE_BEGIN (receiver_collection, receiver) {
    /* If the object has receiver collection configure do not consider it as a receiver, avoiding
     * dependency cycles. */
    if (receiver->light_linking.receiver_collection == nullptr) {
      Emitters &emitters_of_receiver = emitters_of_receivers_.lookup_or_add_cb(
          receiver, []() { return Emitters(); });
      emitters_of_receiver.append(emitter);
    }
  }
  FOREACH_COLLECTION_OBJECT_RECURSIVE_END;
}

Span<const Object *> Cache::get_original_emitters_for_receiver(const Object *receiver) const
{
  /* TODO(sergey): Fix the const correctness of the depsgraph query API and avoid the cons cast. */
  const Object *receiver_orig = DEG_get_original_object(const_cast<Object *>(receiver));

  return emitters_of_receivers_.lookup_default(receiver_orig, Span<const Object *>{});
}

uint64_t Cache::get_emitter_mask(const Object *emitter) const
{
  const LightLinking &light_linking = emitter->light_linking;

  const Collection *receiver_collection = light_linking.receiver_collection;
  if (!receiver_collection) {
    return 0;
  }

  /* TODO(sergey): Fix the const correctness of the depsgraph query API and avoid the cons cast. */
  const Collection *receiver_collection_orig = reinterpret_cast<const Collection *>(
      DEG_get_original_id(const_cast<ID *>(&receiver_collection->id)));

  return receiver_collection_emitter_mask_map_.lookup_default(receiver_collection_orig, 0);
}

}  // namespace blender::deg::light_linking

/* \} */
