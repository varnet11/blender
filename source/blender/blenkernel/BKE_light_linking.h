/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2001-2023 Blender Foundation. */

#pragma once

/** \file
 * \ingroup bke
 *
 * API to manage light linking.
 */

#ifdef __cplusplus
extern "C" {
#endif

struct ID;
struct Main;
struct Object;
struct Collection;
struct ReportList;
struct Scene;
struct ViewLayer;

/* Create new collection and assign it as a receiver collection for the light linking configuration
 * of the given object.
 *
 * The collection is created outside of the view layer collections.
 * If the object has already receiver collection set up it is unreferenced from the object.
 *
 * Returns the newly created collection. */
struct Collection *BKE_light_linking_receiver_collection_new(struct Main *bmain,
                                                             struct Object *object);

/* Remove the given ID from the receiver collection of the given object.
 * The ID is expected to either be collection or an object.
 *
 * Returns true if the ID was unlinked from the receiver collection, false otherwise. The unlinking
 * will be unsuccessful if, for example, the receiver collection is a linked data-block.
 *
 * The optional reports argument is used to provide human-readable details about why unlinking was
 * not successful. */
bool BKE_light_linking_unlink_id_from_receiver_collection(struct Main *bmain,
                                                          struct Object *object,
                                                          struct ID *id,
                                                          struct ReportList *reports);

/* Select all objects which receive the light from the given emitter via the light linking
 * configuration. */
void BKE_light_linking_select_receivers_of_emitter(struct Scene *scene,
                                                   struct ViewLayer *view_layer,
                                                   struct Object *emitter);

/* Link receiver object to the given emitter.
 * if the emitter already has receiver collection specified the object is added to that collection.
 * Otherwise, first a new collection is created and assigned as a receiver collection, and the
 * receiver is added to it. */
void BKE_light_linking_receiver_to_emitter(struct Main *bmain,
                                           struct Object *emitter,
                                           struct Object *receiver);

/* Add receiver to the given receiver collection.
 * The ID is expected to either be collection or an object. */
void BKE_light_linking_receiver_to_collection(struct Main *bmain,
                                              struct Collection *receiver_collection,
                                              struct ID *receiver);

#ifdef __cplusplus
}
#endif
