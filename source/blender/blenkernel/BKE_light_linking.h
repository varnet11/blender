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

struct Main;
struct Object;
struct Collection;

/* Create new collection and assign it as a receiver collection for the light linking configuration
 * of the given object.
 *
 * The collection is created outside of the view layer collections.
 * If the object has already receiver collection set up it is unreferenced from the object.
 *
 * Returns the newly created collection. */
struct Collection *BKE_light_linking_receiver_collection_new(struct Main *bmain,
                                                             struct Object *object);

#ifdef __cplusplus
}
#endif
