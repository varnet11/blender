/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2023 Blender Foundation. */

/** \file
 * \ingroup depsgraph
 */

#pragma once

#include <cstdint>

#include "BLI_span.hh"

struct Object;
struct Depsgraph;

namespace blender::deg::light_linking {

/* Get emitters which has light linking configured to affect the given receiver.
 * The receiver can either be evaluated or original. The function returns original objects. */
Span<const Object *> get_original_emitters_for_receiver(const ::Depsgraph *depsgraph,
                                                        const Object *receiver);

/* Get pre-calculated emission mask of the given emitter.
 *
 * The mask is a bit-field with a single bit set, which corresponds to an identifier of the
 * emitter within the evaluation context.
 *
 * If the given emitter does not have light linking configured (the receiver collection is nullptr)
 * the function returns 0.
 *
 * The emitter can either be evaluated or original. The function returns mask which is only valid
 * within the given dependency graph. */
uint64_t get_emitter_mask(const ::Depsgraph *depsgraph, const Object *emitter);

}  // namespace blender::deg::light_linking
