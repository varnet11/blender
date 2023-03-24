/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2023 Blender Foundation */

#pragma once

#include <cstdint>

#include "MEM_guardedalloc.h"
#include "RNA_blender_cpp.h"

CCL_NAMESPACE_BEGIN

uint64_t light_linking_get_set_membership(const BL::Object &object);
uint light_linking_get_receiver_set(const BL::Object &object);

CCL_NAMESPACE_END
