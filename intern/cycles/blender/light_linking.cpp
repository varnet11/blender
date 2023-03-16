/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2023 Blender Foundation */

#include "blender/light_linking.h"

#include "DNA_object_types.h"

CCL_NAMESPACE_BEGIN

static const ::Object *get_blender_object(const BL::Object &object)
{
  return reinterpret_cast<::Object *>(object.ptr.data);
}

uint64_t light_linking_get_emitter_mask(const BL::Object &object)
{
  const ::Object *blender_object = get_blender_object(object);
  return blender_object->light_linking.runtime.emitter_mask;
}

uint64_t light_linking_get_receiver_mask(const BL::Object &object)
{
  const ::Object *blender_object = get_blender_object(object);
  return blender_object->light_linking.runtime.receiver_mask;
}

CCL_NAMESPACE_END
