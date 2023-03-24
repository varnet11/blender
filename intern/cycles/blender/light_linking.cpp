/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2023 Blender Foundation */

#include "blender/light_linking.h"

#include "DNA_object_types.h"

CCL_NAMESPACE_BEGIN

static const ::Object *get_blender_object(const BL::Object &object)
{
  return reinterpret_cast<::Object *>(object.ptr.data);
}

uint64_t light_linking_get_set_membership(const BL::Object &object)
{
  const ::Object *blender_object = get_blender_object(object);
  return blender_object->light_linking.runtime.set_membership;
}

uint light_linking_get_receiver_set(const BL::Object &object)
{
  const ::Object *blender_object = get_blender_object(object);
  return blender_object->light_linking.runtime.receiver_set;
}

CCL_NAMESPACE_END
