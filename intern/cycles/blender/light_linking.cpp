/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2023 Blender Foundation */

#include "blender/light_linking.h"

#include "DNA_object_types.h"

CCL_NAMESPACE_BEGIN

static const ::Object *get_blender_object(const BL::Object &object)
{
  return reinterpret_cast<::Object *>(object.ptr.data);
}

uint64_t BlenderLightLink::get_light_set_membership(const BL::Object &object)
{
  const ::Object *blender_object = get_blender_object(object);
  return blender_object->light_linking.runtime.light_set_membership;
}

uint BlenderLightLink::get_receiver_light_set(const BL::Object &object)
{
  const ::Object *blender_object = get_blender_object(object);
  return blender_object->light_linking.runtime.receiver_light_set;
}

uint64_t BlenderLightLink::get_shadow_set_membership(const BL::Object &object)
{
  const ::Object *blender_object = get_blender_object(object);
  return blender_object->light_linking.runtime.shadow_set_membership;
}

uint BlenderLightLink::get_blocker_shadow_set(const BL::Object &object)
{
  const ::Object *blender_object = get_blender_object(object);
  return blender_object->light_linking.runtime.blocker_shadow_set;
}

CCL_NAMESPACE_END
