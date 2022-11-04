/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. */

/** \file
 * \ingroup draw
 */

#pragma once

#include "BKE_gpencil.h"
#include "BKE_image.h"
#include "DRW_gpu_wrapper.hh"
#include "DRW_render.h"

#include "draw_manager.hh"
#include "draw_pass.hh"

namespace blender::gpencil {

using namespace draw;

class LayerModule {
 private:
  /** Contains all Objects in the scene. Indexed by gpObject.layer_offset + layer_id. */
  StorageVectorBuffer<gpLayer> layers_buf_ = "gp_layers_buf";

 public:
  void begin_sync()
  {
    layers_buf_.clear();
  }

  void sync(const Object *object, const bGPDlayer *gpl, bool &do_layer_blending)
  {
    UNUSED_VARS(object, gpl);
    /* TODO(fclem): All of this is placeholder. */
    gpLayer layer;
    layer.vertex_color_opacity = 0.0f;
    layer.opacity = 1.0f;
    layer.thickness_offset = 0.0f;
    layer.tint = float4(1.0f, 1.0f, 1.0f, 0.0f);

    layers_buf_.append(layer);

    do_layer_blending = false;
  }

  void end_sync()
  {
    layers_buf_.push_update();
  }

  void bind_resources(PassMain::Sub &sub)
  {
    sub.bind_ssbo(GPENCIL_LAYER_SLOT, &layers_buf_);
  }

  uint object_offset_get() const
  {
    return layers_buf_.size();
  }
};

}  // namespace blender::gpencil
