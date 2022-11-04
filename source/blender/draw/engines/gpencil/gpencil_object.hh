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

#include "gpencil_layer.hh"
#include "gpencil_material.hh"
#include "gpencil_shader.hh"

namespace blender::gpencil {

using namespace draw;

class ObjectModule {
 private:
  /* TODO(fclem): This is a workaround to the current GPencil data structure.
   * Rendering is much easier if we have frame containing layers. */
  struct LayerData {
    /* Layer / Frame tuple representing a layer inside a frame. */
    const bGPDlayer *gpl;
    const bGPDframe *gpf;
  };
  struct FrameData {
    /* First layer frame to access frame members. */
    const bGPDframe *gpf;
    Vector<LayerData, 8> layers;
  };

  LayerModule &layers_;
  MaterialModule &materials_;
  ShaderModule &shaders_;

  /** Contains all Objects in the scene. Indexed by drw_ResourceID. */
  StorageArrayBuffer<gpObject> objects_buf_ = "gp_objects_buf";

  /** Contains all gpencil objects in the scene as well as their effect sub-passes. */
  PassSortable main_ps_ = {"gp_main_ps"};

  /** Contains all composited GPencil layers from one object if is uses VFX. */
  TextureFromPool object_color_tx_ = {"gp_color_object_tx"};
  TextureFromPool object_reveal_tx_ = {"gp_reveal_object_tx"};
  Framebuffer object_fb_ = {"gp_object_fb"};
  bool is_object_fb_needed_ = false;

  /** Contains all strokes from one layer if is uses blending. (also used as target for VFX) */
  TextureFromPool layer_color_tx_ = {"gp_color_layer_tx"};
  TextureFromPool layer_reveal_tx_ = {"gp_reveal_layer_tx"};
  Framebuffer layer_fb_ = {"gp_layer_fb"};
  bool is_layer_fb_needed_ = false;

  bool use_onion_ = true;
  bool use_stroke_fill_ = true;
  bool use_vfx_ = true;
  bool is_render_ = true;
  /** Forward vector used to sort gpencil objects. */
  float3 camera_forward_;
  /** Scene current frame. */
  float current_frame_ = 0;

  /** \note Needs not to be temporary variable since it is dereferenced later. */
  std::array<float4, 2> clear_colors_ = {float4(0.0f, 0.0f, 0.0f, 0.0f),
                                         float4(1.0f, 1.0f, 1.0f, 1.0f)};

 public:
  ObjectModule(LayerModule &layers, MaterialModule &materials, ShaderModule &shaders)
      : layers_(layers), materials_(materials), shaders_(shaders){};

  void init(const View3D *v3d, const Scene *scene)
  {
    const bool is_viewport = (v3d != nullptr);

    if (is_viewport) {
      /* TODO(fclem): Avoid access to global DRW. */
      const struct bContext *evil_C = DRW_context_state_get()->evil_C;
      const bool playing = (evil_C != nullptr) ?
                               ED_screen_animation_playing(CTX_wm_manager(evil_C)) != nullptr :
                               false;
      const bool hide_overlay = ((v3d->flag2 & V3D_HIDE_OVERLAYS) != 0);
      const bool show_onion = ((v3d->gp_flag & V3D_GP_SHOW_ONION_SKIN) != 0);
      use_onion_ = show_onion && !hide_overlay && !playing;
      use_stroke_fill_ = GPENCIL_SIMPLIFY_FILL(scene, playing);
      use_vfx_ = GPENCIL_SIMPLIFY_FX(scene, playing);
      is_render_ = false;
    }
    else {
      use_stroke_fill_ = GPENCIL_SIMPLIFY_FILL(scene, false);
      use_vfx_ = GPENCIL_SIMPLIFY_FX(scene, false);
    }
  }

  void begin_sync(Depsgraph *depsgraph, const View &main_view)
  {
    camera_forward_ = float3(main_view.viewinv()[2]);
    current_frame_ = DEG_get_ctime(depsgraph);

    is_object_fb_needed_ = false;
    is_layer_fb_needed_ = false;

    /* TODO(fclem): Shrink buffer. */
    // objects_buf_.shrink();
  }

  void sync_gpencil(Manager &manager,
                    ObjectRef &object_ref,
                    Framebuffer &main_fb,
                    PassSortable &main_ps)
  {
    Object *object = object_ref.object;
    bGPdata *gpd = static_cast<bGPdata *>(object->data);
    ListBaseWrapper<const bGPDlayer> layers(&gpd->layers);

    if (BLI_listbase_is_empty(&gpd->layers)) {
      return;
    }

    const bool is_stroke_order_3d = (gpd->draw_mode == GP_DRAWMODE_3D);
    bool do_material_holdout = false;
    bool do_layer_blending = false;
    bool object_has_vfx = false;  // TODO vfx.object_has_vfx(gpd);

    uint material_offset = materials_.object_offset_get();
    for (auto i : IndexRange(BKE_object_material_count_eval(object))) {
      materials_.sync(object, i, do_material_holdout);
    }

    uint layer_offset = layers_.object_offset_get();
    for (const bGPDlayer *layer : layers) {
      layers_.sync(object, layer, do_layer_blending);
    }

    /* Order rendering using camera Z distance. */
    float3 position = float3(object->object_to_world[3]);
    float camera_z = math::dot(position, camera_forward_);

    PassMain::Sub &object_subpass = main_ps.sub("GPObject", camera_z);
    object_subpass.framebuffer_set((object_has_vfx) ? &object_fb_ : &main_fb);
    object_subpass.clear_depth(is_stroke_order_3d ? 1.0f : 0.0f);
    if (object_has_vfx) {
      object_subpass.clear_multi(clear_colors_);
    }

    DRWState state = DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_BLEND_ALPHA_PREMUL;
    /* For 2D mode, we render all strokes with uniform depth (increasing with stroke id). */
    state |= (is_stroke_order_3d) ? DRW_STATE_DEPTH_LESS_EQUAL : DRW_STATE_DEPTH_GREATER;
    /* Always write stencil. Only used as optimization for blending. */
    state |= DRW_STATE_WRITE_STENCIL | DRW_STATE_STENCIL_ALWAYS;

    object_subpass.state_set(state);
    object_subpass.shader_set(shaders_.static_shader_get(GREASE_PENCIL));

    Vector<FrameData, 5> frames;
    displayed_frame_select(frames, layers);

    for (const FrameData &frame : frames) {
      /* TODO(fclem): Pass per frame object matrix here. */
      ResourceHandle handle = manager.resource_handle(object_ref);
      gpObject &ob = objects_buf_.get_or_resize(handle.resource_index());
      ob.is_shadeless = false;
      ob.stroke_order3d = false;
      ob.tint = frame_tint_get(gpd, frame.gpf, current_frame_);
      ob.layer_offset = layer_offset;
      ob.material_offset = material_offset;

      GPUVertBuf *position_tx = DRW_cache_gpencil_position_buffer_get(object, frame.gpf->framenum);
      GPUVertBuf *color_tx = DRW_cache_gpencil_color_buffer_get(object, frame.gpf->framenum);
      GPUBatch *geom = DRW_cache_gpencil_get(object, frame.gpf->framenum);

      if (do_layer_blending) {
        for (const LayerData &layer : frame.layers) {
          UNUSED_VARS(layer);
          // if (has_blending(layer)) {
          //   object_subpass.framebuffer_set(*vfx_fb.current());
          // }

          /* TODO(fclem): Only draw subrange of geometry for this layer. */
          object_subpass.draw(geom, handle);

          // if (has_blending(layer)) {
          //   layer_blend_sync(object_ref, object_subpass);
          // }
        }
      }
      else {
        /* Fast path. */
        object_subpass.bind_texture("gp_pos_tx", position_tx);
        object_subpass.bind_texture("gp_col_tx", color_tx);
        object_subpass.draw(geom, handle);
      }
    }

#if 0
    if (object_has_vfx) {
      VfxContext vfx_ctx(object_subpass,
                         layer_fb_,
                         object_fb_,
                         object_color_tx_,
                         layer_color_tx_,
                         object_reveal_tx_,
                         layer_reveal_tx_,
                         is_render_);

      /* \note Update this boolean as the actual number of vfx drawn might differ. */
      object_has_vfx = vfx.object_sync(main_fb_, object_ref, vfx_ctx, do_material_holdout);

      if (object_has_vfx || do_layer_blending) {
        is_layer_fb_needed_ = true;
      }
    }
#endif
  }

  void end_sync()
  {
    objects_buf_.push_update();
  }

  void bind_resources(PassMain::Sub &sub)
  {
    sub.bind_ssbo(GPENCIL_OBJECT_SLOT, &objects_buf_);
  }

  void acquire_temporary_buffers(int2 render_size, eGPUTextureFormat format)
  {
    object_color_tx_.acquire(render_size, format);
    object_reveal_tx_.acquire(render_size, format);
    object_fb_.ensure(GPU_ATTACHMENT_NONE,
                      GPU_ATTACHMENT_TEXTURE(object_color_tx_),
                      GPU_ATTACHMENT_TEXTURE(object_reveal_tx_));
    if (is_layer_fb_needed_) {
      layer_color_tx_.acquire(render_size, format);
      layer_reveal_tx_.acquire(render_size, format);
      layer_fb_.ensure(GPU_ATTACHMENT_NONE,
                       GPU_ATTACHMENT_TEXTURE(layer_color_tx_),
                       GPU_ATTACHMENT_TEXTURE(layer_reveal_tx_));
    }
  }

  void release_temporary_buffers()
  {
    object_color_tx_.release();
    object_reveal_tx_.release();

    layer_color_tx_.release();
    layer_reveal_tx_.release();
  }

  bool scene_has_visible_gpencil_object() const
  {
    return objects_buf_.size() > 0;
  }

 private:
  static float4 frame_tint_get(const bGPdata *gpd, const bGPDframe *gpf, int /* current_frame */)
  {
    /* TODO(fclem): Onion color should rely on time and or frame id and not on runtime.onion_id.
     * This should be evaluated at draw time as it is just a way of displaying the data. */
    const bool use_onion_custom_col = (gpd->onion_flag & GP_ONION_GHOST_PREVCOL) != 0;
    const bool use_onion_fade = (gpd->onion_flag & GP_ONION_FADE) != 0;
    const bool use_next_col = gpf->runtime.onion_id > 0.0f;

    const float *onion_col_custom = (use_onion_custom_col) ?
                                        (use_next_col ? gpd->gcolor_next : gpd->gcolor_prev) :
                                        U.gpencil_new_layer_col;

    float4 tint = {UNPACK3(onion_col_custom), 1.0f};

    tint[3] = use_onion_fade ? (1.0f / abs(gpf->runtime.onion_id)) : 0.5f;
    tint[3] *= gpd->onion_factor;
    tint[3] = (gpd->onion_factor > 0.0f) ? clamp_f(tint[3], 0.1f, 1.0f) :
                                           clamp_f(tint[3], 0.01f, 1.0f);
    return tint;
  }

  void displayed_frame_select(Vector<FrameData, 5> &frames,
                              ListBaseWrapper<const bGPDlayer> layers)
  {
    /* TODO(fclem): Select onion skin frames. */
    /** \note Change data layout to be Frame major instead of Layer major.
     * Hopefully the GPencil data layout will be closer to that in the future. */
    FrameData frame_data;
    frame_data.gpf = layers.get(0)->actframe;
    for (const bGPDlayer *layer : layers) {
      LayerData layer_data;
      layer_data.gpf = layer->actframe;
      layer_data.gpl = layer;
      frame_data.layers.append(layer_data);
    }
    frames.append(frame_data);
  }
};

}  // namespace blender::gpencil
