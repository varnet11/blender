/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2018 Blender Foundation. All rights reserved. */

#include "node_shader_util.hh"

#include "UI_interface.h"
#include "UI_resources.h"

namespace blender::nodes::node_shader_bsdf_hair_microfacet_cc {

/* Color, melanin and absorption coefficient default to approximately same brownish hair. */
static void node_declare(NodeDeclarationBuilder &b)
{
  b.add_input<decl::Color>(N_("Color")).default_value({0.017513f, 0.005763f, 0.002059f, 1.0f});
  b.add_input<decl::Float>(N_("Melanin"))
      .default_value(0.8f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR);
  b.add_input<decl::Float>(N_("Melanin Redness"))
      .default_value(1.0f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR);
  b.add_input<decl::Color>(N_("Tint")).default_value({1.0f, 1.0f, 1.0f, 1.0f});
  b.add_input<decl::Vector>(N_("Absorption Coefficient"))
      .default_value({0.245531f, 0.52f, 1.365f})
      .min(0.0f)
      .max(1000.0f);
  b.add_input<decl::Float>(N_("Eccentricity"))
      .default_value(0.85f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR);
  b.add_input<decl::Float>(N_("Random Axis"))
      .default_value(0.0f)
      .min(-1.0f)
      .max(1.0f);
  b.add_input<decl::Float>(N_("Twist Rate"))
      .default_value(0.0f)
      .min(0.0f)
      .max(2.0f)
      .subtype(PROP_FACTOR);
  b.add_input<decl::Float>(N_("Roughness"))
      .default_value(0.3f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR);
  b.add_input<decl::Float>(N_("R lobe"))
      .default_value(1.0f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR);
  b.add_input<decl::Float>(N_("TT lobe"))
      .default_value(1.0f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR);
  b.add_input<decl::Float>(N_("TRT lobe"))
      .default_value(1.0f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR);
  b.add_input<decl::Float>(N_("IOR")).default_value(1.55f).min(0.0f).max(1000.0f);
  b.add_input<decl::Float>(N_("Offset"))
      .default_value(2.0f * ((float)M_PI) / 180.0f)
      .min(-M_PI_2)
      .max(M_PI_2)
      .subtype(PROP_ANGLE);
  b.add_input<decl::Float>(N_("Blur"))
      .default_value(1.0f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR);
  b.add_input<decl::Float>(N_("Random Color"))
      .default_value(0.0f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR);
  b.add_input<decl::Float>(N_("Random Roughness"))
      .default_value(0.0f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR);
  b.add_input<decl::Float>(N_("Random")).hide_value();
  b.add_input<decl::Float>(N_("Weight")).unavailable();
  b.add_output<decl::Shader>(N_("BSDF"));
}

static void node_shader_buts_microfacet_hair(uiLayout *layout, bContext * /*C*/, PointerRNA *ptr)
{
  uiItemR(layout, ptr, "parametrization", UI_ITEM_R_SPLIT_EMPTY_NAME, "", ICON_NONE);
  uiItemR(layout, ptr, "model_type", UI_ITEM_R_SPLIT_EMPTY_NAME, "", ICON_NONE);
}

/* Initialize the custom Parametrization property to Color. */
static void node_shader_init_hair_microfacet(bNodeTree * /*ntree*/, bNode *node)
{
  node->custom1 = SHD_MICROFACET_HAIR_REFLECTANCE;
  node->custom2 = SHD_MICROFACET_HAIR_CIRCULAR_GGX;
}

/* Triggers (in)visibility of some sockets when changing Parametrization. */
static void node_shader_update_hair_microfacet(bNodeTree *ntree, bNode *node)
{
  int parametrization = node->custom1;
  int model_type = node->custom2;

  bool circular = (model_type == SHD_MICROFACET_HAIR_CIRCULAR_GGX ||
                   model_type == SHD_MICROFACET_HAIR_CIRCULAR_GGX_ANALYTIC ||
                   model_type == SHD_MICROFACET_HAIR_CIRCULAR_BECKMANN);
  bool elliptical = (model_type == SHD_MICROFACET_HAIR_ELLIPTIC_GGX ||
                     model_type == SHD_MICROFACET_HAIR_ELLIPTIC_BECKMANN);

  LISTBASE_FOREACH (bNodeSocket *, sock, &node->inputs) {
    if (STREQ(sock->name, "Color")) {
      nodeSetSocketAvailability(ntree, sock, parametrization == SHD_MICROFACET_HAIR_REFLECTANCE);
    }
    else if (STREQ(sock->name, "Melanin")) {
      nodeSetSocketAvailability(
          ntree, sock, parametrization == SHD_MICROFACET_HAIR_PIGMENT_CONCENTRATION);
    }
    else if (STREQ(sock->name, "Melanin Redness")) {
      nodeSetSocketAvailability(
          ntree, sock, parametrization == SHD_MICROFACET_HAIR_PIGMENT_CONCENTRATION);
    }
    else if (STREQ(sock->name, "Tint")) {
      nodeSetSocketAvailability(
          ntree, sock, parametrization == SHD_MICROFACET_HAIR_PIGMENT_CONCENTRATION);
    }
    else if (STREQ(sock->name, "Absorption Coefficient")) {
      nodeSetSocketAvailability(
          ntree, sock, parametrization == SHD_MICROFACET_HAIR_DIRECT_ABSORPTION);
    }
    else if (STREQ(sock->name, "Random Color")) {
      nodeSetSocketAvailability(
          ntree, sock, parametrization == SHD_MICROFACET_HAIR_PIGMENT_CONCENTRATION);
    }
    else if (STREQ(sock->name, "Eccentricity")) {
      nodeSetSocketAvailability(ntree, sock, elliptical);
    }
    else if (STREQ(sock->name, "Random Axis")) {
      nodeSetSocketAvailability(ntree, sock, elliptical);
    }
    else if (STREQ(sock->name, "Twist Rate")) {
      nodeSetSocketAvailability(ntree, sock, elliptical);
    }
  }
}

static int node_shader_gpu_hair_microfacet(GPUMaterial *mat,
                                           bNode *node,
                                           bNodeExecData * /*execdata*/,
                                           GPUNodeStack *in,
                                           GPUNodeStack *out)
{
  return GPU_stack_link(mat, node, "node_bsdf_hair_microfacet", in, out);
}

}  // namespace blender::nodes::node_shader_bsdf_hair_microfacet_cc

/* node type definition */
void register_node_type_sh_bsdf_hair_microfacet()
{
  namespace file_ns = blender::nodes::node_shader_bsdf_hair_microfacet_cc;

  static bNodeType ntype;

  sh_node_type_base(
      &ntype, SH_NODE_BSDF_HAIR_MICROFACET, "Microfacet Hair BSDF", NODE_CLASS_SHADER);
  ntype.declare = file_ns::node_declare;
  ntype.draw_buttons = file_ns::node_shader_buts_microfacet_hair;
  node_type_size_preset(&ntype, NODE_SIZE_LARGE);
  ntype.initfunc = file_ns::node_shader_init_hair_microfacet;
  ntype.updatefunc = file_ns::node_shader_update_hair_microfacet;
  ntype.gpu_fn = file_ns::node_shader_gpu_hair_microfacet;

  nodeRegisterType(&ntype);
}
