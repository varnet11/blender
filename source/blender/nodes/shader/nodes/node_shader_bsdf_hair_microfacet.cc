/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2018 Blender Foundation. All rights reserved. */

#include "node_shader_util.hh"

#include "UI_interface.h"
#include "UI_resources.h"

namespace blender::nodes::node_shader_bsdf_hair_microfacet_cc {

/* Color, melanin and absorption coefficient default to approximately same brownish hair. */
static void node_declare(NodeDeclarationBuilder &b)
{
  b.add_input<decl::Color>(N_("Color"))
      .default_value({0.017513f, 0.005763f, 0.002059f, 1.0f})
      .description("The RGB color of the strand. Only used in Direct Coloring");
  b.add_input<decl::Float>(N_("Melanin"))
      .default_value(0.8f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR)
      .description("Hair pigment. Specify its absolute quantity between 0 and 1");
  b.add_input<decl::Float>(N_("Melanin Redness"))
      .default_value(1.0f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR)
      .description(
          "Fraction of pheomelanin in melanin, gives yellowish to reddish color, as opposed to "
          "the brownish to black color of eumelanin");
  b.add_input<decl::Color>(N_("Tint"))
      .default_value({1.0f, 1.0f, 1.0f, 1.0f})
      .description("Additional color used for dyeing the hair.");
  b.add_input<decl::Vector>(N_("Absorption Coefficient"))
      .default_value({0.245531f, 0.52f, 1.365f})
      .min(0.0f)
      .max(1000.0f);
  b.add_input<decl::Float>(N_("Aspect Ratio"))
      .default_value(0.85f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR)
      .description(
          "For elliptical hair cross-section, the aspect ratio is the ratio of the minor axis to "
          "the major axis. Recommended values are 0.8~1 for Asian hair, 0.65~0.9 for Caucasian "
          "hair, 0.5~0.65 for African hair");
  b.add_input<decl::Float>(N_("Roughness"))
      .default_value(0.3f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR)
      .description("Microfacet roughness");
  b.add_input<decl::Float>(N_("Reflection"))
      .default_value(1.0f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR)
      .description(
          "The first light bounce off the hair surface. The color of this component is always "
          "white");
  b.add_input<decl::Float>(N_("Transmission"))
      .default_value(1.0f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR)
      .description(
          "The component that is transmitted through the hair. Picks up the color of the pigment "
          "inside the hair");
  b.add_input<decl::Float>(N_("Secondary Reflection"))
      .default_value(1.0f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR)
      .description(
          "The component that is transmitted into the hair, reflected off the backside of the "
          "hair and then transmitted out of the hair, oriented approximately around the incoming "
          "direction. Picks up the color of the pigment inside the hair");
  b.add_input<decl::Float>(N_("IOR")).default_value(1.55f).min(0.0f).max(1000.0f).description(
      "Index of refraction determines how much the ray is bent. At 1.0 rays pass straight through "
      "like in a transparent material; higher values cause larger deflection in angle. Default "
      "value is 1.55 (the IOR of keratin)");
  b.add_input<decl::Float>(N_("Offset"))
      .default_value(2.0f * ((float)M_PI) / 180.0f)
      .min(-M_PI_2)
      .max(M_PI_2)
      .subtype(PROP_ANGLE)
      .description(
          "The tilt angle of the cuticle scales (the outermost part of the hair). They are always "
          "tilted towards the hair root. The value is usually between 2 and 4 for human hair");
  b.add_input<decl::Float>(N_("Blur"))
      .default_value(1.0f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR);
  b.add_input<decl::Float>(N_("Random Color"))
      .default_value(0.0f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR)
      .description("Vary the melanin concentration for each strand");
  b.add_input<decl::Float>(N_("Random Roughness"))
      .default_value(0.0f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR)
      .description("Vary roughness values for each strand");
  b.add_input<decl::Float>(N_("Random")).hide_value();
  b.add_input<decl::Float>(N_("Weight")).unavailable();
  b.add_output<decl::Shader>(N_("BSDF"));
}

static void node_shader_buts_microfacet_hair(uiLayout *layout, bContext * /*C*/, PointerRNA *ptr)
{
  uiItemR(layout, ptr, "parametrization", UI_ITEM_R_SPLIT_EMPTY_NAME, "", ICON_NONE);
  uiItemR(layout, ptr, "cross_section_type", UI_ITEM_R_SPLIT_EMPTY_NAME, "", ICON_NONE);
  uiItemR(layout, ptr, "distribution_type", UI_ITEM_R_SPLIT_EMPTY_NAME, "", ICON_NONE);
}

/* Initialize the custom Parametrization property to Color. */
static void node_shader_init_hair_microfacet(bNodeTree * /*ntree*/, bNode *node)
{
  NodeShaderHairMicrofacet *data = MEM_cnew<NodeShaderHairMicrofacet>(__func__);

  data->parametrization = SHD_MICROFACET_HAIR_REFLECTANCE;
  data->cross_section = SHD_MICROFACET_HAIR_CIRCULAR;
  data->distribution = SHD_MICROFACET_HAIR_GGX;

  node->storage = data;
}

/* Triggers (in)visibility of some sockets when changing Parametrization. */
static void node_shader_update_hair_microfacet(bNodeTree *ntree, bNode *node)
{
  NodeShaderHairMicrofacet *data = static_cast<NodeShaderHairMicrofacet *>(node->storage);

  int parametrization = data->parametrization;
  bool elliptical = (data->cross_section == SHD_MICROFACET_HAIR_ELLIPTIC);

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
    else if (STREQ(sock->name, "Aspect Ratio")) {
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
  node_type_storage(
      &ntype, "NodeShaderHairMicrofacet", node_free_standard_storage, node_copy_standard_storage);

  nodeRegisterType(&ntype);
}
