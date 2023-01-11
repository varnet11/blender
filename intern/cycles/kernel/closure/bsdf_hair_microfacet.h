/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2018-2022 Blender Foundation */

#pragma once

#ifndef __KERNEL_GPU__
#  include <fenv.h>
#endif

#include "kernel/util/color.h"

CCL_NAMESPACE_BEGIN

typedef struct MicrofacetHairExtra {
  float R;
  float TT;
  float TRT;

  float aspect_ratio;

  /* Geometry data. */
  float4 geom;
} MicrofacetHairExtra;

typedef struct MicrofacetHairBSDF {
  SHADER_CLOSURE_BASE;

  /* Absorption coefficient. */
  Spectrum sigma;
  /* Microfacet distribution roughness. */
  float roughness;
  /* Cuticle tilt angle. */
  float alpha;
  /* IOR. */
  float eta;

  /* Blur. */
  float blur;

  /* GGX/Beckmann. */
  int distribution_type;

  /* Circular/Elliptical */
  int cross_section_type;

  /* Extra closure. */
  ccl_private MicrofacetHairExtra *extra;
} MicrofacetHairBSDF;

static_assert(sizeof(ShaderClosure) >= sizeof(MicrofacetHairBSDF),
              "MicrofacetHairBSDF is too large!");
static_assert(sizeof(ShaderClosure) >= sizeof(MicrofacetHairExtra),
              "MicrofacetHairExtra is too large!");

#ifdef __HAIR__
/* Set up the hair closure. */
ccl_device int bsdf_microfacet_hair_setup(ccl_private ShaderData *sd,
                                          ccl_private MicrofacetHairBSDF *bsdf)
{
  bsdf->type = CLOSURE_BSDF_HAIR_MICROFACET_ID;

  bsdf->roughness = clamp(bsdf->roughness, 0.001f, 1.0f);

  /* Compute local frame, aligned to curve tangent and ray direction. */
  float3 Y = safe_normalize(sd->dPdu);
  float3 X = safe_normalize(cross(Y, sd->I));

  /* h -1..0..1 means the rays goes from grazing the hair, to hitting it at
   * the center, to grazing the other edge. This is the sine of the angle
   * between sd->Ng and Z, as seen from the tangent X. */

  float h = (sd->type & PRIMITIVE_CURVE_RIBBON) ? -sd->v : -dot(X, sd->N);

  kernel_assert(fabsf(h) < 1.0f + 1e-4f);
  kernel_assert(isfinite_safe(X));
  kernel_assert(isfinite_safe(h));

  if (bsdf->cross_section_type == NODE_MICROFACET_HAIR_ELLIPTIC) {
    /* Local frame is independent of the ray direction for elliptical hairs. */
    bsdf->extra->geom.w = h;
  }
  else {
    bsdf->extra->geom = make_float4(X.x, X.y, X.z, h);
  }

  return SD_BSDF | SD_BSDF_HAS_EVAL | SD_BSDF_NEEDS_LCG | SD_BSDF_HAS_TRANSMISSION;
}

#endif /* __HAIR__ */

/* -------------------------------------------------------------------- */
/** \name Hair coordinate system utils.
 * \{ */

/* Returns sin(theta) of the given direction. */
ccl_device_inline float sin_theta(const float3 w)
{
  return w.y;
}

/* Returns cos(theta) of the given direction. */
ccl_device_inline float cos_theta(const float3 w)
{
  return safe_sqrtf(sqr(w.x) + sqr(w.z));
}

/* Returns tan(theta) of the given direction. */
ccl_device_inline float tan_theta(const float3 w)
{
  return sin_theta(w) / cos_theta(w);
}

/* Returns sin(phi) and cos(phi) of the given direction. */
ccl_device float sin_phi(const float3 w)
{
  return w.x / cos_theta(w);
}

ccl_device float2 sincos_phi(const float3 w)
{
  float c = cos_theta(w);
  return make_float2(w.x / c, w.z / c);
}

/* Extract the theta coordinate from the given direction.
 * -pi < theta < pi */
ccl_device_inline float dir_theta(const float3 w)
{
  return atan2f(sin_theta(w), cos_theta(w));
}

/* Extract the phi coordinate from the given direction, assuming phi(wi) = 0.
 * -pi < phi < pi */
ccl_device_inline float dir_phi(const float3 w)
{
  return atan2f(w.x, w.z);
}

/* Extract theta and phi coordinates from the given direction, assuming phi(wi) = 0.
 * -pi/2 < theta < pi/2, -pi < phi < pi */
ccl_device_inline float2 dir_sph(const float3 w)
{
  return make_float2(dir_theta(w), dir_phi(w));
}

/* Compute the vector direction given spherical coordinates */
ccl_device_inline float3 sph_dir(float theta, float gamma)
{
  float sin_theta, cos_theta, sin_gamma, cos_gamma;
  fast_sincosf(theta, &sin_theta, &cos_theta);
  fast_sincosf(gamma, &sin_gamma, &cos_gamma);
  return make_float3(sin_gamma * cos_theta, sin_theta, cos_gamma * cos_theta);
}

/* Utility functions for elliptical cross-sections. */

/* Conversion between gamma and phi. Notations see Figure 5 in the paper. */
ccl_device float to_phi(float gamma, float a, float b)
{
  float sin_gamma, cos_gamma;
  fast_sincosf(gamma, &sin_gamma, &cos_gamma);
  return atan2f(b * sin_gamma, a * cos_gamma);
}

ccl_device float to_gamma(float phi, float a, float b)
{
  float sin_phi, cos_phi;
  fast_sincosf(phi, &sin_phi, &cos_phi);
  return atan2f(a * sin_phi, b * cos_phi);
}

/* Compute the coordinate on the ellipse, given gamma, the semi-major and semi-minor axes. */
ccl_device float2 to_point(float gamma, float a, float b)
{
  float sin_gamma, cos_gamma;
  fast_sincosf(gamma, &sin_gamma, &cos_gamma);
  return make_float2(a * sin_gamma, b * cos_gamma);
}

/* Compute the vector direction given by theta and gamma. */
ccl_device float3 sphg_dir(float theta, float gamma, float a, float b)
{
  float sin_theta, cos_theta, sin_gamma, cos_gamma;
  fast_sincosf(theta, &sin_theta, &cos_theta);
  fast_sincosf(gamma, &sin_gamma, &cos_gamma);
  float tan_gamma = sin_gamma / cos_gamma;
  float tan_phi = b / a * tan_gamma;
  float cos_phi = signf(cos_gamma) / sqrtf(sqr(tan_phi) + 1.0f);
  float sin_phi = cos_phi * tan_phi;
  return make_float3(sin_phi * cos_theta, sin_theta, cos_phi * cos_theta);
}

/** \} */

/* sample microfacets from a tilted mesonormal */
ccl_device_inline float3 sample_wh(KernelGlobals kg,
                                   const bool beckmann,
                                   const float roughness,
                                   const float3 wi,
                                   const float3 wm,
                                   const float2 rand)
{
  /* Coordinate transformation for microfacet sampling */
  float3 s, t;
  const float3 n = wm;
  make_orthonormals(n, &s, &t);

  const float3 wi_wm = make_float3(dot(wi, s), dot(wi, t), dot(wi, n));

  float G1o;
  const float3 wh_wm = microfacet_sample_stretched(
      kg, wi_wm, roughness, roughness, rand.x, rand.y, beckmann, &G1o);

  const float3 wh = wh_wm.x * s + wh_wm.y * t + wh_wm.z * n;
  return wh;
}

/* Check micronormal/mesonormal direct visiblity from v */
ccl_device_inline bool microfacet_visible(const float3 v, const float3 m, const float3 h)
{
  return (dot(v, h) > 0.0f && dot(v, m) > 0.0f);
}

/* Check micronormal/mesonormal direct visiblity from wi and wo */
ccl_device_inline bool microfacet_visible(const float3 wi,
                                          const float3 wo,
                                          const float3 m,
                                          const float3 h)
{
  return microfacet_visible(wi, m, h) && microfacet_visible(wo, m, h);
}

/* Check micronormal/mesonormal statistical visiblity from v: Smith's separable shadowing/masking
 * term */
ccl_device_inline float smith_g1(
    const bool beckmann, const float roughness, const float3 v, const float3 m, const float3 h)
{
  /* Assume consistent orientation (can't see the back of the microfacet from the front and vice
   * versa) */
  float cos_vm = dot(v, m);
  if (dot(v, h) <= 0.0f || cos_vm <= 0.0f)
    return 0.0f;

  return beckmann ? bsdf_beckmann_G1(roughness, cos_vm) :
                    2.0f / (1.0f + sqrtf(1.0f + sqr(roughness) * (1.0f / sqr(cos_vm) - 1.0f)));
}

/* Smith's separable shadowing-masking approximation */
ccl_device_inline float G(const bool beckmann,
                          const float roughness,
                          const float3 wi,
                          const float3 wo,
                          const float3 m,
                          const float3 h)
{
  return smith_g1(beckmann, roughness, wi, m, h) * smith_g1(beckmann, roughness, wo, m, h);
}

/* Normal Distribution Function */
ccl_device float D(const bool beckmann, const float roughness, const float3 m, const float3 h)
{
  float cos_theta = dot(h, m);

  float result;
  const float roughness2 = sqr(roughness);
  const float cos_theta2 = sqr(cos_theta);

  if (beckmann) {
    result = expf((1.0f - 1.0f / cos_theta2) / roughness2) /
             (M_PI_F * roughness2 * sqr(cos_theta2));
  }
  else { /* GGX */
    result = roughness2 / (M_PI_F * sqr(1.0f + (roughness2 - 1.0f) * cos_theta2));
  }

  /* Prevent potential numerical issues in other stages of the model */
  return (result * cos_theta > 1e-20f) ? result : 0.0f;
}

/* Compute fresnel reflection. Also return the dot product of the refracted ray and the normal as
 * `cos_theta_t`, as it is used when computing the direction of the refracted ray. */
ccl_device float fresnel(float cos_theta_i, float eta, ccl_private float *cos_theta_t)
{
  kernel_assert(!isnan_safe(cos_theta_i));

  /* Special cases. */
  if (eta == 1.0f) {
    return 0.0f;
  }
  if (cos_theta_i == 0.0f) {
    return 1.0f;
  }

  cos_theta_i = fabsf(cos_theta_i);

  /* Using Snell's law, calculate the squared cosine of the angle between the surface normal and
   * the transmitted ray. */
  float cos_theta_t_sqr = 1.0f - (1.0f - cos_theta_i * cos_theta_i) / (eta * eta);
  *cos_theta_t = safe_sqrtf(cos_theta_t_sqr);

  if (cos_theta_t_sqr <= 0) {
    /* Total internal reflection. */
    return 1.0f;
  }

  /* Amplitudes of reflected waves. */
  float a_s = (cos_theta_i - eta * (*cos_theta_t)) / (cos_theta_i + eta * (*cos_theta_t));
  float a_p = (*cos_theta_t - eta * cos_theta_i) / (*cos_theta_t + eta * cos_theta_i);

  float r = 0.5f * (sqr(a_s) + sqr(a_p));

  /* Adjust the sign of the transmitted direction to be relative to the surface normal. */
  *cos_theta_t = -(*cos_theta_t);

  return r;
}

ccl_device_inline float3 refract_angle(const float3 incident,
                                       const float3 normal,
                                       const float cos_theta_t,
                                       const float inv_eta)
{
  return inv_eta * incident - (inv_eta * dot(normal, incident) + cos_theta_t) * normal;
}

ccl_device float3 bsdf_microfacet_hair_eval_r_circular(ccl_private const ShaderClosure *sc,
                                                       const float3 wi,
                                                       const float3 wo)
{
  ccl_private MicrofacetHairBSDF *bsdf = (ccl_private MicrofacetHairBSDF *)sc;
  const float tilt = -bsdf->alpha;
  const float roughness = bsdf->roughness;
  const float eta = bsdf->eta;
  const bool beckmann = (bsdf->distribution_type == NODE_MICROFACET_HAIR_BECKMANN);

  float3 R = zero_float3();
  if (bsdf->extra->R <= 0.0f)
    return R;

  const float3 wh = normalize(wi + wo);
  const float phi_o = dir_phi(wo);

  /* dot(wi, wmi) > 0 */
  const float tan_tilt = tanf(tilt);
  float phi_m_max1 = acosf(fmaxf(-tan_tilt * tan_theta(wi), 0.0f));
  if (isnan_safe(phi_m_max1))
    return R;
  const float phi_m_min1 = -phi_m_max1;

  /* dot(wo, wmi) > 0 */
  const float phi_m_max2 = acosf(fmaxf(-tan_tilt * tan_theta(wo), 0.0f)) + phi_o;
  if (isnan_safe(phi_m_max2))
    return R;
  const float phi_m_min2 = -phi_m_max2 + 2.0f * phi_o;

  const float phi_m_min = fmaxf(phi_m_min1, phi_m_min2) + 1e-5f;
  const float phi_m_max = fminf(phi_m_max1, phi_m_max2) - 1e-5f;
  if (phi_m_min > phi_m_max)
    return R;

  float integral = 0.0f;

  /* Maximal sample resolution. */
  float res = roughness * 0.7f;
  /* Number of intervals should be even. */
  const size_t intervals = 2 * (size_t)ceilf((phi_m_max - phi_m_min) / res * 0.5f);

  /* Modified resolution based on numbers of intervals. */
  res = (phi_m_max - phi_m_min) / float(intervals);

  /* Integrate using Simpson's rule. */
  for (size_t i = 0; i <= intervals; i++) {

    const float phi_m = phi_m_min + i * res;
    const float3 wm = sph_dir(tilt, phi_m);

    if (microfacet_visible(wi, wo, make_float3(wm.x, 0.0f, wm.z), wh)) {
      const float weight = (i == 0 || i == intervals) ? 0.5f : (i % 2 + 1);
      integral += weight * D(beckmann, roughness, wm, wh) * G(beckmann, roughness, wi, wo, wm, wh);
    }
  }
  integral *= (2.0f / 3.0f * res);

  const float F = fresnel_dielectric_cos(dot(wi, wh), eta);

  R = make_float3(bsdf->extra->R * 0.125f * F * fmaxf(0.0f, integral));
  return R;
}

ccl_device float3 bsdf_microfacet_hair_eval_tt_trt_circular(KernelGlobals kg,
                                                            ccl_private const ShaderClosure *sc,
                                                            const float3 wi,
                                                            const float3 wo,
                                                            uint rng_quadrature)
{
  ccl_private MicrofacetHairBSDF *bsdf = (ccl_private MicrofacetHairBSDF *)sc;
  const float tilt = -bsdf->alpha;
  const float roughness = bsdf->roughness;
  const float eta = bsdf->eta;
  const bool beckmann = (bsdf->distribution_type == NODE_MICROFACET_HAIR_BECKMANN);

  if (bsdf->extra->TT <= 0.0f && bsdf->extra->TRT <= 0.0f)
    return zero_float3();

  /* dot(wi, wmi) > 0 */
  const float tan_tilt = tanf(tilt);
  const float phi_m_max = acosf(fmaxf(-tan_tilt * tan_theta(wi), 0.0f)) - 1e-3f;
  if (isnan_safe(phi_m_max))
    return zero_float3();
  const float phi_m_min = -phi_m_max;

  /* dot(wo, wmo) < 0 */
  float tmp1 = acosf(fminf(tan_tilt * tan_theta(wo), 0.0f));
  if (isnan_safe(tmp1))
    return zero_float3();

  const float3 mu_a = bsdf->sigma;
  const float inv_eta = 1.0f / eta;

  float res = roughness * 0.8f;
  const size_t intervals = 2 * (size_t)ceilf((phi_m_max - phi_m_min) / res * 0.5f);
  res = (phi_m_max - phi_m_min) / intervals;

  float3 S_tt = zero_float3();
  float3 S_trt = zero_float3();
  for (size_t i = 0; i <= intervals; i++) {

    const float phi_mi = phi_m_min + i * res;
    const float3 wmi = sph_dir(tilt, phi_mi);

    /* sample wh1 */
    const float2 sample1 = make_float2(lcg_step_float(&rng_quadrature),
                                       lcg_step_float(&rng_quadrature));

    const float3 wh1 = sample_wh(kg, beckmann, roughness, wi, wmi, sample1);
    const float dot_wi_wh1 = dot(wi, wh1);
    if (!(dot_wi_wh1 > 0)) {
      continue;
    }

    float cos_theta_t1;
    const float T1 = 1.0f - fresnel(dot_wi_wh1, eta, &cos_theta_t1);

    /* refraction at the first interface */
    const float3 wt = -refract_angle(wi, wh1, cos_theta_t1, inv_eta);
    const float phi_t = dir_phi(wt);
    const float phi_mt = 2.0f * phi_t - phi_mi;
    const float3 wmt = sph_dir(-tilt, phi_mt);

    const float G1 = G(beckmann, roughness, wi, -wt, wmi, wh1);
    if (G1 == 0.0f || !microfacet_visible(wi, -wt, make_float3(wmi.x, 0.0f, wmi.z), wh1)) {
      continue;
    }

    /* Simpson's rule weight */
    float weight = (i == 0 || i == intervals) ? 0.5f : (i % 2 + 1);

    float3 A_t = exp(mu_a * 2.0f * cosf(phi_t - phi_mi) / cos_theta(wt));

    /* TT */
    if (bsdf->extra->TT > 0.0f) {

      /* total internal reflection */
      if (dot(wo, wt) >= inv_eta - 1e-5f) {

        /* microfacet visiblity from macronormal */
        float3 wh2 = -wt + inv_eta * wo;
        if (dot(wmt, wh2) >= 0.0f) {

          const float rcp_norm_wh2 = 1.0f / len(wh2);
          wh2 *= rcp_norm_wh2;

          float dot_wt_wh2 = dot(-wt, wh2);

          const float T2 = 1.0f - fresnel_dielectric_cos(dot_wt_wh2, inv_eta);
          float D2 = D(beckmann, roughness, wh2, wmt) * G(beckmann, roughness, -wt, -wo, wmt, wh2);

          const float3 result = T1 * T2 * D2 * A_t * dot_wt_wh2 * dot(wo, wh2) *
                                sqr(rcp_norm_wh2) / dot(wt, wmi) * weight *
                                smith_g1(beckmann, roughness, -wt, wmi, wh1) * dot(wi, wmi);

          if (isfinite_safe(result))
            S_tt += bsdf->extra->TT * result;
        }
      }
    }

    /* TRT */
    if (bsdf->extra->TRT > 0.0f) {

      /* sample wh2 */
      const float2 sample2 = make_float2(lcg_step_float(&rng_quadrature),
                                         lcg_step_float(&rng_quadrature));

      float3 wh2 = sample_wh(kg, beckmann, roughness, -wt, wmt, sample2);

      const float cos_th2 = dot(-wt, wh2);
      if (!(cos_th2 > 0)) {
        continue;
      }

      const float R2 = fresnel_dielectric_cos(cos_th2, inv_eta);
      float3 wtr = -reflect(wt, wh2);

      float G2 = G(beckmann, roughness, -wt, -wtr, wmt, wh2);
      if (G2 == 0.0f || !microfacet_visible(-wt, -wtr, make_float3(wmt.x, 0.0f, wmt.z), wh2)) {
        continue;
      }

      /* total internal reflection */
      if (dot(-wtr, wo) < inv_eta - 1e-5f) {
        continue;
      }

      float phi_tr = dir_phi(wtr);
      float phi_mtr = phi_mi - 2.0f * (phi_t - phi_tr) + M_PI_F;
      const float3 wmtr = sph_dir(-tilt, phi_mtr);

      float3 wh3 = wtr + inv_eta * wo;
      float G3 = G(beckmann, roughness, wtr, -wo, wmtr, wh3);
      if (dot(wmtr, wh3) < 0.0f || G3 == 0.0f ||
          !microfacet_visible(wtr, -wo, make_float3(wmtr.x, 0.0f, wmtr.z), wh3)) {
        continue;
      }

      float rcp_norm_wh3 = 1.0f / len(wh3);
      wh3 *= rcp_norm_wh3;
      const float cos_trh3 = dot(wh3, wtr);

      const float T3 = 1.0f - fresnel_dielectric_cos(cos_trh3, inv_eta);
      const float D3 = D(beckmann, roughness, wh3, wmtr) * G3;

      const float3 A_tr = exp(mu_a * 2.0f * cosf(phi_tr - phi_mt) / cos_theta(wtr));

      const float3 result = T1 * R2 * T3 * D3 * cos_trh3 * dot(wh3, wo) * sqr(rcp_norm_wh3) * A_t *
                            A_tr * weight / (dot(wt, wmi) * dot(wtr, wmt)) *
                            smith_g1(beckmann, roughness, -wt, wmi, wh1) *
                            smith_g1(beckmann, roughness, -wtr, wmt, wh2) * dot(wi, wmi) *
                            dot(wt, wmt);

      if (isfinite_safe(result)) {
        S_trt += bsdf->extra->TRT * result;
      }
    }
  }

  return (S_tt + S_trt) / 3.0f * res * sqr(inv_eta);
}

ccl_device Spectrum bsdf_microfacet_hair_eval_circular(KernelGlobals kg,
                                                       ccl_private const ShaderData *sd,
                                                       ccl_private const ShaderClosure *sc,
                                                       const float3 omega_in,
                                                       ccl_private float *pdf)
{
  /* Define local coordinate system:
   * . X along the hair fiber width & perpendicular to the incoming ray (sd->I).
   * . Y along the hair fiber tangent.
   * . Z = X ^ Y (facing the incoming ray (sd->I) as much as possible). */
  ccl_private MicrofacetHairBSDF *bsdf = (ccl_private MicrofacetHairBSDF *)sc;
  const float3 X = float4_to_float3(bsdf->extra->geom);
  const float3 Y = safe_normalize(sd->dPdu);
  const float3 Z = safe_normalize(cross(X, Y));

  /* Get local wi/wo (convention is reversed from other hair bcsdfs). */
  const float3 wi = make_float3(dot(sd->I, X), dot(sd->I, Y), dot(sd->I, Z));
  const float3 wo = make_float3(dot(omega_in, X), dot(omega_in, Y), dot(omega_in, Z));

  /* Evaluate R, TT, TRT terms. */
  const float3 R = bsdf_microfacet_hair_eval_r_circular(sc, wi, wo) +
                   bsdf_microfacet_hair_eval_tt_trt_circular(kg, sc, wi, wo, sd->lcg_state);

  /* TODO: better estimation of the pdf */
  *pdf = 1.0f;

  return rgb_to_spectrum(R / cos_theta(wi));
}

ccl_device int bsdf_microfacet_hair_sample_circular(const KernelGlobals kg,
                                                    ccl_private const ShaderClosure *sc,
                                                    ccl_private ShaderData *sd,
                                                    float randu,
                                                    float randv,
                                                    ccl_private Spectrum *eval,
                                                    ccl_private float3 *omega_in,
                                                    ccl_private float *pdf,
                                                    ccl_private float2 *sampled_roughness,
                                                    ccl_private float *eta)
{
  ccl_private MicrofacetHairBSDF *bsdf = (ccl_private MicrofacetHairBSDF *)sc;
  *sampled_roughness = make_float2(bsdf->roughness, bsdf->roughness);
  *eta = bsdf->eta;
  const float inv_eta = 1.0f / *eta;

  if (bsdf->extra->R <= 0.0f && bsdf->extra->TT <= 0.0f && bsdf->extra->TRT <= 0.0f) {
    /* early out for inactive lobe */
    *pdf = 0.0f;
    return LABEL_NONE;
  }

  /* Define local coordinate system:
   * . X along the hair fiber width & perpendicular to the incoming ray (sd->I).
   * . Y along the hair fiber tangent.
   * . Z = X ^ Y (facing the incoming ray (sd->I) as much as possible). */
  const float3 X = float4_to_float3(bsdf->extra->geom);
  const float3 Y = safe_normalize(sd->dPdu);
  const float3 Z = safe_normalize(cross(X, Y));

  /* Get local wi (convention is reversed from other hair bcsdfs). */
  const float3 wi = make_float3(dot(sd->I, X), dot(sd->I, Y), dot(sd->I, Z));

  const float tilt = -bsdf->alpha;
  const float roughness = bsdf->roughness;
  const bool beckmann = (bsdf->distribution_type == NODE_MICROFACET_HAIR_BECKMANN);

  /* generate sample */
  float sample_lobe = randu;
  const float sample_h = randv;
  const float2 sample_h1 = make_float2(lcg_step_float(&sd->lcg_state),
                                       lcg_step_float(&sd->lcg_state));
  const float2 sample_h2 = make_float2(lcg_step_float(&sd->lcg_state),
                                       lcg_step_float(&sd->lcg_state));
  const float2 sample_h3 = make_float2(lcg_step_float(&sd->lcg_state),
                                       lcg_step_float(&sd->lcg_state));

  /* sample offset h = -sin(phi_m) */
  const float sin_phi_mi = sample_h * 2.0f - 1.0f;
  const float cos_phi_mi = safe_sqrtf(1.0f - sqr(sin_phi_mi));

  float st, ct;
  fast_sincosf(tilt, &st, &ct);

  const float3 wmi = make_float3(sin_phi_mi * ct, st, cos_phi_mi * ct); /* mesonormal */
  const float3 wmi_ = make_float3(sin_phi_mi, 0.0f, cos_phi_mi);        /* macronormal */

  if (dot(wmi, wi) < 0.0f || dot(wmi_, wi) < 0.0f) {
    /* macro/mesonormal invisible */
    *pdf = 0.0f;
    return LABEL_NONE;
  }

  /* sample R lobe */
  const float3 wh1 = sample_wh(kg, beckmann, roughness, wi, wmi, sample_h1);
  const float3 wr = -reflect(wi, wh1);

  /* ensure that this is a valid sample */
  if (!(dot(wr, wh1) > 0.0f) || !(dot(wr, wmi) > 0.0f) || !microfacet_visible(wi, wr, wmi_, wh1)) {
    *pdf = 0.0f;
    return LABEL_NONE;
  }

  float3 TT = zero_float3();
  float3 TRT = zero_float3();

  float cos_theta_t1;
  float R1 = fresnel(dot(wi, wh1), *eta, &cos_theta_t1);
  float3 R = make_float3(bsdf->extra->R * R1);

  /* sample TT lobe */
  const float3 wt = -refract_angle(wi, wh1, cos_theta_t1, inv_eta);
  const float phi_t = dir_phi(wt);

  float phi_mi = atan2f(sin_phi_mi, cos_phi_mi);
  float phi_mt = 2.0f * phi_t - phi_mi;

  const float3 wmt = sph_dir(-tilt, phi_mt);
  const float3 wmt_ = sph_dir(0.0f, phi_mt);

  const float3 wh2 = sample_wh(kg, beckmann, roughness, -wt, wmt, sample_h2);
  const float3 wtr = -reflect(wt, wh2);

  float3 wh3;
  float3 wtt, wtrt;
  float3 wmtr, wmtr_;

  if (dot(wt, wh2) < 0.0f && dot(wmt, wt) < 0.0f &&
      microfacet_visible(-wt, -wtr, make_float3(wmt.x, 0.0f, wmt.z), wh2)) {

    const float3 mu_a = bsdf->sigma;
    const float cos_gamma_t = -cosf(phi_t - phi_mi);
    const float cos_theta_wt = sqrtf(1.0f - sqr(wt.y));
    const float3 A_t = exp(-mu_a * (2.0f * cos_gamma_t / cos_theta_wt));

    float cos_theta_t2;
    const float R2 = fresnel(dot(-wt, wh2), inv_eta, &cos_theta_t2);
    const float3 T1 = make_float3(1.0f - R1);
    const float3 T2 = make_float3(1.0f - R2);

    wtt = -refract_angle(-wt, wh2, cos_theta_t2, *eta);

    if (dot(wtt, wmt) < 0.0f && cos_theta_t2 != 0.0f) /* total internal reflection */
      TT = bsdf->extra->TT * T1 * A_t * T2;

    /* sample TRT lobe */
    const float phi_tr = dir_phi(wtr);
    const float phi_mtr = phi_mi - 2.0f * (phi_t - phi_tr) + M_PI_F;
    wmtr = sph_dir(-tilt, phi_mtr);
    wmtr_ = sph_dir(0.0f, phi_mtr);

    wh3 = sample_wh(kg, beckmann, roughness, wtr, wmtr, sample_h3);

    float cos_theta_t3;
    const float R3 = fresnel(dot(wtr, wh3), inv_eta, &cos_theta_t3);

    wtrt = -refract_angle(wtr, wh3, cos_theta_t3, *eta);

    if (cos_theta_t3 != 0.0f && dot(wtr, wh3) > 0.0f && dot(wmtr, wtr) > 0.0f &&
        dot(wtrt, wmtr) < 0.0f &&
        microfacet_visible(wtr, -wtrt, make_float3(wmtr.x, 0.0f, wmtr.z), wh3)) {

      const float3 T3 = make_float3(1.0f - R3);
      const float cos_gamma_t2 = -cos(phi_tr - phi_mt);
      const float cos_theta_wtr = sqrtf(1.0f - sqr(wtr.y));
      const float3 A_tr = exp(-mu_a * (2.0f * cos_gamma_t2 / cos_theta_wtr));

      TRT = bsdf->extra->TRT * T1 * R2 * T3 * A_t * A_tr;
    }
  }

  /* select lobe based on energy */
  const float r = average(R);
  const float tt = average(TT);
  const float trt = average(TRT);
  const float total_energy = r + tt + trt;

  if (total_energy <= 0.0f) {
    *pdf = 0.0f;
    return LABEL_NONE;
  }

  float3 wo;
  float visibility = 0.0f;
  int label = LABEL_GLOSSY;

  sample_lobe *= total_energy;
  if (sample_lobe < r) {
    wo = wr;
    *eval = rgb_to_spectrum(R / r * total_energy);

    if (microfacet_visible(wi, wr, wmi_, wh1)) {
      visibility = smith_g1(beckmann, roughness, wr, wmi, wh1);
    }

    label |= LABEL_REFLECT;
  }
  else if (sample_lobe < (r + tt)) {
    wo = wtt;
    *eval = rgb_to_spectrum(TT / tt * total_energy);

    if (microfacet_visible(wi, -wt, wmi_, wh1) && microfacet_visible(-wt, -wtt, wmt_, wh2)) {
      visibility = smith_g1(beckmann, roughness, -wt, wmi, wh1) *
                   smith_g1(beckmann, roughness, -wtt, wmt, wh2);
    }

    label |= LABEL_TRANSMIT;
  }
  else { /* if (sample_lobe >= (r + tt)) */
    wo = wtrt;
    *eval = rgb_to_spectrum(TRT / trt * total_energy);

    if (microfacet_visible(wi, -wt, wmi_, wh1) && microfacet_visible(-wt, -wtr, wmt_, wh2) &&
        microfacet_visible(wtr, -wtrt, wmtr_, wh3)) {
      visibility = smith_g1(beckmann, roughness, -wt, wmi, wh1) *
                   smith_g1(beckmann, roughness, -wtr, wmt, wh2) *
                   smith_g1(beckmann, roughness, -wtrt, wmtr, wh3);
    }

    label |= LABEL_TRANSMIT;
  }

  *eval *= visibility;

  *omega_in = wo.x * X + wo.y * Y + wo.z * Z;

  /* correction of the cosine foreshortening term
   *eval *= dot(wi, wmi) / dot(wi, wmi_); */

  /* ensure the same pdf is returned for BSDF and emitter sampling. The importance sampling pdf is
   * already factored in the value so this value is only used for mis */
  *pdf = 1.0f;

  return label;
}

/* Elliptic specific */

ccl_device float3 bsdf_microfacet_hair_eval_r_elliptic(ccl_private const ShaderClosure *sc,
                                                       const float3 wi_,
                                                       const float3 wo_)
{
  ccl_private MicrofacetHairBSDF *bsdf = (ccl_private MicrofacetHairBSDF *)sc;
  const float tilt = -bsdf->alpha;
  const float roughness = bsdf->roughness;
  const float eta = bsdf->eta;
  const bool beckmann = (bsdf->distribution_type == NODE_MICROFACET_HAIR_BECKMANN);

  float3 R = zero_float3();
  if (bsdf->extra->R <= 0.0f) {
    return R;
  }

  /* get elliptical cross section characteristic */
  const float a = 1.0f;
  /* TODO: rename this as aspect ratio? e is in [0, 0.85], b is in [0.52, 1]. */
  const float b = bsdf->extra->aspect_ratio;
  const float e2 = 1.0f - sqr(b / a);

  /* this follows blender's convention (unlike the circular case?) */
  const float3 wo = wi_;
  const float3 wi = wo_;

  const float phi_i = dir_phi(wi);
  const float phi_o = dir_phi(wo);
  const float3 wh = normalize(wi + wo);

  /* dot(wi, wmi) > 0 */
  const float tan_tilt = tanf(tilt);
  float phi_m_max1 = acosf(fmaxf(-tan_tilt * tan_theta(wi), 0.0f)) + phi_i;
  if (isnan_safe(phi_m_max1))
    return R;
  float phi_m_min1 = -phi_m_max1 + 2.0f * phi_i;

  /* dot(wo, wmi) > 0 */
  float phi_m_max2 = acosf(fmaxf(-tan_tilt * tan_theta(wo), 0.0f)) + phi_o;
  if (isnan_safe(phi_m_max2)) {
    return R;
  }
  float phi_m_min2 = -phi_m_max2 + 2.0f * phi_o;

  /* try to wrap range */
  if ((phi_m_max2 - phi_m_min1) > M_2PI_F) {
    phi_m_min2 -= M_2PI_F;
    phi_m_max2 -= M_2PI_F;
  }
  if ((phi_m_max1 - phi_m_min2) > M_2PI_F) {
    phi_m_min1 -= M_2PI_F;
    phi_m_max1 -= M_2PI_F;
  }

  const float phi_m_min = fmaxf(phi_m_min1, phi_m_min2) + 1e-3f;
  const float phi_m_max = fminf(phi_m_max1, phi_m_max2) - 1e-3f;
  if (phi_m_min > phi_m_max) {
    return R;
  }

  const float gamma_m_min = to_gamma(phi_m_min, a, b);
  float gamma_m_max = to_gamma(phi_m_max, a, b);
  if (gamma_m_max < gamma_m_min) {
    gamma_m_max += M_2PI_F;
  }

  /* Maximal sample resolution. */
  float res = roughness * .7f;
  /* Number of intervals should be even. */
  const size_t intervals = 2 * (size_t)ceilf((gamma_m_max - gamma_m_min) / res * 0.5f);

  /* Modified resolution based on numbers of intervals. */
  res = (gamma_m_max - gamma_m_min) / float(intervals);

  /* Integrate using Simpson's rule. */
  float integral = 0.0f;
  for (size_t i = 0; i <= intervals; i++) {

    const float gamma_m = gamma_m_min + i * res;
    const float3 wm = sphg_dir(tilt, gamma_m, a, b);

    if (microfacet_visible(wi, wo, make_float3(wm.x, 0.0f, wm.z), wh)) {
      const float weight = (i == 0 || i == intervals) ? 0.5f : (i % 2 + 1);
      const float arc_length = sqrtf(1.0f - e2 * sqr(sinf(gamma_m)));

      integral += weight * D(beckmann, roughness, wm, wh) *
                  G(beckmann, roughness, wi, wo, wm, wh) * arc_length;
    }
  }

  integral *= (2.0f / 3.0f * res);

  const float F = fresnel_dielectric_cos(dot(wi, wh), eta);
  const float d_o_inv = 1.0f / sqrtf(1.0f - e2 * sqr(sinf(phi_o)));

  R = make_float3(bsdf->extra->R * 0.125f * F * integral * d_o_inv);
  return R;
}

ccl_device float3 bsdf_microfacet_hair_eval_tt_trt_elliptic(KernelGlobals kg,
                                                            ccl_private const ShaderClosure *sc,
                                                            const float3 wi_,
                                                            const float3 wo_,
                                                            uint rng_quadrature)
{
  ccl_private MicrofacetHairBSDF *bsdf = (ccl_private MicrofacetHairBSDF *)sc;
  const float tilt = -bsdf->alpha;
  const float roughness = bsdf->roughness;
  const float eta = bsdf->eta;
  const bool beckmann = (bsdf->distribution_type == NODE_MICROFACET_HAIR_BECKMANN);

  if (bsdf->extra->TT <= 0.0f && bsdf->extra->TRT <= 0.0f) {
    return zero_float3();
  }

  /* this follows blender's convention (unlike the circular case?) */
  const float3 wo = wi_;
  const float3 wi = wo_;

  const float phi_i = dir_phi(wi);
  const float phi_o = dir_phi(wo);

  /* dot(wi, wmi) > 0 */
  const float tan_tilt = tanf(tilt);
  float phi_m_max = acosf(fmaxf(-tan_tilt * tan_theta(wi), 0.0f)) + phi_i;
  if (isnan_safe(phi_m_max)) {
    return zero_float3();
  }
  float phi_m_min = -phi_m_max + 2.0f * phi_i;

  /* dot(wo, wmo) < 0 */
  float tmp1 = acosf(fminf(tan_tilt * tan_theta(wo), 0.0f));
  if (isnan_safe(tmp1)) {
    return zero_float3();
  }

  const float3 mu_a = bsdf->sigma;
  const float inv_eta = 1.0f / eta;

  /* get elliptical cross section characteristic */
  const float a = 1.0f;
  const float b = bsdf->extra->aspect_ratio;
  const float e2 = 1.0f - sqr(b / a); /* Squared Eccentricity. */

  float gamma_m_min = to_gamma(phi_m_min, a, b) + 1e-3f;
  float gamma_m_max = to_gamma(phi_m_max, a, b) - 1e-3f;
  if (gamma_m_max < gamma_m_min) {
    gamma_m_max += M_2PI_F;
  }

  float res = roughness * 0.8f;
  const size_t intervals = 2 * (size_t)ceilf((gamma_m_max - gamma_m_min) / res * 0.5f);
  res = (gamma_m_max - gamma_m_min) / intervals;

  float3 S_tt = zero_float3();
  float3 S_trt = zero_float3();
  for (size_t i = 0; i <= intervals; i++) {

    const float gamma_mi = gamma_m_min + i * res;

    const float3 wmi = sphg_dir(tilt, gamma_mi, a, b);
    const float3 wmi_ = sphg_dir(0.0f, gamma_mi, a, b);

    /* sample wh1 */
    const float2 sample1 = make_float2(lcg_step_float(&rng_quadrature),
                                       lcg_step_float(&rng_quadrature));

    const float3 wh1 = sample_wh(kg, beckmann, roughness, wi, wmi, sample1);
    const float dot_wi_wh1 = dot(wi, wh1);
    if (!(dot_wi_wh1 > 0)) {
      continue;
    }

    float cos_theta_t1;
    const float T1 = 1.0f - fresnel(dot_wi_wh1, eta, &cos_theta_t1);

    /* refraction at the first interface */
    const float3 wt = -refract_angle(wi, wh1, cos_theta_t1, inv_eta);
    const float phi_t = dir_phi(wt);
    const float gamma_mt = 2.0f * to_phi(phi_t, a, b) - gamma_mi;
    const float3 wmt = sphg_dir(-tilt, gamma_mt, a, b);
    const float3 wmt_ = sphg_dir(0.0f, gamma_mt, a, b);

    const float G1 = G(beckmann, roughness, wi, -wt, wmi, wh1);
    if (G1 == 0.0f || !microfacet_visible(wi, -wt, wmi_, wh1)) {
      continue;
    }

    /* Simpson's rule weight */
    const float weight = (i == 0 || i == intervals) ? 0.5f : (i % 2 + 1);

    const float2 pi = to_point(gamma_mi, a, b);
    const float2 pt = to_point(gamma_mt + M_PI_F, a, b);
    /* TODO: check the definition of mu_a. */
    const float3 A_t = exp(-mu_a * len(pi - pt) / cos_theta(wt));

    /* TT */
    if (bsdf->extra->TT > 0.0f) {

      /* Total internal reflection otherwise. */
      if (dot(wo, wt) >= inv_eta - 1e-5f) {

        /* microfacet visiblity from macronormal */
        float3 wh2 = -wt + inv_eta * wo;
        if (dot(wmt, wh2) >= 0.0f) {

          const float rcp_norm_wh2 = 1.0f / len(wh2);
          wh2 *= rcp_norm_wh2;

          const float dot_wt_wh2 = dot(-wt, wh2);

          const float T2 = 1.0f - fresnel_dielectric_cos(dot_wt_wh2, inv_eta);
          const float D2 = D(beckmann, roughness, wh2, wmt) *
                           G(beckmann, roughness, -wt, -wo, wmt, wh2);

          const float3 result = T1 * T2 * D2 * A_t * dot_wt_wh2 * dot(wo, wh2) *
                                sqr(rcp_norm_wh2) / dot(wt, wmi) * weight *
                                smith_g1(beckmann, roughness, -wt, wmi, wh1) * dot(wi, wmi);

          if (isfinite_safe(result)) {
            const float arc_length = sqrtf(1.0f - e2 * sqr(sinf(gamma_mt)));
            S_tt += bsdf->extra->TT * result * arc_length;
          }
        }
      }
    }

    /* TRT */
    if (bsdf->extra->TRT > 0.0f) {

      /* sample wh2 */
      const float2 sample2 = make_float2(lcg_step_float(&rng_quadrature),
                                         lcg_step_float(&rng_quadrature));

      const float3 wh2 = sample_wh(kg, beckmann, roughness, -wt, wmt, sample2);

      const float cos_th2 = dot(-wt, wh2);
      if (!(cos_th2 > 0)) {
        continue;
      }

      const float R2 = fresnel_dielectric_cos(cos_th2, inv_eta);
      const float3 wtr = -reflect(wt, wh2);

      const float G2 = G(beckmann, roughness, -wt, -wtr, wmt, wh2);
      if (G2 == 0.0f || !microfacet_visible(-wt, -wtr, wmt_, wh2)) {
        continue;
      }

      /* total internal reflection */
      if (dot(-wtr, wo) < inv_eta - 1e-5f) {
        continue;
      }

      const float phi_tr = dir_phi(wtr);
      const float gamma_mtr = gamma_mi - 2.0f * (to_phi(phi_t, a, b) - to_phi(phi_tr, a, b)) +
                              M_PI_F;
      const float3 wmtr = sphg_dir(-tilt, gamma_mtr, a, b);
      const float3 wmtr_ = sphg_dir(0.0f, gamma_mtr, a, b);

      float3 wh3 = wtr + inv_eta * wo;
      const float G3 = G(beckmann, roughness, wtr, -wo, wmtr, wh3);
      if (dot(wmtr, wh3) < 0.0f || G3 == 0.0f || !microfacet_visible(wtr, -wo, wmtr_, wh3)) {
        continue;
      }

      const float rcp_norm_wh3 = 1.0f / len(wh3);
      wh3 *= rcp_norm_wh3;

      const float cos_trh3 = dot(wh3, wtr);

      const float T3 = 1.0f - fresnel_dielectric_cos(cos_trh3, inv_eta);

      const float D3 = D(beckmann, roughness, wh3, wmtr) * G3;

      const float2 ptr = to_point(gamma_mtr + M_PI_F, a, b);
      const float3 A_tr = exp(-mu_a * len(pt - ptr) / cos_theta(wtr));

      const float3 result = T1 * R2 * T3 * D3 * cos_trh3 * dot(wh3, wo) * sqr(rcp_norm_wh3) * A_t *
                            A_tr * weight / (dot(wt, wmi) * dot(wtr, wmt)) *
                            smith_g1(beckmann, roughness, -wt, wmi, wh1) *
                            smith_g1(beckmann, roughness, -wtr, wmt, wh2) * dot(wi, wmi) *
                            dot(wt, wmt);

      if (isfinite_safe(result)) {
        const float arc_length = sqrtf(1.0f - e2 * sqr(sin(gamma_mtr)));
        S_trt += bsdf->extra->TRT * result * arc_length;
      }
    }
  }

  const float d_o_inv = 1.0f / sqrtf(1.0f - e2 * sqr(sin(phi_o)));
  return (S_tt + S_trt) / 3.0f * res * sqr(inv_eta) * d_o_inv;
}

ccl_device Spectrum bsdf_microfacet_hair_eval_elliptic(KernelGlobals kg,
                                                       ccl_private const ShaderData *sd,
                                                       ccl_private const ShaderClosure *sc,
                                                       const float3 omega_in,
                                                       ccl_private float *pdf)
{
  /* get local wi(convention is reversed from other hair bcsdfs) */
  ccl_private MicrofacetHairBSDF *bsdf = (ccl_private MicrofacetHairBSDF *)sc;
  const float3 X = float4_to_float3(bsdf->extra->geom);
  float3 Y = sd->dPdu;
  const float3 Z = safe_normalize(cross(X, Y));
  Y = safe_normalize(cross(Z, X));

  const float3 wi = make_float3(dot(sd->I, X), dot(sd->I, Y), dot(sd->I, Z));
  const float3 wo = make_float3(dot(omega_in, X), dot(omega_in, Y), dot(omega_in, Z));

  /* Treat as transparent material if intersection lies outside of the projected radius. */
  const float e2 = 1.0f - sqr(bsdf->extra->aspect_ratio);
  const float radius = sqrtf(1.0f - e2 * sqr(sin_phi(wi)));
  if (fabsf(bsdf->extra->geom.w) > radius) {
    *pdf = 0.0f;
    return zero_spectrum();
  }

  /* evaluate */
  const float3 R = bsdf_microfacet_hair_eval_r_elliptic(sc, wi, wo) +
                   bsdf_microfacet_hair_eval_tt_trt_elliptic(kg, sc, wi, wo, sd->lcg_state);

  *pdf = 1.0f;

  return rgb_to_spectrum(R / cos_theta(wi));
}

ccl_device int bsdf_microfacet_hair_sample_elliptic(const KernelGlobals kg,
                                                    ccl_private const ShaderClosure *sc,
                                                    ccl_private ShaderData *sd,
                                                    float randu,
                                                    float randv,
                                                    ccl_private Spectrum *eval,
                                                    ccl_private float3 *omega_in,
                                                    ccl_private float *pdf,
                                                    ccl_private float2 *sampled_roughness,
                                                    ccl_private float *eta)
{
  ccl_private MicrofacetHairBSDF *bsdf = (ccl_private MicrofacetHairBSDF *)sc;
  *sampled_roughness = make_float2(bsdf->roughness, bsdf->roughness);
  *eta = bsdf->eta;
  const float inv_eta = 1.0f / *eta;

  /* get local wi (convention is reversed from other hair bcsdfs) */
  const float3 X = float4_to_float3(bsdf->extra->geom);
  float3 Y = sd->dPdu;
  const float3 Z = safe_normalize(cross(X, Y));
  Y = safe_normalize(cross(Z, X));

  const float3 wi = make_float3(dot(sd->I, X), dot(sd->I, Y), dot(sd->I, Z));

  /* get elliptical cross section characteristic */
  const float a = 1.0f;
  const float b = bsdf->extra->aspect_ratio;
  const float e2 = 1.0f - sqr(b / a);

  /* macronormal */
  const float2 sincos_phi_i = sincos_phi(wi);
  const float sin_phi_i = sincos_phi_i.x;
  const float cos_phi_i = sincos_phi_i.y;
  const float d_i = sqrtf(1.0f - e2 * sqr(sin_phi_i));

  /* Treat as transparent material if intersection lies outside of the projected radius. */
  if (fabsf(bsdf->extra->geom.w) > d_i) {
    *omega_in = -sd->I;
    *pdf = 1;
    *eval = one_spectrum();
    return LABEL_TRANSMIT | LABEL_TRANSPARENT;
  }

  if (bsdf->extra->R <= 0.0f && bsdf->extra->TT <= 0.0f && bsdf->extra->TRT <= 0.0f) {
    /* early out for inactive lobe */
    *pdf = 0.0f;
    return LABEL_NONE;
  }

  const float tilt = -bsdf->alpha;
  const float roughness = bsdf->roughness;
  const bool beckmann = (bsdf->distribution_type == NODE_MICROFACET_HAIR_BECKMANN);

  /* generate sample */
  float sample_lobe = randu;
  const float sample_h = randv;
  const float2 sample_h1 = make_float2(lcg_step_float(&sd->lcg_state),
                                       lcg_step_float(&sd->lcg_state));
  const float2 sample_h2 = make_float2(lcg_step_float(&sd->lcg_state),
                                       lcg_step_float(&sd->lcg_state));
  const float2 sample_h3 = make_float2(lcg_step_float(&sd->lcg_state),
                                       lcg_step_float(&sd->lcg_state));

  const float h = d_i * (sample_h * 2.0f - 1.0f);
  const float gamma_mi = atan2f(cos_phi_i, -b / a * sin_phi_i) -
                         acosf(h / sqrtf(sqr(cos_phi_i) + sqr(b / a * sin_phi_i)));
  float sin_gamma_mi, cos_gamma_mi;
  fast_sincosf(gamma_mi, &sin_gamma_mi, &cos_gamma_mi);
  const float3 wmi_ = normalize(make_float3(b * sin_gamma_mi, 0.0f, a * cos_gamma_mi));

  /* mesonormal */
  float st, ct;
  fast_sincosf(tilt, &st, &ct);
  const float3 wmi = make_float3(wmi_.x * ct, st, wmi_.z * ct);

  if (dot(wmi, wi) < 0.0f || dot(wmi_, wi) < 0.0f) {
    /* macro/mesonormal invisible */
    *pdf = 0.0f;
    return LABEL_NONE;
  }

  /* sample R lobe */
  const float3 wh1 = sample_wh(kg, beckmann, roughness, wi, wmi, sample_h1);

  const float3 wr = -reflect(wi, wh1);

  /* ensure that this is a valid sample */
  if (!(dot(wr, wh1) > 0.0f) || !(dot(wr, wmi) > 0.0f) || !microfacet_visible(wi, wr, wmi_, wh1)) {
    *pdf = 0.0f;
    return LABEL_NONE;
  }

  float3 TT = zero_float3();
  float3 TRT = zero_float3();

  float cos_theta_t1;
  const float R1 = fresnel(dot(wi, wh1), *eta, &cos_theta_t1);
  float3 R = make_float3(bsdf->extra->R * R1);

  /* sample TT lobe */
  const float3 wt = -refract_angle(wi, wh1, cos_theta_t1, inv_eta);
  const float phi_t = dir_phi(wt);

  const float gamma_mt = 2.0f * to_phi(phi_t, a, b) - gamma_mi;
  const float3 wmt = sphg_dir(-tilt, gamma_mt, a, b);
  const float3 wmt_ = sphg_dir(0.0f, gamma_mt, a, b);

  const float3 wh2 = sample_wh(kg, beckmann, roughness, -wt, wmt, sample_h2);

  const float3 wtr = -reflect(wt, wh2);

  float3 wh3;
  float3 wtt, wtrt;
  float3 wmtr, wmtr_;

  if (dot(wt, wh2) < 0.0f && dot(wmt, wt) < 0.0f &&
      microfacet_visible(-wt, -wtr, make_float3(wmt.x, 0.0f, wmt.z), wh2)) {

    const float3 mu_a = bsdf->sigma;
    const float2 pi = to_point(gamma_mi, a, b);
    const float2 pt = to_point(gamma_mt + M_PI_F, a, b);
    const float3 A_t = exp(-mu_a * len(pi - pt) / cos_theta(wt));

    float cos_theta_t2;
    const float R2 = fresnel(dot(-wt, wh2), inv_eta, &cos_theta_t2);
    const float3 T1 = make_float3(1.0f - R1);
    const float3 T2 = make_float3(1.0f - R2);

    wtt = -refract_angle(-wt, wh2, cos_theta_t2, *eta);

    if (dot(wtt, wmt) < 0.0f && cos_theta_t2 != 0.0f) {
      TT = bsdf->extra->TT * T1 * A_t * T2;
    }

    /* sample TRT lobe */
    const float phi_tr = dir_phi(wtr);
    float gamma_mtr = gamma_mi - 2.0f * (to_phi(phi_t, a, b) - to_phi(phi_tr, a, b)) + M_PI_F;
    wmtr = sphg_dir(-tilt, gamma_mtr, a, b);
    wmtr_ = sphg_dir(0.0f, gamma_mtr, a, b);

    wh3 = sample_wh(kg, beckmann, roughness, wtr, wmtr, sample_h3);

    float cos_theta_t3;
    const float R3 = fresnel(dot(wtr, wh3), inv_eta, &cos_theta_t3);

    wtrt = -refract_angle(wtr, wh3, cos_theta_t3, *eta);

    if (cos_theta_t3 != 0.0f && dot(wtr, wh3) > 0.0f && dot(wmtr, wtr) > 0.0f &&
        dot(wtrt, wmtr) < 0.0f &&
        microfacet_visible(wtr, -wtrt, make_float3(wmtr.x, 0.0f, wmtr.z), wh3)) {

      const float3 T3 = make_float3(1.0f - R3);

      const float2 ptr = to_point(gamma_mtr + M_PI_F, a, b);
      const float3 A_tr = exp(-mu_a * len(pt - ptr) / cos_theta(wtr));

      TRT = bsdf->extra->TRT * T1 * R2 * T3 * A_t * A_tr;
    }
  }

  /* select lobe based on energy */
  const float r = average(R);
  const float tt = average(TT);
  const float trt = average(TRT);
  const float total_energy = r + tt + trt;

  if (total_energy == 0.0f) {
    *pdf = 0.0f;
    return LABEL_NONE;
  }

  float3 wo;
  float visibility = 0.0f;
  int label = LABEL_GLOSSY;

  sample_lobe *= total_energy;
  if (sample_lobe < r) {
    wo = wr;
    *eval = rgb_to_spectrum(R / r * total_energy);

    if (microfacet_visible(wi, wr, wmi_, wh1)) {
      visibility = smith_g1(beckmann, roughness, wr, wmi, wh1);
    }

    label |= LABEL_REFLECT;
  }
  else if (sample_lobe < (r + tt)) {
    wo = wtt;
    *eval = rgb_to_spectrum(TT / tt * total_energy);

    if (microfacet_visible(wi, -wt, wmi_, wh1) && microfacet_visible(-wt, -wtt, wmt_, wh2)) {
      visibility = smith_g1(beckmann, roughness, -wt, wmi, wh1) *
                   smith_g1(beckmann, roughness, -wtt, wmt, wh2);
    }

    label |= LABEL_TRANSMIT;
  }
  else { /* if (sample_lobe >= (r + tt)) */
    wo = wtrt;
    *eval = rgb_to_spectrum(TRT / trt * total_energy);

    if (microfacet_visible(wi, -wt, wmi_, wh1) && microfacet_visible(-wt, -wtr, wmt_, wh2) &&
        microfacet_visible(wtr, -wtrt, wmtr_, wh3)) {
      visibility = smith_g1(beckmann, roughness, -wt, wmi, wh1) *
                   smith_g1(beckmann, roughness, -wtr, wmt, wh2) *
                   smith_g1(beckmann, roughness, -wtrt, wmtr, wh3);
    }

    label |= LABEL_TRANSMIT;
  }

  *eval *= visibility;
  *omega_in = wo.x * X + wo.y * Y + wo.z * Z;

  /* correction of the cosine foreshortening term
   *eval *= dot(wi, wmi) / dot(wi, wmi_); */

  /* ensure the same pdf is returned for BSDF and emitter sampling */
  *pdf = 1.0f;

  return label;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main sample and eval functions selecting model
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ccl_device Spectrum bsdf_microfacet_hair_eval(KernelGlobals kg,
                                              ccl_private const ShaderData *sd,
                                              ccl_private const ShaderClosure *sc,
                                              const float3 omega_in,
                                              ccl_private float *pdf)
{
  ccl_private MicrofacetHairBSDF *bsdf = (ccl_private MicrofacetHairBSDF *)sc;

  if (bsdf->cross_section_type == NODE_MICROFACET_HAIR_CIRCULAR ||
      bsdf->extra->aspect_ratio == 1.0f) {
    return bsdf_microfacet_hair_eval_circular(kg, sd, sc, omega_in, pdf);
  }

  return bsdf_microfacet_hair_eval_elliptic(kg, sd, sc, omega_in, pdf);
}

ccl_device int bsdf_microfacet_hair_sample(KernelGlobals kg,
                                           ccl_private const ShaderClosure *sc,
                                           ccl_private ShaderData *sd,
                                           float randu,
                                           float randv,
                                           ccl_private Spectrum *eval,
                                           ccl_private float3 *omega_in,
                                           ccl_private float *pdf,
                                           ccl_private float2 *sampled_roughness,
                                           ccl_private float *eta)
{
  ccl_private MicrofacetHairBSDF *bsdf = (ccl_private MicrofacetHairBSDF *)sc;

  if (bsdf->cross_section_type == NODE_MICROFACET_HAIR_CIRCULAR ||
      bsdf->extra->aspect_ratio == 1.0f) {
    return bsdf_microfacet_hair_sample_circular(
        kg, sc, sd, randu, randv, eval, omega_in, pdf, sampled_roughness, eta);
  }

  return bsdf_microfacet_hair_sample_elliptic(
      kg, sc, sd, randu, randv, eval, omega_in, pdf, sampled_roughness, eta);
}

/* Implements Filter Glossy by capping the effective roughness. */
ccl_device void bsdf_microfacet_hair_blur(ccl_private ShaderClosure *sc, float roughness)
{
  ccl_private MicrofacetHairBSDF *bsdf = (ccl_private MicrofacetHairBSDF *)sc;

  if (bsdf->blur > 0) {
    bsdf->roughness = fmaxf(roughness, bsdf->roughness);
  }
}

/* Hair Albedo */

ccl_device_inline float bsdf_microfacet_hair_albedo_roughness_scale(
    const float azimuthal_roughness)
{
  const float x = azimuthal_roughness;
  return (((((0.245f * x) + 5.574f) * x - 10.73f) * x + 2.532f) * x - 0.215f) * x + 5.969f;
}

ccl_device float3 bsdf_microfacet_hair_albedo(ccl_private const ShaderClosure *sc)
{
  ccl_private MicrofacetHairBSDF *bsdf = (ccl_private MicrofacetHairBSDF *)sc;
  return exp(-sqrt(bsdf->sigma) * bsdf_microfacet_hair_albedo_roughness_scale(bsdf->roughness));
}

ccl_device_inline float3
bsdf_microfacet_hair_sigma_from_reflectance(const float3 color, const float azimuthal_roughness)
{
  const float3 sigma = log(color) /
                       bsdf_microfacet_hair_albedo_roughness_scale(azimuthal_roughness);
  return sigma * sigma;
}

ccl_device_inline float3 bsdf_microfacet_hair_sigma_from_concentration(const float eumelanin,
                                                                       const float pheomelanin)
{
  return eumelanin * make_float3(0.506f, 0.841f, 1.653f) +
         pheomelanin * make_float3(0.343f, 0.733f, 1.924f);
}

CCL_NAMESPACE_END
