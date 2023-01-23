/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2018-2022 Blender Foundation */

/* This code implements the paper [A Microfacet-based Hair Scattering
 * Model](https://onlinelibrary.wiley.com/doi/full/10.1111/cgf.14588) by Weizhen Huang, Matthias B.
 * Hullin and Johannes Hanika. */

#pragma once

#ifndef __KERNEL_GPU__
#  include <fenv.h>
#endif

#include "kernel/util/color.h"

CCL_NAMESPACE_BEGIN

typedef struct MicrofacetHairExtra {
  /* TODO: is this necessary? */
  float R;
  float TT;
  float TRT;

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
  float tilt;
  /* IOR. */
  float eta;

  /* GGX/Beckmann. */
  int distribution_type;

  /* The ratio of the minor axis to the major axis. */
  float aspect_ratio;

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

  /* Compute local frame. The Y axis is aligned with the curve tangent; the X axis is perpendicular
   to the ray direction for circular cross-sections, or aligned with the major axis for elliptical
   cross-sections. */
  const float3 Y = safe_normalize(sd->dPdu);
  const float3 X = safe_normalize(cross(Y, sd->wi));

  /* h -1..0..1 means the rays goes from grazing the hair, to hitting it at the center, to grazing
   * the other edge. This is the cosine of the angle between sd->N and X. */
  const float h = (sd->type & PRIMITIVE_CURVE_RIBBON) ? -sd->v : -dot(X, sd->N);

  kernel_assert(fabsf(h) < 1.0f + 1e-4f);
  kernel_assert(isfinite_safe(X));
  kernel_assert(isfinite_safe(h));

  if (bsdf->aspect_ratio != 1.0f) {
    if (bsdf->aspect_ratio > 1.0f) {
      bsdf->aspect_ratio = 1.0f / bsdf->aspect_ratio;

      /* Switch major and minor axis. */
      const float3 minor_axis = safe_normalize(cross(
          sd->dPdu, make_float3(bsdf->extra->geom.x, bsdf->extra->geom.y, bsdf->extra->geom.z)));
      const float3 major_axis = safe_normalize(cross(minor_axis, sd->dPdu));

      bsdf->extra->geom = make_float4(major_axis.x, major_axis.y, major_axis.z, h);
    }
    else {
      bsdf->extra->geom.w = h;
    }
  }
  else {
    /* Align local frame with the ray direction so that `phi_i == 0`. */
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

/* Conversion between gamma and phi. Notations see Figure 5 in the paper. */
ccl_device_inline float to_phi(float gamma, float b)
{
  if (b == 1.0f) {
    return gamma;
  }
  float sin_gamma, cos_gamma;
  fast_sincosf(gamma, &sin_gamma, &cos_gamma);
  return atan2f(b * sin_gamma, cos_gamma);
}

ccl_device_inline float to_gamma(float phi, float b)
{
  if (b == 1.0f) {
    return phi;
  }
  float sin_phi, cos_phi;
  fast_sincosf(phi, &sin_phi, &cos_phi);
  return atan2f(sin_phi, b * cos_phi);
}

/* Compute the coordinate on the ellipse, given gamma and the aspect ratio between the minor axis
 * and the major axis. */
ccl_device_inline float2 to_point(float gamma, float b)
{
  float sin_gamma, cos_gamma;
  fast_sincosf(gamma, &sin_gamma, &cos_gamma);
  return make_float2(sin_gamma, b * cos_gamma);
}

/* Compute the vector direction given by theta and gamma. */
ccl_device_inline float3 sphg_dir(float theta, float gamma, float b)
{
  float sin_theta, cos_theta, sin_gamma, cos_gamma, sin_phi, cos_phi;

  fast_sincosf(theta, &sin_theta, &cos_theta);
  fast_sincosf(gamma, &sin_gamma, &cos_gamma);

  if (b == 1.0f) {
    sin_phi = sin_gamma;
    cos_phi = cos_gamma;
  }
  else {
    float tan_gamma = sin_gamma / cos_gamma;
    float tan_phi = b * tan_gamma;
    cos_phi = signf(cos_gamma) / sqrtf(sqr(tan_phi) + 1.0f);
    sin_phi = cos_phi * tan_phi;
  }
  return make_float3(sin_phi * cos_theta, sin_theta, cos_phi * cos_theta);
}

ccl_device_inline float arc_length(float e2, float gamma)
{
  return e2 == 0 ? 1.0f : sqrtf(1.0f - e2 * sqr(sinf(gamma)));
}

ccl_device_inline float projected_radius(float e2, float phi)
{
  return e2 == 0 ? 1.0f : sqrtf(1.0f - e2 * sqr(sinf(phi)));
}

/** \} */

/* Sample microfacets from a tilted mesonormal. */
template<MicrofacetType m_type>
ccl_device_inline float3 sample_wh(
    KernelGlobals kg, const float roughness, const float3 wi, const float3 wm, const float2 rand)
{
  /* Coordinate transformation for microfacet sampling. */
  float3 s, t;
  const float3 n = wm;
  make_orthonormals(n, &s, &t);

  const float3 wi_wm = make_float3(dot(wi, s), dot(wi, t), dot(wi, n));

  float discard;
  const float3 wh_wm = microfacet_sample_stretched<m_type>(
      kg, wi_wm, roughness, roughness, rand.x, rand.y, &discard);

  const float3 wh = wh_wm.x * s + wh_wm.y * t + wh_wm.z * n;
  return wh;
}

/* Check micronormal/mesonormal direct visiblity from direction v. */
ccl_device_inline bool microfacet_visible(const float3 v, const float3 m, const float3 h)
{
  return (dot(v, h) > 0.0f && dot(v, m) > 0.0f);
}

/* Check micronormal/mesonormal direct visiblity from directinos wi and wo. */
ccl_device_inline bool microfacet_visible(const float3 wi,
                                          const float3 wo,
                                          const float3 m,
                                          const float3 h)
{
  return microfacet_visible(wi, m, h) && microfacet_visible(wo, m, h);
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
  const float cos_theta_t_sqr = 1.0f - (1.0f - cos_theta_i * cos_theta_i) / (eta * eta);
  *cos_theta_t = safe_sqrtf(cos_theta_t_sqr);

  if (cos_theta_t_sqr <= 0) {
    /* Total internal reflection. */
    return 1.0f;
  }

  /* Amplitudes of reflected waves. */
  const float a_s = (cos_theta_i - eta * (*cos_theta_t)) / (cos_theta_i + eta * (*cos_theta_t));
  const float a_p = (*cos_theta_t - eta * cos_theta_i) / (*cos_theta_t + eta * cos_theta_i);

  /* Adjust the sign of the transmitted direction to be relative to the surface normal. */
  *cos_theta_t = -(*cos_theta_t);

  return 0.5f * (sqr(a_s) + sqr(a_p));
}

/* Refract the incident ray, given the cosine of the refraction angle and the inverse IOR. */
ccl_device_inline float3 refract_angle(const float3 incident,
                                       const float3 normal,
                                       const float cos_theta_t,
                                       const float inv_eta)
{
  return inv_eta * incident - (inv_eta * dot(normal, incident) + cos_theta_t) * normal;
}

template<MicrofacetType m_type>
ccl_device float3 bsdf_microfacet_hair_eval_r(ccl_private const ShaderClosure *sc,
                                              const float3 wi,
                                              const float3 wo)
{
  ccl_private MicrofacetHairBSDF *bsdf = (ccl_private MicrofacetHairBSDF *)sc;
  const float tilt = -bsdf->tilt;
  const float roughness = bsdf->roughness;
  const float roughness2 = sqr(roughness);
  const float eta = bsdf->eta;

  if (bsdf->extra->R <= 0.0f) {
    return zero_float3();
  }

  /* Get elliptical cross section characteristic. Assuming major axis is 1. */
  const float b = bsdf->aspect_ratio;
  const float e2 = 1.0f - sqr(b); /* Squared Eccentricity. */
  const bool is_circular = (b == 1.0f);

  const float phi_i = is_circular ? 0.0f : dir_phi(wi);
  const float phi_o = dir_phi(wo);
  const float3 wh = normalize(wi + wo);

  /* dot(wi, wmi) > 0 */
  const float tan_tilt = tanf(tilt);
  float phi_m_max1 = acosf(fmaxf(-tan_tilt * tan_theta(wi), 0.0f)) + phi_i;
  if (isnan_safe(phi_m_max1)) {
    return zero_float3();
  }
  float phi_m_min1 = -phi_m_max1 + 2.0f * phi_i;

  /* dot(wo, wmi) > 0 */
  float phi_m_max2 = acosf(fmaxf(-tan_tilt * tan_theta(wo), 0.0f)) + phi_o;
  if (isnan_safe(phi_m_max2)) {
    return zero_float3();
  }
  float phi_m_min2 = -phi_m_max2 + 2.0f * phi_o;

  if (!is_circular) {
    /* Try to wrap range. */
    if ((phi_m_max2 - phi_m_min1) > M_2PI_F) {
      phi_m_min2 -= M_2PI_F;
      phi_m_max2 -= M_2PI_F;
    }
    if ((phi_m_max1 - phi_m_min2) > M_2PI_F) {
      phi_m_min1 -= M_2PI_F;
      phi_m_max1 -= M_2PI_F;
    }
  }

  const float phi_m_min = fmaxf(phi_m_min1, phi_m_min2) + 1e-3f;
  const float phi_m_max = fminf(phi_m_max1, phi_m_max2) - 1e-3f;
  if (phi_m_min > phi_m_max) {
    return zero_float3();
  }

  const float gamma_m_min = to_gamma(phi_m_min, b);
  float gamma_m_max = to_gamma(phi_m_max, b);
  if (gamma_m_max < gamma_m_min) {
    gamma_m_max += M_2PI_F;
  }

  /* Maximal sample resolution. */
  float res = roughness * 0.7f;
  /* Number of intervals should be even. */
  const size_t intervals = 2 * (size_t)ceilf((gamma_m_max - gamma_m_min) / res * 0.5f);

  /* Modified resolution based on numbers of intervals. */
  res = (gamma_m_max - gamma_m_min) / float(intervals);

  /* Integrate using Composite Simpson's 1/3 rule. */
  float integral = 0.0f;
  for (size_t i = 0; i <= intervals; i++) {

    const float gamma_m = gamma_m_min + i * res;
    const float3 wm = sphg_dir(tilt, gamma_m, b);

    if (microfacet_visible(wi, wo, make_float3(wm.x, 0.0f, wm.z), wh)) {
      const float weight = (i == 0 || i == intervals) ? 0.5f : (i % 2 + 1);
      integral += weight * bsdf_D<m_type>(roughness2, dot(wm, wh)) *
                  bsdf_G<m_type>(roughness2, dot(wi, wm), dot(wo, wm)) * arc_length(e2, gamma_m);
    }
  }

  integral *= (2.0f / 3.0f * res);

  const float F = fresnel_dielectric_cos(dot(wi, wh), eta);

  return make_float3(bsdf->extra->R * 0.125f * F * integral / projected_radius(e2, phi_i));
}

template<MicrofacetType m_type>
ccl_device float3 bsdf_microfacet_hair_eval_tt_trt(KernelGlobals kg,
                                                   ccl_private const ShaderClosure *sc,
                                                   const float3 wi,
                                                   const float3 wo,
                                                   uint rng_quadrature)
{
  ccl_private MicrofacetHairBSDF *bsdf = (ccl_private MicrofacetHairBSDF *)sc;
  const float tilt = -bsdf->tilt;
  const float roughness = bsdf->roughness;
  const float roughness2 = sqr(roughness);
  const float eta = bsdf->eta;

  if (bsdf->extra->TT <= 0.0f && bsdf->extra->TRT <= 0.0f) {
    return zero_float3();
  }

  /* Get elliptical cross section characteristic. Assuming major axis is 1. */
  const float b = bsdf->aspect_ratio;
  const float e2 = 1.0f - sqr(b); /* Squared Eccentricity. */
  const bool is_circular = (b == 1.0f);

  const float phi_i = is_circular ? 0.0f : dir_phi(wi);

  const float tan_tilt = tanf(tilt);
  const float phi_m_max = acosf(fmaxf(-tan_tilt * tan_theta(wi), 0.0f)) + phi_i;
  if (isnan_safe(phi_m_max)) {
    /* Early detection of dot(wi, wmi) < 0. */
    return zero_float3();
  }
  const float phi_m_min = -phi_m_max + 2.0f * phi_i;

  if (tan_tilt * tan_theta(wo) < -1.0f) {
    /* Early detection of dot(wo, wmo) < 0. */
    return zero_float3();
  }

  const float3 mu_a = bsdf->sigma;
  const float inv_eta = 1.0f / eta;

  const float gamma_m_min = to_gamma(phi_m_min, b) + 1e-3f;
  float gamma_m_max = to_gamma(phi_m_max, b) - 1e-3f;
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

    const float3 wmi = sphg_dir(tilt, gamma_mi, b);
    const float3 wmi_ = sphg_dir(0.0f, gamma_mi, b);

    /* Sample wh1. */
    const float2 sample1 = make_float2(lcg_step_float(&rng_quadrature),
                                       lcg_step_float(&rng_quadrature));

    const float3 wh1 = sample_wh<m_type>(kg, roughness, wi, wmi, sample1);
    const float cos_hi1 = dot(wi, wh1);
    if (!(cos_hi1 > 0)) {
      continue;
    }

    float cos_theta_t1;
    const float T1 = 1.0f - fresnel(cos_hi1, eta, &cos_theta_t1);

    /* Refraction at the first interface. */
    const float3 wt = -refract_angle(wi, wh1, cos_theta_t1, inv_eta);
    const float phi_t = dir_phi(wt);
    const float gamma_mt = 2.0f * to_phi(phi_t, b) - gamma_mi;
    const float3 wmt = sphg_dir(-tilt, gamma_mt, b);
    const float3 wmt_ = sphg_dir(0.0f, gamma_mt, b);

    const float cos_mi1 = dot(wi, wmi);
    const float cos_mo1 = dot(-wt, wmi);
    const float cos_mi2 = dot(-wt, wmt);
    const float G1 = bsdf_G<m_type>(roughness2, cos_mi1, cos_mo1);
    if (G1 == 0.0f || !microfacet_visible(wi, -wt, wmi_, wh1)) {
      continue;
    }

    const float weight = (i == 0 || i == intervals) ? 0.5f : (i % 2 + 1);

    const float3 A_t = exp(mu_a / cos_theta(wt) *
                           (is_circular ?
                                2.0f * cosf(gamma_mi - phi_t) :
                                -len(to_point(gamma_mi, b) - to_point(gamma_mt + M_PI_F, b))));

    /* TT */
    if (bsdf->extra->TT > 0.0f) {
      if (dot(wo, wt) >= inv_eta - 1e-5f) { /* Total internal reflection otherwise. */
        float3 wh2 = -wt + inv_eta * wo;
        const float rcp_norm_wh2 = 1.0f / len(wh2);
        wh2 *= rcp_norm_wh2;
        const float cos_mh2 = dot(wmt, wh2);
        if (cos_mh2 >= 0.0f) { /* Microfacet visiblity from macronormal. */
          const float cos_hi2 = dot(-wt, wh2);
          const float cos_ho2 = dot(-wo, wh2);
          const float cos_mo2 = dot(-wo, wmt);

          const float T2 = 1.0f - fresnel_dielectric_cos(cos_hi2, inv_eta);
          const float D2 = bsdf_D<m_type>(roughness2, cos_mh2);
          const float G2 = bsdf_G<m_type>(roughness2, cos_mi2, cos_mo2);

          const float3 result = weight * T1 * T2 * D2 * G2 * A_t / cos_mo1 * cos_mi1 * cos_hi2 *
                                cos_ho2 * sqr(rcp_norm_wh2) * bsdf_G1<m_type>(roughness2, cos_mo1);

          if (isfinite_safe(result)) {
            S_tt += bsdf->extra->TT * result * arc_length(e2, gamma_mt);
          }
        }
      }
    }

    /* TRT */
    if (bsdf->extra->TRT > 0.0f) {
      /* Sample wh2. */
      const float2 sample2 = make_float2(lcg_step_float(&rng_quadrature),
                                         lcg_step_float(&rng_quadrature));
      const float3 wh2 = sample_wh<m_type>(kg, roughness, -wt, wmt, sample2);
      const float cos_hi2 = dot(-wt, wh2);
      if (!(cos_hi2 > 0)) {
        continue;
      }
      const float R2 = fresnel_dielectric_cos(cos_hi2, inv_eta);

      const float3 wtr = -reflect(wt, wh2);
      if (dot(-wtr, wo) < inv_eta - 1e-5f) {
        /* Total internal reflection. */
        continue;
      }
      const float cos_mo2 = dot(-wtr, wmt);
      const float G2 = bsdf_G<m_type>(roughness2, cos_mi2, cos_mo2);
      if (G2 == 0.0f || !microfacet_visible(-wt, -wtr, wmt_, wh2)) {
        continue;
      }

      const float phi_tr = dir_phi(wtr);
      const float gamma_mtr = gamma_mi - 2.0f * (to_phi(phi_t, b) - to_phi(phi_tr, b)) + M_PI_F;
      const float3 wmtr = sphg_dir(-tilt, gamma_mtr, b);
      const float3 wmtr_ = sphg_dir(0.0f, gamma_mtr, b);

      const float G3 = bsdf_G<m_type>(roughness2, dot(wtr, wmtr), dot(-wo, wmtr));
      float3 wh3 = wtr + inv_eta * wo;
      const float rcp_norm_wh3 = 1.0f / len(wh3);
      wh3 *= rcp_norm_wh3;
      const float cos_mh3 = dot(wmtr, wh3);
      if (cos_mh3 < 0.0f || G3 == 0.0f || !microfacet_visible(wtr, -wo, wmtr_, wh3)) {
        continue;
      }

      const float cos_hi3 = dot(wh3, wtr);
      const float cos_ho3 = dot(wh3, -wo);

      const float T3 = 1.0f - fresnel_dielectric_cos(cos_hi3, inv_eta);
      const float D3 = bsdf_D<m_type>(roughness2, cos_mh3);

      const float3 A_tr = exp(mu_a / cos_theta(wtr) *
                              (is_circular ?
                                   2.0f * cosf(phi_tr - gamma_mt) :
                                   -len(to_point(gamma_mtr, b) - to_point(gamma_mt, b))));

      const float3 result = weight * T1 * R2 * T3 * D3 * G3 * A_t * A_tr / (cos_mo1 * cos_mo2) *
                            cos_mi1 * cos_mi2 * cos_hi3 * cos_ho3 * sqr(rcp_norm_wh3) *
                            bsdf_G<m_type>(roughness2, cos_mo1, cos_mo2);

      if (isfinite_safe(result)) {
        S_trt += bsdf->extra->TRT * result * arc_length(e2, gamma_mtr);
      }
    }
  }

  return (S_tt + S_trt) / 3.0f * res * sqr(inv_eta) * projected_radius(e2, phi_i);
}

template<MicrofacetType m_type>
ccl_device int bsdf_microfacet_hair_sample(const KernelGlobals kg,
                                           ccl_private const ShaderClosure *sc,
                                           ccl_private ShaderData *sd,
                                           float randu,
                                           float randv,
                                           ccl_private Spectrum *eval,
                                           ccl_private float3 *wo,
                                           ccl_private float *pdf,
                                           ccl_private float2 *sampled_roughness,
                                           ccl_private float *eta)
{
  ccl_private MicrofacetHairBSDF *bsdf = (ccl_private MicrofacetHairBSDF *)sc;

  *sampled_roughness = make_float2(bsdf->roughness, bsdf->roughness);
  *eta = bsdf->eta;
  const float inv_eta = 1.0f / *eta;

  if (bsdf->extra->R <= 0.0f && bsdf->extra->TT <= 0.0f && bsdf->extra->TRT <= 0.0f) {
    /* Early out for inactive lobe. */
    *pdf = 0.0f;
    return LABEL_NONE;
  }

  /* Get local coordinate system:
   * . X major axis.
   * . Y along the fiber tangent.
   * . Z minor axis. */
  const float3 X = float4_to_float3(bsdf->extra->geom);
  const float3 Z = safe_normalize(cross(X, sd->dPdu));
  const float3 Y = safe_normalize(cross(Z, X));

  /* Transform wi from global coordinate system to local. */
  const float3 wi = make_float3(dot(sd->wi, X), dot(sd->wi, Y), dot(sd->wi, Z));

  /* Get elliptical cross section characteristic. Assuming major axis is 1. */
  const float b = bsdf->aspect_ratio;
  const float e2 = 1.0f - sqr(b); /* Squared Eccentricity. */
  const bool is_circular = (b == 1.0f);

  /* Macronormal. */
  const float2 sincos_phi_i = sincos_phi(wi);
  const float sin_phi_i = sincos_phi_i.x;
  const float cos_phi_i = sincos_phi_i.y;
  const float d_i = projected_radius(e2, sin_phi_i);

  /* Treat as transparent material if intersection lies outside of the projected radius. */
  if (fabsf(bsdf->extra->geom.w) > d_i) {
    *wo = -sd->wi;
    *pdf = 1;
    *eval = one_spectrum();
    return LABEL_TRANSMIT | LABEL_TRANSPARENT;
  }

  const float tilt = -bsdf->tilt;
  const float roughness = bsdf->roughness;
  const float roughness2 = sqr(roughness);

  /* Generate samples. */
  float sample_lobe = randu;
  const float sample_h = randv;
  const float2 sample_h1 = make_float2(lcg_step_float(&sd->lcg_state),
                                       lcg_step_float(&sd->lcg_state));
  const float2 sample_h2 = make_float2(lcg_step_float(&sd->lcg_state),
                                       lcg_step_float(&sd->lcg_state));
  const float2 sample_h3 = make_float2(lcg_step_float(&sd->lcg_state),
                                       lcg_step_float(&sd->lcg_state));

  const float h = sample_h * 2.0f - 1.0f;
  const float gamma_mi = is_circular ?
                             asinf(h) :
                             atan2f(cos_phi_i, -b * sin_phi_i) -
                                 acosf(h * d_i / sqrtf(sqr(cos_phi_i) + sqr(b * sin_phi_i)));

  const float3 wmi_ = sphg_dir(0, gamma_mi, b); /* Macronormal. */

  /* Mesonormal. */
  float st, ct;
  fast_sincosf(tilt, &st, &ct);
  const float3 wmi = make_float3(wmi_.x * ct, st, wmi_.z * ct);

  if (dot(wmi, wi) < 0.0f || dot(wmi_, wi) < 0.0f) {
    /* Macro/mesonormal invisible. */
    *pdf = 0.0f;
    return LABEL_NONE;
  }

  /* Sample R lobe. */
  const float3 wh1 = sample_wh<m_type>(kg, roughness, wi, wmi, sample_h1);
  const float3 wr = -reflect(wi, wh1);

  /* Ensure that this is a valid sample. */
  if (!(dot(wr, wh1) > 0.0f) || !(dot(wr, wmi) > 0.0f) || !microfacet_visible(wi, wr, wmi_, wh1)) {
    *pdf = 0.0f;
    return LABEL_NONE;
  }

  float cos_theta_t1;
  const float R1 = fresnel(dot(wi, wh1), *eta, &cos_theta_t1);
  const float3 R = make_float3(bsdf->extra->R * R1);

  /* Sample TT lobe. */
  const float3 wt = -refract_angle(wi, wh1, cos_theta_t1, inv_eta);
  const float phi_t = dir_phi(wt);

  const float gamma_mt = 2.0f * to_phi(phi_t, b) - gamma_mi;
  const float3 wmt = sphg_dir(-tilt, gamma_mt, b);
  const float3 wmt_ = sphg_dir(0.0f, gamma_mt, b);

  const float3 wh2 = sample_wh<m_type>(kg, roughness, -wt, wmt, sample_h2);

  const float3 wtr = -reflect(wt, wh2);

  float3 wh3;
  float3 wtt, wtrt;
  float3 wmtr, wmtr_;
  float3 TT = zero_float3();
  float3 TRT = zero_float3();

  if (dot(wt, wh2) < 0.0f && dot(wmt, wt) < 0.0f && microfacet_visible(-wt, -wtr, wmt_, wh2)) {
    const float3 mu_a = bsdf->sigma;
    const float3 A_t = exp(mu_a / cos_theta(wt) *
                           (is_circular ?
                                2.0f * cosf(phi_t - gamma_mi) :
                                -len(to_point(gamma_mi, b) - to_point(gamma_mt + M_PI_F, b))));

    float cos_theta_t2;
    const float R2 = fresnel(dot(-wt, wh2), inv_eta, &cos_theta_t2);
    const float3 T1 = make_float3(1.0f - R1);
    const float3 T2 = make_float3(1.0f - R2);

    wtt = -refract_angle(-wt, wh2, cos_theta_t2, *eta);

    if (dot(wtt, wmt) < 0.0f && cos_theta_t2 != 0.0f) {
      TT = bsdf->extra->TT * T1 * A_t * T2;
    }

    /* Sample TRT lobe. */
    const float phi_tr = dir_phi(wtr);
    const float gamma_mtr = gamma_mi - 2.0f * (to_phi(phi_t, b) - to_phi(phi_tr, b)) + M_PI_F;
    wmtr = sphg_dir(-tilt, gamma_mtr, b);

    wh3 = sample_wh<m_type>(kg, roughness, wtr, wmtr, sample_h3);

    float cos_theta_t3;
    const float R3 = fresnel(dot(wtr, wh3), inv_eta, &cos_theta_t3);

    wtrt = -refract_angle(wtr, wh3, cos_theta_t3, *eta);

    if (cos_theta_t3 != 0.0f && dot(wtr, wh3) > 0.0f && dot(wmtr, wtr) > 0.0f &&
        dot(wtrt, wmtr) < 0.0f &&
        microfacet_visible(wtr, -wtrt, make_float3(wmtr.x, 0.0f, wmtr.z), wh3)) {
      const float3 T3 = make_float3(1.0f - R3);

      const float3 A_tr = exp(mu_a / cos_theta(wtr) *
                              (is_circular ?
                                   2.0f * cos(phi_tr - gamma_mt) :
                                   -len(to_point(gamma_mt, b) - to_point(gamma_mtr, b))));

      TRT = bsdf->extra->TRT * T1 * R2 * T3 * A_t * A_tr;
    }
  }

  /* Select lobe based on energy. */
  const float r = average(R);
  const float tt = average(TT);
  const float trt = average(TRT);
  const float total_energy = r + tt + trt;

  if (total_energy == 0.0f) {
    *pdf = 0.0f;
    return LABEL_NONE;
  }

  float3 local_O;
  float visibility = 0.0f;
  int label = LABEL_GLOSSY;

  sample_lobe *= total_energy;
  if (sample_lobe < r) {
    local_O = wr;
    *eval = rgb_to_spectrum(R / r * total_energy);

    if (microfacet_visible(wi, wr, wmi_, wh1)) {
      visibility = bsdf_G1<m_type>(roughness2, dot(wr, wmi));
    }

    label |= LABEL_REFLECT;
  }
  else if (sample_lobe < (r + tt)) {
    local_O = wtt;
    *eval = rgb_to_spectrum(TT / tt * total_energy);

    if (microfacet_visible(wi, -wt, wmi_, wh1) && microfacet_visible(-wt, -wtt, wmt_, wh2)) {
      visibility = bsdf_G<m_type>(roughness2, dot(-wt, wmi), dot(-wtt, wmt));
    }

    label |= LABEL_TRANSMIT;
  }
  else { /* if (sample_lobe >= (r + tt)) */
    local_O = wtrt;
    *eval = rgb_to_spectrum(TRT / trt * total_energy);

    if (microfacet_visible(wi, -wt, wmi_, wh1)) {
      visibility = bsdf_G<m_type>(roughness2, dot(-wt, wmi), dot(-wtr, wmt)) *
                   bsdf_G1<m_type>(roughness2, dot(-wtrt, wmtr));
    }

    label |= LABEL_TRANSMIT;
  }

  *eval *= visibility;
  *wo = local_O.x * X + local_O.y * Y + local_O.z * Z;

  /* Ensure the same pdf is returned for BSDF and emitter sampling. The importance sampling pdf is
   * already factored in the value so this value is only used for MIS. */
  *pdf = 1.0f;

  return label;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main sample and eval functions selecting model
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ccl_device Spectrum bsdf_microfacet_hair_eval(KernelGlobals kg,
                                              ccl_private const ShaderData *sd,
                                              ccl_private const ShaderClosure *sc,
                                              const float3 wo,
                                              ccl_private float *pdf)
{
  ccl_private MicrofacetHairBSDF *bsdf = (ccl_private MicrofacetHairBSDF *)sc;

  /* Get local coordinate system:
   * . X major axis.
   * . Y along the fiber tangent.
   * . Z minor axis. */
  const float3 X = float4_to_float3(bsdf->extra->geom);
  const float3 Z = safe_normalize(cross(X, sd->dPdu));
  const float3 Y = safe_normalize(cross(Z, X));

  /* Transform wi/wo from global coordinate system to local. */
  const float3 local_I = make_float3(dot(sd->wi, X), dot(sd->wi, Y), dot(sd->wi, Z));
  const float3 local_O = make_float3(dot(wo, X), dot(wo, Y), dot(wo, Z));

  /* Treat as transparent material if intersection lies outside of the projected radius. */
  const float e2 = 1.0f - sqr(bsdf->aspect_ratio);
  if (fabsf(bsdf->extra->geom.w) > projected_radius(e2, dir_phi(local_I))) {
    *pdf = 0.0f;
    return zero_spectrum();
  }

  /* Evaluate. */
  float3 R;
  if (bsdf->distribution_type == NODE_MICROFACET_HAIR_BECKMANN) {
    R = bsdf_microfacet_hair_eval_r<MicrofacetType::BECKMANN>(sc, local_I, local_O) +
        bsdf_microfacet_hair_eval_tt_trt<MicrofacetType::BECKMANN>(
            kg, sc, local_I, local_O, sd->lcg_state);
  }
  else {
    R = bsdf_microfacet_hair_eval_r<MicrofacetType::GGX>(sc, local_I, local_O) +
        bsdf_microfacet_hair_eval_tt_trt<MicrofacetType::GGX>(
            kg, sc, local_I, local_O, sd->lcg_state);
  }

  /* TODO: better estimation of the pdf */
  *pdf = 1.0f;

  return rgb_to_spectrum(R / cos_theta(local_I));
}

ccl_device int bsdf_microfacet_hair_sample(KernelGlobals kg,
                                           ccl_private const ShaderClosure *sc,
                                           ccl_private ShaderData *sd,
                                           float randu,
                                           float randv,
                                           ccl_private Spectrum *eval,
                                           ccl_private float3 *wo,
                                           ccl_private float *pdf,
                                           ccl_private float2 *sampled_roughness,
                                           ccl_private float *eta)
{
  ccl_private MicrofacetHairBSDF *bsdf = (ccl_private MicrofacetHairBSDF *)sc;

  if (bsdf->distribution_type == NODE_MICROFACET_HAIR_BECKMANN) {
    return bsdf_microfacet_hair_sample<MicrofacetType::BECKMANN>(
        kg, sc, sd, randu, randv, eval, wo, pdf, sampled_roughness, eta);
  }
  return bsdf_microfacet_hair_sample<MicrofacetType::GGX>(
      kg, sc, sd, randu, randv, eval, wo, pdf, sampled_roughness, eta);
}

/* Implements Filter Glossy by capping the effective roughness. */
ccl_device void bsdf_microfacet_hair_blur(ccl_private ShaderClosure *sc, float roughness)
{
  ccl_private MicrofacetHairBSDF *bsdf = (ccl_private MicrofacetHairBSDF *)sc;

  bsdf->roughness = fmaxf(roughness, bsdf->roughness);
}

/* Hair Albedo. */
ccl_device float3 bsdf_microfacet_hair_albedo(ccl_private const ShaderClosure *sc)
{
  ccl_private MicrofacetHairBSDF *bsdf = (ccl_private MicrofacetHairBSDF *)sc;
  return exp(-sqrt(bsdf->sigma) * bsdf_hair_albedo_roughness_scale(bsdf->roughness));
}

CCL_NAMESPACE_END
