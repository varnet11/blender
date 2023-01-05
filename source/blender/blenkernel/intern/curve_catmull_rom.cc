/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 */

#include "BKE_attribute_math.hh"
#include "BKE_curves.hh"
#include "BLI_math_rotation_legacy.hh"
#include "BLI_math_vector.hh"

namespace blender::bke::curves::catmull_rom {

int calculate_evaluated_num(const int points_num, const bool cyclic, const int resolution)
{
  const int eval_num = resolution * segments_num(points_num, cyclic);
  if (cyclic) {
    /* Make sure there is a single evaluated point for the single-point curve case. */
    return std::max(eval_num, 1);
  }
  /* If the curve isn't cyclic, one last point is added to the final point. */
  return eval_num + 1;
}

/* Adapted from Cycles #catmull_rom_basis_eval function. */
void calculate_basis(const float parameter, float4 &r_weights)
{
  const float t = parameter;
  const float s = 1.0f - parameter;
  r_weights[0] = -t * s * s;
  r_weights[1] = 2.0f + t * t * (3.0f * t - 5.0f);
  r_weights[2] = 2.0f + s * s * (3.0f * s - 5.0f);
  r_weights[3] = -s * t * t;
}

/* Adapted from Cycles #catmull_rom_basis_derivative function. */
void calculate_basis_derivative(const float parameter, float4 &r_weights)
{
  const float t = parameter;
  const float s = 1.0f - parameter;
  r_weights[0] = s * (2.0f * t - s);
  r_weights[1] = t * (9.0f * t - 10.0f);
  r_weights[2] = -s * (9.0f * s - 10.0f);
  r_weights[3] = t * (t - 2.0f * s);
}

/* Adapted from Cycles #catmull_rom_basis_derivative2 function. */
void calculate_basis_derivative2(const float parameter, float4 &r_weights)
{
  const float t = parameter;
  const float s = 1.0f - parameter;
  r_weights[0] = -6.0f * t + 4.0f;
  r_weights[1] = 18.0f * t - 10.0f;
  r_weights[2] = 18.0f * s - 10.0f;
  r_weights[3] = 6.0f * t - 2.0f;
}

template<typename T>
static void evaluate_segment(
    const int order, const T &a, const T &b, const T &c, const T &d, MutableSpan<T> dst)
{
  const float step = 1.0f / dst.size();
  dst.first() = b;
  IndexRange index_range = (order == 0) ? dst.index_range().drop_front(1) : dst.index_range();
  for (const int i : index_range) {
    dst[i] = interpolate(order, a, b, c, d, i * step);
  }
}

/**
 * \param range_fn: Returns an index range describing where in the #dst span each segment should be
 * evaluated to, and how many points to add to it. This is used to avoid the need to allocate an
 * actual offsets array in typical evaluation use cases where the resolution is per-curve.
 * \param order: The order of the derivative. order == 0 is the default.
 */
template<typename T, typename RangeForSegmentFn>
static void interpolate_to_evaluated(const int order,
                                     const Span<T> src,
                                     const bool cyclic,
                                     const RangeForSegmentFn &range_fn,
                                     MutableSpan<T> dst)

{
  /* - First deal with one and two point curves need special attention.
   * - Then evaluate the first and last segment(s) whose control points need to wrap around
   *   to the other side of the source array.
   * - Finally evaluate all of the segments in the middle in parallel. */

  if (src.size() == 1) {
    switch (order) {
      case 0:
        dst.first() = src.first();
        break;
      default:
        dst.first() = T(0);
        break;
    }
    return;
  }

  const IndexRange first = range_fn(0);

  if (src.size() == 2) {
    evaluate_segment(order, src.first(), src.first(), src.last(), src.last(), dst.slice(first));
    if (cyclic) {
      const IndexRange last = range_fn(1);
      evaluate_segment(order, src.last(), src.last(), src.first(), src.first(), dst.slice(last));
    }
    else {
      dst.last() = interpolate(order, src.first(), src.first(), src.last(), src.last(), 1);
    }
    return;
  }

  const IndexRange second_to_last = range_fn(src.index_range().last(1));
  const IndexRange last = range_fn(src.index_range().last());
  if (cyclic) {
    evaluate_segment(order, src.last(), src[0], src[1], src[2], dst.slice(first));
    evaluate_segment(
        order, src.last(2), src.last(1), src.last(), src.first(), dst.slice(second_to_last));
    evaluate_segment(order, src.last(1), src.last(), src[0], src[1], dst.slice(last));
  }
  else {
    evaluate_segment(order, src[0], src[0], src[1], src[2], dst.slice(first));
    evaluate_segment(
        order, src.last(2), src.last(1), src.last(), src.last(), dst.slice(second_to_last));
    /* For non-cyclic curves, the last segment should always just have a single point. We could
     * assert that the size of the provided range is 1 here, but that would require specializing
     * the #range_fn implementation for the last point, which may have a performance cost. */
    dst.last() = interpolate(order, src.last(2), src.last(1), src.last(), src.last(), 1);
  }

  /* Evaluate every segment that isn't the first or last. */
  const IndexRange inner_range = src.index_range().drop_back(2).drop_front(1);
  threading::parallel_for(inner_range, 512, [&](IndexRange range) {
    for (const int i : range) {
      const IndexRange segment = range_fn(i);
      evaluate_segment(order, src[i - 1], src[i], src[i + 1], src[i + 2], dst.slice(segment));
    }
  });
}

template<typename T>
static void interpolate_to_evaluated(const int order,
                                     const Span<T> src,
                                     const bool cyclic,
                                     const int resolution,
                                     MutableSpan<T> dst)

{
  BLI_assert(dst.size() == calculate_evaluated_num(src.size(), cyclic, resolution));
  interpolate_to_evaluated(
      order,
      src,
      cyclic,
      [resolution](const int segment_i) -> IndexRange {
        return {segment_i * resolution, resolution};
      },
      dst);
}

template<typename T>
static void interpolate_to_evaluated(const Span<T> src,
                                     const bool cyclic,
                                     const Span<int> evaluated_offsets,
                                     MutableSpan<T> dst)

{
  interpolate_to_evaluated(
      0,
      src,
      cyclic,
      [evaluated_offsets](const int segment_i) -> IndexRange {
        return bke::offsets_to_range(evaluated_offsets, segment_i);
      },
      dst);
}

void interpolate_to_evaluated(const GSpan src,
                              const bool cyclic,
                              const int resolution,
                              GMutableSpan dst)
{
  attribute_math::convert_to_static_type(src.type(), [&](auto dummy) {
    using T = decltype(dummy);
    interpolate_to_evaluated(0, src.typed<T>(), cyclic, resolution, dst.typed<T>());
  });
}

void interpolate_to_evaluated(const GSpan src,
                              const bool cyclic,
                              const Span<int> evaluated_offsets,
                              GMutableSpan dst)
{
  attribute_math::convert_to_static_type(src.type(), [&](auto dummy) {
    using T = decltype(dummy);
    interpolate_to_evaluated(src.typed<T>(), cyclic, evaluated_offsets, dst.typed<T>());
  });
}

/* Calculate analytical tangents from the first derivative. */
void calculate_tangents(const Span<float3> positions,
                        const bool is_cyclic,
                        const int resolution,
                        MutableSpan<float3> evaluated_tangents)
{
  if (evaluated_tangents.size() == 1) {
    evaluated_tangents[0] = float3(0.0f, 0.0f, 1.0f);
    return;
  }

  interpolate_to_evaluated(1, positions, is_cyclic, resolution, evaluated_tangents);

  for (const int i : evaluated_tangents.index_range()) {
    evaluated_tangents[i] = math::normalize(evaluated_tangents[i]);
  }
}

/* Use a combination of Newton iterations and bisection to find the root in the monotonic interval
 * [x_begin, x_end], given precision. Idea taken from the paper High-Performance Polynomial Root
 * Finding for Graphics by Cem Yuksel, HPG 2022. */
/* TODO: move this function elsewhere as a utility function? Also test corner cases. */
static float find_root_newton_bisection(float x_begin,
                                        float x_end,
                                        std::function<float(float)> eval,
                                        std::function<float(float)> eval_derivative,
                                        const float precision)
{
  float x_mid = 0.5f * (x_begin + x_end);
  float y_begin = eval(x_begin);

  BLI_assert(signf(y_begin) != signf(eval(x_end)));

  while (x_mid - x_begin > precision) {
    const float y_mid = eval(x_mid);
    if (signf(y_begin) == signf(y_mid)) {
      x_begin = x_mid;
      y_begin = y_mid;
    }
    else {
      x_end = x_mid;
    }
    const float x_newton = x_mid - y_mid / eval_derivative(x_mid);
    x_mid = (x_newton > x_begin && x_newton < x_end) ? x_newton : 0.5f * (x_begin + x_end);
  }
  return x_mid;
}

/* Calculate normals by minimizing the potential energy due to twist and bending. Global
 * estimation involves integration and is thus too costly. Instead, we start with the first
 * curvature vector and propogate it along the curve. */
void calculate_normals(const Span<float3> positions,
                       const bool is_cyclic,
                       const int resolution,
                       MutableSpan<float3> evaluated_normals)
{
  if (evaluated_normals.size() == 1) {
    evaluated_normals[0] = float3(1.0f, 0.0f, 0.0f);
    return;
  }

  /* TODO: check if derivatives == 0.*/
  Vector<float3> first_derivatives_(evaluated_normals.size());
  MutableSpan<float3> first_derivatives = first_derivatives_;
  interpolate_to_evaluated(1, positions, is_cyclic, resolution, first_derivatives);

  Vector<float3> second_derivatives_(evaluated_normals.size());
  MutableSpan<float3> second_derivatives = second_derivatives_;
  interpolate_to_evaluated(2, positions, is_cyclic, resolution, second_derivatives);

  /* TODO: below are hair-specific values, assuming elliptical cross-section. Maybe make them more
   * general by specifying user-defined weight? */
  /* Values taken from Table 9.5 in book Chemical and Physical Behavior of Human Hair by Clarence
   * R. Robbins, 5th edition. Unit: GPa. */
  const float youngs_modulus = 3.89f;
  const float shear_modulus = 0.89f;
  /* Assuming major axis a = 1. */
  const float aspect_ratio = 0.5f; /* Minimal realistic value. */
  const float aspect_ratio_squared = aspect_ratio * aspect_ratio;
  const float torsion_constant = M_PI * aspect_ratio_squared * aspect_ratio /
                                 (1.0f + aspect_ratio_squared);
  const float moment_of_inertia_coefficient = M_PI * aspect_ratio * 0.25f *
                                              (1.0f - aspect_ratio_squared);

  const float dt = 1.0f / resolution;
  /* Compute evaluated_normals. */
  for (const int i : evaluated_normals.index_range()) {
    /* B = r' x r'' */
    const float3 binormal = math::cross(first_derivatives[i], second_derivatives[i]);
    /* N = B x T */
    /* TODO: Catmull-Rom is not C2 continuous. Some smoothening maybe? */
    float3 curvature_vector = math::normalize(math::cross(binormal, first_derivatives[i]));

    if (i == 0) {
      /* TODO: Does it make sense to propogate from where the curvature is the largest? */
      evaluated_normals[i] = curvature_vector;
      continue;
    }

    const float first_derivative_norm = math::length(first_derivatives[i]);
    /* Arc length depends on the unit, assuming meter. */
    const float arc_length = first_derivative_norm * dt;
    const float curvature = math::length(binormal) / pow3f(first_derivative_norm);

    const float weight_twist = 0.5f * shear_modulus * torsion_constant / arc_length;
    const float weight_bending = 0.5f * arc_length * curvature * curvature * youngs_modulus *
                                 moment_of_inertia_coefficient;

    const float3 last_tangent = math::normalize(first_derivatives[i - 1]);
    const float3 current_tangent = first_derivatives[i] / first_derivative_norm;
    const float angle = angle_normalized_v3v3(last_tangent, current_tangent);
    const float3 last_normal = evaluated_normals[i - 1];
    const float3 minimal_twist_normal =
        angle > 0 ?
            math::rotate_direction_around_axis(
                last_normal, math::normalize(math::cross(last_tangent, current_tangent)), angle) :
            last_normal;

    /* Angle between the curvature vector and the minimal twist normal. */
    float theta_t = angle_normalized_v3v3(curvature_vector, minimal_twist_normal);
    if (theta_t > M_PI_2) {
      curvature_vector = -curvature_vector;
      theta_t = M_PI - theta_t;
    }

    /* Minimizing the total potential energy U(θ) = wt * (θ - θt)^2 + wb * sin^2(θ) by solving the
     * equation: 2 * wt * (θ - θt)  + wb * sin(2θ) = 0, with θ being the angle between the computed
     * normal and the curvature vector. */
    const float theta_begin = 0.0f;
    const float theta_end = weight_twist + weight_bending * cosf(2.0f * theta_t) < 0 ?
                                0.5f * acosf(-weight_twist / weight_bending) :
                                theta_t;

    const float theta = find_root_newton_bisection(
        theta_begin,
        theta_end,
        [weight_bending, weight_twist, theta_t](float t) -> float {
          return 2.0f * weight_twist * (t - theta_t) + weight_bending * sinf(2.0f * t);
        },
        [weight_bending, weight_twist](float t) -> float {
          return 2.0f * (weight_twist + weight_bending * cosf(2.0f * t));
        },
        1e-5f);

    const float3 axis = math::dot(current_tangent,
                                  math::cross(curvature_vector, minimal_twist_normal)) > 0 ?
                            current_tangent :
                            -current_tangent;
    evaluated_normals[i] = math::rotate_direction_around_axis(curvature_vector, axis, theta);
  }
}

}  // namespace blender::bke::curves::catmull_rom
