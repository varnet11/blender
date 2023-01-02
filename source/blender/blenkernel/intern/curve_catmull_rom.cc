/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 */

#include "BKE_attribute_math.hh"
#include "BKE_curves.hh"

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
      case 1:
        dst.first() = float3(0.0f, 0.0f, 1.0f);
        break;
      case 2:
        dst.first() = float3(1.0f, 0.0f, 0.0f);
        break;
      default:
        dst.first() = src.first();
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
static void interpolate_to_evaluated(const Span<T> src,
                                     const bool cyclic,
                                     const int resolution,
                                     MutableSpan<T> dst)

{
  BLI_assert(dst.size() == calculate_evaluated_num(src.size(), cyclic, resolution));
  interpolate_to_evaluated(
      0,
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
    interpolate_to_evaluated(src.typed<T>(), cyclic, resolution, dst.typed<T>());
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
  interpolate_to_evaluated(
      1,
      positions,
      is_cyclic,
      [resolution](const int segment_i) -> IndexRange {
        return {segment_i * resolution, resolution};
      },
      evaluated_tangents);

  for (const int i : evaluated_tangents.index_range()) {
    evaluated_tangents[i] = math::normalize(evaluated_tangents[i]);
  }
}

/* Calculate analytical normals from the first and the second derivative. */
void calculate_normals(const Span<float3> positions,
                       const bool is_cyclic,
                       const int resolution,
                       const Span<float3> evaluated_tangents,
                       MutableSpan<float3> evaluated_normals)
{
  /* Compute the second derivative (r'') of the control points, evaluated at t = 0. */
  interpolate_to_evaluated(
      2,
      positions,
      is_cyclic,
      [resolution](const int segment_i) -> IndexRange {
        return {segment_i * resolution, 1};
      },
      evaluated_normals);

  /* The normals along the segment is interpolated between the control points. */
  const float step = 1.0f / resolution;
  for (const int i : positions.index_range()) {
    if (i == positions.size() - 1 && !is_cyclic) {
      break;
    }
    const float3 last_normal = evaluated_normals[i * resolution];
    float3 next_normal = (i == positions.size() - 1) ? evaluated_normals[0] :
                                                       evaluated_normals[(i + 1) * resolution];
    if (!is_cyclic && math::dot(last_normal, next_normal) < 0) {
      /* Prevent sudden normal changes. */
      next_normal = -next_normal;
      evaluated_normals[(i + 1) * resolution] = next_normal;
    }
    for (int j = 1; j < resolution; j++) {
      const float t = j * step;
      evaluated_normals[i * resolution + j] = (1.0f - t) * last_normal + t * next_normal;
    }
  }

  MutableSpan<float3> binormals(evaluated_normals);

  /* Compute evaluated_normals of the control points. */
  for (const int i : evaluated_normals.index_range()) {
    /* B = normalize(r' x r'') */
    binormals[i] = math::cross(evaluated_tangents[i], evaluated_normals[i]);
    /* N = B x T */
    evaluated_normals[i] = math::normalize(math::cross(binormals[i], evaluated_tangents[i]));
  }

  if (positions.size() == 1) {
    return;
  }
}

}  // namespace blender::bke::curves::catmull_rom
