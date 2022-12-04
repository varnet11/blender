/* SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

#include <map>

#include "BLI_compute_context.hh"
#include "BLI_map.hh"

#include "BKE_geometry_set.hh"

namespace blender::bke::sim {

/* TODO: Take into account subframes. But do we make caches for subframes? Probably not. */
struct TimePoint {
  int frame;
  float time;

  uint64_t hash()
  {
    return get_default_hash(this->frame);
  }

  friend bool operator==(const TimePoint &a, const TimePoint &b)
  {
    return a.frame == b.frame;
  }
};

/* TODO: Clear cache when editing nodes? Only sometimes, when persistent caching is turned off. */
class SimulationCache {

  struct CacheValues {
    /* TODO: This will need to be a generic value at some point. */
    /* Map from simulation time index (see #SimulationCache::times) to value. */
    Vector<GeometrySet> persistent_cache;

    std::optional<GeometrySet> non_persistent_value;

    GeometrySet lookup_index(const int index) const
    {
      if (!persistent_cache.index_range().contains(index)) {
        return persistent_cache[index];
      }
      return {};
    }
  };

  /* Map from cache data identifier (socket name) to values at stored times. */
  Map<std::string, CacheValues> caches_;

  /* Ordered list of cached simulation frames. */
  Vector<TimePoint> persistent_times_;

  std::optional<TimePoint> start_time_;

  std::optional<TimePoint> last_run_time_;

 public:
  std::optional<TimePoint> start_time() const
  {
    return start_time_;
  }
  std::optional<TimePoint> last_run_time() const
  {
    return last_run_time_;
  }

  std::optional<GeometrySet> value_at_or_before_time(StringRef data_name, TimePoint time)
  {
    const CacheValues *values = caches_.lookup_ptr(data_name);
    if (!values) {
      return std::nullopt;
    }
    if (last_run_time_->time < time.time) {
      if (values->non_persistent_value) {
        return std::move(values->non_persistent_value);
      }
    }
    const int index = this->index_at_or_before_time(time);
    if (!values->persistent_cache.index_range().contains(index)) {
      return std::nullopt;
    }
    return values->persistent_cache[index];
  }

  std::optional<GeometrySet> value_before_time(StringRef data_name, TimePoint time)
  {
    const CacheValues *values = caches_.lookup_ptr(data_name);
    if (!values) {
      return std::nullopt;
    }
    if (last_run_time_->time < time.time) {
      if (values->non_persistent_value) {
        return std::move(values->non_persistent_value);
      }
    }
    const int index = this->index_before_time(time);
    if (!values->persistent_cache.index_range().contains(index)) {
      return std::nullopt;
    }
    return values->persistent_cache[index];
  }

  std::optional<GeometrySet> value_at_time(StringRef data_name, TimePoint time)
  {
    const CacheValues *values = caches_.lookup_ptr(data_name);
    if (!values) {
      return std::nullopt;
    }
    const std::optional<int> index = this->index_at_time(time);
    if (!index) {
      return std::nullopt;
    }
    return values->persistent_cache[*index];
  }

  void store_temporary(const StringRef data_name, const TimePoint time, GeometrySet value)
  {
    last_run_time_.emplace(time);
    if (!start_time_) {
      start_time_.emplace(time);
    }
    CacheValues &values = caches_.lookup_or_add_default_as(data_name);
    values.non_persistent_value.emplace(std::move(value));
  }

  void store_persistent(const StringRef data_name, const TimePoint time, GeometrySet value)
  {
    last_run_time_.emplace(time);
    if (!start_time_) {
      start_time_.emplace(time);
    }
    const int index = this->index_before_time(time);
    persistent_times_.resize(index);
    persistent_times_[index] = time;
    CacheValues &values = caches_.lookup_or_add_default_as(data_name);
    values.persistent_cache.resize(index);
    values.persistent_cache[index] = std::move(value);
  }

  void clear()
  {
    persistent_times_.clear();
    caches_.clear();
  }

  bool is_empty() const
  {
    return persistent_times_.is_empty();
  }

 private:
  int index_at_or_before_time(const TimePoint time) const
  {
    if (persistent_times_.is_empty()) {
      return 0;
    }
    int insert_index = 0;
    for (const int i : persistent_times_.index_range()) {
      if (persistent_times_[i].frame <= time.frame) {
        break;
      }
      insert_index++;
    }
    return insert_index;
  }
  int index_before_time(const TimePoint time) const
  {
    if (persistent_times_.is_empty()) {
      return 0;
    }
    int insert_index = 0;
    for (const int i : persistent_times_.index_range()) {
      if (persistent_times_[i].frame < time.frame) {
        break;
      }
      insert_index++;
    }
    return insert_index;
  }
  std::optional<int> index_at_time(const TimePoint time) const
  {
    for (const int i : persistent_times_.index_range()) {
      if (persistent_times_[i].frame == time.frame) {
        return i;
      }
    }
    return std::nullopt;
  }
};

struct ComputeCaches {
 private:
  mutable std::mutex mutex;
  Map<ComputeContextHash, SimulationCache> cache_per_context;

 public:
  ComputeCaches() = default;
  ComputeCaches(const ComputeCaches &other)
  {
    cache_per_context = other.cache_per_context;
  }

  SimulationCache *lookup_context(const ComputeContextHash &context_hash)
  {
    std::scoped_lock lock{mutex};
    return cache_per_context.lookup_ptr(context_hash);
  }
  const SimulationCache *lookup_context(const ComputeContextHash &context_hash) const
  {
    std::scoped_lock lock{mutex};
    return cache_per_context.lookup_ptr(context_hash);
  }

  SimulationCache &ensure_for_context(const ComputeContextHash &context_hash)
  {
    std::scoped_lock lock{mutex};
    return cache_per_context.lookup_or_add_default(context_hash);
  }

  bool is_empty() const
  {
    return cache_per_context.is_empty();
  }
};

}  // namespace blender::bke::sim
