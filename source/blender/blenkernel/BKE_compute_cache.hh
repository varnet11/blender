/* SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

#include <map>

#include "BLI_bit_vector.hh"
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

  /* TODO: These will need to be a generic value at some point. */
  struct CacheValues {
    Vector<GeometrySet> persistent_cache;
    BitVector<> persistent_caches_filled;

    std::optional<GeometrySet> non_persistent_value;

    void ensure_size(const int size)
    {
      persistent_cache.resize(size);
      persistent_caches_filled.resize(size);
    }
  };

  /* Map from cache data identifier (socket name) to values at stored times. */
  Map<std::string, CacheValues> caches_;

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
    if (std::optional<GeometrySet> value = this->value_at_time(data_name, time)) {
      return value;
    }
    if (std::optional<GeometrySet> value = this->value_before_time(data_name, time)) {
      return value;
    }
    return std::nullopt;
  }

  std::optional<GeometrySet> value_before_time(StringRef data_name, TimePoint time)
  {
    const CacheValues *values = caches_.lookup_ptr(data_name);
    if (!values) {
      return std::nullopt;
    }
    if (values->non_persistent_value) {
      if (last_run_time_->time < time.time) {
        return std::move(values->non_persistent_value);
      }
    }
    /* TODO: Maybe separate retrieval of persistent and temporary cache values?
     * Though that doesn't really provide a benefit right now. */
    if (values->persistent_cache.is_empty()) {
      return std::nullopt;
    }
    int index = std::min<int>(values->persistent_cache.index_range().last(),
                              time.frame - start_time_->frame);
    if (index < 0) {
      return std::nullopt;
    }
    while (!values->persistent_caches_filled[index]) {
      index--;
      if (index < 0) {
        return std::nullopt;
      }
    }
    return values->persistent_cache[index];
  }

  std::optional<GeometrySet> value_at_time(StringRef data_name, TimePoint time)
  {
    const CacheValues *values = caches_.lookup_ptr(data_name);
    if (!values) {
      return std::nullopt;
    }
    if (values->non_persistent_value) {
      if (last_run_time_->frame == time.frame) {
        return std::move(values->non_persistent_value);
      }
    }
    const int index = time.frame - start_time_->frame;
    if (!values->persistent_cache.index_range().contains(index)) {
      return std::nullopt;
    }
    if (!values->persistent_caches_filled[index]) {
      return std::nullopt;
    }
    return values->persistent_cache[index];
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
    CacheValues &values = caches_.lookup_or_add_default_as(data_name);
    const int index = time.frame - start_time_->frame;
    values.ensure_size(index + 1);
    values.persistent_cache[index] = std::move(value);
    values.persistent_caches_filled[index].set();
  }

  void clear()
  {
    caches_.clear();
    start_time_.reset();
    last_run_time_.reset();
  }

  bool is_empty() const
  {
    for (const CacheValues &values : caches_.values()) {
      if (values.non_persistent_value) {
        return false;
      }
      if (!values.persistent_cache.is_empty()) {
        return false;
      }
    }
    return true;
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
