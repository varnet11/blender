/* SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

#include "BKE_geometry_set.hh"

#include "BLI_map.hh"

namespace blender::bke::sim {

class SimulationStateItem {
 public:
  virtual ~SimulationStateItem() = default;
};

class GeometrySimulationStateItem : public SimulationStateItem {
 private:
  GeometrySet geometry_;

 public:
  GeometrySimulationStateItem(GeometrySet geometry) : geometry_(std::move(geometry)) {}

  const GeometrySet &geometry() const
  {
    return geometry_;
  }
};

class SimulationZoneState {
 public:
  Vector<std::unique_ptr<SimulationStateItem>> items;
};

struct SimulationZoneID {
  Vector<int> node_ids;

  uint64_t hash() const
  {
    return get_default_hash(this->node_ids);
  }

  friend bool operator==(const SimulationZoneID &a, const SimulationZoneID &b)
  {
    return a.node_ids == b.node_ids;
  }
};

class ModifierSimulationState {
 private:
  mutable std::mutex mutex_;
  Map<SimulationZoneID, std::unique_ptr<SimulationZoneState>> zone_states_;

 public:
  const SimulationZoneState *get_zone_state(const SimulationZoneID &zone_id) const
  {
    std::lock_guard lock{mutex_};
    if (auto *ptr = zone_states_.lookup_ptr(zone_id)) {
      return ptr->get();
    }
    return nullptr;
  }

  SimulationZoneState &get_zone_state_for_write(const SimulationZoneID &zone_id)
  {
    std::lock_guard lock{mutex_};
    return *zone_states_.lookup_or_add_cb(
        zone_id, []() { return std::make_unique<SimulationZoneState>(); });
  }
};

class ModifierSimulationCache {
 private:
  Map<float, std::unique_ptr<ModifierSimulationState>> states_by_time_;
  bool invalid_ = false;

 public:
  bool has_state_at_time(const float time) const
  {
    return states_by_time_.contains(time);
  }

  const ModifierSimulationState *get_state_at_time(const float time) const
  {
    if (auto *ptr = states_by_time_.lookup_ptr(time)) {
      return ptr->get();
    }
    return nullptr;
  }

  ModifierSimulationState &get_state_for_write(const float time)
  {
    return *states_by_time_.lookup_or_add_cb(
        time, []() { return std::make_unique<ModifierSimulationState>(); });
  }

  std::pair<float, const ModifierSimulationState *> try_get_last_state_before(
      const float time) const
  {
    float last_time = -FLT_MAX;
    const ModifierSimulationState *last_state = nullptr;
    for (const auto &item : states_by_time_.items()) {
      if (item.key < time && item.key > last_time) {
        last_time = item.key;
        last_state = item.value.get();
      }
    }
    return {last_time, last_state};
  }

  void invalidate()
  {
    invalid_ = true;
  }

  bool is_invalid() const
  {
    return invalid_;
  }

  void reset()
  {
    states_by_time_.clear();
    invalid_ = false;
  }
};

}  // namespace blender::bke::sim
