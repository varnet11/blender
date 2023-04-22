/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "BKE_simulation_state.hh"
#include "BKE_simulation_state_serialize.hh"

#include "BLI_fileops.hh"
#include "BLI_path_util.h"

namespace blender::bke::sim {

void ModifierSimulationCache::try_discover_bake(const StringRefNull meta_dir,
                                                const StringRefNull bdata_dir)
{
  if (failed_finding_bake_) {
    return;
  }
  if (!BLI_is_dir(meta_dir.c_str()) || !BLI_is_dir(bdata_dir.c_str())) {
    failed_finding_bake_ = true;
    return;
  }

  direntry *dir_entries = nullptr;
  const int dir_entries_num = BLI_filelist_dir_contents(meta_dir.c_str(), &dir_entries);
  BLI_SCOPED_DEFER([&]() { BLI_filelist_free(dir_entries, dir_entries_num); });

  if (dir_entries_num == 0) {
    failed_finding_bake_ = true;
    return;
  }

  this->reset();

  for (const int i : IndexRange(dir_entries_num)) {
    const direntry &dir_entry = dir_entries[i];
    const StringRefNull dir_entry_path = dir_entry.path;
    if (!dir_entry_path.endswith(".json")) {
      continue;
    }
    char modified_file_name[FILENAME_MAX];
    BLI_strncpy(modified_file_name, dir_entry.relname, sizeof(modified_file_name));
    BLI_str_replace_char(modified_file_name, '_', '.');

    const SubFrame frame = std::stof(modified_file_name);

    auto new_state_at_frame = std::make_unique<ModifierSimulationStateAtFrame>();
    new_state_at_frame->frame = frame;
    new_state_at_frame->state.bdata_dir_ = bdata_dir;
    new_state_at_frame->state.meta_path_ = dir_entry.path;
    new_state_at_frame->state.owner_ = this;
    states_at_frames_.append(std::move(new_state_at_frame));
  }

  bdata_sharing_ = std::make_unique<BDataSharing>();

  cache_state_ = CacheState::Baked;
}

void ModifierSimulationState::ensure_bake_loaded() const
{
  std::scoped_lock lock{mutex_};
  if (bake_loaded_) {
    return;
  }
  if (!meta_path_ || !bdata_dir_) {
    return;
  }

  const std::shared_ptr<io::serialize::Value> io_root_value = io::serialize::read_json_file(
      *meta_path_);
  if (!io_root_value) {
    return;
  }
  const DictionaryValue *io_root = io_root_value->as_dictionary_value();
  if (!io_root) {
    return;
  }

  const DiskBDataReader bdata_reader{*bdata_dir_};
  deserialize_modifier_simulation_state(*io_root,
                                        bdata_reader,
                                        *owner_->bdata_sharing_,
                                        const_cast<ModifierSimulationState &>(*this));
  bake_loaded_ = true;
}

void ModifierSimulationCache::reset()
{
  states_at_frames_.clear();
  bdata_sharing_.reset();
  cache_state_ = CacheState::Valid;
}

}  // namespace blender::bke::sim
