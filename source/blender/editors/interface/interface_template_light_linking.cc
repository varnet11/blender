/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2023 Blender Foundation. */

/** \file
 * \ingroup edinterface
 */

#include "UI_interface.h"

#include <cstdio>
#include <memory>

#include "BLT_translation.h"

#include "DNA_collection_types.h"
#include "DNA_object_types.h"

#include "RNA_access.h"
#include "RNA_prototypes.h"

#include "UI_interface.h"
#include "UI_interface.hh"
#include "UI_resources.h"
#include "UI_tree_view.hh"

namespace blender {

namespace {

class ReceiverCollectionViewItem : public ui::BasicTreeViewItem {
 public:
  ReceiverCollectionViewItem(ID &id, const BIFIconID icon)
      : ui::BasicTreeViewItem(id.name + 2, icon), id_(&id)
  {
  }

  void build_row(uiLayout &row) override
  {
    add_label(row);

    if (!is_hovered()) {
      return;
    }

    PointerRNA id_ptr{id_, &RNA_ID, (ID *)id_};
    uiLayoutSetContextPointer(&row, "id", &id_ptr);

    uiItemO(&row, "", ICON_X, "OBJECT_OT_light_linking_unlink_from_receiver_collection");
  }

 private:
  ID *id_{nullptr};
};

class ReceiverCollectionView : public ui::AbstractTreeView {
 public:
  explicit ReceiverCollectionView(Collection &collection) : collection_(&collection)
  {
  }

  void build_tree() override
  {
    LISTBASE_FOREACH (CollectionChild *, collection_child, &collection_->children) {
      Collection *child_collection = collection_child->collection;
      add_tree_item<ReceiverCollectionViewItem>(child_collection->id, ICON_OUTLINER_COLLECTION);
    }

    LISTBASE_FOREACH (CollectionObject *, collection_object, &collection_->gobject) {
      Object *child_object = collection_object->ob;
      add_tree_item<ReceiverCollectionViewItem>(child_object->id, ICON_OBJECT_DATA);
    }
  }

 private:
  Collection *collection_{nullptr};
};

}  // namespace

}  // namespace blender

namespace ui = blender::ui;

void uiTemplateLightLinkingReceiverCollection(struct uiLayout *layout,
                                              struct PointerRNA *ptr,
                                              const char *propname)
{
  if (!ptr->data) {
    return;
  }

  PropertyRNA *prop = RNA_struct_find_property(ptr, propname);
  if (!prop) {
    printf(
        "%s: property not found: %s.%s\n", __func__, RNA_struct_identifier(ptr->type), propname);
    return;
  }

  if (RNA_property_type(prop) != PROP_POINTER) {
    printf("%s: expected pointer property for %s.%s\n",
           __func__,
           RNA_struct_identifier(ptr->type),
           propname);
    return;
  }

  const PointerRNA collection_ptr = RNA_property_pointer_get(ptr, prop);
  if (!collection_ptr.data) {
    return;
  }
  if (collection_ptr.type != &RNA_Collection) {
    printf("%s: expected collection pointer property for %s.%s\n",
           __func__,
           RNA_struct_identifier(ptr->type),
           propname);
    return;
  }

  Collection *collection = static_cast<Collection *>(collection_ptr.data);

  uiBlock *block = uiLayoutGetBlock(layout);

  ui::AbstractTreeView *tree_view = UI_block_add_view(
      *block,
      "Receiver Collection Tree View",
      std::make_unique<blender::ReceiverCollectionView>(*collection));

  ui::TreeViewBuilder builder(*block);
  builder.build_tree_view(*tree_view);
}
