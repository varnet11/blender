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

#include "BKE_context.h"
#include "BKE_light_linking.h"

#include "RNA_access.h"
#include "RNA_prototypes.h"

#include "UI_interface.h"
#include "UI_interface.hh"
#include "UI_resources.h"
#include "UI_tree_view.hh"

#include "WM_api.h"

#include "ED_undo.h"

namespace blender::ui::light_linking {

namespace {

class CollectionDropTarget : public AbstractViewItemDropTarget {
 public:
  CollectionDropTarget(AbstractView &view, Collection &collection)
      : AbstractViewItemDropTarget(view), collection_(collection)
  {
  }

  bool can_drop(const wmDrag &drag, const char **r_disabled_hint) const override
  {
    if (drag.type != WM_DRAG_ID) {
      return false;
    }

    const wmDragID *drag_id = static_cast<wmDragID *>(drag.ids.first);
    if (!drag_id) {
      return false;
    }

    /* The dragged IDs are guaranteed to be the same type, so only check the type of the first one.
     */
    const ID_Type id_type = GS(drag_id->id->name);
    if (!ELEM(id_type, ID_OB, ID_GR)) {
      *r_disabled_hint = "Can only add objects and collections to the light linking collection";
      return false;
    }

    return true;
  }

  std::string drop_tooltip(const wmDrag & /*drag*/) const override
  {
    return TIP_("Add to light linking collection");
  }

  bool on_drop(struct bContext *C, const wmDrag &drag) const override
  {
    Main *bmain = CTX_data_main(C);
    Scene *scene = CTX_data_scene(C);

    LISTBASE_FOREACH (wmDragID *, drag_id, &drag.ids) {
      BKE_light_linking_receiver_to_collection(bmain, &collection_, drag_id->id);
    }

    /* It is possible that the light linking collection is also used by the view layer.
     * For this case send a notifier so that the UI is updated for the changes in the collection
     * content. */
    WM_event_add_notifier(C, NC_SCENE | ND_LAYER_CONTENT, scene);

    ED_undo_push(C, "Add to light linking collection");

    return true;
  }

 private:
  Collection &collection_;
};

class CollectionViewItem : public BasicTreeViewItem {
 public:
  CollectionViewItem(ID &id, const BIFIconID icon) : BasicTreeViewItem(id.name + 2, icon), id_(&id)
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

    uiItemO(&row, "", ICON_X, "OBJECT_OT_light_linking_unlink_from_collection");
  }

 private:
  ID *id_{nullptr};
};

class CollectionView : public AbstractTreeView {
 public:
  explicit CollectionView(Collection &collection) : collection_(collection)
  {
  }

  void build_tree() override
  {
    LISTBASE_FOREACH (CollectionChild *, collection_child, &collection_.children) {
      Collection *child_collection = collection_child->collection;
      add_tree_item<CollectionViewItem>(child_collection->id, ICON_OUTLINER_COLLECTION);
    }

    LISTBASE_FOREACH (CollectionObject *, collection_object, &collection_.gobject) {
      Object *child_object = collection_object->ob;
      add_tree_item<CollectionViewItem>(child_object->id, ICON_OBJECT_DATA);
    }
  }

  std::unique_ptr<AbstractViewDropTarget> create_drop_target() override
  {
    return std::make_unique<CollectionDropTarget>(*this, collection_);
  }

 private:
  Collection &collection_;
};

}  // namespace

}  // namespace blender::ui::light_linking

namespace ui = blender::ui;

void uiTemplateLightLinkingCollection(struct uiLayout *layout,
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
      "Light Linking Collection Tree View",
      std::make_unique<blender::ui::light_linking::CollectionView>(*collection));
  tree_view->set_min_rows(3);

  ui::TreeViewBuilder::build_tree_view(*tree_view, *layout);
}
