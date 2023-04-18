/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2001-2023 Blender Foundation */

/** \file
 * \ingroup edobj
 */

#include "object_intern.h"

#include "BKE_context.h"
#include "BKE_light_linking.h"

#include "ED_object.h"
#include "ED_screen.h"

#include "WM_api.h"
#include "WM_types.h"

#include "RNA_prototypes.h"

#include "DEG_depsgraph.h"

/* -------------------------------------------------------------------- */
/** \name Create New Light Linking Receiver/Blocker Collection Operators
 * \{ */

template<LightLinkingType link_type>
static int light_linking_collection_new_exec(bContext *C, wmOperator * /*op*/)
{
  Main *bmain = CTX_data_main(C);
  Object *object = ED_object_active_context(C);

  BKE_light_linking_collection_new(bmain, object, link_type);

  return OPERATOR_FINISHED;
}

void OBJECT_OT_light_linking_receiver_collection_new(wmOperatorType *ot)
{
  /* identifiers */
  ot->name = "New Light Linking Collection";
  ot->description = "Create new light linking collection used by the active emitter";
  ot->idname = "OBJECT_OT_light_linking_receiver_collection_new";

  /* api callbacks */
  ot->exec = light_linking_collection_new_exec<LIGHT_LINKING_RECEIVER>;
  ot->poll = ED_operator_object_active_editable;

  /* flags */
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

void OBJECT_OT_light_linking_blocker_collection_new(wmOperatorType *ot)
{
  /* identifiers */
  ot->name = "New Light Linking Collection";
  ot->description = "Create new light linking collection used by the active emitter";
  ot->idname = "OBJECT_OT_light_linking_blocker_collection_new";

  /* api callbacks */
  ot->exec = light_linking_collection_new_exec<LIGHT_LINKING_BLOCKER>;
  ot->poll = ED_operator_object_active_editable;

  /* flags */
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Select Light Linking Receivers/Blockers Operator
 * \{ */

template<LightLinkingType link_type>
static int light_linking_select_exec(bContext *C, wmOperator * /*op*/)
{
  Scene *scene = CTX_data_scene(C);
  ViewLayer *view_layer = CTX_data_view_layer(C);
  Object *emitter = ED_object_active_context(C);

  BKE_light_linking_select_receivers_of_emitter(scene, view_layer, emitter, link_type);

  WM_event_add_notifier(C, NC_SCENE | ND_OB_SELECT, scene);

  return OPERATOR_FINISHED;
}

void OBJECT_OT_light_linking_receivers_select(wmOperatorType *ot)
{
  /* identifiers */
  ot->name = "Select Light Linking Receivers";
  ot->description = "Select all objects which receive light from this emitter";
  ot->idname = "OBJECT_OT_light_linking_receivers_select";

  /* api callbacks */
  ot->exec = light_linking_select_exec<LIGHT_LINKING_RECEIVER>;
  ot->poll = ED_operator_object_active;

  /* flags */
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

void OBJECT_OT_light_linking_blockers_select(wmOperatorType *ot)
{
  /* identifiers */
  ot->name = "Select Light Linking Blockers";
  ot->description = "Select all objects which block light from this emitter";
  ot->idname = "OBJECT_OT_light_linking_blockers_select";

  /* api callbacks */
  ot->exec = light_linking_select_exec<LIGHT_LINKING_BLOCKER>;
  ot->poll = ED_operator_object_active;

  /* flags */
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Link Receivers/Blockers to Emitter Operator
 * \{ */

template<LightLinkingType link_type>
static int light_linking_link_exec(bContext *C, wmOperator * /*op*/)
{
  Main *bmain = CTX_data_main(C);
  Scene *scene = CTX_data_scene(C);
  Object *emitter = ED_object_active_context(C);

  CTX_DATA_BEGIN (C, Object *, receiver, selected_objects) {
    if (receiver == emitter) {
      continue;
    }

    BKE_light_linking_link_receiver_to_emitter(bmain, emitter, receiver, link_type);
  }
  CTX_DATA_END;

  /* It is possible that the receiver collection is also used by the view layer.
   * For this case send a notifier so that the UI is updated for the changes in the collection
   * content. */
  WM_event_add_notifier(C, NC_SCENE | ND_LAYER_CONTENT, scene);

  return OPERATOR_FINISHED;
}

void OBJECT_OT_light_linking_receivers_link(wmOperatorType *ot)
{
  /* identifiers */
  ot->name = "Light Link Receivers to Emitter";
  ot->description = "Light link selected receivers to the active emitter object";
  ot->idname = "OBJECT_OT_light_linking_receivers_link";

  /* api callbacks */
  ot->exec = light_linking_link_exec<LIGHT_LINKING_RECEIVER>;
  ot->poll = ED_operator_object_active_editable;

  /* flags */
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

void OBJECT_OT_light_linking_blockers_link(wmOperatorType *ot)
{
  /* identifiers */
  ot->name = "Light Link Blockers to Emitter";
  ot->description = "Light link selected blockers to the active emitter object";
  ot->idname = "OBJECT_OT_light_linking_blockers_link";

  /* api callbacks */
  ot->exec = light_linking_link_exec<LIGHT_LINKING_BLOCKER>;
  ot->poll = ED_operator_object_active_editable;

  /* flags */
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Unlink from the Light Linking Collection Operator
 * \{ */

static int light_linking_unlink_from_collection_exec(bContext *C, wmOperator *op)
{
  Main *bmain = CTX_data_main(C);

  ID *id = static_cast<ID *>(CTX_data_pointer_get_type(C, "id", &RNA_ID).data);
  Collection *collection = static_cast<Collection *>(
      CTX_data_pointer_get_type(C, "collection", &RNA_Collection).data);

  if (!id || !collection) {
    return OPERATOR_PASS_THROUGH;
  }

  if (!BKE_light_linking_unlink_id_from_collection(bmain, collection, id, op->reports)) {
    return OPERATOR_CANCELLED;
  }

  /* Copy notifiers from the Outliner's "Unlink" operation for objects and collections. */
  WM_event_add_notifier(C, NC_SCENE | ND_LAYER, nullptr);
  WM_event_add_notifier(C, NC_ID | NA_EDITED, nullptr);
  WM_event_add_notifier(C, NC_SPACE | ND_SPACE_OUTLINER, nullptr);

  return OPERATOR_FINISHED;
}

void OBJECT_OT_light_linking_unlink_from_collection(struct wmOperatorType *ot)
{
  /* identifiers */
  ot->name = "Remove From Light Linking Collection";
  ot->description = "Remove this object or collection from the light linking collection";
  ot->idname = "OBJECT_OT_light_linking_unlink_from_collection";

  /* api callbacks */
  ot->exec = light_linking_unlink_from_collection_exec;
  ot->poll = ED_operator_object_active_editable;

  /* flags */
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

/** \} */
