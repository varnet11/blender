/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2001-2023 Blender Foundation. */

/** \file
 * \ingroup edobj
 */

#include "object_intern.h"

#include "BKE_context.h"
#include "BKE_light_linking.h"

#include "ED_screen.h"

#include "WM_api.h"
#include "WM_types.h"

/* -------------------------------------------------------------------- */
/** \name Create New Light Linking Collection Operator
 * \{ */

static int light_linking_receiver_collection_new_exec(bContext *C, wmOperator *op)
{
  Main *bmain = CTX_data_main(C);
  Object *object = CTX_data_active_object(C);

  BKE_light_linking_receiver_collection_new(bmain, object);

  return OPERATOR_FINISHED;
}

void OBJECT_OT_light_linking_receiver_collection_new(wmOperatorType *ot)
{
  /* identifiers */
  ot->name = "New Light Linking Collection";
  ot->description = "Create new light linking collection used by the active emitter";
  ot->idname = "OBJECT_OT_light_linking_receiver_collection_new";

  /* api callbacks */
  ot->exec = light_linking_receiver_collection_new_exec;
  ot->poll = ED_operator_object_active_editable;

  /* flags */
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

/** \} */
