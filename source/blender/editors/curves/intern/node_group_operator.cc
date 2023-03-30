/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup edcurves
 */

#include "ED_curves.h"
#include "ED_object.h"
#include "ED_screen.h"
#include "ED_select_utils.h"
#include "ED_view3d.h"

#include "WM_api.h"

#include "BKE_attribute_math.hh"
#include "BKE_context.h"
#include "BKE_curves.hh"
#include "BKE_geometry_set.hh"
#include "BKE_layer.h"
#include "BKE_lib_id.h"
#include "BKE_node_runtime.hh"
#include "BKE_object.h"
#include "BKE_report.h"

#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "DEG_depsgraph.h"
#include "DEG_depsgraph_query.h"

#include "RNA_access.h"
#include "RNA_define.h"
#include "RNA_enum_types.h"
#include "RNA_prototypes.h"

#include "UI_interface.h"
#include "UI_resources.h"

namespace blender::ed::curves {

static int run_node_group_exec(bContext *C, wmOperator *op)
{
  Scene *scene = CTX_data_scene(C);
  ViewLayer *view_layer = CTX_data_view_layer(C);
  Object *object = CTX_data_active_object(C);
  if (!object) {
    return OPERATOR_CANCELLED;
  }
  eObjectMode mode = eObjectMode(object->mode);

  const uint32_t session_uuid = RNA_int_get(op->ptr, "session_uuid");
  const ID *id = BKE_libblock_find_session_uuid(CTX_data_main(C), ID_NT, session_uuid);
  if (!id) {
    return OPERATOR_CANCELLED;
  }
  const bNodeTree &node_tree = reinterpret_cast<const bNodeTree &>(*id);
  const nodes::GeometryNodesLazyFunctionGraphInfo *lf_graph_info =
      nodes::ensure_geometry_nodes_lazy_function_graph(tree);
  if (lf_graph_info == nullptr) {
    BKE_report(op->reports, RPT_ERROR, "Cannot evaluate node group");
    return;
  }

  uint objects_len = 0;
  Object **objects = BKE_view_layer_array_from_objects_in_mode_unique_data(
      scene, view_layer, CTX_wm_view3d(C), &objects_len, mode);

  return OPERATOR_FINISHED;
}

void CURVES_OT_node_group(wmOperatorType *ot)
{
  ot->name = "Run Node Group";
  ot->idname = __func__;
  ot->description = "Dummy";  // TODO: Retrieve from node group.

  ot->exec = node_group::run_node_group_exec;
  ot->poll = editable_curves_poll;

  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;

  PropertyRNA *prop = RNA_def_int(ot->srna,
                                  "session_uuid",
                                  0,
                                  INT32_MIN,
                                  INT32_MAX,
                                  "Session UUID",
                                  "Session UUID of the node group",
                                  INT32_MIN,
                                  INT32_MAX);
  RNA_def_property_flag(prop, PropertyFlag(PROP_SKIP_SAVE | PROP_HIDDEN));
}

}  // namespace blender::ed::curves
