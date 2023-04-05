/* SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

#include "BKE_node.h"

#ifdef __cplusplus
extern "C" {
#endif

extern struct bNodeTreeType *ntreeType_Geometry;

void register_node_type_geo_custom_group(bNodeType *ntype);

/* Simulation Input Node API */

/**
 * Pair a simulation input node with an output node.
 * @return True if pairing the node was successful.
 */
bool node_geo_simulation_input_pair_with_output(const bNodeTree *node_tree,
                                                bNode *simulation_input_node,
                                                const bNode *simulation_output_node);

#ifdef __cplusplus
}
#endif
