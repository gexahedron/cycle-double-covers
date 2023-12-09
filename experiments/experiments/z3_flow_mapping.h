/*
 * File:   z3_flow_mapping.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
 *
 * Created on June 1, 2019
 *
 */

#pragma once

#include "graph.h"

namespace ExpZ3FlowMapping {

void gen_all_z3_flows(Graph& graph);
void find_all_z3_flow_mappings(Graph& graph, const Graph& petersen_graph);

} // ExpZ3FlowMapping
