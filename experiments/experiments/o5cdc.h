/*
 * File:   o5cdc.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
 *
 */

#pragma once

#include "graph.h"

using namespace std;


namespace Exp5cdc {
bool start_or_continue_build_5cdc(int possible_edge_lower_bound, Graph& graph, int cur_cycle_layer, bool only_find);
void find_all_o5cdc(Graph& graph, bool only_find = false);
} // Exp5cdc
