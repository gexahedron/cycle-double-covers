/*
 * File:   flows.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
 *
 */

#pragma once

#include "graph.h"
#include "constants.h"

#include <vector>

using namespace std;


namespace UtilFlows {
Mask build_full_cycle_from_nz5_flow(Graph& graph, const vector<int>& f);
bool check_3flow(Graph& graph, const Mask c);
} // UtilFlows
