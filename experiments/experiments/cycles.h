/*
 * File:   cycles.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
 *
 */

#pragma once

#include "graph.h"

using namespace std;


namespace ExpCycles {

void start_or_continue_build_cycle(int possible_edge_lower_bound, Graph& graph, bool hamiltonian);
void build_cycle(int cur_vertex, int min_possible_edge, Graph& graph, bool hamiltonian);
void find_all_cycles(Graph& graph, bool hamiltonian = false);

} // ExpCycles
