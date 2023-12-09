/*
 * File:   cycles_data.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
 *
 */

#pragma once

#include "constants.h"

#include <unordered_set>
#include <unordered_map>
#include <vector>

using namespace std;


namespace ExpCycles {

struct Data {
  int oddness = 0;
  unordered_set<Mask> all_cycles;
  unordered_set<Mask> all_circuits;
  unordered_set<Mask> all_even_cycles;
  unordered_set<Mask> all_even_v_minus_4_cycles;
  unordered_set<Mask> all_full_cycles;
  unordered_set<Mask> all_dominating_circuits;
  unordered_set<Mask> all_dominating_cycles;
  map<Mask, Mask> ignored_vertices_by_dominating_cycle;
  unordered_set<Mask> all_dominating_vertex_sets;
  unordered_set<Mask> all_nonseparating_cycles;
  unordered_map<Mask, vector<vector<int>>> cycles_as_circuits;
  unordered_map<Mask, vector<Mask>> cycles_as_circuit_masks;
};

} // ExpCycles
