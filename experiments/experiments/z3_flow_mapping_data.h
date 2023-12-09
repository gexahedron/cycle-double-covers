/*
 * File:   z3_flow_mapping_data.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
 *
 * Created on June 1, 2019
 *
 */

#pragma once

#include <vector>
#include <set>

using namespace std;


namespace ExpZ3FlowMapping {

struct Data {
  set<vector<int>> all_z3_flows;
  bool has_z3_flow_mapping = false;
  set<vector<pair<int, int>>> z3_edge_mappings;
};

} // ExpZ3FlowMapping
