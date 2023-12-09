/*
 * File:   tree_cycle_matching_data.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
 *
 */

#pragma once

#include <unordered_set>
#include <vector>

using namespace std;


struct TreeCycleMatching {
  unordered_set<int> tree;
  unordered_set<int> cycle;
  unordered_set<int> matching;
};

namespace ExpTreeCycleMatching {

struct Data {
  vector<TreeCycleMatching> tree_cycle_matchings;
};

} // ExpTreeCycleMatching
