/*
 * File:   o6c4c_data.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
 *
 */

#pragma once

#include "constants.h"

#include <set>
#include <map>

using namespace std;


namespace Exp6c4c {

struct Data {
  set<multiset<Mask>> all_6c4c;
  set<multiset<Mask>> all_o6c4c;
  bool has_nz_mod5_from_o6c4c = false;
  set<string> o6c4c_profiles;
  map<string, set<int>> u6c4c_profiles_to_parities;
  map<string, set<int>> u6c4c_profiles_to_orientations;
  //map<string, size_t> o6c4c_profiles_min_oriented_vertices;
  set<int> possible_oriented_vertices;

  set<pair<int, int>> s0s1_pet;
  set<pair<int, int>> s0s1_par0;

  // 6c4c vs 5cdc
  map<multiset<Mask>, set<set<Mask>>> u244_6c4c_5cdc_pairs;
  map<set<Mask>, set<multiset<Mask>>> u244_5cdc_6c4c_pairs;

  map<vector<int>, set<Mask>> oriented_vertices_by_nz5;
  map<vector<int>, set<Mask>> oriented_vertices_by_nzmod5; // TODO: remove
};

} // Exp6c4c
