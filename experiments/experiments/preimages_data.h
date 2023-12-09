/*
 * File:   preimages_data.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
 *
 */

#pragma once

#include "constants.h"

#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>

using namespace std;


namespace ExpPreimages {

struct Data {
  unordered_set<Mask> full_cycles_from_petersen;
  map<Mask, set<Mask>> all_cycles_from_petersen;
  map<set<Mask>, vector<int>> petersen_5cdc_to_colouring_idx;
  map<set<Mask>, vector<set<Mask>>> petersen_5cdc_to_pet_5cdc;
  map<multiset<Mask>, vector<int>> petersen_6c4c_to_colouring_idx;
  unordered_set<Mask> all_pet_dominating_circuits;

  set<set<Mask>> petersen_5cdc;
  set<multiset<Mask>> petersen_6c4c;

  map<set<Mask>, set<multiset<Mask>>> petersen_5cdc_6c4c_pairs;
  map<multiset<Mask>, set<set<Mask>>> petersen_6c4c_5cdc_pairs;
};

} // ExpPreimages
