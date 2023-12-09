/*
 * File:   o5cdc_data.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
 *
 */

#pragma once

#include "constants.h"

#include <unordered_set>
#include <set>

using namespace std;


namespace Exp5cdc {

struct Data {
  set<set<Mask>> all_5cdc;
  set<set<Mask>> all_5cdc_with_domcyc;
  set<set<Mask>> all_5cdc_with_domcirc;
  set<set<Mask>> all_o5cdc;
  set<set<pair<int, int>>> all_o5cdc_edge_pairs;
  map<set<Mask>, set<set<Mask>>> all_5cdc_by_triples;
  map<Mask, map<int, int>> windings;

  // TODO: these come from o5cdc; how come?
  set<int> dominating_petersens;
  unordered_set<Mask> pet_5cdc_dominating_circuits;
  unordered_set<Mask> all_33pp_dominating_circuits;

  map<pair<Mask, pair<Mask, Mask>>, set<set<Mask>>> u5cdc_from_33pp;
  set<set<Mask>> u5cdc_with_33pp;
  set<set<Mask>> u5cdc_with_domcyc_with_33pp;
  set<set<Mask>> u5cdc_with_domcyc_with_33pp_relaxed;
  set<set<Mask>> u5cdc_with_domcirc_with_33pp;
  set<set<Mask>> u5cdc_with_domcirc_with_33pp_relaxed;

  map<vector<int>, set<Mask>> ignored_vertices_by_nz5;
  map<vector<int>, set<Mask>> ignored_vertices_by_nzmod5; // todo: remove

  map<vector<int>, set<pair<Mask, pair<Mask, Mask>>>> from_nz5_to_33pp;
};

} // Exp5cdc
