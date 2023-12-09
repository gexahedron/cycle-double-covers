/*
 * File:   flow_parity_pairs_data.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
 *
 */

#pragma once

#include "constants.h"

#include <set>

using namespace std;


namespace ExpFlowParityPairs {

struct Data {
  set<Mask> all_33pp_cycles;
  set<Mask> all_33pp_circuits;
  set<Mask> all_33pp_full_cycles;

  set<Mask> all_333pp_cycles;
  set<Mask> all_333pp_even_cycles;
  set<Mask> all_333pp_full_cycles;
};

} // ExpFlowParityPairs
