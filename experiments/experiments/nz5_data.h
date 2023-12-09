/*
 * File:   nz5_data.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
 *
 */

#pragma once

#include <vector>

using namespace std;


namespace ExpNZ5 {

struct Data {
  vector<vector<int>> all_nz5_flows;
  vector<vector<int>> all_nz_mod5_flows;
  vector<vector<int>> all_fancy_4colourings;
  vector<vector<vector<int>>> all_nz_complex4_flows;
};

} // ExpNZ5
