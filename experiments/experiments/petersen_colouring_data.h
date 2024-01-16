/*
 * File:   petersen_colouring_data.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 */

#pragma once

#include <unordered_map>
#include <vector>
#include <string>

using namespace std;


namespace ExpPetersenColouring {

struct Data {
  vector<vector<int>> normal5_colourings; // just colours
  vector<vector<int>> petersen_colourings; // petersen edge numbers
  vector<string> petersen_profiles;
  vector<int> petersen_poor_counts;
  unordered_map<string, int> profiles;
  unordered_map<string, vector<int>> petersen_colourings_by_profiles;
  vector<bool> petersen_is_strong;
  bool has_strong_petersen_colouring = false;
};

} // ExpPetersenColouring
