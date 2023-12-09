/*
 * File:   args.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * Created on 21 August 2019
 *
 * Contains data structure for parsing arguments
 *
 */

#pragma once

#include <string>
#include <unordered_set>

using namespace std;


struct Args {
  int start_idx;
  int finish_idx;
  unordered_set<int> idxs;

  string filetype;

  Args(int argc, char** argv);
};
