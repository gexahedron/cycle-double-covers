/*
 * File:   args.cpp
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * Created on 21 August 2019
 *
 * Contains data structure for parsing arguments
 *
 */

#include "args.h"
#include "constants.h"

#include <iostream>
#include <cstdlib>
#include <fstream>

using namespace std;


Args::Args(int argc, char** argv)
  : start_idx(NONE)
  , finish_idx(NONE)
{
  if (argc < 2 || argc > 4) {
    cerr << "Error: invalid number of arguments" << endl;
    cerr << "Usage: " << argv[0] << " <filetype: mc, g6, adj, bghm> " <<
        "<index_of_start_graph|file> <index_of_finish_graph>" << endl;
    exit(1);
  } else {
    filetype = argv[1];
    if (argc == 2) {
      start_idx = 0;
    } else {
      assert(argc >= 3);
      if (strcmp(argv[2], "file") == 0) { // if strings are equal
        assert(argc >= 4);
        std::ifstream infile(argv[3]);
        int idx;
        while (infile >> idx) {
          idxs.insert(idx);
        }
      } else {
        if (argc >= 3) {
          start_idx = atoi(argv[2]);
        }
        if (argc >= 4) {
          finish_idx = atoi(argv[3]);
        }
      }
    }
  }
}
