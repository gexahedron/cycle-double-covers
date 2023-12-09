/*
 * File:   matching_5cdc_and_6c4c.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
 * Depends on: preimages, o6c4c
 *
 */

#pragma once

#include "graph.h"
#include "common.h"

// TODO: remove unused includes
#include <cassert>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <map>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>

using namespace std;


namespace Exp5cdcVs6c4c {

/*********************************Methods*********************************/

void compare_6c4c_5cdc_pairs(Graph& graph) {
    //cerr << "all 6c4c: " << all_6c4c_solutions << endl;
    cerr << "all 6c4c (in 6c4c-5cdc pairs): " << graph.u244_6c4c_5cdc_pairs.size() << endl;
    cerr << "pet 6c4c: " << graph.petersen_6c4c_5cdc_pairs.size() << endl;
    for (const auto& u6c4c_pairs : graph.u244_6c4c_5cdc_pairs) {
        for (const auto& u5cdc : u6c4c_pairs.second) {
            if (graph.petersen_6c4c_5cdc_pairs.find(u6c4c_pairs.first) != graph.petersen_6c4c_5cdc_pairs.end()) {
                if (graph.petersen_6c4c_5cdc_pairs[u6c4c_pairs.first].find(u5cdc) != graph.petersen_6c4c_5cdc_pairs[u6c4c_pairs.first].end()) {
                    cerr << "both" << endl;
                    //return;
                }
            }
        }
    }
    for (const auto& petersen_6c4c_pairs : graph.petersen_6c4c_5cdc_pairs) {
        cerr << "is petersen 6c4c oriented: " << (graph.all_o6c4c.find(petersen_6c4c_pairs.first) != graph.all_o6c4c.end()) << endl;
        for (const auto& petersen_5cdc : petersen_6c4c_pairs.second) {
            //cerr << "is petersen 5cdc oriented: " << (graph.all_o5cdc.find(petersen_5cdc) != graph.all_o5cdc.end()) << endl;
        }
    }
    for (const auto& petersen_5cdc_pairs : graph.petersen_5cdc_6c4c_pairs) {
        //cerr << "is petersen 5cdc in 6c4c-244-flows-33-pp construction: " << (u244_5cdc_6c4c_pairs.find(petersen_5cdc_pairs.first) != u244_5cdc_6c4c_pairs.end()) << endl;
    }
    for (const auto& petersen_6c4c : graph.petersen_6c4c) {
        //cerr << "is petersen 6c4c oriented: " << (graph.all_o6c4c.find(petersen_6c4c) != graph.all_o6c4c.end()) << endl;
    }
}

} // namespace Exp5cdcVs6c4c
