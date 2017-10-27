/*
 * File:   preimages.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
 *
 */

#pragma once

#include "graph.h"

// TODO: remove unused includes
#include <cassert>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <map>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>

namespace NExpPreimages {

using namespace std;

unordered_set<TMask> full_cycles_from_petersen;

/*********************************Methods*********************************/

void preimage_full_cycles(const TGraph& petersen_graph, TGraph& graph) {
    full_cycles_from_petersen.clear();
    for (int col_idx = 0; col_idx < graph.petersen_colourings.size(); ++col_idx) {
        for (const auto c : petersen_graph.all_full_cycles) {
            TMask cycle_mask = 0;
            for (int e = 0; e < graph.number_of_edges; ++e) {
                int petersen_edge = graph.petersen_colourings[col_idx][e];
                if (c & BIT(petersen_edge)) {
                    cycle_mask += BIT(e);
                }
            }
            //cerr << "from: " << c << " to: " << cycle_mask << endl;
            full_cycles_from_petersen.insert(cycle_mask);
        }
    }
}

void preimage_6c4c_5cdc_cycles(const TGraph& petersen_graph, TGraph& graph) {
    for (int col_idx = 0; col_idx < graph.petersen_colourings.size(); ++col_idx) {
        for (const auto& petersen_6c4c : petersen_graph.all_6c4c) {
            set<TMask> cur_6c4c;
            for (const auto& c : petersen_6c4c) {
                TMask cycle_mask = 0;
                for (int e = 0; e < graph.number_of_edges; ++e) {
                    int petersen_edge = graph.petersen_colourings[col_idx][e];
                    if (c & BIT(petersen_edge)) {
                        cycle_mask += BIT(e);
                    }
                }
                cur_6c4c.insert(cycle_mask);
            }
            graph.petersen_6c4c.insert(cur_6c4c);
            for (const auto& petersen_5cdc : petersen_graph.all_5cdc) {
                set<TMask> cur_5cdc;
                for (const auto& c : petersen_5cdc) {
                    TMask cycle_mask = 0;
                    for (int e = 0; e < graph.number_of_edges; ++e) {
                        int petersen_edge = graph.petersen_colourings[col_idx][e];
                        if (c & BIT(petersen_edge)) {
                            cycle_mask += BIT(e);
                        }
                    }
                    cur_5cdc.insert(cycle_mask);
                }
                graph.petersen_5cdc.insert(cur_5cdc);
                graph.petersen_6c4c_5cdc_pairs[cur_6c4c].insert(cur_5cdc);
                graph.petersen_5cdc_6c4c_pairs[cur_5cdc].insert(cur_6c4c);
            }
        }
    }
}

} // NExpPreimages
