/*
 * File:   experiments.cpp
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * Created on 14 August 2017
 * Here I will study properties and relations between various constructions
 * relevant to snarks: normal (Petersen) colouring, o6c4c, o5cdc, nz5 (nz-mod5, 33-pp)
 *
 */

#pragma once

#include "graph.h"

#include "util/flows.h"

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

namespace NExpNZ5 {

using namespace NUtilFlows;

using namespace std;

int nz5_flow_vals[8] = {-4, -3, -2, -1, 1, 2, 3, 4};
int nz5_edge_flow[REG * MAXN / 2];
int nz5_vertex_flow[MAXN];
set<vector<int>> nz5_orientations;

set<TMask> full_cycles_from_nz5;

/*********************************Methods*********************************/

bool gen_all_nz5_flows(TGraph& graph, int cur_edge_idx) {
    if (cur_edge_idx == graph.number_of_edges) {
        vector<int> flow_vec;
        vector<int> orientations;
        for (int e = 0; e < graph.number_of_edges; ++e) {
            flow_vec.push_back(nz5_edge_flow[e]);
            if (nz5_edge_flow[e] < 0) {
                orientations.push_back(-1);
            } else {
                orientations.push_back(1);
            }
        }
        graph.all_nz5_flows.push_back(flow_vec);
        nz5_orientations.insert(orientations);
        return false;
    }

    int cur_edge = graph.faster_edge_order[cur_edge_idx];
    
    nz5_edge_flow[cur_edge] = 0;
    int v1 = graph.e2v[cur_edge][0];
    int v2 = graph.e2v[cur_edge][1];
    int not_flowed1 = 0;
    for (int j = 0; j < REG; ++j) {
        if (nz5_edge_flow[graph.v2e[v1][j]] == 0) {
            ++not_flowed1;
        }
    }
    int not_flowed2 = 0;
    for (int j = 0; j < REG; ++j) {
        if (nz5_edge_flow[graph.v2e[v2][j]] == 0) {
            ++not_flowed2;
        }
    }

    for (int f_idx = 0; f_idx < 8; ++f_idx) {
        int f = nz5_flow_vals[f_idx];

        if ((not_flowed1 == 1) && nz5_vertex_flow[v1] != f) {
            continue;
        }

        if ((not_flowed2 == 1) && nz5_vertex_flow[v2] != -f) {
            continue;
        }

        nz5_edge_flow[cur_edge] = f;

        nz5_vertex_flow[v1] -= f;
        nz5_vertex_flow[v2] += f;
        if (gen_all_nz5_flows(graph, cur_edge_idx + 1)) {
            return true;
        }
        nz5_vertex_flow[v1] += f;
        nz5_vertex_flow[v2] -= f;

    }
    nz5_edge_flow[cur_edge] = 0;
    return false;
}

void gen_full_cycles_from_nz5_flows(TGraph& graph) {
    full_cycles_from_nz5.clear();
    int good_count = 0;
    int maybe_count = 0;
    int bad_count = 0;

    for (const auto& f : graph.all_nz5_flows) {
        TMask c = build_full_cycle_from_nz5_flow(graph, f);
        if (c == 0) {
            ++bad_count;
        } else if (graph.all_full_cycles.find(c) == graph.all_full_cycles.end()) {
            ++maybe_count;
        } else {
            ++good_count;
            full_cycles_from_nz5.insert(c);
        }
    }
    cerr << "good  nz5: " << (float) good_count / graph.all_nz5_flows.size() << endl;
    cerr << "bad   nz5: " << (float) bad_count / graph.all_nz5_flows.size() << endl;
    cerr << "maybe nz5: " << (float) maybe_count / graph.all_nz5_flows.size() << endl;
    cerr << "full cycles from nz5 vs all: " << full_cycles_from_nz5.size() << " " << graph.all_full_cycles.size() << endl;
}

void find_all_nz5_flows(TGraph& graph) {
    nz5_orientations.clear();
    for (int v = 0; v < graph.number_of_vertices; ++v) {
        nz5_vertex_flow[v] = 0;
    }
    for (int e = 0; e < graph.number_of_edges; ++e) {
        nz5_edge_flow[e] = 0;
    }

    gen_all_nz5_flows(graph, 0);
    gen_full_cycles_from_nz5_flows(graph);

    cerr << "number of nz5 flows: " << graph.all_nz5_flows.size() << endl;
    cerr << "number of nz5 orientations: " << nz5_orientations.size() << endl;
}

} // NExpNZ5

