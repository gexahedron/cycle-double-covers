/*
 * File:   o5cdc.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
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
#include <algorithm>
#include <string>

namespace NExp5cdc {

using namespace std;

set<pair<TMask, TMask>> u33pp_pairs;

set<set<TMask>> all_5cdc;
set<TMask> cur_5cdc;

set<TMask> u333pp_cycles_from_o5cdc;
set<TMask> u333pp_even_cycles_from_o5cdc;

set<TMask> all_cycles_from_5cdc;
set<TMask> all_full_cycles_from_5cdc;
set<TMask> all_even_cycles_from_5cdc;
set<TMask> all_circuits_from_5cdc;

int start_vertex_in_5cdc[5][MAXN];
bool vertex_in_5cdc[5][MAXN];
bool edge_in_5cdc[5][REG * MAXN / 2];
int cycle_length_in_5cdc[5];
int number_of_circuits_in_5cdc[5];
int separate_circuits_in_5cdc[5][MAXN][2 * MAXN];
int separate_circuits_length_in_5cdc[5][MAXN];
int vertex_count_in_5cdc[MAXN / 2];
int edge_count_in_5cdc[REG * MAXN / 2];
int total_edge_count_in_5cdc;
//int edge_pair_count_in_5cdc[REG * MAXN];
TMask u5cdc_cycles[5];

int o5cdc_aggregated_solutions;
int all_o5cdc_solutions;
int all_5cdc_solutions;
int same_cycles_different_orientations_in_5cdc;
vector<vector<int>> all_circuits_in_5cdc;
int edge_orientation_count_in_5cdc[REG * MAXN / 2][2];
int orientations_in_5cdc[5 * MAXN];
int layer_in_5cdc[5 * MAXN];
int layer_flow_in_5cdc[5][REG * MAXN / 2];

int u33pp_indices[15][5] = {{0, 1, 2, 3, 4}, {0, 2, 1, 3, 4}, {0, 3, 1, 2, 4},
                            {0, 1, 2, 4, 3}, {0, 2, 1, 4, 3}, {0, 4, 1, 2, 3},
                            {0, 1, 3, 4, 2}, {0, 3, 1, 4, 2}, {0, 4, 1, 3, 2},
                            {0, 2, 3, 4, 1}, {0, 3, 2, 4, 1}, {0, 4, 2 ,3, 1},
                            {1, 2, 3, 4, 0}, {1, 3, 2, 4, 0}, {1, 4, 2, 3, 0}};

int u33pp_count;

/*********************************Methods*********************************/

void gen_33pp_from_o5cdc(TGraph& graph) {
    for (int i = 0; i < 15; ++i) {
        TMask m1 = 0;
        TMask m2 = 0;
        for (int e = 0; e < graph.number_of_edges; ++e) {
            bool in[5];
            for (int j = 0; j < 5; ++j) {
                in[j] = (BIT(e) & u5cdc_cycles[u33pp_indices[i][j]]) > 0;
            }
            if ((in[0] && in[1]) || (in[2] && in[3]) || in[4]) {
                m1 |= BIT(e);
                m2 |= BIT(e);
            }
            if ((in[0] && in[2]) || (in[1] && in[3])) {
                m1 |= BIT(e);
            }
            if ((in[0] && in[3]) || (in[1] && in[2])) {
                m2 |= BIT(e);
            }
        }
        if (m1 > m2) {
            swap(m1, m2);
        }
        u33pp_pairs.insert(make_pair(m1, m2));
    }
    for (int i = 0; i < 5; ++i) {
        u333pp_cycles_from_o5cdc.insert(u5cdc_cycles[i]);
        if (graph.all_even_cycles.find(u5cdc_cycles[i]) != graph.all_even_cycles.end()) {
            u333pp_even_cycles_from_o5cdc.insert(u5cdc_cycles[i]);
        }
    }
}

bool orient_5cdc(TGraph& graph, int cur_circuit) {
    if (cur_circuit == all_circuits_in_5cdc.size()) {
        ++same_cycles_different_orientations_in_5cdc;
        if (same_cycles_different_orientations_in_5cdc > 1) {
            cerr << "wut" << endl;
        }
        gen_33pp_from_o5cdc(graph);
        return false;//true;
    }
    int max_orientation = 1;
    if (cur_circuit == 0) {
        max_orientation = 0;
    }
    for (int orientation = 0; orientation <= max_orientation; ++orientation) {
        orientations_in_5cdc[cur_circuit] = orientation;
        int vi = 0;
        while (vi < all_circuits_in_5cdc[cur_circuit].size() - 1) {
            int v1 = all_circuits_in_5cdc[cur_circuit][vi];
            int v2 = all_circuits_in_5cdc[cur_circuit][vi + 1];
            int ei = graph.edge_index[v1][v2];
            int cur_edge_orientation = orientation;
            if (v1 > v2) {
                cur_edge_orientation = 1 - cur_edge_orientation;
            }
            if (edge_orientation_count_in_5cdc[ei][cur_edge_orientation] == 1) {
                break;
            }
            ++edge_orientation_count_in_5cdc[ei][cur_edge_orientation];
            ++vi;
        }
        if (vi == all_circuits_in_5cdc[cur_circuit].size() - 1) {
            if (orient_5cdc(graph, cur_circuit + 1)) {
                return true;
            }
        }
        --vi;
        while (vi >= 0) {
            int v1 = all_circuits_in_5cdc[cur_circuit][vi];
            int v2 = all_circuits_in_5cdc[cur_circuit][vi + 1];
            int ei = graph.edge_index[v1][v2];
            int cur_edge_orientation = 0;
            if (v1 > v2) {
                cur_edge_orientation = 1 - cur_edge_orientation;
            }
            if (orientation == 1) {
                cur_edge_orientation = 1 - cur_edge_orientation;
            }
            --edge_orientation_count_in_5cdc[ei][cur_edge_orientation];
            --vi;
        }
    }
    return false;
}

bool check_orientability_5cdc(TGraph& graph) {
    vector<vector<int>> circuits;
    all_circuits_in_5cdc.clear();
    for (int i = 0; i < 5; ++i) {
        circuits = graph.cycles_as_circuits[u5cdc_cycles[i]];
        for (const auto& circuit : circuits) {
            layer_in_5cdc[all_circuits_in_5cdc.size()] = i;
            all_circuits_in_5cdc.push_back(circuit);
        }
    }
    //cerr << "all circuits: " << all_circuits_in_5cdc.size() << endl;
    for (int e = 0; e < graph.number_of_edges; ++e) {
        for (int orientation = 0; orientation < 2; ++orientation) {
            edge_orientation_count_in_5cdc[e][orientation] = 0;
        }
    }
    same_cycles_different_orientations_in_5cdc = 0;
    orient_5cdc(graph, 0);
    if (same_cycles_different_orientations_in_5cdc > 0) {
        ++o5cdc_aggregated_solutions;
        all_o5cdc_solutions += same_cycles_different_orientations_in_5cdc;
        graph.all_o5cdc.insert(cur_5cdc);
    }
    /*if (same_cycles_different_orientations_in_5cdc > 1) {
        cerr << "found " << same_cycles_different_orientations_in_5cdc << " same cycles different orientations for 5cdc" << endl;
    }*/
    return false;
}

void gen_33pp_from_5cdc_with_orientations(TGraph& graph) {
    for (int i = 0; i < 15; ++i) {
        TMask m1 = 0;
        TMask m2 = 0;
        for (int e = 0; e < graph.number_of_edges; ++e) {
            int flows[4];
            for (int j = 0; j < 4; ++j) {
                flows[j] = layer_flow_in_5cdc[u33pp_indices[i][j]][e];
            }
            if ((flows[0] - flows[1]) + (flows[2] - flows[3]) != 0) {
                m1 |= BIT(e);
            }
            if ((flows[0] - flows[1]) - (flows[2] - flows[3]) != 0) {
                m2 |= BIT(e);
            }
        }
        if ((m1 | m2) != BIT(graph.number_of_edges) - 1) {
            continue;
        }
        if (m1 > m2) {
            swap(m1, m2);
        }

        u5cdc_from_33pp[make_pair(m1, m2)].insert(cur_5cdc);
        ++u33pp_count;
    }
}

void gen_33pp_from_5cdc(TGraph& graph, int cur_circuit) {
    if (cur_circuit == all_circuits_in_5cdc.size()) {
        gen_33pp_from_5cdc_with_orientations(graph);
        return;
    }
    int max_orientation = 1;
    if (cur_circuit == 0) {
        max_orientation = 0;
    }
    for (int orientation = 0; orientation <= max_orientation; ++orientation) {
        orientations_in_5cdc[cur_circuit] = orientation;
        int vi = 0;
        while (vi < all_circuits_in_5cdc[cur_circuit].size() - 1) {
            int v1 = all_circuits_in_5cdc[cur_circuit][vi];
            int v2 = all_circuits_in_5cdc[cur_circuit][vi + 1];
            int ei = graph.edge_index[v1][v2];
            int cur_edge_orientation = orientation;
            if (v1 > v2) {
                cur_edge_orientation = 1 - cur_edge_orientation;
            }
            //if (edge_orientation_count_in_5cdc[ei][cur_edge_orientation] == 1) {
                // TODO: record this event as inconsistency between layers
            //}
            if (cur_edge_orientation == 0) {
                layer_flow_in_5cdc[layer_in_5cdc[cur_circuit]][ei] += 1;
            } else {
                layer_flow_in_5cdc[layer_in_5cdc[cur_circuit]][ei] -= 1;
            }

            ++edge_orientation_count_in_5cdc[ei][cur_edge_orientation];
            ++vi;
        }

        gen_33pp_from_5cdc(graph, cur_circuit + 1);

        --vi;
        while (vi >= 0) {
            int v1 = all_circuits_in_5cdc[cur_circuit][vi];
            int v2 = all_circuits_in_5cdc[cur_circuit][vi + 1];
            int ei = graph.edge_index[v1][v2];
            int cur_edge_orientation = 0;
            if (v1 > v2) {
                cur_edge_orientation = 1 - cur_edge_orientation;
            }
            if (orientation == 1) {
                cur_edge_orientation = 1 - cur_edge_orientation;
            }
            if (cur_edge_orientation == 0) {
                layer_flow_in_5cdc[layer_in_5cdc[cur_circuit]][ei] -= 1;
            } else {
                layer_flow_in_5cdc[layer_in_5cdc[cur_circuit]][ei] += 1;
            }
            --edge_orientation_count_in_5cdc[ei][cur_edge_orientation];
            --vi;
        }
    }
    return;
}

void start_or_continue_build_5cdc(int possible_edge_lower_bound, TGraph& graph, int cur_cycle_layer, bool only_find);

void build_5cdc(int cur_vertex, int min_possible_edge, TGraph& graph, int cur_cycle_layer, bool only_find) {
    if (cur_vertex == start_vertex_in_5cdc[cur_cycle_layer][number_of_circuits_in_5cdc[cur_cycle_layer]]) {
        bool has_leaf = false;
        for (int v = 0; v < graph.number_of_vertices; ++v) {
            int deg = REG;
            for (int j = 0; j < REG; ++j) {
                if (edge_count_in_5cdc[graph.v2e[v][j]] == 2) {
                    --deg;
                }
            }
            if (deg == 1) {
                has_leaf = true;
                break;
            }
        }
        if (!has_leaf) {
            TMask bit_cycle = 0;
            for (int e = 0; e < graph.number_of_edges; ++e) {
                if (edge_in_5cdc[cur_cycle_layer][e]) {
                    bit_cycle += BIT(e);
                }
            }
            u5cdc_cycles[cur_cycle_layer] = bit_cycle;
            start_or_continue_build_5cdc(0, graph, cur_cycle_layer + 1, only_find);
        }
    }
    if (cur_vertex == start_vertex_in_5cdc[cur_cycle_layer][number_of_circuits_in_5cdc[cur_cycle_layer]] && separate_circuits_length_in_5cdc[cur_cycle_layer][number_of_circuits_in_5cdc[cur_cycle_layer]] > 0) {
        start_or_continue_build_5cdc(min_possible_edge + 1, graph, cur_cycle_layer, only_find);
        return;
    }

    for (int j = 0; j < REG; ++j) {
        int next_vertex = graph.v2v[cur_vertex][j];
        int ei = graph.v2e[cur_vertex][j];

        // conditions
        if (ei < min_possible_edge || edge_in_5cdc[cur_cycle_layer][ei] || vertex_in_5cdc[cur_cycle_layer][next_vertex]) {
            continue;
        }
        // specific 5cdc conditions
        if (edge_count_in_5cdc[ei] == 2) {
            continue;
        }

        // initialization
        separate_circuits_in_5cdc[cur_cycle_layer][number_of_circuits_in_5cdc[cur_cycle_layer]][separate_circuits_length_in_5cdc[cur_cycle_layer][number_of_circuits_in_5cdc[cur_cycle_layer]]] = next_vertex;
        ++cycle_length_in_5cdc[cur_cycle_layer];
        ++separate_circuits_length_in_5cdc[cur_cycle_layer][number_of_circuits_in_5cdc[cur_cycle_layer]];
        vertex_in_5cdc[cur_cycle_layer][next_vertex] = true;
        edge_in_5cdc[cur_cycle_layer][ei] = true;
        ++vertex_count_in_5cdc[next_vertex];
        ++edge_count_in_5cdc[ei];
        ++total_edge_count_in_5cdc;

        // recursion
        build_5cdc(next_vertex, min_possible_edge, graph, cur_cycle_layer, only_find);

        // undo
        --cycle_length_in_5cdc[cur_cycle_layer];
        --separate_circuits_length_in_5cdc[cur_cycle_layer][number_of_circuits_in_5cdc[cur_cycle_layer]];
        vertex_in_5cdc[cur_cycle_layer][next_vertex] = false;
        edge_in_5cdc[cur_cycle_layer][ei] = false;
        --vertex_count_in_5cdc[next_vertex];
        --edge_count_in_5cdc[ei];
        --total_edge_count_in_5cdc;
    }
}

void start_or_continue_build_5cdc(int possible_edge_lower_bound, TGraph& graph, int cur_cycle_layer, bool only_find) {
    if (cur_cycle_layer == 4) {
        u5cdc_cycles[4] = (BIT(graph.number_of_edges) - 1) * 2 - u5cdc_cycles[0] - u5cdc_cycles[1] - u5cdc_cycles[2] - u5cdc_cycles[3];
        if (graph.all_cycles.find(u5cdc_cycles[4]) == graph.all_cycles.end()) {
            return;
        }
    //if (cur_cycle_layer == 5) {
        /*for (int e = 0; e < graph.number_of_edges; ++e) {
            if (edge_count_in_5cdc[e] != 2) {
                return;
            }
        }*/

        cur_5cdc.clear(); // actually it's more right to choose vector, not set; e.g., for 10c6c
        for (int i = 0; i < 5; ++i) {
            cur_5cdc.insert(u5cdc_cycles[i]);
        }
        if (all_5cdc.find(cur_5cdc) != all_5cdc.end()) {
            return;
        }
        all_5cdc.insert(cur_5cdc);
        graph.all_5cdc.insert(cur_5cdc);
        if (only_find) {
            return;
        }

        vector<vector<int>> circuits;
        all_circuits_in_5cdc.clear();
        for (int i = 0; i < 5; ++i) {
            circuits = graph.cycles_as_circuits[u5cdc_cycles[i]];
            for (const auto& circuit : circuits) {
                layer_in_5cdc[all_circuits_in_5cdc.size()] = i;
                all_circuits_in_5cdc.push_back(circuit);
            }
        }

        for (int e = 0; e < graph.number_of_edges; ++e) {
            for (int i = 0; i < 5; ++i) {
                layer_flow_in_5cdc[i][e] = 0;
            }
        }

        u33pp_count = 0;
        gen_33pp_from_5cdc(graph, 0);
        if (graph.petersen_5cdc.find(cur_5cdc) != graph.petersen_5cdc.end()) {
            //cerr << "pet 5cdc: " << u33pp_count << endl;
        }

        for (int i = 0; i < 5; ++i) {
            all_cycles_from_5cdc.insert(u5cdc_cycles[i]);
            if (graph.all_full_cycles.find(u5cdc_cycles[i]) != graph.all_full_cycles.end()) {
                all_full_cycles_from_5cdc.insert(u5cdc_cycles[i]);
            }
            if (graph.all_even_cycles.find(u5cdc_cycles[i]) != graph.all_even_cycles.end()) {
                all_even_cycles_from_5cdc.insert(u5cdc_cycles[i]);
            }
            if (graph.all_circuits.find(u5cdc_cycles[i]) != graph.all_circuits.end()) {
                all_circuits_from_5cdc.insert(u5cdc_cycles[i]);
            }

        }

        check_orientability_5cdc(graph);
        return;
    }
    if (possible_edge_lower_bound == 0) {
        if (2 * graph.number_of_edges - total_edge_count_in_5cdc > graph.number_of_vertices * (5 - cur_cycle_layer)) {
            return;
        }
        if (cur_cycle_layer == 3) {
            for (int v = 0; v < graph.number_of_vertices; ++v) {
                if (vertex_count_in_5cdc[v] == 0) {
                    return;
                }
            }
        }
    }
    for (int min_possible_edge = possible_edge_lower_bound; min_possible_edge < graph.number_of_edges; ++min_possible_edge) {
        int v1 = graph.e2v[min_possible_edge][0];
        int v2 = graph.e2v[min_possible_edge][1];

        // conditions
        if (edge_in_5cdc[cur_cycle_layer][min_possible_edge] || vertex_in_5cdc[cur_cycle_layer][v1] || vertex_in_5cdc[cur_cycle_layer][v2]) {
            continue;
        }
        // specific 5cdc conditions
        if (edge_count_in_5cdc[min_possible_edge] == 2) {
            continue;
        }

        // initialization
        ++number_of_circuits_in_5cdc[cur_cycle_layer];
        ++cycle_length_in_5cdc[cur_cycle_layer];
        separate_circuits_in_5cdc[cur_cycle_layer][number_of_circuits_in_5cdc[cur_cycle_layer]][0] = v2;
        separate_circuits_length_in_5cdc[cur_cycle_layer][number_of_circuits_in_5cdc[cur_cycle_layer]] = 1;
        start_vertex_in_5cdc[cur_cycle_layer][number_of_circuits_in_5cdc[cur_cycle_layer]] = v1;
        vertex_in_5cdc[cur_cycle_layer][v2] = true;
        edge_in_5cdc[cur_cycle_layer][min_possible_edge] = true;
        ++vertex_count_in_5cdc[v2];
        ++edge_count_in_5cdc[min_possible_edge];
        ++total_edge_count_in_5cdc;

        // recursion
        build_5cdc(v2, min_possible_edge, graph, cur_cycle_layer, only_find);

        // undo
        --number_of_circuits_in_5cdc[cur_cycle_layer];
        --cycle_length_in_5cdc[cur_cycle_layer];
        vertex_in_5cdc[cur_cycle_layer][graph.e2v[min_possible_edge][1]] = false;
        edge_in_5cdc[cur_cycle_layer][min_possible_edge] = false;
        --vertex_count_in_5cdc[v2];
        --edge_count_in_5cdc[min_possible_edge];
        --total_edge_count_in_5cdc;

        if (possible_edge_lower_bound == 0) {
            break;
        }
    }
}

void prepare_build_5cdc(TGraph& graph, bool only_find) {
    for (int layer = 0; layer < 5; ++layer) {
        cycle_length_in_5cdc[layer] = 0;
        separate_circuits_length_in_5cdc[layer][0] = 0;

        for (int v = 0; v < graph.number_of_vertices; ++v) {
            vertex_in_5cdc[layer][v] = false;
        }
        for (int e = 0; e < graph.number_of_edges; ++e) {
            edge_in_5cdc[layer][e] = false;
        }
        number_of_circuits_in_5cdc[layer] = -1;
    }

    for (int v = 0; v < graph.number_of_vertices; ++v) {
        vertex_count_in_5cdc[v] = 0;
    }
    for (int e = 0; e < graph.number_of_edges; ++e) {
        edge_count_in_5cdc[e] = 0;
    }
    total_edge_count_in_5cdc = 0;

    start_or_continue_build_5cdc(0, graph, 0, only_find);
}

void find_all_o5cdc(TGraph& graph, bool only_find = false) {
    u5cdc_from_33pp.clear();
    all_5cdc.clear();
    all_5cdc_solutions = 0;
    all_o5cdc_solutions = 0;
    u33pp_pairs.clear();
    u333pp_cycles_from_o5cdc.clear();
    u333pp_even_cycles_from_o5cdc.clear();

    all_cycles_from_5cdc.clear();
    all_full_cycles_from_5cdc.clear();
    all_even_cycles_from_5cdc.clear();
    all_circuits_from_5cdc.clear();

    prepare_build_5cdc(graph, only_find);

    //cerr << "all 5cdc solutions: " << all_5cdc_solutions << endl;
    //cerr << "all o5cdc solutions: " << all_o5cdc_solutions << endl;
    // FIXME: these numbers are wrong, they are only upper bounds right now
}

} // NExp5cdc
