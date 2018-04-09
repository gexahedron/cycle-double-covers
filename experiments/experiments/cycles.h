/*
 * File:   cycles.h
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

using namespace std;

namespace NExpCycles {

bool vertex_in_cycle[MAXN];
bool edge_in_cycle[REG * MAXN / 2];
int cycle_length;

int start_vertex[MAXN];
int number_of_circuits = 0;
int separate_circuits[MAXN][2 * MAXN];
int separate_circuits_length[MAXN];
int cycle_count_for_edges[REG * MAXN / 2];

/*********************************Methods*********************************/

void start_or_continue_build_cycle(int possible_edge_lower_bound, TGraph& graph);

void build_cycle(int cur_vertex, int min_possible_edge, TGraph& graph) {
    if (cur_vertex == start_vertex[number_of_circuits]) {
        TMask bit_cycle = 0;
        for (int e = 0; e < graph.number_of_edges; ++e) {
            if (edge_in_cycle[e]) {
                bit_cycle += BIT(e);
                ++cycle_count_for_edges[e];
            }
        }
        graph.all_cycles.insert(bit_cycle);
        if (cycle_length == graph.number_of_vertices) {
            graph.all_full_cycles.insert(bit_cycle);
        }
        if (number_of_circuits == 0) {
            graph.all_circuits.insert(bit_cycle);
            bool is_dominating = true;
            for (int e = 0; e < graph.number_of_edges; ++e) {
                if (!vertex_in_cycle[graph.e2v[e][0]] && !vertex_in_cycle[graph.e2v[e][1]]) {
                    is_dominating = false;
                    break;
                }
            }
            if (is_dominating) {
                graph.all_dominating_circuits.insert(bit_cycle);
            }
        }

        int odd_cycle_count = 0;
        for (int i = 0; i <= number_of_circuits; ++i) {
            if (separate_circuits_length[i] % 2 == 1) {
                ++odd_cycle_count;
            }
        }
        if (odd_cycle_count == 0) {
            graph.all_even_cycles.insert(bit_cycle);
            if (cycle_length == graph.number_of_vertices - 4) {
                graph.all_even_v_minus_4_cycles.insert(bit_cycle);
            }
        }
        if (cycle_length == graph.number_of_vertices) {
            graph.oddness = min(graph.oddness, odd_cycle_count);
        }

        vector<vector<int>> circuits;
        for (int i = 0; i <= number_of_circuits; ++i) {
            vector<int> circuit;
            for (int j = 0; j < separate_circuits_length[i]; ++j) {
                circuit.push_back(separate_circuits[i][j]);
            }
            circuit.push_back(separate_circuits[i][0]);
            circuits.push_back(circuit);
        }
        graph.cycles_as_circuits[bit_cycle] = circuits;

        // section for Petersen graph
        if (graph.number_of_vertices == 10) {
            /*if (cycle_length <= 6) {
                check_tcm_decomposition_in_petersen_graph(graph);
            }*/
        }
    }
    if (cur_vertex == start_vertex[number_of_circuits] && separate_circuits_length[number_of_circuits] > 0) {
        start_or_continue_build_cycle(min_possible_edge + 1, graph);
        return;
    }

    for (int j = 0; j < REG; ++j) {
        int next_vertex = graph.v2v[cur_vertex][j];
        int ei = graph.v2e[cur_vertex][j];

        // conditions
        if (ei < min_possible_edge || edge_in_cycle[ei] || vertex_in_cycle[next_vertex]) {
            continue;
        }

        // initialization
        separate_circuits[number_of_circuits][separate_circuits_length[number_of_circuits]] = next_vertex;
        ++cycle_length;
        ++separate_circuits_length[number_of_circuits];
        vertex_in_cycle[next_vertex] = true;
        edge_in_cycle[ei] = true;

        // recursion
        build_cycle(next_vertex, min_possible_edge, graph);

        // undo
        --cycle_length;
        --separate_circuits_length[number_of_circuits];
        vertex_in_cycle[next_vertex] = false;
        edge_in_cycle[ei] = false;
    }
}

void start_or_continue_build_cycle(int possible_edge_lower_bound, TGraph& graph) {
    for (int min_possible_edge = possible_edge_lower_bound; min_possible_edge < graph.number_of_edges; ++min_possible_edge) {
        int v1 = graph.e2v[min_possible_edge][0];
        int v2 = graph.e2v[min_possible_edge][1];

        // conditions
        if (edge_in_cycle[min_possible_edge] || vertex_in_cycle[v1] || vertex_in_cycle[v2]) {
            continue;
        }

        // initialization
        ++number_of_circuits;
        ++cycle_length;
        separate_circuits[number_of_circuits][0] = v2;
        separate_circuits_length[number_of_circuits] = 1;
        start_vertex[number_of_circuits] = v1;
        vertex_in_cycle[v2] = true;
        edge_in_cycle[min_possible_edge] = true;

        // recursion
        build_cycle(v2, min_possible_edge, graph);

        // undo
        --number_of_circuits;
        --cycle_length;
        vertex_in_cycle[graph.e2v[min_possible_edge][1]] = false;
        edge_in_cycle[min_possible_edge] = false;
    }
}

void prepare_build_cycle(TGraph& graph) {
    graph.all_cycles.clear();
    for (int e = 0; e < graph.number_of_edges; ++e) {
        cycle_count_for_edges[e] = 0;
    }
    cycle_length = 0;
    separate_circuits_length[0] = 0;
    for (int v = 0; v < graph.number_of_vertices; ++v) {
        vertex_in_cycle[v] = false;
    }
    for (int e = 0; e < graph.number_of_edges; ++e) {
        edge_in_cycle[e] = false;
    }
    number_of_circuits = -1;

    start_or_continue_build_cycle(0, graph);
}

} // NExpCycles

