/*
 * File:   flow_parity_pairs.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
 * Depends on: cycles
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

namespace NExpFlowParityPairs {

using namespace NUtilFlows;

using namespace std;

int nz5_from_33pp_coefs[4][2] = {{3, 1}, {3, -1}, {1, 3}, {1, -3}};

// TODO: there is duplicate code in o6c4c about 33pp
bool vertex_in_33pp_cycle[MAXN];
bool edge_in_33pp_cycle[REG * MAXN / 2];
int u33pp_partition[REG * MAXN / 2];

set<TMask> all_33pp_cycles;
set<TMask> all_33pp_circuits;
set<TMask> all_33pp_full_cycles;

set<TMask> all_333pp_cycles;
set<TMask> all_333pp_even_cycles;
set<TMask> all_333pp_full_cycles;
/*********************************Methods*********************************/

bool check_33pp(TGraph& graph, const TMask c) {
    int three_flow_count = 0;
    TMask part_mask[3];
    bool has_3flow[3];
    bool had_more_than_one = false;
    bool edge_in_cur_part[3][REG * MAXN / 2];
    bool vertex_in_cur_part[3][MAXN];
    int vertex_colour[3][MAXN];
    int max_colour_pair[3];
    int edge_orientation[3][REG * MAXN / 2];
    int edge_colour[3][REG * MAXN / 2];
    for (int part = 0; part < 3; ++part) {
        part_mask[part] = 0;
        for (int v = 0; v < graph.number_of_vertices; ++v) {
            vertex_colour[part][v] = -1;
            vertex_in_cur_part[part][v] = false;
        }
        for (int e = 0; e < graph.number_of_edges; ++e) {
            if (edge_in_33pp_cycle[e]) {
                edge_in_cur_part[part][e] = true;
            } else {
                edge_in_cur_part[part][e] = (u33pp_partition[e] != part);

                if (edge_in_cur_part[part][e]) {
                    for (int i = 0; i < 2; ++i) {
                       if (vertex_in_33pp_cycle[graph.e2v[e][i]]) {
                            vertex_in_cur_part[part][graph.e2v[e][i]] = true;
                        }
                    }
                }
            }
            if (edge_in_cur_part[part][e]) {
                part_mask[part] += BIT(e);
            }
        }

        bool edge_visited[REG * MAXN / 2];
        for (int e = 0; e < graph.number_of_edges; ++e) {
            edge_visited[e] = !edge_in_cur_part[part][e];
        }
        int queue[MAXN];
        int queue_size;
        has_3flow[part] = true;
        bool has_one = false;
        bool has_more = false;
        int cur_colour = -2;
        max_colour_pair[part] = -1;
        for (int v = 0; v < graph.number_of_vertices; ++v) {
            // FIXME: not all vertices and edges are here
            // 1. some of them could be not connected to cycle, they will have flow 2
            // 2. there could possibly be a circuit with flow 1 without edges of flow 2
            if (vertex_in_cur_part[part][v] && vertex_colour[part][v] == -1) {
                cur_colour += 2;
                ++max_colour_pair[part];
                if (has_one) {
                    has_more = true;
                }
                has_one = true;
                vertex_colour[part][v] = cur_colour;
                queue_size = 1;
                queue[0] = v;
                for (int cur_idx = 0; cur_idx < queue_size; ++cur_idx) {
                    int cur_vertex = queue[cur_idx];
                    for (int j = 0; j < REG; ++j) {
                        int ei_start = graph.v2e[cur_vertex][j];
                        if (edge_visited[ei_start]) {
                            continue;
                        }
                        edge_visited[ei_start] = true;
                        edge_colour[part][ei_start] = max_colour_pair[part];
                        int v1 = cur_vertex;
                        int v2 = graph.v2v[cur_vertex][j];
                        edge_orientation[part][ei_start] = 1;
                        if (v1 > v2) {
                            edge_orientation[part][ei_start] *= -1;
                        }
                        if (vertex_colour[part][cur_vertex] == cur_colour + 1) {
                            edge_orientation[part][ei_start] *= -1;
                        }
                        if (edge_in_33pp_cycle[ei_start]) {
                            edge_orientation[part][ei_start] *= -1;
                        }
                        while (!vertex_in_cur_part[part][v2]) {
                            for (int j2 = 0; j2 < REG; ++j2) {
                                int e = graph.v2e[v2][j2];
                                int v3 = graph.v2v[v2][j2];
                                if (!edge_visited[e] && v3 != v1) {
                                    v1 = v2;
                                    v2 = v3;
                                    edge_visited[e] = true;
                                    edge_colour[part][e] = max_colour_pair[part];
                                    edge_orientation[part][e] = 1;
                                    if (v1 > v2) {
                                        edge_orientation[part][e] *= -1;
                                    }
                                    if (vertex_colour[part][cur_vertex] == cur_colour + 1) {
                                        edge_orientation[part][e] *= -1;
                                    }
                                    if (edge_in_33pp_cycle[e]) {
                                        edge_orientation[part][e] *= -1;
                                    }
                                    break;
                                }
                            }
                        }
                        if (vertex_colour[part][v2] == -1) {
                            vertex_colour[part][v2] = 2 * cur_colour + 1 - vertex_colour[part][cur_vertex];
                            queue[queue_size] = v2;
                            ++queue_size;
                        } else if (vertex_colour[part][v2] == vertex_colour[part][cur_vertex]) {
                            has_3flow[part] = false;
                            break;
                        } else if (vertex_colour[part][v2] < cur_colour) {
                            cerr << "wtf" << endl;
                        }
                    }
                    if (!has_3flow[part]) {
                        break;
                    }
                }
            }
            if (!has_3flow[part]) {
                break;
            }
        }
        if (has_3flow[part]) {
            ++three_flow_count;
            if (has_more) {
                had_more_than_one = true;
            }
        }
    }
 
    if (three_flow_count < 3) {
        return false;
    }
    /*if (had_more_than_one) {
        cerr << "more than 1" << endl;
    }*/
    all_33pp_cycles.insert(c);
    if (graph.all_full_cycles.find(c) != graph.all_full_cycles.end()) {
        all_33pp_full_cycles.insert(c);
    }
    if (graph.all_circuits.find(c) != graph.all_circuits.end()) {
        all_33pp_circuits.insert(c);
    }

    /*if (three_flow_count == 3) {
        cerr << "tf\t";
    }*/

    // now i can find 33-pp pairs and build nz5,
    // using weights 3/2 and +- 1/2, 1/2 and +- 3/2
    // and using arrays: edge_colour, edge_orientation, max_colour_pair
    vector<vector<int>> part_pairs = {{0, 1}, {0, 2}, {1, 2}};
    for (const auto& parts_pair : part_pairs) {
        if (!has_3flow[parts_pair[0]] || !has_3flow[parts_pair[1]]) {
            continue;
        }
        int edge_flows[2][REG * MAXN / 2];
        int masks[2];
        for (int m1 = 0; m1 < 1; ++m1) {//BIT(max_colour_pair[parts_pair[0]]); ++m1) {
            masks[0] = m1;
            for (int m2 = 0; m2 < 1; ++m2) {//BIT(max_colour_pair[parts_pair[1]] + 1); ++m2) {
                masks[1] = m2;
                for (int i = 0; i < 2; ++i) {
                    int cur_part = parts_pair[i];
                    for (int e = 0; e < graph.number_of_edges; ++e) {
                        int flow_val = 0;
                        if (edge_in_33pp_cycle[e]) { // flow = 1
                            flow_val = 1;
                        } else if (edge_in_cur_part[cur_part][e]) { // flow = 2
                            flow_val = 2;
                        } else { // flow = 0
                            edge_flows[i][e] = 0;
                            continue;
                        }
                        flow_val *= edge_orientation[cur_part][e];
                        if ((BIT(edge_colour[cur_part][e]) & masks[i]) > 0) {
                            flow_val *= -1;
                        }
                        edge_flows[i][e] = flow_val;
                    }
                }

                // sanity checks
                bool is_hard = false;
                for (int e = 0; e < graph.number_of_edges; ++e) {
                    int f1 = abs(edge_flows[0][e]);
                    int f2 = abs(edge_flows[1][e]);
                    // FIXME: should fix and remove all is_hard
                    if (f1 == 0 && f2 == 0) {
                        is_hard = true;
                        break;
                    }
                    assert(f1 != 0 || f2 != 0);
                    if (f1 % 2 != f2 % 2) {
                        is_hard = true;
                        break;
                    }
                    assert((f1 % 2) == (f2 % 2));
                }
                if (is_hard) {
                    continue;
                }
                for (int i = 0; i < 2; ++i) {
                    for (int v = 0; v < graph.number_of_vertices; ++v) {
                        int flow = 0;
                        for (int j = 0; j < REG; ++j) {
                            int v2 = graph.v2v[v][j];
                            int edge_flow = edge_flows[i][graph.v2e[v][j]];
                            if (v > v2) {
                                edge_flow *= -1;
                            }
                            flow += edge_flow;
                        }
                        // FIXME
                        if (flow != 0) {
                            //print_graph(graph);
                            /*cerr << "parts: " << parts_pair[0] << " " << parts_pair[1] << endl;
                            for (int v = 0; v < graph.number_of_vertices; ++v) {
                                cerr << v << " ";
                                for (int j = 0; j < REG; ++j) {
                                    int e = graph.v2e[v][j];
                                    cerr << graph.v2v[v][j] << "(" << edge_in_33pp_cycle[e] << ", ";
                                    cerr << edge_in_cur_part[parts_pair[0]][e] << " " << edge_in_cur_part[parts_pair[1]][e];
                                    cerr << ") ";
                                }
                                cerr << endl;
                            }*/
                            is_hard = true;
                            break;
                        }
                        assert(flow == 0);
                        // still fails on 18 vertices
                    }
                    if (is_hard) {
                        break;
                    }
                }
                if (is_hard) {
                    continue;
                }

                set<TMask> full_cycles;
                int bad_count = 0;
                int maybe_count = 0;
                int good_count = 0;
                for (int coef_pair = 0; coef_pair < 4; ++coef_pair) {
                    vector<int> f;
                    for (int e = 0; e < graph.number_of_edges; ++e) {
                        f.push_back((nz5_from_33pp_coefs[coef_pair][0] * edge_flows[0][e] + nz5_from_33pp_coefs[coef_pair][1] * edge_flows[1][e]) / 2);
                        assert(f[f.size() - 1] != 0);
                    }
                    TMask c = build_full_cycle_from_nz5_flow(graph, f);
                    if (graph.all_full_cycles.find(c) != graph.all_full_cycles.end()) {
                        full_cycles.insert(c);
                        ++good_count;
                    } else {
                        if (c == 0) {
                            ++bad_count;
                        } else {
                            ++maybe_count;
                        }
                    }
                }
                if (full_cycles.size() >= 3) {
                    cerr << three_flow_count << "\t" << "yeah!" << "\t" << full_cycles.size();

                    set<set<TMask>> full_cycles_as_sets;
                    if (full_cycles.size() == 3) {
                        full_cycles_as_sets.insert(full_cycles);
                    } else {
                        for (const auto& c : full_cycles) {
                            set<TMask> full_cycles_without_c;
                            for (const auto& c2 : full_cycles) {
                                if (c2 != c) {
                                    full_cycles_without_c.insert(c2);
                                }
                            }
                            full_cycles_as_sets.insert(full_cycles_without_c);
                        }
                    }

                    for (const auto& full_cycles_set : full_cycles_as_sets) {
                        //FIXME: had to temporarily remove this code
                        /*if (all_6c4c_triples.find(full_cycles_set) != all_6c4c_triples.end()) {
                            cerr << "\tfound!";
                        }
                        if (all_o6c4c_triples.find(full_cycles_set) != all_o6c4c_triples.end()) {
                            cerr << "\torientable!";
                        }*/
                    }
                    cerr << endl;
                } /*else if (three_flow_count == 3) {
                    cerr << full_cycles.size() << " (" << good_count << "," << maybe_count << "," << bad_count << ")" << "\t";
                }*/
            }
        }
    }

    /*if (three_flow_count == 3) {
        cerr << endl;
    }*/

    return true;

    if (three_flow_count == 2) {
        return false;
    }
    all_333pp_cycles.insert(c);
    if (graph.all_even_cycles.find(c) != graph.all_even_cycles.end()) {
        all_333pp_even_cycles.insert(c);
    }
    if (graph.all_full_cycles.find(c) != graph.all_full_cycles.end()) {
        all_333pp_full_cycles.insert(c);
    }
    return true;
}

bool build_33pp_partition(TGraph& graph, const TMask c, int e) {
    if (e == graph.number_of_edges) {
        return check_33pp(graph, c);
    }
    if ((BIT(e) & c) > 0) {
        return build_33pp_partition(graph, c, e + 1);
    }
    bool used[3];
    for (int i = 0; i < 3; ++i) {
        used[i] = false;
    }
    for (int i = 0; i < 2; ++i) { // TODO: rewrite with lower_index_neib_edges structure
        for (int j = 0; j < REG; ++j) {
            int cur_part = u33pp_partition[graph.v2e[graph.e2v[e][i]][j]];
            if (cur_part != -1) {
                used[cur_part] = true;
            }
        }
    }
    for (int part = 0; part < 3; ++part) {
        if (used[part]) {
            continue;
        }
        u33pp_partition[e] = part;
        if (build_33pp_partition(graph, c, e + 1)) {
            return true;
        }
    }
    u33pp_partition[e] = -1;
    return false;
}

void find_some_33pp_for_cycle(TGraph& graph, const TMask c) {
    for (int v = 0; v < graph.number_of_vertices; ++v) {
        vertex_in_33pp_cycle[v] = false;
    }
    for (int e = 0; e < graph.number_of_edges; ++e) {
        if ((BIT(e) & c) > 0) {
            edge_in_33pp_cycle[e] = true;
            for (int i = 0; i < 2; ++i) {
                vertex_in_33pp_cycle[graph.e2v[e][i]] = true;
            }
        } else {
            edge_in_33pp_cycle[e] = false;
        }
    }

    for (int e = 0; e < graph.number_of_edges; ++e) {
        u33pp_partition[e] = -1;
    }
    build_33pp_partition(graph, c, 0);
}

void find_all_33pp(TGraph& graph) {
    all_33pp_cycles.clear();
    all_33pp_circuits.clear();
    all_33pp_full_cycles.clear();
    all_333pp_cycles.clear();
    all_333pp_even_cycles.clear();
    all_333pp_full_cycles.clear();
    for (const auto c : graph.all_cycles) {
        find_some_33pp_for_cycle(graph, c);
    }
}

} // NExpFlowParityPairs
