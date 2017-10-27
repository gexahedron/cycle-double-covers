/*
 * File:   o6c4c.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
 *
 */

#pragma once

#include "graph.h"
#include "common.h"

#include "experiments/mnk_flows.h"

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

namespace NExp6c4c {

using namespace NExpMNKFlows;

using namespace std;

//using namespace NMNKFlows;

set<TMask> cur_6c4c;

vector<TMask> bit_cycles;
int o6c4c_aggregated_solutions;
int all_o6c4c_solutions;
int all_6c4c_solutions;
int same_cycles_different_orientations;
int all_nz_mod5_from_o6c4c;
vector<vector<int>> all_circuits_in_6c4c;
int edge_orientation_count[REG * MAXN / 2][2];
int orientations[6 * MAXN];
int layer[6 * MAXN];
int cur_flow[REG * MAXN / 2];
int layer_flow[6][REG * MAXN / 2];

unordered_map<TMask, int> full_cycle_count;

bool has_dominating_circuit = false;
map<pair<int, int>, int> edge_pairs;
set<set<int>> all_oriented_vertices;
set<int> oriented_vertices;
bool o6c4c_edge_is_poor[MAXN * REG / 2];
TMask poor_mask = 0;


bool o6c4c_always_1[REG * MAXN / 2];
bool o6c4c_always_2[REG * MAXN / 2];
bool o6c4c_always_3[REG * MAXN / 2];
bool o6c4c_always_4[REG * MAXN / 2];
bool o6c4c_never_1[REG * MAXN / 2];
bool o6c4c_never_2[REG * MAXN / 2];
bool o6c4c_never_3[REG * MAXN / 2];
bool o6c4c_never_4[REG * MAXN / 2];

set<pair<int, int>> u6c4c_edge_pair_counts[REG * MAXN / 2];

int u33pp_solutions;
int all_33pp_solutions;

set<TMask> u333pp_cycles_from_o6c4c;

bool has_33pp_from_3pm;
bool has_333pp_from_3pm;

set<vector<int>> common_triples;
set<vector<int>> all_33pp_triples;

bool u6c4c_always_poor[REG * MAXN / 2];
bool u6c4c_always_rich[REG * MAXN / 2];
int u6c4c_min_poor;
int u6c4c_max_poor;

bool has_or_comb = false;
bool has_un_comb = false;

set<set<TMask>> all_6c4c_triples;
set<set<TMask>> all_o6c4c_triples;

/*********************************Methods*********************************/

/*bool find_nz_mod5_from_o6c4c(TGraph& graph, int cur_layer) {
    if (cur_layer == 6) {
        for (int e = 0; e < graph.number_of_edges; ++e) {
            if (cur_flow[e] % 5 == 0) {
                return false;
            }
        }
        ++all_nz_mod5_from_o6c4c;
        //for (int v = 0; v < graph.number_of_vertices; ++v) {
        //    int vertex_flow = 0;
        //    for (int j = 0; j < REG; ++j) {
        //        if (v < graph.v2v[v][j]) {
        //            vertex_flow += cur_flow[graph.v2e[v][j]];
        //        } else {
        //            vertex_flow -= cur_flow[graph.v2e[v][j]];
        //        }
        //    }
        //    cerr << vertex_flow % 5 << " ";
        //}
        //cerr << endl;
        return true;
    }

    for (int w = 0; w < 5; ++w) {
        if (w > 0) {
            for (int e = 0; e < graph.number_of_edges; ++e) {
                cur_flow[e] += layer_flow[cur_layer][e];
            }
        }
        if (find_nz_mod5_from_o6c4c(graph, cur_layer + 1)) {
            return true;
        }
    }
    for (int e = 0; e < graph.number_of_edges; ++e) {
        cur_flow[e] -= layer_flow[cur_layer][e] * 4;
    }
    return false;
}*/

// checks whether there exists dominating circuit, which doesn't go through oriented vertices
bool has_compatible_dominating_circuit(TGraph& graph) {
    /*unordered_set<int> ignored_vertices;
    for (int v = 0; v < graph.number_of_vertices; ++v) {
        bool all_oriented = true;
        for (int j = 0; j < REG; ++j) {
            if (oriented_vertices.find(graph.v2v[v][j]) == oriented_vertices.end()) {
                all_oriented = false;
                break;
            }
        }
        if (all_oriented && oriented_vertices.find(v) != oriented_vertices.end()) {
            cerr << "wut" << endl;
        } else if (all_oriented) {
            ignored_vertices.insert(v);
        }
    }*/
    int count = 0;
    for (const auto& c : graph.all_dominating_circuits) {     
        bool has_all_oriented_vertices = true;
        for (const auto& v : oriented_vertices) {
            bool has_edge = false;
            for (int j = 0; j < REG; ++j) {
                if ((BIT(graph.v2e[v][j]) & c) > 0) {
                    has_edge = true;
                    break;
                }
            }
            if (!has_edge) {
                has_all_oriented_vertices = false;
                break;
            }
        }
        if (!has_all_oriented_vertices) {
            continue;
        }

        /*bool has_ignored_vertex = false;
        for (const auto& v : ignored_vertices) {
            bool has_edge = false;
            for (int j = 0; j < REG; ++j) {
                if ((BIT(graph.v2e[v][j]) & c) > 0) {
                    has_edge = true;
                    break;
                }
            }
            if (has_edge) {
                has_ignored_vertex = true;
                break;
            }
        }
        if (has_ignored_vertex) {
            continue;
        }*/

        ++count;
        has_dominating_circuit = true;
        //return true;
    }
    //cerr << "dominat: " << count << endl;
    return false;
}

bool orient_6c4c(TGraph& graph, int cur_circuit) {
    if (cur_circuit == all_circuits_in_6c4c.size()) {
        ++same_cycles_different_orientations;

        for (int e = 0; e < graph.number_of_edges; ++e) {
            cur_flow[e] = 0;
            for (int i = 0; i < 6; ++i) {
                layer_flow[i][e] = 0;
            }
        }

        edge_pairs.clear();

        for (int c = 0; c < all_circuits_in_6c4c.size(); ++c) {
            for (int vi = 0; vi < all_circuits_in_6c4c[c].size() - 1; ++vi) {
                int v1 = all_circuits_in_6c4c[c][vi];
                int v2 = all_circuits_in_6c4c[c][vi + 1];
                int ei = graph.edge_index[v1][v2];
                int cur_edge_orientation = orientations[c];
                    
                int e1 = -1;
                int e2 = -1;
                int v0 = -1;
                if (vi > 0) {
                    v0 = all_circuits_in_6c4c[c][vi - 1];
                } else {
                    v0 = all_circuits_in_6c4c[c][all_circuits_in_6c4c[c].size() - 2];
                }
                int prev_ei = graph.edge_index[v0][v1];
                if (orientations[c] == 0) {
                    e1 = prev_ei;
                    e2 = ei;
                } else {
                    e1 = ei;
                    e2 = prev_ei;
                }
                ++edge_pairs[make_pair(e1, e2)];

                if (v1 > v2) {
                    cur_edge_orientation = 1 - cur_edge_orientation;
                }
                if (cur_edge_orientation == 0) {
                    layer_flow[layer[c]][ei] += 1;
                } else {
                    layer_flow[layer[c]][ei] -= 1;
                }
            }
        }

        oriented_vertices.clear();
        for (int v = 0; v < graph.number_of_vertices; ++v) {
            int e1 = graph.v2e[v][0];
            int e2 = graph.v2e[v][1];
            if (edge_pairs[make_pair(e1, e2)] != 1) {
                oriented_vertices.insert(v);
            }
        }
        //all_oriented_vertices.insert(oriented_vertices);

        //cerr << "   o6c4c edges: ";
        for (int e = 0; e < graph.number_of_edges; ++e) {
            if (u6c4c_edge_pair_counts[e].size() != 2 && u6c4c_edge_pair_counts[e].size() != 4) {
                cerr << "wut" << endl;
                continue;
            }
            bool v1_is_oriented = (oriented_vertices.find(graph.e2v[e][0]) != oriented_vertices.end());
            bool v2_is_oriented = (oriented_vertices.find(graph.e2v[e][1]) != oriented_vertices.end());
            int oriented_count = 0;
            if (v1_is_oriented) {
                ++oriented_count;
            }
            if (v2_is_oriented) {
                ++oriented_count;
            }
            if (u6c4c_edge_pair_counts[e].size() == 2) {
                o6c4c_always_3[e] = false;
                o6c4c_always_4[e] = false;
                o6c4c_edge_is_poor[e] = true;
                if (oriented_count == 2) {
                    o6c4c_always_2[e] = false;
                    o6c4c_never_1[e] = false;
                    //cerr << "1 ";
                } else if (oriented_count == 0) {
                    o6c4c_always_1[e] = false;
                    o6c4c_never_2[e] = false;
                    //cerr << "2 ";
                } else {
                    cerr << "wat1" << endl;
                }
            } else {
                o6c4c_edge_is_poor[e] = false;
                o6c4c_always_1[e] = false;
                o6c4c_always_2[e] = false;
                if (oriented_count == 1) {
                    o6c4c_always_4[e] = false;
                    o6c4c_never_3[e] = false;
                    //cerr << "3 ";
                } else if (oriented_count == 0) {
                    o6c4c_always_3[e] = false;
                    o6c4c_never_4[e] = false;
                    //cerr << "4 ";
                } else {
                    cerr << "wat2" << endl;
                }
            }
        }

        //cerr << endl;

        has_compatible_dominating_circuit(graph);
        /*if (!has_dominating_circuit) {
            cerr << "wow" << endl;
        }
        has_dominating_circuit = false;*/

        return false;
        //find_nz_mod5_from_o6c4c(graph, 1); // skip first layer - for canonicity - it will always have weight 0
        return true;
        return false;//true;
    }
    int max_orientation = 1;
    if (cur_circuit == 0) {
        max_orientation = 0;
    }
    for (int orientation = 0; orientation <= max_orientation; ++orientation) {
        orientations[cur_circuit] = orientation;
        int vi = 0;
        while (vi < all_circuits_in_6c4c[cur_circuit].size() - 1) {
            int v1 = all_circuits_in_6c4c[cur_circuit][vi];
            int v2 = all_circuits_in_6c4c[cur_circuit][vi + 1];
            int ei = graph.edge_index[v1][v2];
            int cur_edge_orientation = orientation;
            if (v1 > v2) {
                cur_edge_orientation = 1 - cur_edge_orientation;
            }
            if (edge_orientation_count[ei][cur_edge_orientation] == 2) {
                break;
            }
            ++edge_orientation_count[ei][cur_edge_orientation];
            ++vi;
        }
        if (vi == all_circuits_in_6c4c[cur_circuit].size() - 1) {
            if (orient_6c4c(graph, cur_circuit + 1)) {
                return true;
            }
        }
        --vi;
        while (vi >= 0) {
            int v1 = all_circuits_in_6c4c[cur_circuit][vi];
            int v2 = all_circuits_in_6c4c[cur_circuit][vi + 1];
            int ei = graph.edge_index[v1][v2];
            int cur_edge_orientation = 0;
            if (v1 > v2) {
                cur_edge_orientation = 1 - cur_edge_orientation;
            }
            if (orientation == 1) {
                cur_edge_orientation = 1 - cur_edge_orientation;
            }
            --edge_orientation_count[ei][cur_edge_orientation];
            --vi;
        }
    }
    return false;
}

bool check_orientability_6c4c(TGraph& graph) {
    has_dominating_circuit = false;
    /*cerr << "check orientability for cycles ";
    for (int i = 0; i < 6; ++i) {
        cerr << u6c4c_cycles[i] << " ";
    }
    cerr << endl;*/
    vector<vector<int>> circuits;
    all_circuits_in_6c4c.clear();
    for (int i = 0; i < 6; ++i) {
        circuits = graph.cycles_as_circuits[u6c4c_cycles[i]];
        for (const auto& circuit : circuits) {
            layer[all_circuits_in_6c4c.size()] = i;
            all_circuits_in_6c4c.push_back(circuit);
        }
    }
    //cerr << "all circuits: " << all_circuits_in_6c4c.size() << endl;
    for (int e = 0; e < graph.number_of_edges; ++e) {
        for (int orientation = 0; orientation < 2; ++orientation) {
            edge_orientation_count[e][orientation] = 0;
        }
    }

    for (int e = 0; e < graph.number_of_edges; ++e) {
        u6c4c_edge_pair_counts[e].clear();
    }

    for (int c = 0; c < all_circuits_in_6c4c.size(); ++c) {
        for (int vi = 0; vi < all_circuits_in_6c4c[c].size() - 1; ++vi) {
            int v0 = -1;
            if (vi > 0) {
                v0 = all_circuits_in_6c4c[c][vi - 1];
            } else {
                v0 = all_circuits_in_6c4c[c][all_circuits_in_6c4c[c].size() - 2];
            }
            int v1 = all_circuits_in_6c4c[c][vi];
            int v2 = all_circuits_in_6c4c[c][vi + 1];
            int v3 = -1;
            if (vi < all_circuits_in_6c4c[c].size() - 2) {
                v3 = all_circuits_in_6c4c[c][vi + 2];
            } else {
                v3 = all_circuits_in_6c4c[c][1];
            }
            if (v0 > v3) {
                swap(v0, v3);
            }
            int ei = graph.edge_index[v1][v2];
            u6c4c_edge_pair_counts[ei].insert(make_pair(v0, v3));
        }
    }

    same_cycles_different_orientations = 0;
    orient_6c4c(graph, 0);
    if (same_cycles_different_orientations > 0) {
        ++o6c4c_aggregated_solutions;
        all_o6c4c_solutions += same_cycles_different_orientations;
        for (int i = 0; i < 6; ++i) {
            ++full_cycle_count[u6c4c_cycles[i]];
        }
        graph.all_o6c4c.insert(cur_6c4c);
    }
    /*if (same_cycles_different_orientations > 1) {
        cerr << "found " << same_cycles_different_orientations << " same cycles different orientations for o6c4c" << endl;
    }*/
    return false;
}

bool find_33pp_from_3pm(TGraph& graph) { // TODO: rewrite using check_33pp
    bool vertex_in_33pp_cycle[MAXN];
    bool edge_in_33pp_cycle[REG * MAXN / 2];
    for (int v = 0; v < graph.number_of_vertices; ++v) {
        vertex_in_33pp_cycle[v] = false;
    }
    TMask cycle_mask = 0;
    for (int e = 0; e < graph.number_of_edges; ++e) {
        int edge_count = 0;
        for (int i = 0; i < 3; ++i) {
            if ((BIT(e) & u3_inv_pm[i]) > 0) {
                ++edge_count;
            }
        }
        if (edge_count == 1 || edge_count == 3) {
            edge_in_33pp_cycle[e] = true;
            cycle_mask += BIT(e);
            for (int i = 0; i < 2; ++i) {
                vertex_in_33pp_cycle[graph.e2v[e][i]] = true;
            }
        } else {
            edge_in_33pp_cycle[e] = false;
        }
    }

    /*if (u333pp_cycles_from_o5cdc.find(cycle_mask) == u333pp_cycles_from_o5cdc.end()) { // has exceptions - 28.05g1422
        return false;
    }*/

    /*if ((poor_mask & cycle_mask) > 0) { // has exceptions - 22.05g3, g8
        return false;
    }*/

    /*bool has_some_oriented_vertices = false; // has exceptions - 10.05g1
    for (const auto& vs : all_oriented_vertices) {
        bool has_this = true;
        for (const auto& v : vs) {
            if (!vertex_in_33pp_cycle[v]) {
                has_this = false;
                break;
            }
        }
        if (has_this) {
            has_some_oriented_vertices = true;
            break;
        }
    }
    if (!has_some_oriented_vertices) {
        return false;
    }*/

    int three_flow_count = 0;
    TMask part_mask[3];
    bool has_3flow[3];
    for (int part = 0; part < 3; ++part) {
        part_mask[part] = 0;
        bool edge_in_cur_part[REG * MAXN / 2];
        bool vertex_in_cur_part[MAXN];
        int vertex_colour[MAXN];
        for (int v = 0; v < graph.number_of_vertices; ++v) {
            vertex_colour[v] = -1;
            vertex_in_cur_part[v] = false;
        }
        for (int e = 0; e < graph.number_of_edges; ++e) {
            if (edge_in_33pp_cycle[e]) {
                edge_in_cur_part[e] = true;
            } else {
                edge_in_cur_part[e] = ((BIT(e) & u3_inv_pm[part]) > 0);

                if (edge_in_cur_part[e]) {
                    for (int i = 0; i < 2; ++i) {
                       if (vertex_in_33pp_cycle[graph.e2v[e][i]]) {
                            vertex_in_cur_part[graph.e2v[e][i]] = true;
                        }
                    }
                }
            }
            if (edge_in_cur_part[e]) {
                part_mask[part] += BIT(e);
            }
        }

        bool edge_visited[REG * MAXN / 2];
        for (int e = 0; e < graph.number_of_edges; ++e) {
            edge_visited[e] = !edge_in_cur_part[e];
        }
        int queue[MAXN];
        int queue_size;
        has_3flow[part] = true;
        for (int v = 0; v < graph.number_of_vertices; ++v) {
            if (vertex_in_cur_part[v] && vertex_colour[v] == -1) {
                vertex_colour[v] = 0;
                queue_size = 1;
                queue[0] = v;
                for (int cur_idx = 0; cur_idx < queue_size; ++cur_idx) {
                    int cur_vertex = queue[cur_idx];
                    for (int j = 0; j < REG; ++j) {
                        if (edge_visited[graph.v2e[cur_vertex][j]]) {
                            continue;
                        }
                        edge_visited[graph.v2e[cur_vertex][j]] = true;
                        int v1 = cur_vertex;
                        int v2 = graph.v2v[cur_vertex][j];
                        while (!vertex_in_cur_part[v2]) {
                            for (int j2 = 0; j2 < REG; ++j2) {
                                int e = graph.v2e[v2][j2];
                                int v3 = graph.v2v[v2][j2];
                                if (!edge_visited[e] && v3 != v1) {
                                    v1 = v2;
                                    v2 = v3;
                                    edge_visited[e] = true;
                                    break;
                                }
                            }
                        }
                        if (vertex_colour[v2] == -1) {
                            vertex_colour[v2] = 1 - vertex_colour[cur_vertex];
                            queue[queue_size] = v2;
                            ++queue_size;
                        } else if (vertex_colour[v2] == vertex_colour[cur_vertex]) {
                            has_3flow[part] = false;
                            break;
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
        }
    }
    if (three_flow_count < 2) {
        return false;
    }

    has_33pp_from_3pm = true;
    if (three_flow_count == 3) {
        u333pp_cycles_from_o6c4c.insert(cycle_mask);
    }

    for (int part = 0; part < 3; ++part) {
        int i = 0;
        int j = 2;
        if (part == 0) {
            i = 1;
        }
        if (part == 2) {
            j = 1;
        }
        if (!has_3flow[i] || !has_3flow[j])
            continue;
        TMask m1 = part_mask[i];
        TMask m2 = part_mask[j];
        if (m1 > m2) {
            swap(m1, m2);
        }
        for (const auto& u5cdc : u5cdc_from_33pp[make_pair(m1, m2)]) {
            graph.u244_6c4c_5cdc_pairs[cur_6c4c].insert(u5cdc);
            graph.u244_5cdc_6c4c_pairs[u5cdc].insert(cur_6c4c);
            /*cerr << "6c4c-5cdc:\t";
            for (const auto& c : cur_6c4c) {
                cerr << c << " ";
            }
            cerr << "\t-\t";
            for (const auto& c : u5cdc) {
                cerr << c << " ";
            }
            cerr << endl;*/
        }
        /*if (u33pp_pairs.find(make_pair(m1, m2)) == u33pp_pairs.end()) {
            print_graph(graph);
            for (int e = 0; e < graph.number_of_edges; ++e) {
                cerr << e;
                if (edge_in_33pp_cycle[e]) {
                    cerr << ": in cycle";
                } else {
                    for (int part = 0; part < 3; ++part) {
                        if ((BIT(e) & u3_inv_pm[part]) == 0) {
                            cerr << ": p" << part;
                        }
                    }
                }
                cerr << endl;
            }

            //u33pp_solutions += 1;
            //return true;
            break;
        }*/
    }
    return false;
    //return false;
    /*
    if (three_flow_count == 2) {
        u33pp_solutions += 1;
    } else if (three_flow_count == 3) {
        u33pp_solutions += 3;
    }
    return false;*/
}

bool find_33pp_from_6c4c(TGraph& graph) {
    has_33pp_from_3pm = false;
    has_333pp_from_3pm = false;
    u33pp_solutions = 0;
    int min_edge_count = graph.number_of_edges;
    set<vector<int>> triples;
    for (int i = 0; i < 6; ++i) {
        u3_inv_pm[0] = u6c4c_cycles[i];
        for (int j = i + 1; j < 6; ++j) {
            u3_inv_pm[1] = u6c4c_cycles[j];
            for (int k = j + 1; k < 6; ++k) {
                u3_inv_pm[2] = u6c4c_cycles[k];
                int edge_count = 0;
                for (int e = 0; e < graph.number_of_edges; ++e) {
                    if ((BIT(e) & u3_inv_pm[0] & u3_inv_pm[1] & u3_inv_pm[2]) > 0) {
                        ++edge_count;
                    }
                }
                if (edge_count < min_edge_count) {
                    triples.clear();
                    min_edge_count = edge_count;
                }
                //if (edge_count == min_edge_count) {
                    vector<int> anti_triple;
                    for (int m = 0; m < 6; ++m) {
                        if (m != i && m != j && m != k) {
                            anti_triple.push_back(m);
                        }
                    }
                    vector<int> triple = {i, j, k};
                    triples.insert(anti_triple);
                    triples.insert(triple);
                //}
            }
        }
    }
    for (const auto& triple: triples) {
        u3_inv_pm[0] = u6c4c_cycles[triple[0]];
        u3_inv_pm[1] = u6c4c_cycles[triple[1]];
        u3_inv_pm[2] = u6c4c_cycles[triple[2]];
        find_33pp_from_3pm(graph);
        if (has_33pp_from_3pm) {
            has_33pp_from_3pm = false;
            all_33pp_triples.insert(triple);
            common_triples.insert(triple);
        }
        if (u33pp_solutions > 0) {
            return true;
        }
    }
    return false;
}

/*void analyze_edges_of_o6c4c_solution(TGraph& graph) {
    set<pair<int, int>> counts[REG * MAXN / 2];
    for (int e = 0; e < graph.number_of_edges; ++e) {
        counts[e].clear();
    }

    for (int c = 0; c < all_circuits_in_6c4c.size(); ++c) {
        for (int vi = 0; vi < all_circuits_in_6c4c[c].size() - 1; ++vi) {
            int v0 = -1;
            if (vi > 0) {
                v0 = all_circuits_in_6c4c[c][vi - 1];
            } else {
                v0 = all_circuits_in_6c4c[c][all_circuits_in_6c4c[c].size() - 2];
            }
            int v1 = all_circuits_in_6c4c[c][vi];
            int v2 = all_circuits_in_6c4c[c][vi + 1];
            int v3 = -1;
            if (vi < all_circuits_in_6c4c[c].size() - 2) {
                v3 = all_circuits_in_6c4c[c][vi + 2];
            } else {
                v3 = all_circuits_in_6c4c[c][1];
            }
            if (v0 > v3) {
                swap(v0, v3);
            }
            int ei = graph.edge_index[v1][v2];
            counts[ei].insert(make_pair(v0, v3));
        }
    }
    int poor_count = 0;
    for (int e = 0; e < graph.number_of_edges; ++e) {
        if (counts[e].size() != 2 && counts[e].size() != 4) {
            cerr << "wut" << endl;
            continue;
        }
        if (counts[e].size() == 2) {
            ++poor_count;
            u6c4c_always_rich[e] = false;
        } else {
            u6c4c_always_poor[e] = false;
        }
    }
    u6c4c_min_poor = min(u6c4c_min_poor, poor_count);
    u6c4c_max_poor = max(u6c4c_max_poor, poor_count);
}*/

bool gen_o6c4c(TGraph& graph, int cur_cycle_layer, int min_cycle_idx, bool only_find) {
    if (cur_cycle_layer == 6) {
        int edge_count[REG * MAXN / 2];
        for (int e = 0; e < graph.number_of_edges; ++e) {
            edge_count[e] = 0;
        }
        for (int c = 0; c < 6; ++c) {
            for (int e = 0; e < graph.number_of_edges; ++e) {
                if (BIT(e) & u6c4c_cycles[c]) {
                    ++edge_count[e];
                }
            }
        }
        for (int e = 0; e < graph.number_of_edges; ++e) {
            if (edge_count[e] != 4) {
                return false;
            }
        }

        ++all_6c4c_solutions;
        cur_6c4c.clear();
        for (int c = 0; c < 6; ++c) {
            cur_6c4c.insert(u6c4c_cycles[c]);
        }

        graph.all_6c4c.insert(cur_6c4c);

        if (only_find) {
            return false;
        }
        if (check_orientability_6c4c(graph)) {
            cerr << "found o6c4c" << endl;
            return true;
        }
        bool has_o6c4c = (same_cycles_different_orientations > 0);
        //return true;

        for (int i = 0; i < 6; ++i) {
            for (int j = i + 1; j < 6; ++j) {
                for (int k = j + 1; k < 6; ++k) {
                    set<TMask> triple;
                    triple.insert(u6c4c_cycles[i]);
                    triple.insert(u6c4c_cycles[j]);
                    triple.insert(u6c4c_cycles[k]);
                    all_6c4c_triples.insert(triple);
                    if (has_o6c4c) {
                        all_o6c4c_triples.insert(triple);
                    }
                }
            }
        }
        return false;

        all_33pp_triples.clear();
        find_33pp_from_6c4c(graph);
        has_33pp_from_3pm = all_33pp_triples.size() > 0;
        if (!has_33pp_from_3pm) {
            return false;
        }

        //all_oriented_vertices.clear();
        if (!has_o6c4c) {
            return false;
        }

        poor_mask = 0;
        for (int e = 0; e < graph.number_of_edges; ++e) {
            if (o6c4c_edge_is_poor[e]) {
                poor_mask += BIT(e);
            }
        }
        //return false;

        o244_triples.clear();
        find_o244_flows_from_6c4c(all_33pp_triples, graph);
        has_o244_flows = o244_triples.size() > 0;
        if (!has_o244_flows) {
            return false;
        }

        bool have_same_triple = false;
        common_triples.clear();
        for (const auto& triple : all_33pp_triples) {
            if (o244_triples.find(triple) != o244_triples.end()) {
                have_same_triple = true;
                common_triples.insert(triple);
            }
        }
        if (!have_same_triple) {
            return false;
        }
        //has_o244_flows = true;
        //bool have_same_triple = true;

        find_333flows_from_6c4c(common_triples, graph);
        if (!has_all_3flows) {
            return false;
        }

        if (has_o6c4c && has_33pp_from_3pm && has_all_3flows && has_o244_flows && have_same_triple) {// && has_dominating_circuit) {
            has_or_comb = true;
        }
        if (has_or_comb) {
            cerr << "has oriented combination" << endl;
            return true;
        }
        /*
        if (has_33pp_from_3pm && has_333pp_from_3pm && has_all_3flows && has_o244_flows && have_same_triple) {
            has_un_comb = true;
        }*/

        /*if (has_or_comb && has_un_comb) {
            cerr << "has both" << endl;
            return true;
        }*/

        //if ((same_cycles_different_orientations > 0) && has_33pp_from_3pm && has_all_3flows) {
        //cerr << "6c4c:\t" << (same_cycles_different_orientations > 0) << "\t" << has_33pp_from_3pm << "\t" << has_333pp_from_3pm << "\t" << has_all_3flows << "\t" << has_o244_flows << "\t" << have_same_triple << endl;
            //return true;
        //}

        if (same_cycles_different_orientations > 0) {
            all_33pp_solutions += u33pp_solutions;
            if (all_33pp_solutions > 0)
                return true;
        }
        return false;
    }
    for (int i = min_cycle_idx; i < bit_cycles.size(); ++i) {
        u6c4c_cycles[cur_cycle_layer] = bit_cycles[i];
        bool compat = true;
        for (int c1 = 0; c1 < cur_cycle_layer; ++c1) {
            for (int c2 = c1 + 1; c2 < cur_cycle_layer; ++c2) {
                if (inv(graph, u6c4c_cycles[c1]) & inv(graph, u6c4c_cycles[c2]) & inv(graph, u6c4c_cycles[cur_cycle_layer])) {
                    compat = false;
                    break;
                }
            }
            if (!compat) {
                break;
            }
        }
        if (compat) {
            if (gen_o6c4c(graph, cur_cycle_layer + 1, i + 1, only_find)) {
                return true;
            }
        }
    }
    return false;
}

/*
void find_o6c4c_compatible_with_preimages(TGraph& graph) {
    cerr << "overall full cycles: " << graph.all_full_cycles.size() << endl;
    cerr << "full cycles from petersen: " << full_cycles_from_petersen.size() << endl;
    bit_cycles.clear();
    for (const auto c : graph.all_full_cycles) {
        if (full_cycles_from_petersen.find(c) != full_cycles_from_petersen.end()) {
            bit_cycles.push_back(c);
        }
    }
    cerr << "left: " << bit_cycles.size() << endl;
    if (!gen_o6c4c(graph, 0, 0)) {
        cerr << "didn't find o6c4c" << endl;
    }
}
*/

void find_all_o6c4c(TGraph& graph, bool only_find = false) {
    all_6c4c_triples.clear();
    all_o6c4c_triples.clear();
    o6c4c_aggregated_solutions = 0;
    all_o6c4c_solutions = 0;
    all_6c4c_solutions = 0;
    all_nz_mod5_from_o6c4c = 0;
    all_33pp_solutions = 0;
    full_cycle_count.clear();
    bit_cycles.clear();
    u333pp_cycles_from_o6c4c.clear();
    for (const auto c : graph.all_full_cycles) {
        bit_cycles.push_back(c);
        full_cycle_count[c] = 0;
    }
    has_or_comb = false;
    has_un_comb = false;
    for (int e = 0; e < graph.number_of_edges; ++e) {
        u6c4c_always_poor[e] = true;
        u6c4c_always_rich[e] = true;

        o6c4c_always_1[e] = true;
        o6c4c_always_2[e] = true;
        o6c4c_always_3[e] = true;
        o6c4c_always_4[e] = true;

        o6c4c_never_1[e] = true;
        o6c4c_never_2[e] = true;
        o6c4c_never_3[e] = true;
        o6c4c_never_4[e] = true;
    }
    u6c4c_min_poor = graph.number_of_edges;
    u6c4c_max_poor = 0;

    gen_o6c4c(graph, 0, 0, only_find);
    if (!only_find && !has_or_comb) {
        cerr << "no oriented solution found!" << endl;
    }

    //for (const auto c : full_cycle_count) {
    //    cerr << c.first << "\t" << c.second * 1.0 / o6c4c_aggregated_solutions << endl;
    //}
    //cerr << "ratio: " << graph.petersen_colourings.size() * 1.0 / o6c4c_aggregated_solutions << endl;
    //cerr << "ratio2: " << graph.profiles.size() * 1.0 / o6c4c_aggregated_solutions << endl;
    //cerr << "ratio3: " << o6c4c_aggregated_solutions * 1.0 / all_o6c4c_solutions << endl;

    //cerr << "all 6c4c solutions: " << all_6c4c_solutions << endl;
    /*cerr << "all aggregated o6c4c solutions: " << o6c4c_aggregated_solutions << endl;
    cerr << "all o6c4c solutions: " << all_o6c4c_solutions << endl;
    cerr << "all nz-mod5 from o6c4c combinations: " << all_nz_mod5_from_o6c4c << endl;
    cerr << "all 33pp solutions from 6c4c: " << all_33pp_solutions << endl;
    cerr << endl;*/
}

} // NExp6c4c
