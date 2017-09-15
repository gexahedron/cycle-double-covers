/*
 * File:   relations.cpp
 * Author: Nikolay Ulyanov
 *
 * Created on 14 August 2017
 * Here I will study relations between normal (Petersen colouring), o6c4c, o5cdc, nz5 (nz-mod5, 33-pp)
 * Milestones:
 * added TGraph struct
 * added fast Petersen colouring
 * added cc-mapping
 * added cycles
 * added Hoffman-Ostenhof decomposition
 * added check for preimage Hoffman-Ostenhof solutions (fail)
 * added new code for o6c4c
 * added new code for o5cdc
 * added new code for a new construction, which includes 6c4c (244-flows) and 33-pp (which includes nz5 and 5cdc)
 * added code for 33-pp and 333-pp flows
 * added counts for circuits (in the list of cycles, full cycles, even cycles)
 * compared 33-pp and 5cdc
 * added code for 333-flows from 6c4c
 * added code for o244-flows from 6c4c
 * added check for same triple 33-pp and o244-flows
 * added code for dominating circuits
 * added code for poor-rich edges from 6c4c
 * added code for oriented vertices from o6c4c
 * added code for poor-rich edges from o6c4c (4 different types of edges)
 * compared o6c4c and petersen colouring
 * added code for (6c4c, 5cdc) pairs from 6c4c-244-flows-33-pp construction
 * added code for (6c4c, 5cdc) pairs from petersen colouring
 *
 * TODO:
 * rewrite 5cdc so that i won't have duplicates
 * combine code for 5cdc and 6c4c and add support for 9c6c
 * add graph6 format support
 * ...
 */

#include "util.h"

#include <unordered_set>
#include <unordered_map>
#include <set>
#include <map>
#include <climits>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>

using namespace std;

TMask number_of_graphs_read = 0;

const bool PRINT_SOLUTIONS = false;

TGraph petersen_graph;

// additional structures
// petersen colouring
bool edge_is_poor[REG * MAXN / 2];
bool vertex_is_poor[10]; // vertex is poor when it is connected to itself
int edge_colour[REG * MAXN / 2];
int edge_order[REG * MAXN / 2];
int width[REG * MAXN / 2];
unordered_set<int> visited_edges;
int not_coloured_edges_count_near_vertex[MAXN];
bool had_four = false;

// cc-mapping
unordered_map<int, int> ve_to_petersen_edge;
unordered_map<int, int> mask_to_petersen_vertex;

unordered_set<TMask> full_cycles_from_petersen;
unordered_map<TMask, int> full_cycle_count;

set<set<TMask>> all_5cdc;
set<TMask> cur_5cdc;
map<pair<TMask, TMask>, set<set<TMask>>> u5cdc_from_33pp;

set<TMask> cur_6c4c;

map<set<TMask>, set<set<TMask>>> u244_6c4c_5cdc_pairs;
map<set<TMask>, set<set<TMask>>> petersen_6c4c_5cdc_pairs;

map<set<TMask>, set<set<TMask>>> u244_5cdc_6c4c_pairs;
map<set<TMask>, set<set<TMask>>> petersen_5cdc_6c4c_pairs;

/*********************************Methods**************************************/

/******************petersen colouring**********************/
int combine_hash_from_mask_and_edge_colour(int mask, int edge_colour) {
    return mask + edge_colour * BIT(6);
}

bool check_edge(int e, const TGraph& graph) {
    int v1 = graph.e2v[e][0];
    int v2 = graph.e2v[e][1];
    if (not_coloured_edges_count_near_vertex[v1] != 0 || not_coloured_edges_count_near_vertex[v2] != 0) {
        return true;
    }
    unordered_set<int> colours;
    for (int j = 0; j < REG; ++j) {
        int e2 = graph.v2e[v1][j];
        colours.insert(edge_colour[e2]);
    }
    for (int j = 0; j < REG; ++j) {
        int e2 = graph.v2e[v2][j];
        colours.insert(edge_colour[e2]);
    }
    if (colours.size() == 3) {
        edge_is_poor[e] = true;
        return true;
    } else if (colours.size() == 5) {
        edge_is_poor[e] = false;
        return true;
    }
    return false;
}

bool petersen_always_poor[REG * MAXN / 2];
bool petersen_always_rich[REG * MAXN / 2];
int petersen_min_poor;
int petersen_max_poor;

bool normal_colour_edges(int cur_edge_idx, TGraph& graph) {
    if (cur_edge_idx == graph.number_of_edges) {
        vector<int> colours;
        for (int e = 0; e < graph.number_of_edges; ++e) {
            colours.push_back(edge_colour[e]);
        }
        graph.normal5_colourings.push_back(colours);
        string prof;
        int poor_count = 0;
        unordered_map<int, pair<int, int>> poor_masks;
        vector<int> petersen_edges;
        for (int v = 0; v < 10; ++v)
            vertex_is_poor[v] = false;
        //cerr << "petersen edges: ";
        for (int e = 0; e < graph.number_of_edges; ++e) {
            if (edge_is_poor[e]) {
                ++poor_count;
                //cerr << "POOR: ";
                //cerr << "0 ";
                prof += "p";
                petersen_always_rich[e] = false;
            } else {
                //cerr << "rich: ";
                //cerr << "1 ";
                prof += "r";
                petersen_always_poor[e] = false;
            }

            //cerr << graph.e2v[e][0] << " (";
            int mask1 = 0;
            for (int j = 0; j < REG; ++j) {
                //cerr << edge_colour[graph.v2e[graph.e2v[e][0]][j]];
                mask1 += BIT(edge_colour[graph.v2e[graph.e2v[e][0]][j]]);
            }
            //cerr << ") -" << edge_colour[e] << "> " << graph.e2v[e][1] << " (";
            int mask2 = 0;
            for (int j = 0; j < REG; ++j) {
                //cerr << edge_colour[graph.v2e[graph.e2v[e][1]][j]];
                mask2 += BIT(edge_colour[graph.v2e[graph.e2v[e][1]][j]]);
            }
            if (mask1 == mask2) {
                vertex_is_poor[mask_to_petersen_vertex[mask1]] = true;
            }
            //cerr << "); ";
            int petersen_edge = ve_to_petersen_edge[combine_hash_from_mask_and_edge_colour(mask1, edge_colour[e])];
            if (edge_is_poor[e]) {
                poor_masks[petersen_edge] = make_pair(mask1, edge_colour[e]);
            }
            petersen_edges.push_back(petersen_edge);
        }
        //cerr << endl;
        petersen_min_poor = min(petersen_min_poor, poor_count);
        petersen_max_poor = max(petersen_max_poor, poor_count);

        bool dominating = false;
        for (int v = 0; v < 10; ++v) {
            if (!vertex_is_poor[v]) {
                //cerr << "can dominate with " << v << endl;
                dominating = true;
                break;
            }
        }
        if (!dominating) {
            //cerr << "no domination!" << endl;
        }
        graph.petersen_colourings.push_back(petersen_edges);
        ++graph.profiles[prof];
        //cerr << endl;
        //cerr << "#poor: " << poor_count << "; " << "unique poor: " << poor_masks.size() << endl;
        for (const auto& pm_pair : poor_masks) {
            //cerr << "(";
            for (int i = 1; i <= 5; ++i) {
                if (BIT(i) & pm_pair.second.first) {
                    //cerr << i;
                }
            }
            //cerr << " " << pm_pair.second.second << "), ";
        }
        //cerr << endl;
        //cerr << endl;
        return false;
    }

    int cur_edge = edge_order[cur_edge_idx];
    bool had_colour = (edge_colour[cur_edge] != 0);
    int low_colour = 1;
    int high_colour = 5;
    if (had_colour) {
        low_colour = high_colour = edge_colour[cur_edge];
    }
    int v1 = graph.e2v[cur_edge][0];
    int v2 = graph.e2v[cur_edge][1];

    // initialization
    if (!had_colour) {
        --not_coloured_edges_count_near_vertex[v1];
        --not_coloured_edges_count_near_vertex[v2];
    }
    for (int i = low_colour; i <= high_colour; ++i) {
        // more initialization
        bool had_four_backup = had_four;
        if (i == 4)
            had_four = true;

        edge_colour[cur_edge] = i;

        // conditions
        bool checks_passed = true;
        for (int j = 0; j < REG; ++j) {
            int ei = graph.v2e[v1][j];
            if ((ei != cur_edge && edge_colour[ei] == edge_colour[cur_edge]) || (not_coloured_edges_count_near_vertex[v1] == 0 && !check_edge(ei, graph))) {
                checks_passed = false;
                break;
            }
        }
        if (checks_passed) {
            for (int j = 0; j < REG; ++j) {
                int ei = graph.v2e[v2][j];
                if (ei != cur_edge && (edge_colour[ei] == edge_colour[cur_edge] || (not_coloured_edges_count_near_vertex[v2] == 0 && !check_edge(ei, graph)))) {
                    checks_passed = false;
                    break;
                }
            }
        }

        // recursion
        if (checks_passed && normal_colour_edges(cur_edge_idx + 1, graph)) {
            return true;
        }
        
        // undo
        had_four = had_four_backup;
        if (!had_four_backup && i == 4) {
            break;
        }
    }
    // more undo
    if (!had_colour) {
        ++not_coloured_edges_count_near_vertex[v1];
        ++not_coloured_edges_count_near_vertex[v2];
        edge_colour[cur_edge] = 0;
    }
    return false;
}

int find_edge_order(int start_edge, const TGraph& graph) {
    visited_edges.clear();
    edge_order[0] = start_edge;
    width[0] = 0;
    visited_edges.insert(0);
    for (int i = 0; i < graph.number_of_edges; ++i) {
        int v1 = graph.e2v[edge_order[i]][0];
        int v2 = graph.e2v[edge_order[i]][1];
        if (v1 > v2) {
            swap(v1, v2);
        }
        for (int j = 0; j < REG; ++j) {
            int ei = graph.v2e[v1][j];
            if (visited_edges.find(ei) == visited_edges.end()) {
                edge_order[visited_edges.size()] = ei;
                width[visited_edges.size()] = width[i] + 1;
                visited_edges.insert(ei);
            }
        }
        for (int j = 0; j < REG; ++j) {
            int ei = graph.v2e[v2][j];
            if (visited_edges.find(ei) == visited_edges.end()) {
                edge_order[visited_edges.size()] = ei;
                width[visited_edges.size()] = width[i] + 1;
                visited_edges.insert(ei);
            }
        }
    }
    return width[graph.number_of_edges - 1];
}

void find_all_petersen_colourings(TGraph& graph) {
    int best_edge_ans = 0;
    int best_edge_idx = 0;
    for (int ei = 0; ei < graph.number_of_edges; ++ei) {
        int cur_ans = find_edge_order(ei, graph);
        if (cur_ans > best_edge_ans) {
            best_edge_ans = cur_ans;
            best_edge_idx = ei;
        }
    }
    find_edge_order(best_edge_idx, graph);

    graph.profiles.clear();

    for (int e = 0; e < graph.number_of_edges; ++e) {
        edge_colour[e] = 0;
        edge_is_poor[e] = false;
        petersen_always_poor[e] = true;
        petersen_always_rich[e] = true;
    }
    petersen_min_poor = graph.number_of_edges;
    petersen_max_poor = 0;

    for (int v = 0; v < graph.number_of_vertices; ++v) {
        not_coloured_edges_count_near_vertex[v] = REG;
    }

    // colour edges around vertex 0
    for (int j = 0; j < REG; ++j) {
        edge_colour[graph.v2e[0][j]] = j + 1;
        --not_coloured_edges_count_near_vertex[0];
        --not_coloured_edges_count_near_vertex[graph.v2v[0][j]];
    }
    had_four = false;
    normal_colour_edges(0, graph);
    cerr << "number of petersen colourings: " << graph.petersen_colourings.size() << "; number of profiles: " << graph.profiles.size() << endl;
}

/******************cc-mapping**********************/
void create_cc_mapping() {
    for (int v = 0; v < petersen_graph.number_of_vertices; ++v) {
        int mask = 0;
        for (int j = 0; j < REG; ++j) {
            int e = petersen_graph.v2e[v][j];
            mask += BIT(petersen_graph.normal5_colourings[0][e]);
        }
        mask_to_petersen_vertex[mask] = v;
        for (int j = 0; j < REG; ++j) {
            int e = petersen_graph.v2e[v][j];
            ve_to_petersen_edge[combine_hash_from_mask_and_edge_colour(mask, petersen_graph.normal5_colourings[0][e])] = e;
        }
    }
    cerr << "number of mask/edge_colour pairs: " << ve_to_petersen_edge.size() << endl;
}

/**********************hoffman-ostenhof decomposition************************/
bool vertex_in_cycle[MAXN];
bool edge_in_cycle[REG * MAXN / 2];
int cycle_length;

void check_tree_in_hoffman_ostenhof(TGraph& graph) {
    if (graph.number_of_vertices != 10) {
        return;
    }
    bool vertex_in_tree[MAXN];
    for (int v = 0; v < graph.number_of_vertices; ++v) {
        vertex_in_tree[v] = false;
    }

    vertex_in_tree[0] = true;
    int queue_size = 1;
    int queue[MAXN];
    queue[0] = 0;

    for (int cur_idx = 0; cur_idx < queue_size; ++cur_idx) {
        int cur_vertex = queue[cur_idx];
        for (int j = 0; j < REG; ++j) {
            int next_vertex = graph.v2v[cur_vertex][j];
            if (!edge_in_cycle[graph.v2e[cur_vertex][j]] && !vertex_in_tree[next_vertex]) {
                vertex_in_tree[next_vertex] = true;
                queue[queue_size] = next_vertex;
                ++queue_size;
            }
        }
    }
    for (int v = 0; v < graph.number_of_vertices; ++v) {
        if (!vertex_in_tree[v]) {
            return;
        }
    }

    if (cycle_length == 6) {
        THoffmanOstenhofDecomposition solution;
        for (int e = 0; e < graph.number_of_edges; ++e) {
            if (!edge_in_cycle[e]) {
                solution.tree.insert(e);
            } else {
                solution.cycle.insert(e);
            }
        }
        graph.hoffman_ostenhof_solutions.push_back(solution);
    } else {
        for (int pe = 0; pe < graph.number_of_edges; ++pe) {
            if (!vertex_in_cycle[graph.e2v[pe][0]] && !vertex_in_cycle[graph.e2v[pe][1]]) {
                THoffmanOstenhofDecomposition solution;
                solution.matching.insert(pe);
                for (int e = 0; e < graph.number_of_edges; ++e) {
                    if (e != pe) {
                        if (!edge_in_cycle[e]) {
                            solution.tree.insert(e);
                        } else {
                            solution.cycle.insert(e);
                        }
                    }
                }
                graph.hoffman_ostenhof_solutions.push_back(solution);
            }
        }
    }
}

void preimage_hoffman_ostenhof(TGraph& graph) {
    int ho_sols = 0;
    // there's no need to check matching - in preimage we'll always get a matching
    // same goes for the cycle
    // and even more - the edges of the tree are spanning
    // so we only need to check the number of edges in the tree and connectedness
    for (int col_idx = 0; col_idx < graph.petersen_colourings.size(); ++col_idx) {
        for (int ho_idx = 0; ho_idx < petersen_graph.hoffman_ostenhof_solutions.size(); ++ho_idx) {
            /*bool vertex_in_tree[MAXN];
            for (int v = 0; v < graph.number_of_vertices; ++v) {
                vertex_in_tree[v] = false;
            }*/
            bool edge_in_tree[REG * MAXN / 2];
            int tree_edge_count = 0;
            for (int e = 0; e < graph.number_of_edges; ++e) {
                edge_in_tree[e] = false;
                int petersen_edge = graph.petersen_colourings[col_idx][e];
                if (petersen_graph.hoffman_ostenhof_solutions[ho_idx].tree.find(petersen_edge) != petersen_graph.hoffman_ostenhof_solutions[ho_idx].tree.end()) {
                    ++tree_edge_count;
                    edge_in_tree[e] = true;
                    /*for (int j = 0; j < 2; ++j) {
                        vertex_in_tree[graph.e2v[e][j]] = true;
                    }*/
                }
            }
            if (tree_edge_count == graph.number_of_vertices - 1) {

                // copypasta
                bool vertex_in_tree[MAXN];
                for (int v = 0; v < graph.number_of_vertices; ++v) {
                    vertex_in_tree[v] = false;
                }

                vertex_in_tree[0] = true;
                int queue_size = 1;
                int queue[MAXN];
                queue[0] = 0;

                for (int cur_idx = 0; cur_idx < queue_size; ++cur_idx) {
                    int cur_vertex = queue[cur_idx];
                    for (int j = 0; j < REG; ++j) {
                        int next_vertex = graph.v2v[cur_vertex][j];
                        if (edge_in_tree[graph.v2e[cur_vertex][j]] && !vertex_in_tree[next_vertex]) {
                            vertex_in_tree[next_vertex] = true;
                            queue[queue_size] = next_vertex;
                            ++queue_size;
                        }
                    }
                }
                bool has_sol = true;
                for (int v = 0; v < graph.number_of_vertices; ++v) {
                    if (!vertex_in_tree[v]) {
                        has_sol = false;
                        break;
                    }
                }
                if (has_sol) {
                    ++ho_sols;
                }
            }
        }
    }
    cerr << "number of hoffman-ostenhof solutions: " << ho_sols << endl;
}

/******************cycles*********************/
int start_vertex[MAXN];
int number_of_circuits = 0;
int separate_circuits[MAXN][2 * MAXN];
int separate_circuits_length[MAXN];
int cycle_count_for_edges[REG * MAXN / 2];


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

        bool all_even = true;
        for (int i = 0; i <= number_of_circuits; ++i) {
            if (separate_circuits_length[i] % 2 == 1) {
                all_even = false;
                break;
            }
        }
        if (all_even) {
            graph.all_even_cycles.insert(bit_cycle);
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
                check_tree_in_hoffman_ostenhof(graph);
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


void preimage_full_cycles(TGraph& graph) {
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

void preimage_6c4c_5cdc_cycles(TGraph& graph) {
    petersen_6c4c_5cdc_pairs.clear();
    petersen_5cdc_6c4c_pairs.clear();
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
                petersen_6c4c_5cdc_pairs[cur_6c4c].insert(cur_5cdc);
                petersen_5cdc_6c4c_pairs[cur_5cdc].insert(cur_6c4c);
            }
        }
    }
}

TMask inv(TGraph& graph, TMask mask) {
    return BIT(graph.number_of_edges) - 1 - mask;
}

int u6by3_shuffles[10][6] = {
    {0, 1, 2, 3, 4, 5},
    {0, 1, 3, 2, 4, 5},
    {0, 1, 4, 2, 3, 5},
    {0, 1, 5, 2, 3, 4},
    {0, 2, 3, 1, 4, 5},
    {0, 2, 4, 1, 3, 5},
    {0, 2, 5, 1, 3, 4},
    {0, 3, 4, 1, 2, 5},
    {0, 3, 5, 1, 2, 4},
    {0, 4, 5, 1, 2, 3}
};
vector<vector<int>> all_333flows_combinations;
set<string> all_333flows_combinations_str;

void gen_333flows_combinations() {
    int idx[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    int edge_count[6][6];
    do {
        if (idx[0] < idx[1] && idx[2] < idx[3] && idx[4] < idx[5] && idx[0] < idx[2] && idx[2] < idx[4]) {
            string mask;
            for (int i = 0; i < 6; ++i) {
                mask += to_string(idx[i]);
            }
            if (all_333flows_combinations_str.find(mask) == all_333flows_combinations_str.end()) {
                for (int i = 0; i < 6; ++i) {
                    for (int j = i + 1; j < 6; ++j) {
                        edge_count[i][j] = 0;
                        edge_count[j][i] = 0;
                    }
                }
                for (int i = 0; i < 6; i += 2) {
                    set<pair<int, int>> cur_edges;
                    for (int j = 0; j < 2; ++j) {
                        for (int k = 0; k < 6; k += 3) {
                            int v1 = u6by3_shuffles[idx[i + j]][k];
                            int v2 = u6by3_shuffles[idx[i + j]][k + 1];
                            int v3 = u6by3_shuffles[idx[i + j]][k + 2];
                            cur_edges.insert(make_pair(v1, v2));
                            cur_edges.insert(make_pair(v2, v1));
                            cur_edges.insert(make_pair(v1, v3));
                            cur_edges.insert(make_pair(v3, v1));
                            cur_edges.insert(make_pair(v2, v3));
                            cur_edges.insert(make_pair(v3, v2));
                        }
                    }
                    for (const auto& p : cur_edges) {
                        ++edge_count[p.first][p.second];
                    }
                }
                bool all2 = true;
                for (int i = 0; i < 6; ++i) {
                    for (int j = i + 1; j < 6; ++j) {
                        if (edge_count[i][j] != 2) {
                            all2 = false;
                            break;
                        }
                    }
                    if (!all2) {
                        break;
                    }
                }
                if (all2) {
                    all_333flows_combinations_str.insert(mask);
                    vector<int> comb;
                    for (int i = 0; i < 6; ++i) {
                        comb.push_back(idx[i]);
                    }
                    all_333flows_combinations.push_back(comb);
                }
            }
        }
    } while (next_permutation(idx, idx + 10));
    cerr << "combinations for 333-flows: " << all_333flows_combinations.size() << endl;
}

vector<TMask> bit_cycles;
TMask u6c4c_cycles[6];
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

bool find_nz_mod5_from_o6c4c(TGraph& graph, int cur_layer) {
    if (cur_layer == 6) {
        for (int e = 0; e < graph.number_of_edges; ++e) {
            if (cur_flow[e] % 5 == 0) {
                return false;
            }
        }
        ++all_nz_mod5_from_o6c4c;
        /*for (int v = 0; v < graph.number_of_vertices; ++v) {
            int vertex_flow = 0;
            for (int j = 0; j < REG; ++j) {
                if (v < graph.v2v[v][j]) {
                    vertex_flow += cur_flow[graph.v2e[v][j]];
                } else {
                    vertex_flow -= cur_flow[graph.v2e[v][j]];
                }
            }
            cerr << vertex_flow % 5 << " ";
        }
        cerr << endl;*/
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
}

bool has_dominating_circuit = false;
map<pair<int, int>, int> edge_pairs;
set<set<int>> all_oriented_vertices;
set<int> oriented_vertices;
bool o6c4c_edge_is_poor[MAXN * REG / 2];
TMask poor_mask = 0;

bool check_dominating_circuit(TGraph& graph) {
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

bool o6c4c_always_1[REG * MAXN / 2];
bool o6c4c_always_2[REG * MAXN / 2];
bool o6c4c_always_3[REG * MAXN / 2];
bool o6c4c_always_4[REG * MAXN / 2];
bool o6c4c_never_1[REG * MAXN / 2];
bool o6c4c_never_2[REG * MAXN / 2];
bool o6c4c_never_3[REG * MAXN / 2];
bool o6c4c_never_4[REG * MAXN / 2];

set<pair<int, int>> u6c4c_edge_pair_counts[REG * MAXN / 2];

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

        check_dominating_circuit(graph);
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

int u33pp_solutions;
int all_33pp_solutions;
TMask u3_inv_pm[3];
bool vertex_in_33pp_cycle[REG * MAXN / 2];
bool edge_in_33pp_cycle[REG * MAXN / 2];
set<pair<TMask, TMask>> u33pp_pairs;
set<TMask> u333pp_cycles_from_o5cdc;
set<TMask> u333pp_even_cycles_from_o5cdc;

set<TMask> u333pp_cycles_from_o6c4c;

set<TMask> all_cycles_from_5cdc;
set<TMask> all_full_cycles_from_5cdc;
set<TMask> all_even_cycles_from_5cdc;
set<TMask> all_circuits_from_5cdc;
bool has_33pp_from_3pm;
bool has_333pp_from_3pm;

set<TMask> all_33pp_cycles;
set<TMask> all_33pp_circuits;
set<TMask> all_33pp_full_cycles;
set<TMask> all_333pp_cycles;
set<TMask> all_333pp_even_cycles;
set<TMask> all_333pp_full_cycles;
int u33pp_partition[REG * MAXN / 2];

bool find_33pp_from_3pm(TGraph& graph) { // TODO: rewrite using check_33pp
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
            u244_6c4c_5cdc_pairs[cur_6c4c].insert(u5cdc);
            u244_5cdc_6c4c_pairs[u5cdc].insert(cur_6c4c);
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

set<vector<int>> all_33pp_triples;

bool find_33pp_from_6c4c(TGraph& graph) {
    has_33pp_from_3pm = false;
    has_333pp_from_3pm = false;
    u33pp_solutions = 0;
    for (int i = 0; i < 6; ++i) {
        u3_inv_pm[0] = u6c4c_cycles[i];
        for (int j = i + 1; j < 6; ++j) {
            u3_inv_pm[1] = u6c4c_cycles[j];
            for (int k = j + 1; k < 6; ++k) {
                u3_inv_pm[2] = u6c4c_cycles[k];
                find_33pp_from_3pm(graph);
                if (has_33pp_from_3pm) {
                    has_33pp_from_3pm = false;
                    vector<int> triple = {i, j, k};
                    all_33pp_triples.insert(triple);
                }
                if (u33pp_solutions > 0) {
                    return true;
                }
            }
        }
    }
    /*if (u33pp_solutions == 0) {
        cout << "nah" << endl;
    }*/
    return false;
}

bool check_3flow(TGraph& graph, const TMask c) {
    int deg[MAXN];
    for (int v = 0; v < graph.number_of_vertices; ++v) {
        deg[v] = 0;
    }
    for (int e = 0; e < graph.number_of_edges; ++e) {
        if ((BIT(e) & c) > 0) {
            for (int j = 0; j < 2; ++j) {
                ++deg[graph.e2v[e][j]];
            }
        }
    }

    bool vertex_in_cur_part[MAXN];
    bool edge_in_cur_part[REG * MAXN / 2];
    int vertex_colour[MAXN];
    for (int v = 0; v < graph.number_of_vertices; ++v) {
        vertex_colour[v] = -1;
        vertex_in_cur_part[v] = false;
        if (deg[v] == 3) {
            vertex_in_cur_part[v] = true;
        }
    }

    for (int e = 0; e < graph.number_of_edges; ++e) {
        edge_in_cur_part[e] = ((BIT(e) & c) > 0);
    }

    bool edge_visited[REG * MAXN / 2];
    for (int e = 0; e < graph.number_of_edges; ++e) {
        edge_visited[e] = !edge_in_cur_part[e];
    }

    int queue[MAXN];
    int queue_size;
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
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

bool has_all_3flows = true;
bool find_333flows_from_6c4c(TGraph& graph) {
    for (const auto& comb : all_333flows_combinations) {
        has_all_3flows = true;
        for (int i = 0; i < 6; i += 2) {
            TMask c = 0;

            set<pair<int, int>> cur_edges;
            for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 6; k += 3) {
                    int v1 = u6by3_shuffles[comb[i + j]][k];
                    int v2 = u6by3_shuffles[comb[i + j]][k + 1];
                    int v3 = u6by3_shuffles[comb[i + j]][k + 2];
                    cur_edges.insert(make_pair(v1, v2));
                    cur_edges.insert(make_pair(v1, v3));
                    cur_edges.insert(make_pair(v2, v3));
                }
            }

            for (int e = 0; e < graph.number_of_edges; ++e) {
                for (const auto& p : cur_edges) {
                    if (((BIT(e) & inv(graph, u6c4c_cycles[p.first])) > 0) && ((BIT(e) & inv(graph, u6c4c_cycles[p.second])) > 0)) {
                        c += BIT(e);
                        break;
                    }
                }
            }
            if (!check_3flow(graph, c)) {
                has_all_3flows = false;
                break;
            }
        }
        if (has_all_3flows) {
            break;
       }
    }
    return has_all_3flows;
}


bool has_o244_flows = false;
int partition_num[REG * MAXN / 2];

int partition_degs[MAXN][REG];
int number_of_flowed_edges_at_vertex[MAXN][REG];
int edge_flows[REG * MAXN / 2][REG];
int vertex_flows[MAXN][REG];

bool is_oriented_edge_covered[REG * MAXN / 2][2];

bool check_nowhere_zeroness(TGraph& graph, int partition);

bool build_4flow(TGraph& graph, int partition, int edge_index) {
    if (edge_index >= graph.number_of_edges) {
        return check_nowhere_zeroness(graph, partition + 1);
    }
    if (partition_num[edge_index] == partition) {
        return build_4flow(graph, partition, edge_index + 1);
    }

    int right_bound = 3;
    int left_bound = -right_bound;

    int v1 = graph.e2v[edge_index][0];
    int v2 = graph.e2v[edge_index][1];

    if (number_of_flowed_edges_at_vertex[v1][partition] == partition_degs[v1][partition] - 1) {
        right_bound = left_bound = vertex_flows[v1][partition];
    }
    if (number_of_flowed_edges_at_vertex[v2][partition] == partition_degs[v2][partition] - 1) {
        right_bound = left_bound = -vertex_flows[v2][partition];
    }

    for (int flow = left_bound; flow <= right_bound; ++flow) {
        if (flow == 0)
            continue;
        int orientation = 0;
        if (flow < 0) {
            orientation = 1;
        }
        if (is_oriented_edge_covered[edge_index][orientation]) {
            continue;
        }
        is_oriented_edge_covered[edge_index][orientation] = true;
        edge_flows[edge_index][partition] = flow;

        vertex_flows[v1][partition] -= flow;
        vertex_flows[v2][partition] += flow;
        ++number_of_flowed_edges_at_vertex[v1][partition];
        ++number_of_flowed_edges_at_vertex[v2][partition];

        bool all_is_good = (
            (vertex_flows[v1][partition] == 0 || number_of_flowed_edges_at_vertex[v1][partition] != partition_degs[v1][partition]) &&
            (vertex_flows[v2][partition] == 0 || number_of_flowed_edges_at_vertex[v2][partition] != partition_degs[v2][partition]) &&
            (abs(vertex_flows[v1][partition]) < 4 || number_of_flowed_edges_at_vertex[v1][partition] != partition_degs[v1][partition] - 1) &&
            (abs(vertex_flows[v2][partition]) < 4 || number_of_flowed_edges_at_vertex[v2][partition] != partition_degs[v2][partition] - 1)
        );
        if (all_is_good && build_4flow(graph, partition, edge_index + 1)) {
            return true;
        }

        // undoing
        is_oriented_edge_covered[edge_index][orientation] = false;
        vertex_flows[v1][partition] += flow;
        vertex_flows[v2][partition] -= flow;
        --number_of_flowed_edges_at_vertex[v1][partition];
        --number_of_flowed_edges_at_vertex[v2][partition];
    }
    return false;
}

bool check_nowhere_zeroness(TGraph& graph, int partition) {
    if (partition == 3) {
        return true;
    }
    for (int e = 0; e < graph.number_of_edges; ++e) {
        edge_flows[e][partition] = 0;
    }
    for (int v = 0; v < graph.number_of_vertices; ++v) {
        vertex_flows[v][partition] = 0;
        number_of_flowed_edges_at_vertex[v][partition] = 0;
    }
    return build_4flow(graph, partition, 0);
}

bool find_o244_flows_from_3pm(TGraph& graph) {
    for (int v = 0; v < graph.number_of_vertices; ++v) {
        partition_degs[v][1] = REG;
        partition_degs[v][2] = REG;
    }
    for (int e = 0; e < graph.number_of_edges; ++e) {
        int edge_count = 0;
        for (int i = 0; i < 3; ++i) {
            if ((BIT(e) & u3_inv_pm[i]) > 0) {
                ++edge_count;
            }
        }
        if (edge_count == 1) {
            partition_num[e] = 1;
        } else if (edge_count == 3) {
            partition_num[e] = 2;
        } else {
            partition_num[e] = 0;
        }
        for (int j = 0; j < 2; ++j) {
            --partition_degs[graph.e2v[e][j]][partition_num[e]];
        }
    }
    for (int e = 0; e < graph.number_of_edges; ++e) {
        for (int partition = 0; partition < 3; ++partition) {
            edge_flows[e][partition] = 0; // TODO: remove
        }
        is_oriented_edge_covered[e][0] = false;
        is_oriented_edge_covered[e][1] = false;
    }
    if (check_nowhere_zeroness(graph, 1)) {
        has_o244_flows = true;
    }
    return false;
}

set<vector<int>> o244_triples;

bool find_o244_flows_from_6c4c(TGraph& graph) {
    has_o244_flows = false;
    for (int i = 0; i < 6; ++i) {
        u3_inv_pm[0] = u6c4c_cycles[i];
        for (int j = i + 1; j < 6; ++j) {
            u3_inv_pm[1] = u6c4c_cycles[j];
            for (int k = j + 1; k < 6; ++k) {
                u3_inv_pm[2] = u6c4c_cycles[k];
                vector<int> triple = {i, j, k};
                if (all_33pp_triples.find(triple) == all_33pp_triples.end()) {
                    continue;
                }
                find_o244_flows_from_3pm(graph);
                if (has_o244_flows) {
                    o244_triples.insert(triple);
                    return true;
                    has_o244_flows = false;
                    //return true;
                }
            }
        }
    }
    return false;
}

bool u6c4c_always_poor[REG * MAXN / 2];
bool u6c4c_always_rich[REG * MAXN / 2];
int u6c4c_min_poor;
int u6c4c_max_poor;

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

bool has_or_comb = false;
bool has_un_comb = false;

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

        all_33pp_triples.clear();
        find_33pp_from_6c4c(graph);
        has_33pp_from_3pm = all_33pp_triples.size() > 0;
        if (!has_33pp_from_3pm) {
            return false;
        }

        //all_oriented_vertices.clear();
        bool has_o6c4c = (same_cycles_different_orientations > 0);
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

        find_333flows_from_6c4c(graph);
        if (!has_all_3flows) {
            return false;
        }

        o244_triples.clear();
        find_o244_flows_from_6c4c(graph);
        has_o244_flows = o244_triples.size() > 0;
        if (!has_o244_flows) {
            return false;
        }

        bool have_same_triple = false;
        for (const auto& triple : all_33pp_triples) {
            if (o244_triples.find(triple) != o244_triples.end()) {
                have_same_triple = true;
                break;
            }
        }
        if (!have_same_triple) {
            return false;
        }

        if (has_o6c4c && has_33pp_from_3pm && has_all_3flows && has_o244_flows && have_same_triple && has_dominating_circuit) {
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

/******************5cdc*********************/

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
            if (edge_orientation_count[ei][cur_edge_orientation] == 1) {
                break;
            }
            ++edge_orientation_count[ei][cur_edge_orientation];
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
            --edge_orientation_count[ei][cur_edge_orientation];
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
        cerr << "found " << same_cycles_different_orientations_in_5cdc << " same cycles different orientations for o6c4c" << endl;
    }*/
    return false;
}

int u33pp_count;

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
            //if (edge_orientation_count[ei][cur_edge_orientation] == 1) {
                // TODO: record this event as inconsistency between layers
            //}
            if (cur_edge_orientation == 0) {
                layer_flow_in_5cdc[layer_in_5cdc[cur_circuit]][ei] += 1;
            } else {
                layer_flow_in_5cdc[layer_in_5cdc[cur_circuit]][ei] -= 1;
            }

            ++edge_orientation_count[ei][cur_edge_orientation];
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
            --edge_orientation_count[ei][cur_edge_orientation];
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

void find_all_o6c4c(TGraph& graph, bool only_find = false) {
    u244_6c4c_5cdc_pairs.clear();
    u244_5cdc_6c4c_pairs.clear();
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
    if (!has_or_comb) {
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

bool check_33pp(TGraph& graph, const TMask c) {
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
                edge_in_cur_part[e] = (u33pp_partition[e] != part);

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
    all_33pp_cycles.insert(c);
    if (graph.all_full_cycles.find(c) != graph.all_full_cycles.end()) {
        all_33pp_full_cycles.insert(c);
    }
    if (graph.all_circuits.find(c) != graph.all_circuits.end()) {
        all_33pp_circuits.insert(c);
    }

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

/*
void create_graph_from_partition(TGraph& input_graph, const vector<int>& partition, int part, TGraph& output_graph) {
    output_graph.number_of_vertices = input_graph.number_of_vertices;
    for (int v = 0; v < output_graph.number_of_vertices; ++v) {
        output_graph.deg[v] = 0;
    }
    output_graph.number_of_edges = 0;
    for (int e = 0; e < input_graph.number_of_edges; ++e) {
        if (partition[e] != part) { // TODO: rewrite with masks
            add_edge(output_graph, output_graph.number_of_edges, input_graph.e2v[e][0], input_graph.e2v[e][1]);
            ++output_graph.number_of_edges;
        }
    }
}


void find_construction(TGraph& graph) {
    vector<int> four_partition;
    for (int e = 0; e < graph.number_of_edges; ++e) {
        four_partition.push_back(0);
    }
    for (const auto& cycle_and_circuits : graph.cycles_as_circuits) {
        // 1. search for all 244-flows

        // 1.1. check that all circuits are of even length
        bool all_even = true;
        for (const auto& c : cycle_and_circuits.second) {
            if ((c.size() - 1) % 2) {
                all_even = false;
                break;
            }
        }
        if (!all_even) {
            continue;
        }

        // 1.2 now we need to partition edges in this circuit, for 244-flows
        for (int e = 0; e < graph.number_of_edges; ++e) {
            if ((BIT(e) & cycle_and_circuits.first)) {
                four_partition[e] = 0; // 0 means it's not in 2-flow, it's in both 4-flows // TODO: rewrite with masks
            }
        }

        for (TMask four_partition_mask = 0; four_partition_mask < BIT(cycle_and_circuits.second.size() - 1); ++four_partition_mask) {
            int idx = 0;
            for (const auto& c : cycle_and_circuits.second) {
                bool reverse = (BIT(idx) & four_partition_mask) > 0;
                for (int vi = 0; vi < c.size() - 1; ++vi) {
                    int e = graph.edge_index[c[vi]][c[vi + 1]];
                    four_partition[e] = (vi % 2) + 1; // TODO: rewrite with masks
                    if (reverse) {
                        four_partition[e] = 3 - four_partition[e];
                    }
                }
                ++idx;
            }
        }

        // 1.3 create graphs
        TGraph g1, g2;
        create_graph_from_partition(graph, four_partition, 1, g1);
        create_graph_from_partition(graph, four_partition, 2, g2);

        // 1.4 check 4-flows
    }
}*/

void compare_6c4c_5cdc_pairs(TGraph& graph) {
    cerr << "all 6c4c: " << all_6c4c_solutions << endl;
    cerr << "all 6c4c (in 6c4c-5cdc pairs): " << u244_6c4c_5cdc_pairs.size() << endl;
    cerr << "pet 6c4c: " << petersen_6c4c_5cdc_pairs.size() << endl;
    for (const auto& u6c4c_pairs : u244_6c4c_5cdc_pairs) {
        for (const auto& u5cdc : u6c4c_pairs.second) {
            if (petersen_6c4c_5cdc_pairs.find(u6c4c_pairs.first) != petersen_6c4c_5cdc_pairs.end()) {
                if (petersen_6c4c_5cdc_pairs[u6c4c_pairs.first].find(u5cdc) != petersen_6c4c_5cdc_pairs[u6c4c_pairs.first].end()) {
                    cerr << "both" << endl;
                    //return;
                }
            }
        }
    }
    for (const auto& petersen_6c4c_pairs : petersen_6c4c_5cdc_pairs) {
        cerr << "is petersen 6c4c oriented: " << (graph.all_o6c4c.find(petersen_6c4c_pairs.first) != graph.all_o6c4c.end()) << endl;
        for (const auto& petersen_5cdc : petersen_6c4c_pairs.second) {
            //cerr << "is petersen 5cdc oriented: " << (graph.all_o5cdc.find(petersen_5cdc) != graph.all_o5cdc.end()) << endl;
        }
    }
    for (const auto& petersen_5cdc_pairs : petersen_5cdc_6c4c_pairs) {
        //cerr << "is petersen 5cdc in 6c4c-244-flows-33-pp construction: " << (u244_5cdc_6c4c_pairs.find(petersen_5cdc_pairs.first) != u244_5cdc_6c4c_pairs.end()) << endl;
    }
    for (const auto& petersen_6c4c : graph.petersen_6c4c) {
        //cerr << "is petersen 6c4c oriented: " << (graph.all_o6c4c.find(petersen_6c4c) != graph.all_o6c4c.end()) << endl;
    }
}

int main(int argc, char** argv) {
    srand(time(NULL));
    int to_skip = 0;
    int number_of_vertices;
    int number_of_edges;
    if (argc < 3 || argc > 4) {
        cerr << "Error: invalid number of arguments" << endl;
        cerr << "Usage: " << argv[0] << " <path_to_petersen_graph> <number_of_vertices>" << endl;
        exit(1);
    } else {
        number_of_vertices = atoi(argv[2]);
        number_of_edges = REG * number_of_vertices / 2;
        if (number_of_vertices > MAXN) {
            cerr << "Number of vertices is too big (limit is " << MAXN << ")" << endl;
            exit(0);
        }
        if (argc >= 4) {
            to_skip = atoi(argv[3]);
        }
    }

    // initialize petersen graph
    petersen_graph.number_of_vertices = 10;
    petersen_graph.number_of_edges = petersen_graph.number_of_vertices * REG / 2;
    int petersen_codelength = petersen_graph.number_of_vertices * REG / 2 + petersen_graph.number_of_vertices;
    unsigned char petersen_code[petersen_codelength];
    FILE* petersen_file;
    petersen_file = fopen(argv[1], "rb");
    fread(petersen_code, sizeof(unsigned char), petersen_codelength, petersen_file);
    decode_multicode(petersen_code, petersen_codelength, petersen_graph);
    print_graph(petersen_graph);

    find_all_petersen_colourings(petersen_graph);
    create_cc_mapping();
    prepare_build_cycle(petersen_graph);
    cerr << "petersen cycles: " << petersen_graph.all_cycles.size() << endl;
    cerr << "petersen perfect matchings: " << petersen_graph.all_full_cycles.size() << endl;
    //cerr << "petersen hoffman-ostenhof solutions: " << petersen_graph.hoffman_ostenhof_solutions.size() << endl;

    find_all_o5cdc(petersen_graph, true);
    find_all_o6c4c(petersen_graph, true);
    cerr << "petersen 5cdc: " << petersen_graph.all_5cdc.size() << endl;
    cerr << "petersen 6c4c: " << petersen_graph.all_6c4c.size() << endl;

    gen_333flows_combinations();

    int codelength = number_of_vertices * REG / 2 + number_of_vertices;
    unsigned char code[codelength];
    TMask number_of_graphs_without_solution = 0;
    while (fread(code, sizeof(unsigned char), codelength, stdin)) {
        TGraph ref_graph;
        ref_graph.number_of_vertices = number_of_vertices;
        ref_graph.number_of_edges = number_of_edges;
        decode_multicode(code, codelength, ref_graph);
        ++number_of_graphs_read;
        if (to_skip >= number_of_graphs_read)
            continue;
        cerr << "g" << number_of_graphs_read << "\t" << endl << flush;
        cout << "g" << number_of_graphs_read << "\t" << endl << flush;

        prepare_build_cycle(ref_graph);
        //cerr << "graph cycles: " << ref_graph.all_cycles.size() << endl;
        //cerr << "graph not full cycles: " << ref_graph.all_cycles.size() - ref_graph.all_full_cycles.size() << endl;
        //cerr << "graph circuits: " << ref_graph.all_circuits.size() << endl;
        //cerr << "graph dominating circuits: " << ref_graph.all_dominating_circuits.size() << endl;
        //cerr << "graph even cycles: " << ref_graph.all_even_cycles.size() << endl;
        //cerr << "graph perfect matchings: " << ref_graph.all_full_cycles.size() << endl;
        //cerr << endl;
        find_all_petersen_colourings(ref_graph);
        //cerr << "petersen poor: " << petersen_min_poor << "-" << petersen_max_poor << endl;
        /*cerr << "petersen always poor: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << petersen_always_poor[e] << " ";
        }
        cerr << endl;
        cerr << "petersen always rich: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << petersen_always_rich[e] << " ";
        }
        cerr << endl;*/

        //preimage_full_cycles(ref_graph);
        preimage_6c4c_5cdc_cycles(ref_graph);
        //find_o6c4c_compatible_with_preimages(ref_graph);
        //find_all_33pp(ref_graph);
        
        //cerr << "all 33pp cycles: " << all_33pp_cycles.size() << endl;
        //cerr << "all 33pp not full cycles: " << all_33pp_cycles.size() - all_33pp_full_cycles.size() << endl;
        //cerr << "all 33pp circuits: " << all_33pp_circuits.size() << endl;
        //cerr << "all 33pp full cycles: " << all_33pp_full_cycles.size() << endl;
        //cerr << "all 333pp cycles: " << all_333pp_cycles.size() << endl;
        //cerr << "all 333pp even cycles: " << all_333pp_even_cycles.size() << endl;
        //cerr << "all 333pp full cycles: " << all_333pp_full_cycles.size() << endl;
        //cerr << endl;

        find_all_o5cdc(ref_graph);
        //cerr << "all 5cdc cycles: " << all_cycles_from_5cdc.size() << endl;
        //cerr << "all 5cdc not full cycles: " << all_cycles_from_5cdc.size() - all_full_cycles_from_5cdc.size() << endl;
        //cerr << "all 5cdc circuits: " << all_circuits_from_5cdc.size() << endl;
        //cerr << "all 5cdc full cycles: " << all_full_cycles_from_5cdc.size() << endl;
        //cerr << "all even cycles from 5cdc: " << all_even_cycles_from_5cdc.size() << endl;

        cerr << "searching for combination: ";
        find_all_o6c4c(ref_graph);
        //cerr << "all 6c4c poor: " << u6c4c_min_poor << "-" << u6c4c_max_poor << endl;
        /*cerr << "all 6c4c always poor: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << u6c4c_always_poor[e] << " ";
        }
        cerr << endl;
        cerr << "all 6c4c always rich: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << u6c4c_always_rich[e] << " ";
        }*/

        /*cerr << "   o6c4c always    1: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << o6c4c_always_1[e] << " ";
        }
        cerr << endl;
        cerr << "   o6c4c always    2: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << o6c4c_always_2[e] << " ";
        }
        cerr << endl;
        cerr << "   o6c4c always    3: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << o6c4c_always_3[e] << " ";
        }
        cerr << endl;
        cerr << "   o6c4c always    4: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << o6c4c_always_4[e] << " ";
        }
        cerr << endl;
        cerr << endl;

        cerr << "   o6c4c never     1: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << o6c4c_never_1[e] << " ";
        }
        cerr << endl;
        cerr << "   o6c4c never     2: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << o6c4c_never_2[e] << " ";
        }
        cerr << endl;
        cerr << "   o6c4c never     3: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << o6c4c_never_3[e] << " ";
        }
        cerr << endl;
        cerr << "   o6c4c never     4: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << o6c4c_never_4[e] << " ";
        }
        cerr << endl;*/

        //cerr << "333pp cycles from o5cdc: " << u333pp_cycles_from_o5cdc.size() << endl;
        //cerr << "333pp even cycles from o5cdc: " << u333pp_even_cycles_from_o5cdc.size() << endl;
        compare_6c4c_5cdc_pairs(ref_graph);
        cerr << endl;

        //cerr << "333pp even cycles from 6c4c: " << u333pp_cycles_from_o6c4c.size() << endl;
        //cerr << endl;
        //cerr << "33pp pairs: " << u33pp_pairs.size() << endl;
        //cerr << endl;
        //preimage_hoffman_ostenhof(ref_graph);
        //break;

        //find_construction(ref_graph);
    }
    cerr << "fin" << endl;
    return(EXIT_SUCCESS);
}
