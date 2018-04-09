/*
 * File:   util.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * Contains some useful functions about graph parsing and construction and some auxilary functions
 *
 */

#pragma once

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>

using namespace std;

// FIXME: fix number of stars here
/**********************************Defines*************************************/
using TMask = unsigned long long int;

#define BIT(i) (1ULL << (i))

constexpr int MAXN = 42; // Because edge labels are sometimes saved in a bitvector - 42 * 3 / 2 < 64
// TODO: remove MAXN, replace with MAX_VERTEX
constexpr int REG = 3;
constexpr int MAX_VERTEX = 42;
constexpr int MAX_EDGE = REG * MAXN / 2;
// TODO: replace all REG * MAXN / 2 with MAX_EDGE
constexpr int NONE = -1;

struct TTreeCycleMatching {
    unordered_set<int> tree;
    unordered_set<int> cycle;
    unordered_set<int> matching;
};

struct TGraph {
    int number_of_vertices = 0;
    int number_of_edges = 0;
    int oddness = 0;
    unsigned int deg[MAXN];
    unsigned int v2v[MAXN][REG];
    unsigned int v2e[MAXN][REG];
    unsigned int e2v[REG * MAXN / 2][2];

    unsigned int lexi_deg[MAXN];
    unsigned int lexi_v2v[MAXN][REG];

    unsigned int edge_index[MAXN][MAXN];
    // TODO: add edge neighbourhood
    // edges, that are adjacent to current edge and which have lower index

    // petersen colouring
    vector<vector<int>> normal5_colourings; // just colours
    vector<vector<int>> petersen_colourings; // petersen edge numbers
    unordered_map<string, int> profiles;

    // cycles
    unordered_set<TMask> all_cycles;
    unordered_set<TMask> all_circuits;
    unordered_set<TMask> all_even_cycles;
    unordered_set<TMask> all_even_v_minus_4_cycles;
    unordered_set<TMask> all_full_cycles;
    unordered_set<TMask> all_dominating_circuits;
    unordered_map<TMask, vector<vector<int>>> cycles_as_circuits;

    unordered_set<TMask> all_vertex_neib_masks;

    // flows
    vector<vector<int>> all_nz5_flows;

    // cycle covers
    set<set<TMask>> all_5cdc;
    set<set<TMask>> all_o5cdc;
    set<set<TMask>> all_6c4c;
    set<set<TMask>> all_o6c4c;

    set<set<TMask>> petersen_5cdc;
    set<set<TMask>> petersen_6c4c;

    // tree-cycle-matching decompositions
    vector<TTreeCycleMatching> tree_cycle_matchings;

    // preimages
    map<set<TMask>, set<set<TMask>>> petersen_5cdc_6c4c_pairs;
    map<set<TMask>, set<set<TMask>>> petersen_6c4c_5cdc_pairs;

    // o6c4c
    map<set<TMask>, set<set<TMask>>> u244_6c4c_5cdc_pairs;
    map<set<TMask>, set<set<TMask>>> u244_5cdc_6c4c_pairs;

    // misc
    vector<int> faster_edge_order;

    TGraph(unsigned int n = 0);
    void add_edge(unsigned int v, unsigned int w);
    bool has_edge(unsigned int v, unsigned int w) const;
    void finish_init();
    void print() const;
    void find_faster_edge_order();
};

/*********************************Methods*********************************/

// Decodes the code (which is in multicode format) of a graph.
bool decode_multicode(FILE* input, TGraph& graph);

void print_graph(const TGraph& graph);

inline int combine_hash_from_mask_and_edge_colour(int mask, int edge_colour) {
    return mask + edge_colour * BIT(6);
}

inline TMask inv(TGraph& graph, TMask mask) {
    return BIT(graph.number_of_edges) - 1 - mask;
}

void find_faster_edge_order(TGraph& graph);
