/*
 * File:   petersen_colouring.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
 *
 */

#pragma once

#include "graph.h"

//TODO: remove some includes
#include <assert.h>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <map>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>

namespace NExpPetersenColouring {

using namespace std;

bool edge_is_poor[REG * MAXN / 2];
bool petersen_vertex_is_poor[10]; // vertex is poor when at least one of images of petersen vertices is connected to itself
int edge_colour[REG * MAXN / 2];
int not_coloured_edges_count_near_vertex[MAXN];
bool had_four = false;

// cc-mapping
unordered_map<int, int> ve_to_petersen_edge;
unordered_map<int, int> mask_to_petersen_vertex;

bool petersen_always_poor[REG * MAXN / 2];
bool petersen_always_rich[REG * MAXN / 2];
int petersen_min_poor;
int petersen_max_poor;


/*********************************Methods*********************************/

bool petersen_colouring_check_edge(int e, const TGraph& graph) {
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
            petersen_vertex_is_poor[v] = false;
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
                petersen_vertex_is_poor[mask_to_petersen_vertex[mask1]] = true;
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

        // TODO: remove this?
        bool dominating = false;
        for (int v = 0; v < 10; ++v) {
            if (!petersen_vertex_is_poor[v]) {
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

    int cur_edge = graph.faster_edge_order[cur_edge_idx];
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
            if ((ei != cur_edge && edge_colour[ei] == edge_colour[cur_edge]) || (not_coloured_edges_count_near_vertex[v1] == 0 && !petersen_colouring_check_edge(ei, graph))) {
                checks_passed = false;
                break;
            }
        }
        if (checks_passed) {
            for (int j = 0; j < REG; ++j) {
                int ei = graph.v2e[v2][j];
                if (ei != cur_edge && (edge_colour[ei] == edge_colour[cur_edge] || (not_coloured_edges_count_near_vertex[v2] == 0 && !petersen_colouring_check_edge(ei, graph)))) {
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

void find_all_petersen_colourings(TGraph& graph) {
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

// FIXME: fix number of stars here
/********************************cc-mapping********************************/

void create_cc_mapping(const TGraph& petersen_graph) {
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

} // NExpPetersenColouring
