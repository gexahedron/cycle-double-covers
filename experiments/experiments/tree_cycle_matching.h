/*
 * File:   tree_cycle_matching.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
 */

#pragma once

#include "graph.h"

// TODO: remove unused includes
#include <iostream>
#include <cstdlib>

using namespace std;


namespace ExpTreeCycleMatching {

bool vertex_in_tcm_cycle[MAX_VERTEX];
bool edge_in_tcm_cycle[MAX_EDGE];
int tcm_cycle_length;

/*********************************Methods*********************************/

void check_tcm_decomposition_in_petersen_graph(Graph& graph) {
    if (graph.number_of_vertices != 10) {
        return;
    }
    bool vertex_in_tree[MAX_VERTEX];
    for (int v = 0; v < graph.number_of_vertices; ++v) {
        vertex_in_tree[v] = false;
    }

    vertex_in_tree[0] = true;
    int queue_size = 1;
    int queue[MAX_VERTEX];
    queue[0] = 0;

    for (int cur_idx = 0; cur_idx < queue_size; ++cur_idx) {
        int cur_vertex = queue[cur_idx];
        for (int j = 0; j < MAX_DEG; ++j) {
            int next_vertex = graph.v2v[cur_vertex][j];
            if (!edge_in_tcm_cycle[graph.v2e[cur_vertex][j]] && !vertex_in_tree[next_vertex]) {
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

    if (tcm_cycle_length == 6) {
        TreeCycleMatching solution;
        for (int e = 0; e < graph.number_of_edges; ++e) {
            if (!edge_in_tcm_cycle[e]) {
                solution.tree.insert(e);
            } else {
                solution.cycle.insert(e);
            }
        }
        graph.tree_cycle_matchings.push_back(solution);
    } else {
        for (int pe = 0; pe < graph.number_of_edges; ++pe) {
            if (!vertex_in_tcm_cycle[graph.e2v[pe][0]] && !vertex_in_tcm_cycle[graph.e2v[pe][1]]) {
                TreeCycleMatching solution;
                solution.matching.insert(pe);
                for (int e = 0; e < graph.number_of_edges; ++e) {
                    if (e != pe) {
                        if (!edge_in_tcm_cycle[e]) {
                            solution.tree.insert(e);
                        } else {
                            solution.cycle.insert(e);
                        }
                    }
                }
                graph.tree_cycle_matchings.push_back(solution);
            }
        }
    }
}

void preimage_tree_cycle_matchings(const Graph& petersen_graph, Graph& graph) {
    int tcm_count = 0;
    // there's no need to check matching - in preimage we'll always get a matching
    // same goes for the cycle
    // and even more - the edges of the tree are spanning
    // so we only need to check the number of edges in the tree and connectedness
    for (int col_idx = 0; col_idx < graph.petersen_colourings.size(); ++col_idx) {
        for (int ho_idx = 0; ho_idx < petersen_graph.tree_cycle_matchings.size(); ++ho_idx) {
            /*bool vertex_in_tree[MAX_VERTEX];
            for (int v = 0; v < graph.number_of_vertices; ++v) {
                vertex_in_tree[v] = false;
            }*/
            bool edge_in_tree[MAX_EDGE];
            int tree_edge_count = 0;
            for (int e = 0; e < graph.number_of_edges; ++e) {
                edge_in_tree[e] = false;
                int petersen_edge = graph.petersen_colourings[col_idx][e];
                if (petersen_graph.tree_cycle_matchings[ho_idx].tree.find(petersen_edge) !=
                        petersen_graph.tree_cycle_matchings[ho_idx].tree.end())
                {
                    ++tree_edge_count;
                    edge_in_tree[e] = true;
                    /*for (int j = 0; j < 2; ++j) {
                        vertex_in_tree[graph.e2v[e][j]] = true;
                    }*/
                }
            }
            if (tree_edge_count == graph.number_of_vertices - 1) {

                // copypasta
                bool vertex_in_tree[MAX_VERTEX];
                for (int v = 0; v < graph.number_of_vertices; ++v) {
                    vertex_in_tree[v] = false;
                }

                vertex_in_tree[0] = true;
                int queue_size = 1;
                int queue[MAX_VERTEX];
                queue[0] = 0;

                for (int cur_idx = 0; cur_idx < queue_size; ++cur_idx) {
                    int cur_vertex = queue[cur_idx];
                    for (int j = 0; j < MAX_DEG; ++j) {
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
                    ++tcm_count;
                }
            }
        }
    }
    cerr << "number of tree-cycle-matching solutions: " << tcm_count << endl;
}

} // ExpTreeCycleMatching
