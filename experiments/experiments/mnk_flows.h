/*
 * File:   mnk_flows.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * code for double covers with subgraphs with m-, n- and k- nowhere-zero flows
 * 
 * TODO: split into .h and .cpp files
 *
 */

#pragma once

#include "graph.h"
#include "common.h"

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
#include <string>

using namespace std;


namespace Exp6c4c::ExpMNKFlows {

using namespace UtilFlows;

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

bool has_all_3flows = true;

bool has_o244_flows = false; // FIXME: rewrite this flag
int partition_num[MAX_EDGE];
int partition_degs[MAX_VERTEX][MAX_DEG];
int number_of_flowed_edges_at_vertex[MAX_VERTEX][MAX_DEG];
int edge_flows[MAX_EDGE][MAX_DEG];
int vertex_flows[MAX_VERTEX][MAX_DEG];

bool is_oriented_edge_covered[MAX_EDGE][2];

bool check_nowhere_zeroness(Graph& graph, int partition);

set<vector<int>> o244_triples;

/*********************************Methods*********************************/

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

    // TODO: under flag
    cerr << "combinations for 333-flows: " << all_333flows_combinations.size() << endl;
}

bool find_333flows_from_6c4c(const set<vector<int>>& common_triples, Graph& graph) {
    for (const auto& comb : all_333flows_combinations) {
        // check that we have compatible common triple
        bool has_compat_triple = false;
        for (const auto& triple : common_triples) {
            map<int, int> counts;
            bool all_compat = true;
            for (int i = 0; i < 6; ++i) {
                // second triple: u6by3_shuffles[comb[i]][:]
                vector<int> diffs(6);
                vector<int>::iterator it;
                it = set_symmetric_difference(triple.begin(), triple.end(), u6by3_shuffles[comb[i]], u6by3_shuffles[comb[i]] + 3, diffs.begin());
                diffs.resize(it - diffs.begin());
                if (diffs.size() != 2 && diffs.size() != 4) {
                    all_compat = false;
                    break;
                }
                if (diffs.size() == 4) {
                    it = set_symmetric_difference(triple.begin(), triple.end(), u6by3_shuffles[comb[i]] + 3, u6by3_shuffles[comb[i]] + 6, diffs.begin());
                    diffs.resize(it - diffs.begin());
                }
                if (diffs.size() != 2) {
                    cerr << "wat" << endl;
                }
                for (const int& d : diffs) {
                    ++counts[d];
                    if (counts[d] > 2) {
                        all_compat = false;
                        break;
                    }
                }
            }
            if (all_compat) {
                has_compat_triple = true;
                break;
            }
        }
        if (!has_compat_triple) {
            continue;
        }
        has_all_3flows = true;
        for (int i = 0; i < 6; i += 2) {
            Mask c = 0;

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

bool build_4flow(Graph& graph, int partition, int edge_index) {
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

bool check_nowhere_zeroness(Graph& graph, int partition) {
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

bool find_o244_flows_from_3pm(Graph& graph) {
    for (int v = 0; v < graph.number_of_vertices; ++v) {
        partition_degs[v][1] = MAX_DEG;
        partition_degs[v][2] = MAX_DEG;
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

bool find_o244_flows_from_6c4c(const set<vector<int>>& all_33pp_triples, Graph& graph) {
    has_o244_flows = false;
    // int orientation_fail_count = 0;
    for (int i = 0; i < 6; ++i) {
        u3_inv_pm[0] = u6c4c_cycles[i];
        for (int j = i + 1; j < 6; ++j) {
            u3_inv_pm[1] = u6c4c_cycles[j];
            for (int k = j + 1; k < 6; ++k) {
                u3_inv_pm[2] = u6c4c_cycles[k];
                vector<int> triple = {i, j, k};
                if (all_33pp_triples.find(triple) == all_33pp_triples.end()) { // TODO: remove from here
                    continue;
                }
                find_o244_flows_from_3pm(graph);
                if (has_o244_flows) {
                    o244_triples.insert(triple);
                    //return true;
                    has_o244_flows = false;
                    //return true;
                // } else {
                //     // FIXME2024
                //     orientation_fail_count += 1;
                //     if (orientation_fail_count > 1) {
                //         return false;
                //     }
                }
            }
        }
    }
    return false;
}

} // Exp6c4c::ExpMNKFlows
