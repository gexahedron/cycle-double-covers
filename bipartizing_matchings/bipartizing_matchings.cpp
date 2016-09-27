/*
 * File:   bipartizing_matchings.cpp
 * Author: Nikolay Ulyanov
 *
 * Created on 23 March 2016
 */

#include "util.h"

#include <unordered_set>
#include <unordered_map>
#include <set>
#include <climits>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>

using namespace std;

/*************************Defines and global variables*************************/
const bool PRINT_SOLUTIONS = false;

GRAPH graph;
int number_of_vertices;
int number_of_edges;
int number_of_graphs_to_skip = 0;
unsigned int edge_index[MAXN][MAXN];

bool admissible_length[MAXN];
int start_vertex[MAXN];
bool ignored[MAXN];
bool in_circuit[MAXN];
int number_of_ignored = 0;
int circuit_length = 0;
int separate_circuit_lengths[MAXN];
int separate_circuits[MAXN][2 * MAXN];
int degree[MAXN];
int circuit[MAXN];
bool edge_in_circuit[REG * MAXN / 2];

int not_num[REG][2] = {{1, 2}, {2, 0}, {0, 1}};
int permutations[6][REG] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
int just_ignored[MAXN];
int ignored_permutation[MAXN];

int number_of_leftout_edges = 0;
int leftout_edges[MAXN];
int leftout_partition[MAXN];

int number_of_pairs[REG];
bool is_in[MAXN][REG];
int pairs[REG][MAXN][2];
int bipartition[MAXN][REG];

int vertices[REG * MAXN / 2][2];

int number_of_circuits = 0;

// structures for bipartitized graphs
GRAPH new_graph;
int colour[MAXN];
int queue[MAXN];
int queue_size = 0;

int number_of_solutions = 0;

/*********************************Methods**************************************/

bool check_3bms() {
    for (int j = 0; j < REG; ++j) {
        number_of_pairs[j] = 0;
    }
    for (int i = 0; i < number_of_vertices; ++i) {
        for (int j = 0; j < REG; ++j) {
            is_in[i][j] = false;
        }
    }
    for (int i = 0; i < number_of_ignored; ++i) {
        int ignored_vertex = just_ignored[i];
        for (int j = 0; j < REG; ++j) {
            for (int k = 0; k < 2; ++k) {
                int v = graph[ignored_vertex][permutations[ignored_permutation[i]][not_num[j][k]]];
                pairs[j][number_of_pairs[j]][k] = v;
                is_in[v][j] = true;
            }
            ++number_of_pairs[j];
        }
    }
    for (int i = 0; i < number_of_leftout_edges; ++i) {
        int ei = leftout_edges[i];
        for (int j = 0; j < 2; ++j) {
            int part = not_num[leftout_partition[i]][j];
            for (int k = 0; k < 2; ++k) {
                int v = vertices[ei][k];
                pairs[part][number_of_pairs[part]][k] = v;
                is_in[v][part] = true;
            }
            ++number_of_pairs[part];
        }
    }

    // colourize
    bool res[REG];
    for (int part = 0; part < REG; ++part) {
        res[part] = true;
        for (int i = 0; i < number_of_pairs[part]; ++i) {
            new_graph[pairs[part][i][0]][0] = pairs[part][i][1];
            new_graph[pairs[part][i][1]][0] = pairs[part][i][0];
        }
        for (int i = 0; i <= number_of_circuits; ++i) {
            int start_idx = 0;
            for (int j = 0; j < separate_circuit_lengths[i]; ++j) {
                if (is_in[separate_circuits[i][j]][part]) {
                    start_idx = j;
                    break;
                }
            }
            int prev_vertex = separate_circuits[i][start_idx];
            for (int add = 1; add <= separate_circuit_lengths[i]; ++add) {
                int cur_vertex = separate_circuits[i][start_idx + add];
                if (is_in[cur_vertex][part]) {
                    new_graph[prev_vertex][1] = cur_vertex;
                    new_graph[cur_vertex][2] = prev_vertex;
                    prev_vertex = cur_vertex;
                }
            }
        }
        int start_vertex = -1;
        for (int i = 0; i < number_of_vertices; ++i) {
            colour[i] = -1;
            if (is_in[i][part]) {
                start_vertex = i;
            }
        }
        colour[start_vertex] = 0;
        queue_size = 1;
        queue[0] = start_vertex;
        for (int cur_idx = 0; cur_idx < queue_size; ++cur_idx) {
            int cur_vertex = queue[cur_idx];
            for (int j = 0; j < REG; ++j) {
                int neib_vertex = new_graph[cur_vertex][j];
                if (colour[neib_vertex] == -1) {
                    colour[neib_vertex] = 1 - colour[cur_vertex];
                    queue[queue_size] = neib_vertex;
                    ++queue_size;
                } else if (colour[neib_vertex] == colour[cur_vertex]) {
                    res[part] = false;
                    break;
                }
            }
            if (!res[part]) {
                break;
            }
        }
    }

    int res_count = 0;
    for (int i = 0; i < REG; ++i) {
        if (res[i]) {
            ++res_count;
        }
    }

    int threshold = REG;
    if (res_count >= threshold && PRINT_SOLUTIONS) {
        if (number_of_solutions == 0)
            PrintGraph(graph, number_of_vertices);

        cerr << "res: ";
        for (int i = 0; i < REG; ++i)
            cerr << res[i] << " ";
        cerr << endl;

        cerr << "circuit lengths: ";
        for (int i = 0; i <= number_of_circuits; ++i)
            cerr << separate_circuit_lengths[i] << " ";
        cerr << endl;

        cerr << "circuits: " << endl;
        for (int i = 0; i <= number_of_circuits; ++i) {
            separate_circuits[i][separate_circuit_lengths[i]] = separate_circuits[i][0];
            for (int j = 0; j < separate_circuit_lengths[i]; ++j) {
                cerr << separate_circuits[i][j] << " ";
            }
            cerr << endl;
        }

        cerr << "ignored: ";
        for (int i = 0; i < number_of_ignored; ++i) {
            int ignored_vertex = just_ignored[i];
            cerr << ignored_vertex << " (";
            for (int j = 0; j < REG; ++j) {
                cerr << graph[ignored_vertex][permutations[ignored_permutation[i]][j]];
                if (j < REG - 1)
                    cerr << ", ";
            }
            cerr << "); ";
        }
        cerr << endl;
        cerr << number_of_leftout_edges << " leftouts: ";
        for (int i = 0; i < number_of_leftout_edges; ++i) {
            cerr << "(" << vertices[leftout_edges[i]][0] << ", " << vertices[leftout_edges[i]][1] << ") -> " << leftout_partition[i] << "; ";
        }
        cerr << endl;
    }

    if (res_count >= threshold) {
        ++number_of_solutions;
    }
    return res_count >= threshold;
}

bool build_leftout_partitions(int cur_edge_idx) {
    if (cur_edge_idx >= number_of_leftout_edges) {
        return check_3bms();
    }
    for (int i = 0; i < REG; ++i) {
        leftout_partition[cur_edge_idx] = i;
        if (build_leftout_partitions(cur_edge_idx + 1))
            return true;
    }
    return false;
}

bool build_ignored_permutations(int cur_vertex_idx) {
    if (cur_vertex_idx >= number_of_ignored) {
        return build_leftout_partitions(0);
    }
    if (cur_vertex_idx == 0) {
        ignored_permutation[0] = 0;
        return build_ignored_permutations(1);
    }
    for (int i = 0; i < 6; ++i) {
        ignored_permutation[cur_vertex_idx] = i;

        if (build_ignored_permutations(cur_vertex_idx + 1))
            return true;
    }
    return false;
}

bool build_bipartizing_matchings() {
    for (int i = 0; i <= number_of_circuits; ++i) {
        for (int j = 0; j < separate_circuit_lengths[i]; ++j) {
            separate_circuits[i][j + separate_circuit_lengths[i]] = separate_circuits[i][j];
        }
    }

    int cur_count = 0;
    for (int i = 0; i < number_of_vertices; ++i) {
        if (ignored[i]) {
            just_ignored[cur_count] = i;
            ++cur_count;
        }
    }

    number_of_leftout_edges = 0;
    for (int i = 0; i < number_of_edges; ++i) {
        if (!edge_in_circuit[i] && in_circuit[vertices[i][0]] && in_circuit[vertices[i][1]]) {
            leftout_edges[number_of_leftout_edges] = i;
            ++number_of_leftout_edges;
        }
    }

    if (circuit_length + number_of_leftout_edges + REG * number_of_ignored != number_of_edges) {
        PrintGraph(graph, number_of_vertices);
        cerr << "fail: " << circuit_length << " " << number_of_leftout_edges << " " << number_of_ignored << " " << number_of_edges << endl;
        cerr << number_of_ignored + circuit_length << " vs " << number_of_vertices << endl;
        cerr << number_of_circuits << endl;
        for (int i = 0; i <= number_of_circuits; ++i)
            cerr << separate_circuit_lengths[i] << ", ";
        cerr << endl;
        cerr << "circuit: ";
        for (int i = 0; i < circuit_length; ++i)
            cerr << circuit[i] << " ";
        cerr << endl;
        cerr << "ignored: ";
        for (int i = 0; i < number_of_ignored; ++i) {
            int ignored_vertex = just_ignored[i];
            cerr << ignored_vertex << " (";
            for (int j = 0; j < REG; ++j) {
                cerr << graph[ignored_vertex][permutations[ignored_permutation[i]][j]];
                if (j < REG - 1)
                    cerr << ", ";
            }
            cerr << "); ";
        }
        cerr << endl;
        cerr << number_of_leftout_edges << " leftouts: ";
        for (int i = 0; i < number_of_leftout_edges; ++i) {
            cerr << "(" << vertices[leftout_edges[i]][0] << ", " << vertices[leftout_edges[i]][1] << ") -> " << leftout_partition[i] << "; ";
        }
        cerr << endl;
    }
    DEBUGASSERT(circuit_length + number_of_leftout_edges + REG * number_of_ignored == number_of_edges);

    return build_ignored_permutations(0);
}

bool build_dominating_circuit(int cur_vertex) {
    if (cur_vertex == start_vertex[number_of_circuits] && (circuit_length >= number_of_edges / 2) && number_of_ignored + circuit_length == number_of_vertices) {
        return build_bipartizing_matchings();
    }
    if (cur_vertex == start_vertex[number_of_circuits] && separate_circuit_lengths[number_of_circuits] > 0) {
        ++number_of_circuits;
        separate_circuit_lengths[number_of_circuits] = 0;
        for (int i = 0; i < number_of_vertices; ++i) {
            if (!ignored[i] && !in_circuit[i]) {
                // remember start_vertices for each circuit and count circuits
                start_vertex[number_of_circuits] = i;

                if (build_dominating_circuit(start_vertex[number_of_circuits]))
                    return true;

                bool has_ignored_neib = false;
                for (int j = 0; j < REG; ++j) {
                    if (ignored[graph[i][j]]) {
                        has_ignored_neib = true;
                    }
                }
                if (has_ignored_neib)
                    break;

                ignored[i] = true;
                int prev_deg = degree[i];
                degree[i] = 0;
                ++number_of_ignored;
                for (int j = 0; j < REG; ++j) {
                    if (!ignored[graph[i][j]] && !in_circuit[graph[i][j]]) {
                        start_vertex[number_of_circuits] = graph[i][j];
                        if (build_dominating_circuit(start_vertex[number_of_circuits]))
                            return true;
                        break;
                    }
                }
                ignored[i] = false;
                degree[i] = prev_deg;
                --number_of_ignored;
                break;
            }
        }
        --number_of_circuits;
        return false;
    }

    for (int j = 0; j < REG; ++j) {
        int next_vertex = graph[cur_vertex][j];
        if (ignored[next_vertex]) {
            continue;
        }
        if (edge_in_circuit[edge_index[cur_vertex][next_vertex]])
            continue;
        if (!in_circuit[next_vertex]) {
            bool caught_problem = false; // assuming here that we don't have cycles of length 3 in the graph
            if (cur_vertex != start_vertex[number_of_circuits]) {
                for (int k = 0; k < REG; ++k) {
                    int side_vertex = graph[cur_vertex][k];
                    if (k != j && !in_circuit[side_vertex]) {
                        if (degree[side_vertex] == 2) {
                            for (int m = 0; m < REG; ++m) {
                                if (ignored[graph[side_vertex][m]]) {
                                    caught_problem = true;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            if (caught_problem)
                continue;

            in_circuit[next_vertex] = true;
            circuit[circuit_length] = next_vertex;
            separate_circuits[number_of_circuits][separate_circuit_lengths[number_of_circuits]] = next_vertex;
            ++circuit_length;
            ++separate_circuit_lengths[number_of_circuits];
            int ei = edge_index[cur_vertex][next_vertex];
            edge_in_circuit[ei] = true;

            if (cur_vertex != start_vertex[number_of_circuits]) {
                for (int k = 0; k < REG; ++k) {
                    int side_vertex = graph[cur_vertex][k];
                    if (side_vertex == start_vertex[number_of_circuits]) {
                        continue;
                    }
                    if (k != j && !in_circuit[side_vertex]) {
                        --degree[side_vertex];
                        if (degree[side_vertex] == 1) {
                            ignored[side_vertex] = true;
                            ++number_of_ignored;
                        }
                    }
                }
            }

            if (build_dominating_circuit(next_vertex))
                return true;

            // undo
            if (cur_vertex != start_vertex[number_of_circuits]) {
                for (int k = 0; k < REG; ++k) {
                    int side_vertex = graph[cur_vertex][k];
                    if (side_vertex == start_vertex[number_of_circuits]) {
                        continue;
                    }
                    if (k != j && !in_circuit[side_vertex]) {
                        if (degree[side_vertex] == 1) {
                            ignored[side_vertex] = false;
                            --number_of_ignored;
                        }
                        ++degree[side_vertex];
                    }
                }
            }
            edge_in_circuit[ei] = false;
            in_circuit[next_vertex] = false;
            --circuit_length;
            --separate_circuit_lengths[number_of_circuits];
        }
    }
    return false;

}

bool HasThreeBipartiteMatchings() {
    number_of_solutions = 0;
    //Label the edges of the graph
    int edge_label = 0;
    for (int i = 0; i < number_of_vertices; ++i) {
        for (int j = 0; j < REG; ++j) {
            if (i < graph[i][j]) {
                edge_index[i][graph[i][j]] = edge_label;
                edge_index[graph[i][j]][i] = edge_label;
                vertices[edge_label][0] = i;
                vertices[edge_label][1] = graph[i][j];
                ++edge_label;
            }
        }
    }
    if (edge_label != number_of_edges) {
        cerr << "Error: invalid number of edges" << endl;
        exit(1);
    }

    circuit_length = 0;
    separate_circuit_lengths[0] = 0;
    for (int i = 0; i < number_of_vertices; ++i) {
        ignored[i] = false;
        in_circuit[i] = false;
        degree[i] = REG;
    }
    for (int i = 0; i < number_of_edges; ++i) {
        edge_in_circuit[i] = false;
    }
    number_of_ignored = 0;
    number_of_circuits = 0;
    start_vertex[number_of_circuits] = 0;

    if (build_dominating_circuit(start_vertex[number_of_circuits]))
        return true;

    ignored[0] = true;
    degree[0] = 0;
    number_of_ignored = 1;
    start_vertex[number_of_circuits] = graph[0][0];
    return build_dominating_circuit(start_vertex[number_of_circuits]);
}


int main(int argc, char** argv) {
    ParseArgs(argc, argv, number_of_graphs_to_skip);
    ReadGraphs(HasThreeBipartiteMatchings, number_of_graphs_to_skip, number_of_vertices, number_of_edges, graph);
    return(EXIT_SUCCESS);
}
