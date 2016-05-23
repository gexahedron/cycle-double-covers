/*
 * File:   244_flows.cpp
 * Main result: the conjecture of 233-flows is wrong for Petersen graph (and no 234-flows; but it has 333-flows and 244-flows)
 * Program checks existence of 233, 333, 234, 334 and 244 -flows for snarks
 * Actually it looks like 333-flows and 244-flows always exist
 *
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * Created on 10 March 2016, 15:30
 */

#include "util.h"

#include <climits>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <array>

using namespace std;

/********************************** Variables *************************************/

GRAPH graph;
int number_of_vertices;
int number_of_edges;
int number_of_graphs_to_skip = 0;
unsigned int edge_index[MAXN][MAXN];

int vertices[REG * MAXN / 2][2];
int default_trackback_to = REG * MAXN;
int trackback_to;
int vertex_last_edge_index[MAXN];
int partition_num[REG * MAXN / 2];
int partition_degs[MAXN][3];
int number_of_flowed_edges_at_vertex[MAXN];
int edge_flows[REG * MAXN / 2];
int vertex_flows[MAXN];
int stored_edge_flows[MAXN][3];
int max_flow_values[3];

const int number_of_profiles = 5;
int profiles[number_of_profiles][3] = {{2, 3, 3}, {3, 3, 3}, {2, 3, 4}, {3, 3, 4}, {2, 4, 4}};
bool has_profile[number_of_profiles];


/********************************** Functions *************************************/

bool BuildFlow(int partition, int edge_index, int max_flow_value) {
    if (edge_index >= number_of_edges) {
        for (int e = 0; e < number_of_edges; ++e) {
            stored_edge_flows[e][partition] = edge_flows[e];
        }
        return true;
    }
    if (partition_num[edge_index] == partition) {
        return BuildFlow(partition, edge_index + 1, max_flow_value);
    }

    int right_bound = max_flow_value - 1;
    int left_bound = -right_bound;

    int v1 = vertices[edge_index][0];
    int v2 = vertices[edge_index][1];

    if (number_of_flowed_edges_at_vertex[v1] == partition_degs[v1][partition] - 1) {
        right_bound = left_bound = vertex_flows[v1];
    }
    if (number_of_flowed_edges_at_vertex[v2] == partition_degs[v2][partition] - 1) {
        right_bound = left_bound = -vertex_flows[v2];
    }

    for (int flow = left_bound; flow <= right_bound; ++flow) {
        if (flow == 0)
            continue;
        edge_flows[edge_index] = flow;

        vertex_flows[v1] -= flow;
        vertex_flows[v2] += flow;
        ++number_of_flowed_edges_at_vertex[v1];
        ++number_of_flowed_edges_at_vertex[v2];

        bool all_is_good = (
            (vertex_flows[v1] == 0 || number_of_flowed_edges_at_vertex[v1] != partition_degs[v1][partition]) &&
            (vertex_flows[v2] == 0 || number_of_flowed_edges_at_vertex[v2] != partition_degs[v2][partition]) &&
            (abs(vertex_flows[v1]) < max_flow_value || number_of_flowed_edges_at_vertex[v1] != partition_degs[v1][partition] - 1) &&
            (abs(vertex_flows[v2]) < max_flow_value || number_of_flowed_edges_at_vertex[v2] != partition_degs[v2][partition] - 1)
        );
        if (all_is_good && BuildFlow(partition, edge_index + 1, max_flow_value)) {
            return true;
        }

        // undoing
        vertex_flows[v1] += flow;
        vertex_flows[v2] -= flow;
        --number_of_flowed_edges_at_vertex[v1];
        --number_of_flowed_edges_at_vertex[v2];
    }
    return false;
}

bool CheckNowhereZeroness(int partition, int max_flow_value) {
    for (int e = 0; e < number_of_edges; ++e) {
        edge_flows[e] = 0;
    }
    for (int v = 0; v < number_of_vertices; ++v) {
        vertex_flows[v] = 0;
        number_of_flowed_edges_at_vertex[v] = 0;
    }
    return BuildFlow(partition, 0, max_flow_value);
}

bool CheckPartitions() {
    for (int i = 0; i < 3; ++i) {
        for (int v = 0; v < number_of_vertices; ++v) {
            partition_degs[v][i] = 3;
        }
    }
    for (int e = 0; e < number_of_edges; ++e) {
        int partition = partition_num[e];
        --partition_degs[vertices[e][0]][partition];
        --partition_degs[vertices[e][1]][partition];
    }
    for (int i = 0; i < 3; ++i) {
        for (int v = 0; v < number_of_vertices; ++v) {
            if (partition_degs[v][i] == 1) {
                trackback_to = min(trackback_to, vertex_last_edge_index[v]);
            }
        }
    }

    if (max_flow_values[0] == 2) {
        for (int v = 0; v < number_of_vertices; ++v) {
            if (partition_degs[v][0] == 3) {
                trackback_to = min(trackback_to, vertex_last_edge_index[v]);
            }
        }
    }

    if (trackback_to <= number_of_edges) {
        return false;
    }

    if (max_flow_values[0] != 2) {
        return CheckNowhereZeroness(1, max_flow_values[1]) && CheckNowhereZeroness(2, max_flow_values[2]);
    } else {
        return CheckNowhereZeroness(0, max_flow_values[0]) &&
            CheckNowhereZeroness(1, max_flow_values[1]) &&
            CheckNowhereZeroness(2, max_flow_values[2]);
    }
}

bool BuildPartitions(int edge_index) {
    if (edge_index >= number_of_edges)
        return CheckPartitions();
    for (int i = 0; i < 3; ++i) {
        partition_num[edge_index] = i;
        if (BuildPartitions(edge_index + 1)) {
            return true;
        }
        if (trackback_to < edge_index)
            return false;
        trackback_to = default_trackback_to;
    }
    return false;
}

bool HasPartitionedFlows() {
    //Label the edges of the graph
    int edge_label = 0;
    for (int i = 0; i < number_of_vertices; ++i) {
        for (int j = 0; j < REG; ++j) {
            if (i < graph[i][j]) {
                edge_index[i][graph[i][j]] = edge_label;
                edge_index[graph[i][j]][i] = edge_label;
                vertices[edge_label][0] = i;
                vertices[edge_label][1] = graph[i][j];
                vertex_last_edge_index[i] = edge_label;
                vertex_last_edge_index[graph[i][j]] = edge_label;
                ++edge_label;
            }
        }
    }
    if (edge_label != number_of_edges) {
        cerr << "Error: invalid number of edges" << endl;
        exit(1);
    }

    for (int i = 0; i < number_of_profiles; ++i) {
        has_profile[i] = false;
    }

    for (int profile = 0; profile < number_of_profiles; ++profile) {
        for (int i = 0; i < profile; ++i) {
            if (has_profile[i] && profiles[profile][0] >= profiles[i][0] && profiles[profile][1] >= profiles[i][1] && profiles[profile][2] >= profiles[i][2]) {
                has_profile[profile] = true;
                break;
            }
        }
        if (has_profile[profile]) {
            cerr << profiles[profile][0] << profiles[profile][1] << profiles[profile][2] << ":     ; ";
            continue;
        }

        for (int i = 0; i < 3; ++i) {
            max_flow_values[i] = profiles[profile][i];
        }
        trackback_to = default_trackback_to;

        for (int i = 0; i < 3; ++i) {
            for (int v = 0; v < number_of_vertices; ++v) {
                partition_degs[v][i] = 3;
            }
        }

        for (int v = 0; v < number_of_vertices; ++v) {
            vertex_flows[v] = 0;
            number_of_flowed_edges_at_vertex[v] = 0;
        }
        for (int e = 0; e < number_of_edges; ++e) {
            edge_flows[e] = 0;
            partition_num[e] = -1;
        }

        has_profile[profile] = BuildPartitions(0);
        cerr << profiles[profile][0] << profiles[profile][1] << profiles[profile][2] << ": ";
        if (has_profile[profile])
            cerr << "YES!; ";
        else
            cerr << "nope; ";
    }
    cerr << endl;
}

void ParseArgs(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        cerr << "Error: invalid number of arguments" << endl;
        cerr << "Usage: " << argv[0] << " <number_of_vertices>" << endl;
        exit(1);
    } else {
        number_of_vertices = atoi(argv[1]);
        number_of_edges = REG * number_of_vertices / 2;
        if (number_of_vertices > MAXN) {
            cerr << "Number of vertices is too big (limit is " << MAXN << ")" << endl;
            exit(0);
        }

        if (argc >= 3) {
            number_of_graphs_to_skip = atoi(argv[2]);
        }
    }
}

void ReadGraphs(unsigned long long int& number_of_graphs_read) {
    int code_length = number_of_vertices * REG / 2 + number_of_vertices;
    unsigned char code[code_length];
    while (fread(code, sizeof(unsigned char), code_length, stdin)) {
        DecodeMulticode(code, code_length, number_of_vertices, graph);
        ++number_of_graphs_read;
        if (number_of_graphs_to_skip >= number_of_graphs_read) {
            continue;
        }
        cerr << "g" << number_of_graphs_read << "\t";
        HasPartitionedFlows();
    }
}

int main(int argc, char** argv) {
    ParseArgs(argc, argv);
    unsigned long long int number_of_graphs_read = 0;
    ReadGraphs(number_of_graphs_read);
    cerr << "Fin" << endl;
    cerr << "Read " << number_of_graphs_read << " graphs" << endl;
    return(EXIT_SUCCESS);
}
