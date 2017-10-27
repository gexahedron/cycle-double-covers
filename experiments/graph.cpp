/*
 * File:   graph.cpp
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 */

#include "graph.h"
#include <iostream>

using namespace std;

void add_edge(TGraph& graph, unsigned int e, unsigned int v, unsigned int w) {
    if (v > w) {
        swap(v, w);
    }
    graph.v2v[v][graph.deg[v]] = w;
    graph.v2v[w][graph.deg[w]] = v;
    graph.v2e[v][graph.deg[v]] = e;
    graph.v2e[w][graph.deg[w]] = e;
    ++graph.deg[v];
    ++graph.deg[w];

    graph.lexi_v2v[v][graph.lexi_deg[v]] = w;
    ++graph.lexi_deg[v];

    graph.e2v[e][0] = v;
    graph.e2v[e][1] = w;
    graph.edge_index[v][w] = e;
    graph.edge_index[w][v] = e;
}

bool decode_multicode(FILE* input, TGraph& graph) {
    unsigned char code_for_number_of_vertices[1];
    if (!fread(code_for_number_of_vertices, sizeof(unsigned char), 1, input)) {
        return false;
    }

    graph.number_of_vertices = code_for_number_of_vertices[0];
    graph.number_of_edges = graph.number_of_vertices * REG / 2;
    int codelength = graph.number_of_edges + graph.number_of_vertices - 1;
    unsigned char code[codelength];
    if (!fread(code, sizeof(unsigned char), codelength, input)) {
        return false;
    }

    for (int v = 0; v < graph.number_of_vertices; ++v) {
        graph.deg[v] = 0;
        graph.lexi_deg[v] = 0;
    }
    unsigned int edge_number = 0;

    int cur_v = 1;
    // parsing codelength
    for (int i = 0; i < codelength; ++i) {
        if (code[i] == 0) {
            ++cur_v;
        } else {
            add_edge(graph, edge_number, cur_v - 1, code[i] - 1);
            ++edge_number;
        }
    }

    graph.oddness = MAXN;

    for (int v = 0; v < graph.number_of_vertices; ++v) {
        if (graph.deg[v] != REG) {
            if (REG == 3) {
                cerr << "Error: graph is not cubic" << endl;
            } else {
                cerr << "Error: graph is not regular (where vertex degree = " << REG << ")" << endl;
            }
            exit(1);
        }
    }

    if (edge_number != graph.number_of_edges) {
        cerr << "Error: invalid number of edges" << endl;
        exit(1);
    }

	find_faster_edge_order(graph);
    return true;
}

void print_graph(const TGraph& graph) {
    cerr << "Printing graph:" << endl;
    for (int i = 0; i < graph.number_of_vertices; ++i) {
        cerr << i << ":\t";
        for (int j = 0; j < REG; ++j) {
            cerr << graph.v2v[i][j] << "(e" << graph.v2e[i][j] << ")\t";
        }
        cerr << endl;
    }
}

int find_faster_edge_order(int start_edge, TGraph& graph) {
    // TODO: combine visited_edges, faster_edge_order and width into some combined structure
    graph.faster_edge_order.clear();
    graph.faster_edge_order.push_back(start_edge);
    unordered_set<int> visited_edges;
    visited_edges.insert(start_edge);
    int width[REG * MAXN / 2];
    width[0] = 0;
    for (int i = 0; i < graph.number_of_edges; ++i) {
        int v1 = graph.e2v[graph.faster_edge_order[i]][0];
        int v2 = graph.e2v[graph.faster_edge_order[i]][1];
        if (v1 > v2) {
            swap(v1, v2);
        }
        for (int j = 0; j < REG; ++j) {
            int ei = graph.v2e[v1][j];
            if (visited_edges.find(ei) == visited_edges.end()) {
                graph.faster_edge_order.push_back(ei);
                width[visited_edges.size()] = width[i] + 1;
                visited_edges.insert(ei);
            }
        }
        for (int j = 0; j < REG; ++j) {
            int ei = graph.v2e[v2][j];
            if (visited_edges.find(ei) == visited_edges.end()) {
                graph.faster_edge_order.push_back(ei);
                width[visited_edges.size()] = width[i] + 1;
                visited_edges.insert(ei);
            }
        }
    }
    return width[graph.number_of_edges - 1];
}

void find_faster_edge_order(TGraph& graph) {
    int best_edge_ans = 0;
    int best_edge_idx = 0;
    for (int ei = 0; ei < graph.number_of_edges; ++ei) {
        int cur_ans = find_faster_edge_order(ei, graph);
        if (cur_ans > best_edge_ans) {
            best_edge_ans = cur_ans;
            best_edge_idx = ei;
        }
    }
    find_faster_edge_order(best_edge_idx, graph);
}
