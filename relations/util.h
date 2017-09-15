/*
 * File:   util.h
 *
 * Contains some useful defines and auxiliary methods.
 *
 *
 */

#pragma once

#include <iostream>
#include <cstdio>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>

using namespace std;

/**********************************Defines*************************************/
using TMask = unsigned long long int;

#define BIT(i) (1ULL << (i))

//#define MAXN 64 /* Should be 32 on a 32-bit machine */
#define MAXN 42 //Because edge labels are also saved in a bitvector
#define REG 3

#define DEBUGASSERT(assertion) if(!(assertion)) {cerr << __FILE__ << ":" << __LINE__ << " Assertion failed: " << #assertion << endl; exit(1);}

struct THoffmanOstenhofDecomposition {
    unordered_set<int> tree;
    unordered_set<int> cycle;
    unordered_set<int> matching;
};

struct TGraph {
    int number_of_vertices = 0;
    int number_of_edges = 0;
    unsigned int deg[MAXN];
    //unsigned int oriented_deg[MAXN];
    unsigned int v2v[MAXN][REG];
    //unsigned int oriented_v2v[MAXN][REG];
    unsigned int v2e[MAXN][REG];
    unsigned int e2v[REG * MAXN / 2][2];
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
    unordered_set<TMask> all_full_cycles;
    unordered_set<TMask> all_dominating_circuits;
    unordered_map<TMask, vector<vector<int>>> cycles_as_circuits;

    // cycle covers
    set<set<TMask>> all_5cdc;
    set<set<TMask>> all_o5cdc;
    set<set<TMask>> all_6c4c;
    set<set<TMask>> all_o6c4c;

    set<set<TMask>> petersen_5cdc;
    set<set<TMask>> petersen_6c4c;

    // hoffman-ostenhof
    vector<THoffmanOstenhofDecomposition> hoffman_ostenhof_solutions;
};

/***************************Useful auxiliary methods***************************/


void add_edge(TGraph& graph, unsigned int e, unsigned int v, unsigned int w) {
    if (v > w) {
        swap(v, w);
    }
    graph.v2v[v][graph.deg[v]] = w;
    graph.v2v[w][graph.deg[w]] = v;
    //graph.oriented_v2v[v][graph.oriented_deg[v]] = w;
    //++graph.oriented_deg[v];
    graph.v2e[v][graph.deg[v]] = e;
    graph.v2e[w][graph.deg[w]] = e;
    ++graph.deg[v];
    ++graph.deg[w];

    graph.e2v[e][0] = v;
    graph.e2v[e][1] = w;
    graph.edge_index[v][w] = e;
    graph.edge_index[w][v] = e;
}

/**
 * Decodes the code (which is in multicode format) of a graph.
 */
void decode_multicode(unsigned char code[], int codelength, TGraph& graph) {
    if (graph.number_of_vertices != code[0]) {
        cerr << "Error: Wrong number of vertices: expected " << graph.number_of_vertices << " while found " <<  code[0] << " vertices" << endl;
        exit(1);
    }

    for (int v = 0; v < graph.number_of_vertices; ++v) {
        graph.deg[v] = 0;
        //graph.oriented_deg[v] = 0;
    }
    unsigned int edge_number = 0;

    int cur_v = 1;
    // parsing codelength
    for (int i = 1; i < codelength; ++i) {
        if (code[i] == 0) {
            ++cur_v;
        } else {
            add_edge(graph, edge_number, cur_v - 1, code[i] - 1);
            ++edge_number;
        }
    }

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

