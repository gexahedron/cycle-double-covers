/*
 * File:   util.h
 * Contains some useful defines and auxiliary methods.
 * Original author: Jan Goedgebeur (Jan.Goedgebeur@UGent.be)
 * Modified by: Nikolay Ulyanov (ulyanick@gmail.com)
 */

#ifndef _UTIL_H
#define _UTIL_H

#pragma once

#include <iostream>
#include <cstdio>

/********************************** Defines *************************************/

#define BIT(i) (1ULL << (i))

#define DEBUGASSERT(assertion) if(!(assertion)) {std::cerr << __FILE__ << ":" << __LINE__ << " Assertion failed: " << #assertion << std::endl; exit(1);}

//const int MAXN = 64; // Should be 32 on a 32-bit machine
const unsigned int MAXN = 42; // Because edge labels are also saved in a bitvector

const unsigned int REG = 3; // Degree of each vertex

typedef unsigned int GRAPH[MAXN][REG]; // adjacency list


/*************************** Useful auxiliary methods ***************************/

inline void AddEdge(GRAPH graph, unsigned int adj[], unsigned int v, unsigned int w) {
    graph[v][adj[v]] = w;
    graph[w][adj[w]] = v;
    ++adj[v];
    ++adj[w];
}

/**
 * Decodes the code (which is in multicode format) of a graph
 */
inline void DecodeMulticode(unsigned char code[], int code_length, int number_of_vertices, GRAPH graph) {
    unsigned int adj[number_of_vertices];

    if (number_of_vertices != code[0]) {
        std::cerr << "Error: Wrong number of vertices: expected " << number_of_vertices << " while found " << code[0] << " vertices" << std::endl;
        exit(1);
    }

    for (int i = 0; i < number_of_vertices; ++i) {
        adj[i] = 0;
    }

    int current_vertex = 1;
    for (int i = 1; i < code_length; ++i) {
        if (code[i] == 0) {
            ++current_vertex;
        } else {
            AddEdge(graph, adj, current_vertex - 1, code[i] - 1);
        }
    }

    for (int i = 0; i < number_of_vertices; ++i) {
        if (adj[i] != REG) {
            std::cerr << "Error: graph is not cubic" << std::endl;
            exit(1);
        }
    }
}

// TODO: DecodeGraph6(unsigned char code[], int code_length, int number_of_vertices, GRAPH graph) {
//}

inline void PrintGraph(GRAPH graph, int number_of_vertices) {
    std::cerr << "Printing graph:" << std::endl;
    for (int i = 0; i < number_of_vertices; ++i) {
        std::cerr << i << " :";
        for (int j = 0; j < REG; ++j) {
            std::cerr << " " << graph[i][j];
        }
        std::cerr << std::endl;
    }
}

#endif  /* _UTIL_H */
