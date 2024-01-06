/*
 * File:   util.h
 *
 * Contains some useful defines and auxiliary methods.
 *
 * Author: Jan Goedgebeur (Jan.Goedgebeur@UGent.be)
 *
 */

#ifndef _UTIL_H
#define	_UTIL_H

#include <iostream>
#include <cstdio>

/**********************************Defines*************************************/

#define BIT(i) (1ULL << (i))
#define MAX(a, b) (((a) >= (b)) ? (a) : (b))

//#define MAXN 64 /* Should be 32 on a 32-bit machine */
#define MAXN 42 //Because edge labels are also saved in a bitvector

#define REG 3

#define DEBUGASSERT(assertion) if(!(assertion)) {fprintf(stderr, "%s:%u Assertion failed: %s\n", __FILE__, __LINE__, #assertion); exit(1);}

typedef unsigned int GRAPH[MAXN][REG];


int number_of_vertices = 0;
int number_of_edges = 0;


/***************************Useful auxiliary methods***************************/


void add_edge(GRAPH graph, unsigned int adj[], unsigned int v, unsigned int w) {
    graph[v][adj[v]] = w;
    graph[w][adj[w]] = v;
    ++adj[v];
    ++adj[w];
}

/**
 * Decodes the code (which is in multicode format) of a graph.
 */
void decode_multicode(unsigned char code[], int codelength, GRAPH graph) {
    unsigned int adj[number_of_vertices];

    if (number_of_vertices != code[0]) {
        fprintf(stderr, "Error: Wrong number of vertices: expected %d while found %d vertices \n", number_of_vertices, code[0]);
        exit(1);
    }
    //graph[0][0] = number_of_vertices;

    for (int i = 0; i < number_of_vertices; ++i) {
        adj[i] = 0;
    }

    int currentvertex = 1;
    //Codelength
    for (int i = 1; i < codelength; ++i) {
        if (code[i] == 0) {
            ++currentvertex;
        } else {
            add_edge(graph, adj, currentvertex - 1, code[i] - 1);
        }
    }

    for (int i = 0; i < number_of_vertices; ++i) {
        if (adj[i] != REG) {
            fprintf(stderr, "Error: graph is not cubic \n");
            exit(1);
        }
    }

}


void printgraph(GRAPH graph) {
    fprintf(stderr, "Printing graph: \n");
    for (int i = 0; i < number_of_vertices; ++i) {
        fprintf(stderr, "%d :", i);
        for (int j = 0; j < REG; ++j)
            fprintf(stderr, " %d", graph[i][j]);
        fprintf(stderr, "\n");
    }
}


#endif	/* _UTIL_H */
