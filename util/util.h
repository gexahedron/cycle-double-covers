/*
 * File:   util.h
 * Contains some useful defines and auxiliary methods.
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 */

#ifndef _UTIL_H
#define _UTIL_H

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

void AddEdge(GRAPH graph, unsigned int adj[], unsigned int v, unsigned int w);

/**
 * Decodes the code (which is in multicode format) of a graph
 */
void DecodeMulticode(unsigned char code[], int code_length, int number_of_vertices, GRAPH graph);

// TODO: void DecodeGraph6(unsigned char code[], int code_length, int number_of_vertices, GRAPH graph);

void PrintGraph(GRAPH graph, int number_of_vertices);

#endif  /* _UTIL_H */
