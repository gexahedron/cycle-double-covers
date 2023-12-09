/*
 * File:   graph.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * Contains some useful functions about graph parsing and construction and some auxilary functions
 *
 */

#pragma once

#include "constants.h"
#include "experiments/cycles_data.h"
#include "experiments/petersen_coloring_data.h"
#include "experiments/preimages_data.h"
#include "experiments/o5cdc_data.h"
#include "experiments/o6c4c_data.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <set>

using namespace std;


/****************************************** Definitions *******************************************/
#define BIT(i) (1ULL << (i))

constexpr int MAX_DEG = 3;
// why 42: because edge labels are sometimes saved as a bitvector
// we have a restriction on MAX_EDGE < 64
// so we get this number: 42 * 3 / 2 = 63 < 64
constexpr int MAX_VERTEX = 42;
constexpr int MAX_EDGE = MAX_VERTEX * MAX_DEG / 2;

// we always assume that the graph is bridgeless
struct Graph:
  public ExpCycles::Data,
  public ExpPetersenColouring::Data,
  public ExpPreimages::Data,
  public Exp5cdc::Data,
  public Exp6c4c::Data
{
  int number_of_vertices = 0; // TODO: -> vertex_count
  int number_of_edges = 0; // TODO: -> edge_count

  int number = 0;
  bool is_definitely_snark = true; // TODO: it's always true right now

  unsigned int deg[MAX_VERTEX];
  unsigned int v2v[MAX_VERTEX][MAX_DEG];
  unsigned int v2e[MAX_VERTEX][MAX_DEG];
  unsigned int e2v[MAX_EDGE][2];
  unsigned int e2e[MAX_EDGE][2 * MAX_DEG - 2];
  unsigned int vv2e[MAX_VERTEX][MAX_VERTEX];

  unsigned int lexi_deg[MAX_VERTEX];
  unsigned int lexi_v2v[MAX_VERTEX][MAX_DEG];

  bool edges_are_neibs[MAX_EDGE][MAX_EDGE];

  vector<int> faster_edge_order;

  Graph(unsigned int n = 0, bool is_definitely_snark = true);
  void add_edge(unsigned int v, unsigned int w);
  bool has_edge(unsigned int v, unsigned int w) const;
  void finish_init();
  void print() const;
  void print(const Mask& dominating_circuit) const;
  void print(set<int>& oriented_vertices, set<int>& poor_edges) const;
  void find_faster_edge_order();
};

/*********************************Methods*********************************/

// Decodes the code (which is in multicode format) of a graph.
bool decode_multicode(FILE* input, Graph& graph);

// Decodes adjacency list of edges
bool decode_adjacency(istream& input, Graph& graph);

// Decodes format from papers by Brinkmann, Goedgebeur, Hägglund, Markström
bool decode_bghm(istream& input, Graph& graph);

// Decodes graph6 format
bool decode_graph6(istream& input, Graph& graph);

inline int combine_hash_from_mask_and_edge_colour(int mask, int edge_colour) {
  return mask + edge_colour * BIT(6);
}

inline Mask inv(Graph& graph, Mask mask) {
  return BIT(graph.number_of_edges) - 1 - mask;
}

void build_petersen_graph(Graph& graph);
