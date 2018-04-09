/*
 * File:   graph.cpp
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 */

#include "graph.h"
#include <iostream>
#include <cassert>

using namespace std;

TGraph::TGraph(unsigned int n) {
  number_of_vertices = n;
  number_of_edges = 0;
  for (int v = 0; v < number_of_vertices; ++v) {
    deg[v] = 0;
    lexi_deg[v] = 0;
    for (int u = 0; u < number_of_vertices; ++u) {
      edge_index[v][u] = NONE;
    }
  }
}

bool TGraph::has_edge(unsigned int v, unsigned int w) const {
  return (edge_index[v][w] != NONE);
}

void TGraph::add_edge(unsigned int v, unsigned int w) {
  assert(!has_edge(v, w));
  if (v > w) {
    swap(v, w);
  }
  int e = number_of_edges;
  v2v[v][deg[v]] = w;
  v2v[w][deg[w]] = v;
  v2e[v][deg[v]] = e;
  v2e[w][deg[w]] = e;
  ++deg[v];
  ++deg[w];

  lexi_v2v[v][lexi_deg[v]] = w;
  ++lexi_deg[v];

  e2v[e][0] = v;
  e2v[e][1] = w;
  edge_index[v][w] = e;
  edge_index[w][v] = e;
  ++number_of_edges;
}

void TGraph::finish_init() {
  oddness = MAXN;

  for (int v = 0; v < number_of_vertices; ++v) {
    if (deg[v] != REG) {
      if (REG == 3) {
        cerr << "Error: graph is not cubic" << endl;
      } else {
        cerr << "Error: graph is not regular (where vertex degree = " << REG << ")" << endl;
      }
      exit(1);
    }
  }

  if (number_of_edges != number_of_vertices * REG / 2) {
    cerr << "Error: invalid number of edges" << endl;
    exit(1);
  }

  find_faster_edge_order();
}

bool decode_multicode(FILE* input, TGraph& graph) {
  unsigned char code_for_number_of_vertices[1];
  if (!fread(code_for_number_of_vertices, sizeof(unsigned char), 1, input)) {
    return false;
  }

  graph = TGraph(code_for_number_of_vertices[0]);
  const int edge_count = graph.number_of_vertices * REG / 2;
  const int code_length = edge_count + graph.number_of_vertices - 1;
  unsigned char code[code_length];
  if (!fread(code, sizeof(unsigned char), code_length, input)) {
    return false;
  }

  int cur_v = 1;
  // parsing code_length
  for (int i = 0; i < code_length; ++i) {
    if (code[i] == 0) {
      ++cur_v;
    } else {
      graph.add_edge(cur_v - 1, code[i] - 1);
    }
  }

  graph.finish_init();
  return true;
}

void TGraph::print() const {
  cerr << "Printing graph:" << endl;
  for (int i = 0; i < number_of_vertices; ++i) {
    cerr << i << ":\t";
    for (int j = 0; j < REG; ++j) {
      cerr << v2v[i][j] << "(e" << v2e[i][j] << ")\t";
    }
    cerr << endl;
  }
  if (!all_dominating_circuits.empty()) {
    cerr << "dominating cycle: ";
    for (const auto& c : all_dominating_circuits) {
      for (int e = 0; e < number_of_edges; ++e) {
        if (BIT(e) & c) {
          cerr << e << "(" << e2v[e][0] << "," << e2v[e][1] << ") ";
        }
      }
      cerr << endl;
      break;
    }
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

void TGraph::find_faster_edge_order() {
  int best_edge_ans = 0;
  int best_edge_idx = 0;
  for (int ei = 0; ei < number_of_edges; ++ei) {
    int cur_ans = ::find_faster_edge_order(ei, *this);
    if (cur_ans > best_edge_ans) {
      best_edge_ans = cur_ans;
      best_edge_idx = ei;
    }
  }
  ::find_faster_edge_order(best_edge_idx, *this);
}
