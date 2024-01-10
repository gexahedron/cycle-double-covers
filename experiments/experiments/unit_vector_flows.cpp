/*
 * File:   unit_vector_flows.cpp
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
 *
 * Created on August 17, 2019
 *
 */

#include "unit_vector_flows.h"

#include <cassert>
#include <math.h>
#include <map>


namespace ExpUnitVectorFlows {

const double EPS = 1e-8;
const vector<double> zero_vec{0.0, 0.0, 0.0};

int edge_unit_vector_idx[MAX_EDGE];
vector<int> vertex_flowed_edges[MAX_EDGE];

vector<vector<double>> all_unit_vectors;
map<int, int> ops;
map<vector<int>, int> triples;

/*********************************Methods*********************************/

bool same(double f1, double f2) {
  return abs(f1 - f2) < EPS;
}

double dist(const vector<double>& v1, const vector<double>& v2) {
  return sqrt(
      (v1[0] - v2[0]) * (v1[0] - v2[0]) +
      (v1[1] - v2[1]) * (v1[1] - v2[1]) +
      (v1[2] - v2[2]) * (v1[2] - v2[2]));
}

double norm(const vector<double>& v) {
  return dist(v, zero_vec);
}

vector<vector<double>> gen_even_permutations(const vector<vector<double>>& vectors) {
  vector<vector<double>> permutations;
  for (const auto& v : vectors) {
    permutations.push_back(v);
    permutations.push_back({v[1], v[2], v[0]});
    permutations.push_back({v[2], v[0], v[1]});
  }
  return permutations;
}

vector<vector<double>> gen_all_permutations(const vector<vector<double>>& vectors) {
  vector<vector<double>> permutations;
  for (const auto& v : vectors) {
    vector<vector<double>> even = gen_even_permutations({v});
    vector<vector<double>> odd = gen_even_permutations({{v[1], v[0], v[2]}});
    for (const auto& v2 : even) {
      permutations.push_back(v2);
    }
    for (const auto& v2 : odd) {
      permutations.push_back(v2);
    }
  }
  return permutations;
}

vector<vector<double>> gen_plus_minus_variations(const vector<vector<double>>& vectors) {
  vector<vector<double>> permutations;
  for (const auto& v : vectors) {
    vector<vector<int>> multipliers;
    for (const auto& p : v) {
      if (same(p, 0)) {
        multipliers.push_back({1});
      } else {
        multipliers.push_back({1, -1});
      }
    }
    for (const auto& m0 : multipliers[0]) {
      for (const auto& m1 : multipliers[1]) {
        for (const auto& m2 : multipliers[2]) {
          permutations.push_back({v[0] * m0, v[1] * m1, v[2] * m2});
        }
      }
    }
  }
  return permutations;
}

vector<vector<double>> gen_all_plus_minus_permutations(const vector<vector<double>>& vectors) {
  return gen_plus_minus_variations(gen_all_permutations(vectors));
}

vector<vector<double>> gen_unit_vectors() {
  vector<vector<double>> unit_vectors;
  double phi = (1.0 + sqrt(5)) / 2;
  int w = 3;
  for (int a1 = -w; a1 <= w; ++a1) {
    for (int b1 = -w; b1 <= w; ++b1) {
      for (int a2 = a1; a2 <= w; ++a2) {
        for (int b2 = -w; b2 <= w; ++b2) {
          if (a2 == a1 && b2 < b1) {
            continue;
          }
          for (int a3 = a2; a3 <= w; ++a3) {
            for (int b3 = -w; b3 <= w; ++b3) {
              if (a3 == a2 && b3 < b2) {
                continue;
              }
              int s1 = a1 * a1 + b1 * b1 + a2 * a2 + b2 * b2 + a3 * a3 + b3 * b3;
              int s2 = b1 * (2 * a1 + b1) + b2 * (2 * a2 + b2) + b3 * (2 * a3 + b3);
              if (s1 != s2) {
                continue;
              }
              int root = static_cast<int>(round(sqrt(s1)));
              if (root == 0 || (root * root != s1)) {
                continue;
              }
              double p1 = (static_cast<double>(a1) + b1 * phi) / root;
              double p2 = (static_cast<double>(a2) + b2 * phi) / root;
              double p3 = (static_cast<double>(a3) + b3 * phi) / root;
              vector<vector<double>> new_vectors = gen_all_plus_minus_permutations({{p1, p2, p3}});
              for (const auto& v : new_vectors) {
                double v_norm = norm(v);
                unit_vectors.push_back({v[0] / v_norm, v[1] / v_norm, v[2] / v_norm});
              }
            }
          }
        }
      }
    }
  }
  return unit_vectors;
}

void find_vector_relations(const vector<vector<double>>& unit_vectors) {
  vector<vector<double>> tmp_unit_vectors;
  for (size_t i = 0; i < unit_vectors.size(); ++i) {
    bool found_duplicate = false;
    for (size_t j = i + 1; j < unit_vectors.size(); ++j) {
      if (same(dist(unit_vectors[i], unit_vectors[j]), 0)) {
        found_duplicate = true;
        break;
      }
    }
    if (!found_duplicate) {
      tmp_unit_vectors.push_back(unit_vectors[i]);
    }
  }

  vector<bool> has_triple(tmp_unit_vectors.size(), false);
  for (size_t i = 0; i < tmp_unit_vectors.size(); ++i) {
    const auto v0 = tmp_unit_vectors[i];
    for (size_t j = i + 1; j < tmp_unit_vectors.size(); ++j) {
      const auto v1 = tmp_unit_vectors[j];
      for (size_t k = j + 1; k < tmp_unit_vectors.size(); ++k) {
        const auto v2 = tmp_unit_vectors[k];
        if (same(v0[0] + v1[0] + v2[0], 0) &&
            same(v0[1] + v1[1] + v2[1], 0) &&
            same(v0[2] + v1[2] + v2[2], 0)) {
          has_triple[i] = true;
          has_triple[j] = true;
          has_triple[k] = true;
        }
      }
    }
  }
  all_unit_vectors.clear();
  for (size_t i = 0; i < tmp_unit_vectors.size(); ++i) {
    if (has_triple[i]) {
      all_unit_vectors.push_back(tmp_unit_vectors[i]);
    }
  }
  cerr << "unit vectors count: " << all_unit_vectors.size() << endl;

  ops.clear();
  for (int i = 0; i < all_unit_vectors.size(); ++i) {
    if (ops.find(i) != ops.end()) {
      continue;
    }
    const auto v0 = all_unit_vectors[i];
    for (int j = i + 1; j < all_unit_vectors.size(); ++j) {
      const auto v1 = all_unit_vectors[j];
      if (same(v0[0] + v1[0], 0) && same(v0[1] + v1[1], 0) && same(v0[2] + v1[2], 0)) {
        ops[i] = j;
        ops[j] = i;
        break;
      }
    }
    assert(ops.find(i) != ops.end());
  }

  triples.clear();
  for (int i = 0; i < all_unit_vectors.size(); ++i) {
    const auto v0 = all_unit_vectors[i];
    bool has_triple = false;
    for (int j = i + 1; j < all_unit_vectors.size(); ++j) {
      const auto v1 = all_unit_vectors[j];
      for (int k = j + 1; k < all_unit_vectors.size(); ++k) {
        const auto v2 = all_unit_vectors[k];
        if (same(v0[0] + v1[0] + v2[0], 0) &&
            same(v0[1] + v1[1] + v2[1], 0) &&
            same(v0[2] + v1[2] + v2[2], 0)) {
          triples[{i, j}] = k;
          triples[{j, i}] = k;
          triples[{i, k}] = j;
          triples[{k, i}] = j;
          triples[{j, k}] = i;
          triples[{k, j}] = i;
        }
      }
    }
  }
  cerr << "unit vectors triples count: " << triples.size() / 6 << endl;
}

bool gen_unit_vector_flows(Graph& graph, int cur_edge_idx) {
  if (cur_edge_idx == graph.number_of_edges) {
    vector<int> vectors;
    for (int e = 0; e < graph.number_of_edges; ++e) {
      vectors.push_back(edge_unit_vector_idx[e]);
    }
    graph.unit_vector_flows.push_back(vectors);
    return true;
  }

  int cur_edge = graph.faster_edge_order[cur_edge_idx];

  int v1 = graph.e2v[cur_edge][0];
  int v2 = graph.e2v[cur_edge][1];
  int flowed1 = vertex_flowed_edges[v1].size();
  int flowed2 = vertex_flowed_edges[v2].size();
  vector<int> flow_vector_indices;
  if (flowed1 == MAX_DEG - 1) {
    flow_vector_indices.push_back(triples[vertex_flowed_edges[v1]]);
  }
  if (flowed2 == MAX_DEG - 1) {
    int other_flow_idx = ops[triples[vertex_flowed_edges[v2]]];
    if (flowed1 == MAX_DEG - 1) {
      if (other_flow_idx != flow_vector_indices[0]) {
        return false;
      }
    } else {
      flow_vector_indices.push_back(other_flow_idx);
    }
  }
  if (flowed1 != MAX_DEG - 1 && flowed2 != MAX_DEG - 1) {
    for (size_t i = 0; i < all_unit_vectors.size(); ++i) {
      flow_vector_indices.push_back(i);
    }
  }

  for (const auto& v_idx : flow_vector_indices) {
    edge_unit_vector_idx[cur_edge] = v_idx;
    int op_idx = ops[v_idx];
    vertex_flowed_edges[v1].push_back(v_idx);
    vertex_flowed_edges[v2].push_back(op_idx);
    bool can_flow_v1 = (vertex_flowed_edges[v1].size() != 2) ||
        (triples.find(vertex_flowed_edges[v1]) != triples.end());
    bool can_flow_v2 = (vertex_flowed_edges[v2].size() != 2) ||
        (triples.find(vertex_flowed_edges[v2]) != triples.end());
    if (can_flow_v1 && can_flow_v2 && gen_unit_vector_flows(graph, cur_edge_idx + 1)) {
      return true;
    }
    vertex_flowed_edges[v1].pop_back();
    vertex_flowed_edges[v2].pop_back();
  }
  edge_unit_vector_idx[cur_edge] = NONE;
  return false;
}

void find_all_unit_vector_flows(Graph& graph) {
  if (triples.size() == 0) {
    vector<vector<double>> unit_vectors = gen_unit_vectors();
    find_vector_relations(unit_vectors);
  }

  for (int e = 0; e < graph.number_of_edges; ++e) {
    edge_unit_vector_idx[e] = NONE;
  }

  for (int v = 0; v < graph.number_of_vertices; ++v) {
    vertex_flowed_edges[v].clear();
  }

  gen_unit_vector_flows(graph, 0);
  cerr << "number of unit vector flows: " << graph.unit_vector_flows.size() << endl;
}




} // ExpUnitVectorFlows
