/*
 * File:   chord_diagrams.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * Created on May 7, 2018
 *
 * TODO:
 *
 */

#pragma once

#include "graph.h"
#include "common.h"
#include "experiments/o6c4c.h"

// TODO: remove unused includes
#include <cassert>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <map>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <tuple>
#include <random>

using namespace std;


namespace ExpChordDiagrams {

struct GraphTC3 {
  Graph g;
  Graph tc3;
  int special_vertex;
  int cur2more[MAX_VERTEX];
  int more2cur[MAX_VERTEX];

  GraphTC3()
  {}

  GraphTC3(int n)
    : g(n)
    , tc3(n)
    , special_vertex(NONE)
  {
    fill_n(cur2more, MAX_VERTEX, NONE);
    fill_n(more2cur, MAX_VERTEX, NONE);
  }

  void add_edge(int v1, int v2) {
    assert(special_vertex != NONE);
    g.add_edge(v1, v2);
    if (v1 != special_vertex && v2 != special_vertex) {
      tc3.add_edge(v1, v2);
    }
  }

  void add_edge_from_more(int w1, int w2) {
    int v1 = more2cur[w1];
    int v2 = more2cur[w2];
    add_edge(v1, v2);
  }
};

struct Extra {
  vector<GraphTC3> tc3s;
};

bool hamiltonian(Extra& extra) {
  Graph g = extra.tc3s[0].g;
  ExpCycles::prepare_build_cycle(g);
  if (g.all_full_cycles.size() < 3) {
    return false;
  }
  //cerr << "full cycles: " << g.all_full_cycles.size() << endl;
  Exp6c4c::find_all_o6c4c(g, true);
  if (g.all_o6c4c.size() < 2) {
    return false;
  }
  cerr << "o6c4c count: " << g.all_o6c4c.size() << endl;
  return false;
}

bool is_joinable(Extra& extra) {
  const GraphTC3& more_tc3 = extra.tc3s[1];
  const int more_special_vertex = more_tc3.special_vertex;
  for (int jj = 0; jj < MAX_DEG; ++jj) {
    const int v = more_tc3.g.v2v[more_special_vertex][jj];
    // check whether i can off this vertex
    bool can_off = true;
    vector<int> neibs;
    for (int j = 0; j < MAX_DEG - 1; ++j) {
      const int neib = more_tc3.tc3.v2v[v][j];
      if (more_tc3.tc3.deg[neib] != MAX_DEG) {
        can_off = false;
        break;
      } else {
        neibs.push_back(neib);
      }
    }
    if (!can_off) {
      continue;
    }

    // so, now we can off this vertex
    for (int i = 0; i < 2; ++i) {
      const int neib = neibs[i];
      // let's check neighbours of neib, whether they are not connected to each other
      // if so, we can off the neib vertex, add a new edge and recurse
      vector<int> neib_neibs;
      for (int j = 0; j < MAX_DEG; ++j) {
        int n = more_tc3.g.v2v[neib][j];
        if (n != v) {
          neib_neibs.push_back(n);
        }
      }
      const int n1 = neib_neibs[0];
      const int n2 = neib_neibs[1];
      if (more_tc3.g.has_edge(n1, n2)) {
        continue;
      }
      const int offed_edge = more_tc3.tc3.edge_index[n1][n2];

      // create new graph
      GraphTC3& cur_tc3 = extra.tc3s[0];
      cur_tc3 = GraphTC3(more_tc3.g.number_of_vertices - 2);

      int vertex_count = 0;
      for (int w = 0; w < more_tc3.g.number_of_vertices; ++w) {
        if (w != v && w != neib && w != more_tc3.special_vertex) {
          cur_tc3.cur2more[vertex_count] = w;
          cur_tc3.more2cur[w] = vertex_count;
          ++vertex_count;
        }
      }
      assert(vertex_count == cur_tc3.g.number_of_vertices - 1);
      cur_tc3.special_vertex = vertex_count;

      for (int e = 0; e < more_tc3.tc3.number_of_edges; ++e) {
        if (e == offed_edge) {
          continue;
        }
        const int v1 = more_tc3.tc3.e2v[e][0];
        const int v2 = more_tc3.tc3.e2v[e][1];
        if (v1 == v || v1 == neib || v2 == v || v2 == neib) {
          continue;
        }

        cur_tc3.add_edge_from_more(v1, v2);
      }
      cur_tc3.add_edge_from_more(n1, n2);
      for (int v = 0; v < cur_tc3.g.number_of_vertices; ++v) {
        if (cur_tc3.g.deg[v] == MAX_DEG - 1) {
          cur_tc3.add_edge(v, cur_tc3.special_vertex);
        }
      }
      /*cerr << "graph number: " << graph_number << endl;
      for (int v = 0; v < cur_tc3.g.number_of_vertices; ++v) {
        cerr << cur_tc3.g.deg[v] << " ";
      }
      cerr << endl;*/

      // recurse
      if (hamiltonian(extra)) {
        return true;
      }
    }
  }
  return false;
}

void calc_chord_diagram_invariants(Graph& graph) {
  /*Extra extra;
  extra.tc3s.push_back(GraphTC3(graph.number_of_vertices - 2));
  extra.tc3s.push_back(GraphTC3(graph.number_of_vertices));

  for (int v = 0; v < graph.number_of_vertices; ++v) {
    extra.tc3s.back() = GraphTC3(graph.number_of_vertices);
    extra.tc3s.back().special_vertex = v;
    for (int e = 0; e < graph.number_of_edges; ++e) {
      extra.tc3s.back().add_edge(graph.e2v[e][0], graph.e2v[e][1]);
    }

    // now we have TC3 graph
    // let's check joinability
    bool joinable = is_joinable(extra);
    if (joinable) {
      return;
    }
  }
  //cerr << "graph is not joinable" << endl;
  */

  int n = 14;
  int maxn = 30;
  int deg[maxn];
  bool take_pair[maxn * maxn];
  random_device rd;
  //mt19937 gen(42);
  mt19937 gen(rd());

  while (true) {
    vector<pair<int, int>> pii_edges;

    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        pii_edges.push_back(make_pair(i, j));
      }
    }

    bool is_regular = false;
    while (!is_regular) {
      shuffle(pii_edges.begin(), pii_edges.end(), gen);
      for (int i = 0; i < n; ++i) {
        deg[i] = 0;
      }
      for (int i = 0; i < (int) pii_edges.size(); ++i) {
        take_pair[i] = false;
        int v1, v2;
        v1 = pii_edges[i].first;
        v2 = pii_edges[i].second;
        if ((deg[v1] < 3) && (deg[v2] < 3)) {
          take_pair[i] = true;
          ++deg[v1];
          ++deg[v2];
        }
      }
      is_regular = true;
      for (int i = 0; i < n; ++i) {
        is_regular = is_regular && (deg[i] == 3);
        if (!is_regular) {
          break;
        }
      }
    }

    Graph g(n, true);
    for (int i = 0; i < (int) pii_edges.size(); ++i) {
      if (take_pair[i]) {
        int v1, v2;
        v1 = pii_edges[i].first;
        v2 = pii_edges[i].second;
        g.add_edge(v1, v2);
      }
    }

    ExpCycles::prepare_build_cycle(g, true);
    if (g.all_full_cycles.size() < 3) {
      continue;
    }
    /*if (g.all_full_cycles.size() != 18) {
      continue;
    }*/

    Exp6c4c::find_all_o6c4c(g, true);
    if (g.all_o6c4c.size() < 2) {
      continue;
    }

    cerr << "full cycles: " << g.all_full_cycles.size() << endl;
    cerr << "o6c4c count: " << g.all_o6c4c.size() << endl;

    map<Mask, int> intersection_by_cycle;
    map<Mask, vector<int>> pairs_by_cycle;
    for (const auto& c : g.all_full_cycles) {
      int vertex_order[maxn];
      vector<int> pairs;
      int where[maxn];
      for (int i = 0; i < n; ++i) {
        vertex_order[i] = -1;
        pairs.push_back(-1);
        where[i] = -1;
      }
      vertex_order[0] = 0;
      where[0] = 0;
      for (int idx = 0; idx < n - 1; ++idx) {
        int cur = vertex_order[idx];
        vertex_order[idx + 1] = -1;
        int prev = -1;
        if (idx > 0) {
          prev = vertex_order[idx - 1];
        }
        for (int j = 0; j < 3; ++j) {
          int v = g.v2v[cur][j];
          int e = g.v2e[cur][j];
          if (BIT(e) & c) {
            if (v != prev) {
              vertex_order[idx + 1] = v;
              where[v] = idx + 1;
            }
          } else {
            pairs[cur] = v;
            pairs[v] = cur;
          }
        }
      }

      for (int i = 0; i < n; ++i) {
        assert(vertex_order[i] != -1);
        assert(where[i] != -1);
        assert(pairs[i] != -1);
        assert(pairs[i] != i);
      }

      int cur_intersection = 0;
      for (int i = 0; i < n; ++i) {
        int i2 = pairs[i];
        if (where[i2] < where[i]) {
          continue;
        }
        for (int j = i + 1; j < n; ++j) {
          int j2 = pairs[j];
          if (where[j2] < where[j]) {
            continue;
          }

          if (where[i] < where[j] && where[j] < where[i2] && where[i2] < where[j2]) {
            ++cur_intersection;
          } else if (where[j] < where[i] && where[i] < where[j2] && where[j2] < where[i2]) {
            ++cur_intersection;
          }
        }
      }
      intersection_by_cycle[c] = cur_intersection;
      vector<int> pairs_in_cycle;
      for (int i = 0; i < n; ++i) {
        pairs_in_cycle.push_back(-1);
      }
      for (int i = 0; i < n; ++i) {
        pairs_in_cycle[i] = where[pairs[vertex_order[i]]];
      }
      pairs_by_cycle[c] = pairs_in_cycle;
    }

    // compare cycles, just for fun
    // to find some 4-term relations
    int rel0_count = 0;
    int rel2_count = 0;
    for (const auto& c1 : g.all_full_cycles) {
      for (const auto& c2 : g.all_full_cycles) {
        if (c1 > c2) {
          continue;
        }
        int min_add = 0;
        if (c1 == c2) {
          min_add = 1;
        }
        // check all rotations for c2
        for (int add = min_add; add < n - 1; ++add) {
          int rot_pairs[maxn];
          for (int i = 0; i < n; ++i) {
            int rot_i = (i + add) % n;
            int rot_pairs_i = (pairs_by_cycle[c2][i] + add) % n;
            rot_pairs[rot_i] = rot_pairs_i;
          }

          // compare wheres[c1] and rot_where
          int diff_count = 0;
          vector<int> different;
          for (int i = 0; i < n; ++i) {
            if (pairs_by_cycle[c1][i] < i) {
              continue;
            }
            if (pairs_by_cycle[c1][i] != rot_pairs[i]) {
              ++diff_count;
              different.push_back(i);
              different.push_back(pairs_by_cycle[c1][i]);
            }
          }
          if (diff_count > 2) {
            continue;
          }
          if (diff_count == 2) {
            cerr << "differents: " << different[0] << " " << different[1] << "; " << different[2] << " " << different[3] << endl;
            sort(different.begin(), different.end());
            different.push_back(different[0] + n);
            bool good_diff = false;
            for (int idx = 0; idx < different.size() - 1; ++idx) {
              if (different[idx + 1] - different[idx] == 1) {
                good_diff = true;
                break;
              }
            }
            if (!good_diff) {
              continue;
            }
          }
          if (diff_count == 0) {
            ++rel0_count;
            cerr << "rel0: " << c1 << " " << c2 << endl;
          } else {
            cerr << "rel2: " << c1 << " " << c2 << endl;
            ++rel2_count;
          }
        }
      }
    }
    cerr << "rel0 count: " << rel0_count << "; rel2 count: " << rel2_count << endl;

    // count intersection number
    cerr << "intersections: " << endl;
    for (const auto& sol : g.all_o6c4c) {
      int intersection = 0;
      int idx = 0;
      for (const auto& c : sol) {
        intersection += intersection_by_cycle[c];
        cerr << intersection_by_cycle[c];
        if (idx < 5) {
          cerr << " + ";
        }
        ++idx;
      }
      cerr << " = " << intersection << endl;
      cerr << "cycles: ";
      for (const auto& c : sol) {
        cerr << c << " ";
      }
      cerr << endl;
      // TODO: print cycles
    }
    cerr << endl;
    g.print();
    break;
  }
}

} // ExpChordDiagrams
