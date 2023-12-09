/*
 * File:   tc3_joining.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * Created on February 8, 2018
 *
 * TODO:
 *
 */

#pragma once

#include "graph.h"
#include "common.h"

// TODO: remove unused includes
#include <cassert>
#include <set>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <tuple>

namespace ExpTC3Joining {

using namespace std;

constexpr int MAX_GRAPH_COUNT = MAX_VERTEX / 2;

struct MaskedTreeMatching {
  Mask tree_mask;
  Mask matching_mask;

  MaskedTreeMatching()
    : tree_mask(0)
    , matching_mask(0)
  {}

  MaskedTreeMatching(const Mask& t, const Mask& m)
    : tree_mask(t)
    , matching_mask(m)
  {}
};

bool operator<(const MaskedTreeMatching& x, const MaskedTreeMatching& y) {
    return tie(x.tree_mask, x.matching_mask) < tie(y.tree_mask, y.matching_mask);
}

bool operator==(const MaskedTreeMatching& x, const MaskedTreeMatching& y) {
    return (x.tree_mask == y.tree_mask) && (x.matching_mask == y.matching_mask);
}

struct TTreeMatching {
  bool in_tree[MAX_EDGE];
  bool in_matching[MAX_EDGE];
  int tree_deg[MAX_VERTEX];
  int number_of_vertices_;
  int number_of_edges_;

  TTreeMatching()
  {}

  void init(int number_of_vertices, int number_of_edges) {
    assert(number_of_vertices % 2 == 0);
    number_of_vertices_ = number_of_vertices;
    number_of_edges_ = number_of_edges;
    for (int v = 0; v < number_of_vertices; ++v) {
      tree_deg[v] = 0;
    }
    for (int e = 0; e < number_of_edges; ++e) {
      in_tree[e] = false;
      in_matching[e] = false;
    }
  }

  void verify() const {
    int matching_size = 0;
    for (int e = 0; e < number_of_edges_; ++e) {
      if (in_matching[e]) {
        matching_size += 1;
      }
    }
    assert(matching_size == number_of_vertices_ / 2);
  }

  TTreeMatching(int number_of_vertices, int number_of_edges)
  {
    init(number_of_vertices, number_of_edges);
  }

  TTreeMatching(const Graph& graph, const MaskedTreeMatching& m)
  {
    init(graph.number_of_vertices, graph.number_of_edges);
    for (int e = 0; e < number_of_edges_; ++e) {
      if (BIT(e) & m.tree_mask) {
        in_tree[e] = true;
        const int v1 = graph.e2v[e][0];
        const int v2 = graph.e2v[e][1];
        ++tree_deg[v1];
        ++tree_deg[v2];
      } else {
        in_tree[e] = false;
      }
      if (BIT(e) & m.matching_mask) {
        in_matching[e] = true;
      } else {
        in_matching[e] = false;
      }
    }
    verify();
  }

  TTreeMatching(const Graph& graph, const vector<int>& vertices)
  {
    assert(vertices.size() == graph.number_of_vertices);
    init(graph.number_of_vertices, graph.number_of_edges);
    int matching_size = 0;
    for (int i = 0; i < vertices.size() - 1; ++i) {
      const int v = vertices[i];
      const int u = vertices[i + 1];
      const int e = graph.edge_index[v][u];
      in_tree[e] = true;
      if (i % 2 == 0) {
        in_matching[e] = true;
        ++matching_size;
      }
      ++tree_deg[v];
      ++tree_deg[u];
    }
    assert(matching_size == number_of_vertices_ / 2);
  }

  MaskedTreeMatching GenerateMask() const {
    Mask t = 0;
    Mask m = 0;
    int matching_size = 0;
    for (int e = 0; e < number_of_edges_; ++e) {
      if (in_tree[e]) {
        t += BIT(e);
      }
      if (in_matching[e]) {
        m += BIT(e);
        matching_size += 1;
      }
    }
    assert(matching_size == number_of_vertices_ / 2);
    return MaskedTreeMatching{t, m};
  }
};

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

struct TExtra {
  vector<GraphTC3> tc3s;
};

void compare_edges(const GraphTC3& cur_tc3, const GraphTC3& more_tc3,
    vector<int>& cur_edges_also_in_more, vector<int>& more_not_cur_edges,
    int cur2more_edges[MAX_EDGE]) {
  for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
    cur2more_edges[e] = NONE;
    const int v = cur_tc3.g.e2v[e][0];
    const int u = cur_tc3.g.e2v[e][1];
    int more_v = cur_tc3.cur2more[v];
    int more_u = cur_tc3.cur2more[u];
    if (more_v != NONE && more_u != NONE && more_tc3.g.has_edge(more_v, more_u)) {
      cur_edges_also_in_more.push_back(e);
      cur2more_edges[e] = more_tc3.g.edge_index[more_v][more_u];
    }
  }
  assert(cur_edges_also_in_more.size() == cur_tc3.g.number_of_edges - 4);

  for (int e = 0; e < more_tc3.g.number_of_edges; ++e) {
    const int more_v = more_tc3.g.e2v[e][0];
    const int more_u = more_tc3.g.e2v[e][1];
    const int v = cur_tc3.more2cur[more_v];
    const int u = cur_tc3.more2cur[more_u];
    if (v == NONE || u == NONE) {
      more_not_cur_edges.push_back(e);
    }
  }
  assert(more_not_cur_edges.size() == 7);
}

void init_triangle_tree_matchings(const TExtra& extra, vector<TTreeMatching>& tms) {
  // create tm structure for triangle
  // 0-1-2-3, 1-0-2-3,
  // 0-2-1-3, 2-0-1-3,
  // 1-2-0-3, 2-1-0-3
  assert(extra.tc3s[0].special_vertex == 3);
  tms.push_back(TTreeMatching(extra.tc3s[0].g, {0, 1, 2, 3}));
  tms.push_back(TTreeMatching(extra.tc3s[0].g, {1, 0, 2, 3}));
  tms.push_back(TTreeMatching(extra.tc3s[0].g, {0, 2, 1, 3}));
  tms.push_back(TTreeMatching(extra.tc3s[0].g, {2, 0, 1, 3}));
  tms.push_back(TTreeMatching(extra.tc3s[0].g, {1, 2, 0, 3}));
  tms.push_back(TTreeMatching(extra.tc3s[0].g, {2, 1, 0, 3}));
}

void init_triangle_tree_double_covers(const TExtra& extra, const vector<TTreeMatching>& tms,
    set<multiset<MaskedTreeMatching>>& triples) {
  for (int sol1 = 0; sol1 < tms.size(); ++sol1) {
    for (int sol2 = sol1 + 1; sol2 < tms.size(); ++sol2) {
      for (int sol3 = sol2 + 1; sol3 < tms.size(); ++sol3) {
        bool has_tree_double_cover = true;
        for (int e = 0; e < extra.tc3s[0].g.number_of_edges; ++e) {
          const int v1 = extra.tc3s[0].g.e2v[e][0];
          const int v2 = extra.tc3s[0].g.e2v[e][1];
          if (v1 == extra.tc3s[0].special_vertex || v2 == extra.tc3s[0].special_vertex) {
            continue;
          }

          int edge_cover_count = tms[sol1].in_tree[e];
          edge_cover_count += tms[sol2].in_tree[e];
          edge_cover_count += tms[sol3].in_tree[e];

          if (edge_cover_count != 2) {
            has_tree_double_cover = false;
            break;
          }
        }
        if (has_tree_double_cover) {
          multiset<MaskedTreeMatching> triple;
          triple.insert(tms[sol1].GenerateMask());
          triple.insert(tms[sol2].GenerateMask());
          triple.insert(tms[sol3].GenerateMask());
          triples.insert(triple);
        }
      }
    }
  }
}

void init_triangle_tree_6c4cs(const TExtra& extra, const vector<TTreeMatching>& tms,
    set<multiset<MaskedTreeMatching>>& tree_6c4cs) {
  for (int sol1 = 0; sol1 < tms.size(); ++sol1) {
    for (int sol2 = sol1; sol2 < tms.size(); ++sol2) {
      for (int sol3 = sol2 + 1; sol3 < tms.size(); ++sol3) {
        for (int sol4 = sol3; sol4 < tms.size(); ++sol4) {
          for (int sol5 = sol4 + 1; sol5 < tms.size(); ++sol5) {
            for (int sol6 = sol5; sol6 < tms.size(); ++sol6) {

              bool has_tree_6c4c = true;
              for (int e = 0; e < extra.tc3s[0].g.number_of_edges; ++e) {
                const int v1 = extra.tc3s[0].g.e2v[e][0];
                const int v2 = extra.tc3s[0].g.e2v[e][1];

                int edge_matching_cover_count = tms[sol1].in_matching[e];
                edge_matching_cover_count += tms[sol2].in_matching[e];
                edge_matching_cover_count += tms[sol3].in_matching[e];
                edge_matching_cover_count += tms[sol4].in_matching[e];
                edge_matching_cover_count += tms[sol5].in_matching[e];
                edge_matching_cover_count += tms[sol6].in_matching[e];
                if (edge_matching_cover_count != 2) {
                  has_tree_6c4c = false;
                  break;
                }

                if (v1 == extra.tc3s[0].special_vertex || v2 == extra.tc3s[0].special_vertex) {
                  continue;
                }

                int edge_tree_cover_count = tms[sol1].in_tree[e];
                edge_tree_cover_count += tms[sol2].in_tree[e];
                edge_tree_cover_count += tms[sol3].in_tree[e];
                edge_tree_cover_count += tms[sol4].in_tree[e];
                edge_tree_cover_count += tms[sol5].in_tree[e];
                edge_tree_cover_count += tms[sol6].in_tree[e];

                if (edge_tree_cover_count != 4) {
                  has_tree_6c4c = false;
                  break;
                }
              }

              if (has_tree_6c4c) {
                multiset<MaskedTreeMatching> tree_6c4c;
                tree_6c4c.insert(tms[sol1].GenerateMask());
                tree_6c4c.insert(tms[sol2].GenerateMask());
                tree_6c4c.insert(tms[sol3].GenerateMask());
                tree_6c4c.insert(tms[sol4].GenerateMask());
                tree_6c4c.insert(tms[sol5].GenerateMask());
                tree_6c4c.insert(tms[sol6].GenerateMask());
                tree_6c4cs.insert(tree_6c4c);
              }
            }
          }
        }
      }
    }
  }
}


void expand_tree_matching(const TTreeMatching& sol,
      const GraphTC3& cur_tc3, const GraphTC3& more_tc3,
      const vector<int>& cur_edges_also_in_more, const vector<int>& more_not_cur_edges,
      const int cur2more_edges[MAX_EDGE],
      vector<TTreeMatching>& more_tms) {
  // 2.2.1 fill info about common edges
  TTreeMatching more_tm(more_tc3.g.number_of_vertices, more_tc3.g.number_of_edges);
  int tree_edges_left = more_tc3.g.number_of_vertices - 1;
  for (const auto& e : cur_edges_also_in_more) {
    if (sol.in_tree[e]) {
      --tree_edges_left;
      assert(cur2more_edges[e] != NONE);
      const int more_e = cur2more_edges[e];
      more_tm.in_tree[more_e] = true;
      const int more_v1 = more_tc3.g.e2v[more_e][0];
      const int more_v2 = more_tc3.g.e2v[more_e][1];
      ++more_tm.tree_deg[more_v1];
      ++more_tm.tree_deg[more_v2];

      /*if (more_tc3.g.number_of_vertices == 8) {
        cerr << "pre: " << e << "; "  << more_e << ": " << more_v1 << " " << more_v2 << endl;
      }*/
    }
  }

  // 2.2.2 cycle through every mask

  // 2.2.2.0 save pre-mask state
  const int tree_edges_left_seed = tree_edges_left;
  int tree_deg_seed[MAX_VERTEX];
  bool in_tree_seed[MAX_EDGE];
  for (int v = 0; v < more_tc3.g.number_of_vertices; ++v) {
    tree_deg_seed[v] = more_tm.tree_deg[v];
  }
  for (int e = 0; e < more_tc3.g.number_of_edges; ++e) {
    in_tree_seed[e] = more_tm.in_tree[e];
  }
  bool found = false;
  for (Mask mask = 0; mask < BIT(more_not_cur_edges.size()); ++mask) {
    // 2.2.2.0.0 return to pre-mask state
    tree_edges_left = tree_edges_left_seed;
    for (int v = 0; v < more_tc3.g.number_of_vertices; ++v) {
      more_tm.tree_deg[v] = tree_deg_seed[v];
    }
    for (int e = 0; e < more_tc3.g.number_of_edges; ++e) {
      more_tm.in_tree[e] = in_tree_seed[e];
      more_tm.in_matching[e] = false;
    }

    // 2.2.2.1 adding edges from the mask
    for (int i = 0; i < more_not_cur_edges.size(); ++i) {
      const int e = more_not_cur_edges[i];
      if (BIT(i) & mask) {
        more_tm.in_tree[e] = true;
      } else {
        more_tm.in_tree[e] = false;
      }
      if (more_tm.in_tree[e]) {
        --tree_edges_left;
        const int more_v1 = more_tc3.g.e2v[e][0];
        const int more_v2 = more_tc3.g.e2v[e][1];
        ++more_tm.tree_deg[more_v1];
        ++more_tm.tree_deg[more_v2];
      }
    }
    // more_tm is ready

    // 2.2.2.2.0 do i need to check, that tree_deg of special vertex is 1?
    // yes
    if (more_tm.tree_deg[more_tc3.special_vertex] != 1) {
      continue;
    }

    // 2.2.2.2 check that we got a tree (e. g., it has right amount of edges and it's connected)
    if (tree_edges_left != 0) {
      continue;
    }

    vector<int> queue;
    int parent[MAX_VERTEX];
    for (int v = 0; v < more_tc3.g.number_of_vertices; ++v) {
      parent[v] = NONE;
    }
    const int start_vertex = 0;
    queue.push_back(start_vertex);
    parent[start_vertex] = start_vertex;
    int head = 0;
    bool is_a_tree = true;
    while (head < queue.size()) {
      const int v = queue[head];
      ++head;
      for (int j = 0; j < MAX_DEG; ++j) {
        const int e = more_tc3.g.v2e[v][j];
        if (more_tm.in_tree[e]) {
          const int u = more_tc3.g.v2v[v][j];
          if (parent[u] == NONE) {
            parent[u] = v;
            queue.push_back(u);
          } else if (u != parent[v]) {
            is_a_tree = false;
            break;
          }
        }
      }
      if (!is_a_tree) {
        break;
      }
    }
    if (!is_a_tree || queue.size() != more_tc3.g.number_of_vertices) {
      continue;
    }

    // 2.2.2.3 build a matching from this tree
    bool in_tree[MAX_EDGE];
    bool in_matching[MAX_VERTEX];
    int tree_deg[MAX_VERTEX];
    int matching_size = 0;
    queue.clear();
    head = 0;
    bool visited[MAX_VERTEX];
    for (int v = 0; v < more_tc3.g.number_of_vertices; ++v) {
      visited[v] = false;
      tree_deg[v] = more_tm.tree_deg[v];
      in_matching[v] = false;
    }
    for (int e = 0; e < more_tc3.g.number_of_edges; ++e) {
      in_tree[e] = more_tm.in_tree[e];
      const int v1 = more_tc3.g.e2v[e][0];
      const int v2 = more_tc3.g.e2v[e][1];
      if (in_tree[e]) {
        assert(!visited[v1]);
        if (tree_deg[v1] == 1) {
          queue.push_back(v1);
          visited[v1] = true;
        }
        assert(!visited[v2]);
        if (tree_deg[v2] == 1) {
          queue.push_back(v2);
          visited[v2] = true;
        }
      }
    }

    while (head < queue.size()) {
      const int v = queue[head];
      ++head;
      if (tree_deg[v] == 0) {
        continue;
      }
      assert(tree_deg[v] == 1);
      assert(!in_matching[v]);
      int e = NONE;
      for (int j = 0; j < MAX_DEG; ++j) {
        e = more_tc3.g.v2e[v][j];
        if (in_tree[e]) {
          break;
        }
      }
      assert(e != NONE);
      const int u = more_tc3.g.e2v[e][0] + more_tc3.g.e2v[e][1] - v;
      if (in_matching[u]) {
        break; // tree without matching
      }
      in_matching[v] = true;
      in_matching[u] = true;
      ++matching_size;
      more_tm.in_matching[e] = true;

      for (int j = 0; j < MAX_DEG; ++j) {
        const int ee = more_tc3.g.v2e[u][j];
        if (in_tree[ee]) {
          in_tree[ee] = false;
          const int w = more_tc3.g.e2v[ee][0] + more_tc3.g.e2v[ee][1] - u;
          --tree_deg[u];
          --tree_deg[w];
          if (tree_deg[w] == 1 && !visited[w]) {
            queue.push_back(w);
            visited[w] = true;
          }
        }
      }
    }

    if (matching_size == more_tc3.g.number_of_vertices / 2) {
      found = true;
      more_tm.verify();
      more_tms.push_back(more_tm);
    }
  }

  /*if (more_tc3.g.number_of_vertices == 8 && found) {
    cerr << "cool" << endl << endl;
  }*/
  if (!found) {
    /*bool has_tree_double_cover = false;
    for (int sol1 = 0; sol1 < tms.size(); ++sol1) {
      if (sol1 == sol) {
        continue;
      }
      has_tree_double_cover = true;
      vector<int> special_edges;
      Mask mask = 0;
      for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
        int edge_cover_count = tms[sol].in_tree[e];
        edge_cover_count += tms[sol1].in_tree[e];
        const int v1 = cur_tc3.g.e2v[e][0];
        const int v2 = cur_tc3.g.e2v[e][1];
        if (v1 == cur_tc3.special_vertex || v2 == cur_tc3.special_vertex) {
          if (edge_cover_count == 0) {
            special_edges.push_back(e);
          }
          if (edge_cover_count == 2) {
            has_tree_double_cover = false;
            break;
          }
          continue;
        }
        if (edge_cover_count == 0) {
          has_tree_double_cover = false;
          break;
        }
        if (edge_cover_count == 1) {
          mask += BIT(e);
        }
      }
      if (!has_tree_double_cover) {
        continue;
      }
      has_tree_double_cover = false;
      if (special_edges.size() == 1) {
        for (const auto& e : special_edges) {
          const Mask full_mask = mask + BIT(e);
          if (tree_masks.find(full_mask) != tree_masks.end()) {
            has_tree_double_cover = true;
            break;
          }
        }
      }
      if (has_tree_double_cover) {
        break;
      }
    }

    if (has_tree_double_cover) {
      //cerr << "probably found nonextendable tree" << endl;
      ++not_found_count;

      if (print) {
        for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
          if (tms[sol].in_tree[e]) {
            const int v1 = cur_tc3.g.e2v[e][0];
            const int v2 = cur_tc3.g.e2v[e][1];
            cerr << "edge: " << v1 << " " << v2;
            const int more_e = cur2more_edges[e];
            if (more_e != NONE) {
              const int more_v1 = more_tc3.g.e2v[more_e][0];
              const int more_v2 = more_tc3.g.e2v[more_e][1];
              cerr << " vs " << more_v1 << " " << more_v2;
            }
            cerr << endl;
          }
        }
        cerr << "special vertex: " << cur_tc3.special_vertex << " and " << more_tc3.special_vertex << endl;
        for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
          cerr << "e: " << e << "; cur2more_edge: " << cur2more_edges[e] << endl;
        }
        for (int v = 0; v < cur_tc3.g.number_of_vertices; ++v) {
          cerr << "v: " << v << "; cur2more: " << cur_tc3.cur2more[v] << endl;
        }
        more_tc3.g.print();
        cerr << endl << endl;
        return;
      }
    }*/
  }
}

void expand_triple(const multiset<MaskedTreeMatching>& triple,
    const GraphTC3& cur_tc3, const GraphTC3& more_tc3,
    vector<int>& cur_edges_also_in_more, vector<int>& more_not_cur_edges,
    int cur2more_edges[MAX_EDGE], set<multiset<MaskedTreeMatching>>& more_triples) {
  MaskedTreeMatching m1;
  MaskedTreeMatching m2;
  MaskedTreeMatching m3;
  for (const auto& sol : triple) {
    if (!m1.tree_mask) {
      m1 = sol;
    } else if (!m2.tree_mask) {
      m2 = sol;
    } else {
      m3 = sol;
    }
  }
  TTreeMatching tm1{cur_tc3.g, m1};
  TTreeMatching tm2{cur_tc3.g, m2};
  TTreeMatching tm3{cur_tc3.g, m3};

  vector<TTreeMatching> more_tms1;
  vector<TTreeMatching> more_tms2;
  vector<TTreeMatching> more_tms3;

  expand_tree_matching(tm1, cur_tc3, more_tc3, cur_edges_also_in_more, more_not_cur_edges,
      cur2more_edges, more_tms1);
  expand_tree_matching(tm2, cur_tc3, more_tc3, cur_edges_also_in_more, more_not_cur_edges,
      cur2more_edges, more_tms2);
  expand_tree_matching(tm3, cur_tc3, more_tc3, cur_edges_also_in_more, more_not_cur_edges,
      cur2more_edges, more_tms3);

  for (const auto& more_tm1 : more_tms1) {
    for (const auto& more_tm2 : more_tms2) {
      for (const auto& more_tm3 : more_tms3) {
        bool has_tree_double_cover = true;
        vector<int> special_edges;
        for (int e = 0; e < more_tc3.g.number_of_edges; ++e) {
          int edge_cover_count = more_tm1.in_tree[e];
          edge_cover_count += more_tm2.in_tree[e];
          edge_cover_count += more_tm3.in_tree[e];
          const int v1 = more_tc3.g.e2v[e][0];
          const int v2 = more_tc3.g.e2v[e][1];
          if (v1 == more_tc3.special_vertex || v2 == more_tc3.special_vertex) {
            /*if (more_tc3.g.number_of_vertices < 10) {
              if (edge_cover_count != 1) {
                has_tree_double_cover = false;
                break;
              }
            }*/
            continue;
          }
          if (edge_cover_count != 2) {
            has_tree_double_cover = false;
            break;
          }
        }
        if (!has_tree_double_cover) {
          continue;
        }

        if (has_tree_double_cover) {
          multiset<MaskedTreeMatching> more_triple;
          more_triple.insert(more_tm1.GenerateMask());
          more_triple.insert(more_tm2.GenerateMask());
          more_triple.insert(more_tm3.GenerateMask());
          more_triples.insert(more_triple);
        }
      }
    }
  }
}

size_t expand_tree_6c4c(const multiset<MaskedTreeMatching>& tree_6c4c,
    const GraphTC3& cur_tc3, const GraphTC3& more_tc3,
    vector<int>& cur_edges_also_in_more, vector<int>& more_not_cur_edges,
    int cur2more_edges[MAX_EDGE], set<multiset<MaskedTreeMatching>>& more_tree_6c4cs) {
  size_t sol_count = 0;
  MaskedTreeMatching m1;
  MaskedTreeMatching m2;
  MaskedTreeMatching m3;
  MaskedTreeMatching m4;
  MaskedTreeMatching m5;
  MaskedTreeMatching m6;

  for (const auto& sol : tree_6c4c) {
    if (!m1.tree_mask) {
      m1 = sol;
    } else if (!m2.tree_mask) {
      m2 = sol;
    } else if (!m3.tree_mask) {
      m3 = sol;
    } else if (!m4.tree_mask) {
      m4 = sol;
    } else if (!m5.tree_mask) {
      m5 = sol;
    } else {
      m6 = sol;
    }
  }
  TTreeMatching tm1{cur_tc3.g, m1};
  TTreeMatching tm2{cur_tc3.g, m2};
  TTreeMatching tm3{cur_tc3.g, m3};
  TTreeMatching tm4{cur_tc3.g, m4};
  TTreeMatching tm5{cur_tc3.g, m5};
  TTreeMatching tm6{cur_tc3.g, m6};

  vector<TTreeMatching> more_tms1;
  vector<TTreeMatching> more_tms2;
  vector<TTreeMatching> more_tms3;
  vector<TTreeMatching> more_tms4;
  vector<TTreeMatching> more_tms5;
  vector<TTreeMatching> more_tms6;

  expand_tree_matching(tm1, cur_tc3, more_tc3, cur_edges_also_in_more, more_not_cur_edges,
      cur2more_edges, more_tms1);
  expand_tree_matching(tm2, cur_tc3, more_tc3, cur_edges_also_in_more, more_not_cur_edges,
      cur2more_edges, more_tms2);
  expand_tree_matching(tm3, cur_tc3, more_tc3, cur_edges_also_in_more, more_not_cur_edges,
      cur2more_edges, more_tms3);
  expand_tree_matching(tm4, cur_tc3, more_tc3, cur_edges_also_in_more, more_not_cur_edges,
      cur2more_edges, more_tms4);
  expand_tree_matching(tm5, cur_tc3, more_tc3, cur_edges_also_in_more, more_not_cur_edges,
      cur2more_edges, more_tms5);
  expand_tree_matching(tm6, cur_tc3, more_tc3, cur_edges_also_in_more, more_not_cur_edges,
      cur2more_edges, more_tms6);

  for (const auto& more_tm1 : more_tms1) {
    for (const auto& more_tm2 : more_tms2) {
      for (const auto& more_tm3 : more_tms3) {
        for (const auto& more_tm4 : more_tms4) {
          for (const auto& more_tm5 : more_tms5) {
            for (const auto& more_tm6 : more_tms6) {

              bool has_tree_6c4c = true;

              for (int e = 0; e < more_tc3.g.number_of_edges; ++e) {

                const int v1 = more_tc3.g.e2v[e][0];
                const int v2 = more_tc3.g.e2v[e][1];

                int edge_matching_cover_count = more_tm1.in_matching[e];
                edge_matching_cover_count += more_tm2.in_matching[e];
                edge_matching_cover_count += more_tm3.in_matching[e];
                edge_matching_cover_count += more_tm4.in_matching[e];
                edge_matching_cover_count += more_tm5.in_matching[e];
                edge_matching_cover_count += more_tm6.in_matching[e];
                if (edge_matching_cover_count != 2) {
                  has_tree_6c4c = false;
                  break;
                }

                if (v1 == more_tc3.special_vertex || v2 == more_tc3.special_vertex) {
                  continue;
                }

                int edge_tree_cover_count = more_tm1.in_tree[e];
                edge_tree_cover_count += more_tm2.in_tree[e];
                edge_tree_cover_count += more_tm3.in_tree[e];
                edge_tree_cover_count += more_tm4.in_tree[e];
                edge_tree_cover_count += more_tm5.in_tree[e];
                edge_tree_cover_count += more_tm6.in_tree[e];

                if (edge_tree_cover_count != 4) {
                  has_tree_6c4c = false;
                  break;
                }
              }

              if (has_tree_6c4c) {
                ++sol_count;
                multiset<MaskedTreeMatching> more_tree_6c4c;
                more_tree_6c4c.insert(more_tm1.GenerateMask());
                more_tree_6c4c.insert(more_tm2.GenerateMask());
                more_tree_6c4c.insert(more_tm3.GenerateMask());
                more_tree_6c4c.insert(more_tm4.GenerateMask());
                more_tree_6c4c.insert(more_tm5.GenerateMask());
                more_tree_6c4c.insert(more_tm6.GenerateMask());
                more_tree_6c4cs.insert(more_tree_6c4c);
              }
            }
          }
        }
      }
    }
  }
  return sol_count;
}

int build_tree_alternative_6c4c(const TExtra& extra, bool print=true) {
  vector<TTreeMatching> tms;
  init_triangle_tree_matchings(extra, tms);

  set<multiset<MaskedTreeMatching>> tree_6c4cs;
  init_triangle_tree_6c4cs(extra, tms, tree_6c4cs);

  const int number_of_graphs = extra.tc3s.size();
  for (int ii = 0; ii < number_of_graphs - 1; ++ii) {
    const GraphTC3& cur_tc3 = extra.tc3s[ii];
    if (print) {
      cerr << "ii: " << ii << "; number of solutions: " << tree_6c4cs.size() << endl;
      cur_tc3.g.print();
    }
    const GraphTC3& more_tc3 = extra.tc3s[ii + 1];
    vector<int> cur_edges_also_in_more;
    int cur2more_edges[MAX_EDGE];
    vector<int> more_not_cur_edges;
    compare_edges(cur_tc3, more_tc3, cur_edges_also_in_more, more_not_cur_edges, cur2more_edges);


    set<multiset<MaskedTreeMatching>> more_tree_6c4cs;
    for (const auto& tree_6c4c : tree_6c4cs) {
      int sol_count = expand_tree_6c4c(tree_6c4c, cur_tc3, more_tc3, cur_edges_also_in_more,
          more_not_cur_edges, cur2more_edges, more_tree_6c4cs);
      bool expanded = (sol_count > 0);

      /*if (!expanded) {
        cerr << "well, sucks" << endl;
      }*/

      /*bool any_same = false;
      for (const auto& t1 : tp.first) {
        for (const auto& t2 : tp.second) {
          if (t1.matching_mask == t2.matching_mask || t1.tree_mask == t2.tree_mask) {
            any_same = true;
          }
        }
      }*/

      /*if (!expanded && !any_same) {//(tp.first != tp.second)) {


        cerr << endl << endl;
        cerr << "not expanded" << endl;
        cerr << "CUR" << endl;
        cur_tc3.g.print();
        cerr << endl;
        cerr << "MORE" << endl;
        more_tc3.g.print();
        cerr << endl;

        cerr << "CUR2MORE" << endl;
        for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
          const int v1 = cur_tc3.g.e2v[e][0];
          const int v2 = cur_tc3.g.e2v[e][1];
          cerr << "edge: " << v1 << " " << v2;
          const int more_e = cur2more_edges[e];
          if (more_e != NONE) {
            const int more_v1 = more_tc3.g.e2v[more_e][0];
            const int more_v2 = more_tc3.g.e2v[more_e][1];
            cerr << " vs " << more_v1 << " " << more_v2;
          }
          cerr << endl;
        }
        cerr << "special vertex: " << cur_tc3.special_vertex << " and " << more_tc3.special_vertex << endl;
        for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
          cerr << "e: " << e << "; cur2more_edge: " << cur2more_edges[e] << endl;
        }
        for (int v = 0; v < cur_tc3.g.number_of_vertices; ++v) {
          cerr << "v: " << v << "; cur2more: " << cur_tc3.cur2more[v] << endl;
        }
        cerr << endl;

        cerr << "TRIPLE 1" << endl;
        for (const auto& t : tp.first) {
          cerr << t.tree_mask << " " << t.matching_mask << endl;
          const TTreeMatching tm{cur_tc3.g, t};
          cerr << "TREE" << endl;
          for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
            if (tm.in_tree[e]) {
              cerr << e << ": " << cur_tc3.g.e2v[e][0] << " -> " <<  cur_tc3.g.e2v[e][1] << endl;
            }
          }
          cerr << "MATCHING" << endl;
          for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
            if (tm.in_matching[e]) {
              cerr << e << ": " << cur_tc3.g.e2v[e][0] << " -> " <<  cur_tc3.g.e2v[e][1] << endl;
            }
          }
          cerr << endl;
        }
        cerr << endl;
        cerr << "TRIPLE 2" << endl;
        for (const auto& t : tp.second) {
          cerr << t.tree_mask << " " << t.matching_mask << endl;
          const TTreeMatching tm{cur_tc3.g, t};
          cerr << "TREE" << endl;
          for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
            if (tm.in_tree[e]) {
              cerr << e << ": " << cur_tc3.g.e2v[e][0] << " -> " <<  cur_tc3.g.e2v[e][1] << endl;
            }
          }
          cerr << "MATCHING" << endl;
          for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
            if (tm.in_matching[e]) {
              cerr << e << ": " << cur_tc3.g.e2v[e][0] << " -> " <<  cur_tc3.g.e2v[e][1] << endl;
            }
          }
          cerr << endl;
        }
        cerr << endl;

        cerr << endl << endl;
        return;
      }*/
    }
    tree_6c4cs = more_tree_6c4cs;
  }
  if (print) {
    cerr << "final number of solutions: " << tree_6c4cs.size() << endl;
    //extra.tc3s.back().g.print();
  }
  return tree_6c4cs.size();
}

int tree_shuffles[10][6] = {
    {0, 1, 2, 3, 4, 5},
    {0, 1, 3, 2, 4, 5},
    {0, 1, 4, 2, 3, 5},
    {0, 1, 5, 2, 3, 4},
    {0, 2, 3, 1, 4, 5},
    {0, 2, 4, 1, 3, 5},
    {0, 2, 5, 1, 3, 4},
    {0, 3, 4, 1, 2, 5},
    {0, 3, 5, 1, 2, 4},
    {0, 4, 5, 1, 2, 3}
};

int build_tree_second_alternative_6c4c(const TExtra& extra, bool print=true) {
  vector<TTreeMatching> tms;
  init_triangle_tree_matchings(extra, tms);

  set<multiset<MaskedTreeMatching>> triples;
  init_triangle_tree_double_covers(extra, tms, triples);

  set<pair<multiset<MaskedTreeMatching>, multiset<MaskedTreeMatching>>> triple_pairs;
  for (auto it1 = triples.begin(); it1 != triples.end(); ++it1) {
    for (auto it2 = it1; it2 != triples.end(); ++it2) {
      if (it1 == it2) {
        continue;
      }
      int matching_cover_count[MAX_EDGE];
      for (int e = 0; e < extra.tc3s[0].g.number_of_edges; ++e) {
        matching_cover_count[e] = 0;
      }
      for (const auto& t : *it1) {
        for (int e = 0; e < extra.tc3s[0].g.number_of_edges; ++e) {
          if (BIT(e) & t.matching_mask) {
            ++matching_cover_count[e];
          }
        }
      }
      for (const auto& t : *it2) {
        for (int e = 0; e < extra.tc3s[0].g.number_of_edges; ++e) {
          if (BIT(e) & t.matching_mask) {
            ++matching_cover_count[e];
          }
        }
      }
      bool is_6c4c = true;
      for (int e = 0; e < extra.tc3s[0].g.number_of_edges; ++e) {
        if (matching_cover_count[e] != 2) {
          is_6c4c = false;
          break;
        }
      }
      if (is_6c4c) {
        triple_pairs.insert(make_pair(*it1, *it2));
      }
    }
  }

  const int number_of_graphs = extra.tc3s.size();
  for (int ii = 0; ii < number_of_graphs - 1; ++ii) {
    const GraphTC3& cur_tc3 = extra.tc3s[ii];
    if (print) {
      cerr << "ii: " << ii << "; number of solutions: " << triple_pairs.size() << endl;
      cur_tc3.g.print();
    }
    const GraphTC3& more_tc3 = extra.tc3s[ii + 1];
    vector<int> cur_edges_also_in_more;
    int cur2more_edges[MAX_EDGE];
    vector<int> more_not_cur_edges;
    compare_edges(cur_tc3, more_tc3, cur_edges_also_in_more, more_not_cur_edges, cur2more_edges);

    set<pair<multiset<MaskedTreeMatching>, multiset<MaskedTreeMatching>>> more_triple_pairs;
    for (const auto& tp : triple_pairs) {
      bool expanded = false;

      vector<MaskedTreeMatching> mtms;

      for (const auto& sol : tp.first) {
        mtms.push_back(sol);
      }
      for (const auto& sol : tp.second) {
        mtms.push_back(sol);
      }

      for (int jj = 0; jj < 10; ++jj) {
        multiset<MaskedTreeMatching> first_triple;
        for (int k = 0; k < 3; ++k) {
          first_triple.insert(mtms[tree_shuffles[jj][k]]);
        }
        multiset<MaskedTreeMatching> second_triple;
        for (int k = 3; k < 6; ++k) {
          second_triple.insert(mtms[tree_shuffles[jj][k]]);
        }

        set<multiset<MaskedTreeMatching>> more_triples1;
        set<multiset<MaskedTreeMatching>> more_triples2;
        expand_triple(first_triple, cur_tc3, more_tc3, cur_edges_also_in_more, more_not_cur_edges,
            cur2more_edges, more_triples1);
        expand_triple(second_triple, cur_tc3, more_tc3, cur_edges_also_in_more, more_not_cur_edges,
            cur2more_edges, more_triples2);

        for (auto it1 = more_triples1.begin(); it1 != more_triples1.end(); ++it1) {
          for (auto it2 = more_triples2.begin(); it2 != more_triples2.end(); ++it2) {
            int matching_cover_count[MAX_EDGE];
            for (int e = 0; e < more_tc3.g.number_of_edges; ++e) {
              matching_cover_count[e] = 0;
            }
            for (const auto& t : *it1) {
              for (int e = 0; e < more_tc3.g.number_of_edges; ++e) {
                if (BIT(e) & t.matching_mask) {
                  ++matching_cover_count[e];
                }
              }
            }
            for (const auto& t : *it2) {
              for (int e = 0; e < more_tc3.g.number_of_edges; ++e) {
                if (BIT(e) & t.matching_mask) {
                  ++matching_cover_count[e];
                }
              }
            }
            bool is_6c4c = true;
            for (int e = 0; e < more_tc3.g.number_of_edges; ++e) {
              if (matching_cover_count[e] != 2) {
                is_6c4c = false;
                break;
              }
            }
            if (is_6c4c) {
              expanded = true;
              if (*it1 < *it2) {
                more_triple_pairs.insert(make_pair(*it1, *it2));
              } else {
                more_triple_pairs.insert(make_pair(*it2, *it1));
              }
            }
          }
        }
      }

      /*if (!expanded) {
        cerr << "well, sucks" << endl;
      }*/

      bool any_same = false;
      for (const auto& t1 : tp.first) {
        for (const auto& t2 : tp.second) {
          if (t1.matching_mask == t2.matching_mask || t1.tree_mask == t2.tree_mask) {
            any_same = true;
          }
        }
      }

      /*if (!expanded && !any_same) {//(tp.first != tp.second)) {


        cerr << endl << endl;
        cerr << "not expanded" << endl;
        cerr << "CUR" << endl;
        cur_tc3.g.print();
        cerr << endl;
        cerr << "MORE" << endl;
        more_tc3.g.print();
        cerr << endl;

        cerr << "CUR2MORE" << endl;
        for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
          const int v1 = cur_tc3.g.e2v[e][0];
          const int v2 = cur_tc3.g.e2v[e][1];
          cerr << "edge: " << v1 << " " << v2;
          const int more_e = cur2more_edges[e];
          if (more_e != NONE) {
            const int more_v1 = more_tc3.g.e2v[more_e][0];
            const int more_v2 = more_tc3.g.e2v[more_e][1];
            cerr << " vs " << more_v1 << " " << more_v2;
          }
          cerr << endl;
        }
        cerr << "special vertex: " << cur_tc3.special_vertex << " and " << more_tc3.special_vertex << endl;
        for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
          cerr << "e: " << e << "; cur2more_edge: " << cur2more_edges[e] << endl;
        }
        for (int v = 0; v < cur_tc3.g.number_of_vertices; ++v) {
          cerr << "v: " << v << "; cur2more: " << cur_tc3.cur2more[v] << endl;
        }
        cerr << endl;

        cerr << "TRIPLE 1" << endl;
        for (const auto& t : tp.first) {
          cerr << t.tree_mask << " " << t.matching_mask << endl;
          const TTreeMatching tm{cur_tc3.g, t};
          cerr << "TREE" << endl;
          for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
            if (tm.in_tree[e]) {
              cerr << e << ": " << cur_tc3.g.e2v[e][0] << " -> " <<  cur_tc3.g.e2v[e][1] << endl;
            }
          }
          cerr << "MATCHING" << endl;
          for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
            if (tm.in_matching[e]) {
              cerr << e << ": " << cur_tc3.g.e2v[e][0] << " -> " <<  cur_tc3.g.e2v[e][1] << endl;
            }
          }
          cerr << endl;
        }
        cerr << endl;
        cerr << "TRIPLE 2" << endl;
        for (const auto& t : tp.second) {
          cerr << t.tree_mask << " " << t.matching_mask << endl;
          const TTreeMatching tm{cur_tc3.g, t};
          cerr << "TREE" << endl;
          for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
            if (tm.in_tree[e]) {
              cerr << e << ": " << cur_tc3.g.e2v[e][0] << " -> " <<  cur_tc3.g.e2v[e][1] << endl;
            }
          }
          cerr << "MATCHING" << endl;
          for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
            if (tm.in_matching[e]) {
              cerr << e << ": " << cur_tc3.g.e2v[e][0] << " -> " <<  cur_tc3.g.e2v[e][1] << endl;
            }
          }
          cerr << endl;
        }
        cerr << endl;

        cerr << endl << endl;
        return;
      }*/
    }
    triple_pairs = more_triple_pairs;
  }
  if (print) {
    cerr << "final number of solutions: " << triple_pairs.size() << endl;
    //extra.tc3s.back().g.print();
  }
  return triple_pairs.size();
}

int build_tree_6c4c(const TExtra& extra, bool print=true) {
  vector<TTreeMatching> tms;
  init_triangle_tree_matchings(extra, tms);

  set<multiset<MaskedTreeMatching>> triples;
  init_triangle_tree_double_covers(extra, tms, triples);

  set<pair<multiset<MaskedTreeMatching>, multiset<MaskedTreeMatching>>> triple_pairs;
  for (auto it1 = triples.begin(); it1 != triples.end(); ++it1) {
    for (auto it2 = it1; it2 != triples.end(); ++it2) {
      if (it1 == it2) {
        continue;
      }
      int matching_cover_count[MAX_EDGE];
      for (int e = 0; e < extra.tc3s[0].g.number_of_edges; ++e) {
        matching_cover_count[e] = 0;
      }
      for (const auto& t : *it1) {
        for (int e = 0; e < extra.tc3s[0].g.number_of_edges; ++e) {
          if (BIT(e) & t.matching_mask) {
            ++matching_cover_count[e];
          }
        }
      }
      for (const auto& t : *it2) {
        for (int e = 0; e < extra.tc3s[0].g.number_of_edges; ++e) {
          if (BIT(e) & t.matching_mask) {
            ++matching_cover_count[e];
          }
        }
      }
      bool is_6c4c = true;
      for (int e = 0; e < extra.tc3s[0].g.number_of_edges; ++e) {
        if (matching_cover_count[e] != 2) {
          is_6c4c = false;
          break;
        }
      }
      if (is_6c4c) {
        triple_pairs.insert(make_pair(*it1, *it2));
      }
    }
  }

  const int number_of_graphs = extra.tc3s.size();
  for (int ii = 0; ii < number_of_graphs - 1; ++ii) {
    const GraphTC3& cur_tc3 = extra.tc3s[ii];
    if (print) {
      //cerr << "ii: " << ii << "; number of solutions: " << triple_pairs.size() << endl;
      //cur_tc3.g.print();
    }
    const GraphTC3& more_tc3 = extra.tc3s[ii + 1];
    vector<int> cur_edges_also_in_more;
    int cur2more_edges[MAX_EDGE];
    vector<int> more_not_cur_edges;
    compare_edges(cur_tc3, more_tc3, cur_edges_also_in_more, more_not_cur_edges, cur2more_edges);

    set<pair<multiset<MaskedTreeMatching>, multiset<MaskedTreeMatching>>> more_triple_pairs;
    for (const auto& tp : triple_pairs) {
      set<multiset<MaskedTreeMatching>> more_triples1;
      set<multiset<MaskedTreeMatching>> more_triples2;
      expand_triple(tp.first, cur_tc3, more_tc3, cur_edges_also_in_more, more_not_cur_edges,
          cur2more_edges, more_triples1);
      expand_triple(tp.second, cur_tc3, more_tc3, cur_edges_also_in_more, more_not_cur_edges,
          cur2more_edges, more_triples2);
      bool expanded = false;
      for (auto it1 = more_triples1.begin(); it1 != more_triples1.end(); ++it1) {
        for (auto it2 = more_triples2.begin(); it2 != more_triples2.end(); ++it2) {
          int matching_cover_count[MAX_EDGE];
          for (int e = 0; e < more_tc3.g.number_of_edges; ++e) {
            matching_cover_count[e] = 0;
          }
          for (const auto& t : *it1) {
            for (int e = 0; e < more_tc3.g.number_of_edges; ++e) {
              if (BIT(e) & t.matching_mask) {
                ++matching_cover_count[e];
              }
            }
          }
          for (const auto& t : *it2) {
            for (int e = 0; e < more_tc3.g.number_of_edges; ++e) {
              if (BIT(e) & t.matching_mask) {
                ++matching_cover_count[e];
              }
            }
          }
          bool is_6c4c = true;
          for (int e = 0; e < more_tc3.g.number_of_edges; ++e) {
            if (matching_cover_count[e] != 2) {
              is_6c4c = false;
              break;
            }
          }
          if (is_6c4c) {
            expanded = true;
            if (*it1 < *it2) {
              more_triple_pairs.insert(make_pair(*it1, *it2));
            } else {
              more_triple_pairs.insert(make_pair(*it2, *it1));
            }
          }
        }
      }
      /*if (!expanded) {
        cerr << "well, sucks" << endl;
      }*/

      bool any_same = false;
      for (const auto& t1 : tp.first) {
        for (const auto& t2 : tp.second) {
          if (t1.matching_mask == t2.matching_mask || t1.tree_mask == t2.tree_mask) {
            any_same = true;
          }
        }
      }

      /*if (!expanded && !any_same) {//(tp.first != tp.second)) {


        cerr << endl << endl;
        cerr << "not expanded" << endl;
        cerr << "CUR" << endl;
        cur_tc3.g.print();
        cerr << endl;
        cerr << "MORE" << endl;
        more_tc3.g.print();
        cerr << endl;

        cerr << "CUR2MORE" << endl;
        for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
          const int v1 = cur_tc3.g.e2v[e][0];
          const int v2 = cur_tc3.g.e2v[e][1];
          cerr << "edge: " << v1 << " " << v2;
          const int more_e = cur2more_edges[e];
          if (more_e != NONE) {
            const int more_v1 = more_tc3.g.e2v[more_e][0];
            const int more_v2 = more_tc3.g.e2v[more_e][1];
            cerr << " vs " << more_v1 << " " << more_v2;
          }
          cerr << endl;
        }
        cerr << "special vertex: " << cur_tc3.special_vertex << " and " << more_tc3.special_vertex << endl;
        for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
          cerr << "e: " << e << "; cur2more_edge: " << cur2more_edges[e] << endl;
        }
        for (int v = 0; v < cur_tc3.g.number_of_vertices; ++v) {
          cerr << "v: " << v << "; cur2more: " << cur_tc3.cur2more[v] << endl;
        }
        cerr << endl;

        cerr << "TRIPLE 1" << endl;
        for (const auto& t : tp.first) {
          cerr << t.tree_mask << " " << t.matching_mask << endl;
          const TTreeMatching tm{cur_tc3.g, t};
          cerr << "TREE" << endl;
          for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
            if (tm.in_tree[e]) {
              cerr << e << ": " << cur_tc3.g.e2v[e][0] << " -> " <<  cur_tc3.g.e2v[e][1] << endl;
            }
          }
          cerr << "MATCHING" << endl;
          for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
            if (tm.in_matching[e]) {
              cerr << e << ": " << cur_tc3.g.e2v[e][0] << " -> " <<  cur_tc3.g.e2v[e][1] << endl;
            }
          }
          cerr << endl;
        }
        cerr << endl;
        cerr << "TRIPLE 2" << endl;
        for (const auto& t : tp.second) {
          cerr << t.tree_mask << " " << t.matching_mask << endl;
          const TTreeMatching tm{cur_tc3.g, t};
          cerr << "TREE" << endl;
          for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
            if (tm.in_tree[e]) {
              cerr << e << ": " << cur_tc3.g.e2v[e][0] << " -> " <<  cur_tc3.g.e2v[e][1] << endl;
            }
          }
          cerr << "MATCHING" << endl;
          for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
            if (tm.in_matching[e]) {
              cerr << e << ": " << cur_tc3.g.e2v[e][0] << " -> " <<  cur_tc3.g.e2v[e][1] << endl;
            }
          }
          cerr << endl;
        }
        cerr << endl;

        cerr << endl << endl;
        return;
      }*/
    }
    triple_pairs = more_triple_pairs;
  }
  if (print) {
    cerr << "final number of solutions: " << triple_pairs.size() << endl;
    //extra.tc3s.back().g.print();
  }
  return triple_pairs.size();
}

void build_tree_double_covers(const TExtra& extra, bool print=true) {
  vector<TTreeMatching> tms;
  init_triangle_tree_matchings(extra, tms);

  set<multiset<MaskedTreeMatching>> triples;
  init_triangle_tree_double_covers(extra, tms, triples);

  const int number_of_graphs = extra.tc3s.size();
  for (int ii = 0; ii < number_of_graphs - 1; ++ii) {
    const GraphTC3& cur_tc3 = extra.tc3s[ii];
    if (print) {
      cerr << "ii: " << ii << "; number of solutions: " << triples.size() << endl;
      cur_tc3.g.print();
    }
    const GraphTC3& more_tc3 = extra.tc3s[ii + 1];
    vector<int> cur_edges_also_in_more;
    int cur2more_edges[MAX_EDGE];
    vector<int> more_not_cur_edges;
    compare_edges(cur_tc3, more_tc3, cur_edges_also_in_more, more_not_cur_edges, cur2more_edges);

    set<multiset<MaskedTreeMatching>> more_triples;
    for (const auto& triple : triples) {
      expand_triple(triple, cur_tc3, more_tc3, cur_edges_also_in_more, more_not_cur_edges,
          cur2more_edges, more_triples);
    }
    triples = more_triples;
  }
  if (print) {
    cerr << "final number of solutions: " << triples.size() << endl;
    extra.tc3s.back().g.print();
  }
}

void build_tree_matchings(const TExtra& extra, bool print=true) {
  // DONE: testing hypothesis
  // that I can always extend any tms into larger tms
  // (but actually I'm searching for counterexample)
  // regardless of results, this function will then be used
  // to construct solutions for:
  // 1. matchings: form a 6c4c solution
  // 2. trees: go into 2 double coverings of TC3 graph with spanning trees
  // hypothesis is false, unfortunately

  vector<TTreeMatching> tms;
  init_triangle_tree_matchings(extra, tms);

  const int number_of_graphs = extra.tc3s.size();

  // 2. transfer structure further
  for (int ii = 0; ii < number_of_graphs - 1; ++ii) {
    // 2.1 find edges which are common, and edges which are added
    const GraphTC3& cur_tc3 = extra.tc3s[ii];
    if (print) {
      cerr << "ii: " << ii << "; number of solutions: " << tms.size() << endl;
      cur_tc3.g.print();
    }
    const GraphTC3& more_tc3 = extra.tc3s[ii + 1];
    vector<int> cur_edges_also_in_more;
    int cur2more_edges[MAX_EDGE];
    vector<int> more_not_cur_edges;
    compare_edges(cur_tc3, more_tc3, cur_edges_also_in_more, more_not_cur_edges, cur2more_edges);

    // 2.2 iterate through all tms
    // and choose the masks that fit

    // TODO: add add_edge in TTreeMatching
    // TODO: add graph and tc3 into TTreeMatching
    // TODO: maybe convert this code into constructor

    /*int not_found_count = 0;
    set<Mask> tree_masks;
    for (int sol = 0; sol < tms.size(); ++sol) {
      Mask mask = 0;
      for (int e = 0; e < cur_tc3.g.number_of_edges; ++e) {
        if (tms[sol].in_tree[e]) {
          mask += BIT(e);
        }
      }
      tree_masks.insert(mask);
    }*/

    vector<TTreeMatching> more_tms; // tms for extra.tc3s[ii + 1]
    for (int sol = 0; sol < tms.size(); ++sol) {
      expand_tree_matching(tms[sol], cur_tc3, more_tc3, cur_edges_also_in_more, more_not_cur_edges,
          cur2more_edges, more_tms);
    }

    /*if (not_found_count > 0) {
      cerr << "not found: " << not_found_count << endl;
    }*/
    tms = more_tms;
    if (print) {
      cerr << "more_tms: " << tms.size() << endl;
    }
    //break;
  }
}

bool is_joinable(TExtra& extra, int more_graph_number) {
  //cerr << "more graph number: " << more_graph_number << endl;
  if (more_graph_number == 0) {
    int number_of_solutions = build_tree_second_alternative_6c4c(extra);
    return true;//number_of_solutions > 0;//false;// true;
  }

  const GraphTC3& more_tc3 = extra.tc3s[more_graph_number];
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
      const int graph_number = more_graph_number - 1;
      GraphTC3& cur_tc3 = extra.tc3s[graph_number];
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
      if (is_joinable(extra, graph_number)) {
        return true;
      }
    }
  }
  return false;
}

void check_tc3_joinability(Graph& graph) {
  TExtra extra;

  int graph_count = graph.number_of_vertices / 2 - 1;
  for (int i = 0; i < graph_count; ++i) {
    extra.tc3s.push_back(GraphTC3(2 * i + 4));
  }

  for (int v = 0; v < graph.number_of_vertices; ++v) {
    extra.tc3s.back() = GraphTC3(graph.number_of_vertices);
    extra.tc3s.back().special_vertex = v;
    for (int e = 0; e < graph.number_of_edges; ++e) {
      extra.tc3s.back().add_edge(graph.e2v[e][0], graph.e2v[e][1]);
    }

    // now we have TC3 graph
    // let's check joinability
    bool joinable = is_joinable(extra, graph_count - 1);
    if (joinable) {
      return;
    }
  }
  cerr << "graph is not joinable" << endl;
}

} // ExpTC3Joining
