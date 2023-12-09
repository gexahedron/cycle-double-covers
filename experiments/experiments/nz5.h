/*
 * File:   nz5.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * Created around 14 August 2017
 * TODO: add description
 *
 */

#pragma once

#include "graph.h"

#include "util/flows.h"

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

using namespace std;


namespace ExpNZ5 {

using namespace UtilFlows;

vector<int> nz5_flow_vals = {-4, -3, -2, -1, 1, 2, 3, 4};
int nz5_edge_flow[MAX_EDGE];
int nz5_vertex_flow[MAX_VERTEX];
set<vector<int>> nz5_orientations;

set<Mask> full_cycles_from_nz5;

vector<int> nz_mod5_flow_vals = {1, 2, 3, 4};
int nz_mod5_edge_flow[MAX_EDGE];
int nz_mod5_vertex_flow[MAX_VERTEX];
int fancy_4colouring[MAX_VERTEX];
int fancy_single_neibs[5] = {-1, 2, 4, 1, 3};
int fancy_neib_counts[MAX_VERTEX][5];
int fancy_colour_count[5];

// Gaussian integers
vector<int> nz_complex4_flow_vals[12] =
  {{2, 0}, {1, 1}, {0, 2},
   {1, 0}, {0, 1},
   {1, -1}, {-1, 1},
   {0, -1}, {-1, 0},
   {0, -2}, {-1, -1}, {-2, 0}};

// Eisenstein integers
/*vector<int> nz_complex4_flow_vals[18] =
   {{1, 0}, {1, 1}, {0, 1}, {-1, 0}, {-1, -1}, {0, -1},
    {2, 1}, {1, -2}, {-1, 1}, {-2, -1}, {-1, 2}, {1, -1},
    {2, 0}, {2, 2}, {0, 2}, {-2, 0}, {-2, -2}, {0, -1}};*/

 vector<int> nz_complex4_edge_flow[MAX_EDGE];
 vector<int> nz_complex4_vertex_flow[MAX_VERTEX];

/*********************************Methods*********************************/

bool gen_all_nz5_flows(Graph& graph, int cur_edge_idx) {
    if (cur_edge_idx == graph.number_of_edges) {
        vector<int> flow_vec;
        vector<int> orientations;
        for (int e = 0; e < graph.number_of_edges; ++e) {
            flow_vec.push_back(nz5_edge_flow[e]);
            if (nz5_edge_flow[e] < 0) {
                orientations.push_back(-1);
            } else {
                orientations.push_back(1);
            }
        }
        graph.all_nz5_flows.push_back(flow_vec);
        nz5_orientations.insert(orientations);
        return false;
    }

    int cur_edge = graph.faster_edge_order[cur_edge_idx];

    nz5_edge_flow[cur_edge] = 0;
    int v1 = graph.e2v[cur_edge][0];
    int v2 = graph.e2v[cur_edge][1];
    int not_flowed1 = 0;
    for (int j = 0; j < MAX_DEG; ++j) {
        if (nz5_edge_flow[graph.v2e[v1][j]] == 0) {
            ++not_flowed1;
        }
    }
    int not_flowed2 = 0;
    for (int j = 0; j < MAX_DEG; ++j) {
        if (nz5_edge_flow[graph.v2e[v2][j]] == 0) {
            ++not_flowed2;
        }
    }

    for (size_t f_idx = 0; f_idx < nz5_flow_vals.size(); ++f_idx) {
        int f = nz5_flow_vals[f_idx];

        if ((not_flowed1 == 1) && nz5_vertex_flow[v1] != f) {
            continue;
        }

        if ((not_flowed2 == 1) && nz5_vertex_flow[v2] != -f) {
            continue;
        }

        nz5_edge_flow[cur_edge] = f;

        nz5_vertex_flow[v1] -= f;
        nz5_vertex_flow[v2] += f;
        if (gen_all_nz5_flows(graph, cur_edge_idx + 1)) {
            return true;
        }
        nz5_vertex_flow[v1] += f;
        nz5_vertex_flow[v2] -= f;

    }
    nz5_edge_flow[cur_edge] = 0;
    return false;
}

void gen_full_cycles_from_nz5_flows(Graph& graph) {
    full_cycles_from_nz5.clear();
    int good_count = 0;
    int maybe_count = 0;
    int bad_count = 0;

    bool printed = false;
    for (const auto& f : graph.all_nz5_flows) {
        if (!printed) {
            printed = true;
            graph.print();
            for (int v = 0; v < graph.number_of_vertices; ++v) {
              cerr << v << ": ";
              for (int j = 0; j < MAX_DEG; ++j) {
                int flow = f[graph.v2e[v][j]];
                if (flow < 0)
                  flow = 5 + flow;
                if (graph.v2v[v][j] < v) {
                  flow = 5 - flow;
                }
                cerr << flow << " ";
              }
              cerr << endl;
            }
        }
        Mask c = build_full_cycle_from_nz5_flow(graph, f);
        if (c == 0) {
            ++bad_count;
        } else if (graph.all_full_cycles.find(c) == graph.all_full_cycles.end()) {
            ++maybe_count;
        } else {
            ++good_count;
            full_cycles_from_nz5.insert(c);
        }
    }
    cerr << "good  nz5: " << (float) good_count / graph.all_nz5_flows.size() << endl;
    cerr << "bad   nz5: " << (float) bad_count / graph.all_nz5_flows.size() << endl;
    cerr << "maybe nz5: " << (float) maybe_count / graph.all_nz5_flows.size() << endl;
    cerr << "full cycles from nz5 vs all: " << full_cycles_from_nz5.size() << " " << graph.all_full_cycles.size() << endl;
}

void find_all_nz5_flows(Graph& graph) {
    nz5_orientations.clear();
    for (int v = 0; v < graph.number_of_vertices; ++v) {
        nz5_vertex_flow[v] = 0;
    }
    for (int e = 0; e < graph.number_of_edges; ++e) {
        nz5_edge_flow[e] = 0;
    }

    gen_all_nz5_flows(graph, 0);
    //gen_full_cycles_from_nz5_flows(graph); // i don't like this construction

    cerr << "number of nz5 flows: " << graph.all_nz5_flows.size() << endl;
    cerr << "number of nz5 orientations: " << nz5_orientations.size() << endl;
}

bool gen_all_nz_mod5_flows(Graph& graph, int cur_edge_idx) {
    if (cur_edge_idx == graph.number_of_edges) {
        vector<int> flow_vec;
        for (int e = 0; e < graph.number_of_edges; ++e) {
            flow_vec.push_back(nz_mod5_edge_flow[e]);
        }
        graph.all_nz_mod5_flows.push_back(flow_vec);
        return false;
    }

    int cur_edge = graph.faster_edge_order[cur_edge_idx];

    nz_mod5_edge_flow[cur_edge] = 0;
    int v1 = graph.e2v[cur_edge][0];
    int v2 = graph.e2v[cur_edge][1];
    int not_flowed1 = 0;
    for (int j = 0; j < MAX_DEG; ++j) {
        if (nz_mod5_edge_flow[graph.v2e[v1][j]] == 0) {
            ++not_flowed1;
        }
    }
    int not_flowed2 = 0;
    for (int j = 0; j < MAX_DEG; ++j) {
        if (nz_mod5_edge_flow[graph.v2e[v2][j]] == 0) {
            ++not_flowed2;
        }
    }

    const int mod = nz_mod5_flow_vals.size() + 1;
    for (size_t f_idx = 0; f_idx < nz_mod5_flow_vals.size(); ++f_idx) {
        int f = nz_mod5_flow_vals[f_idx];

        if ((not_flowed1 == 1) && (nz_mod5_vertex_flow[v1] + mod - f) % mod != 0) {
            continue;
        }

        if ((not_flowed2 == 1) && (nz_mod5_vertex_flow[v2] + f) % mod != 0) {
            continue;
        }

        nz_mod5_edge_flow[cur_edge] = f;

        nz_mod5_vertex_flow[v1] += mod - f;
        nz_mod5_vertex_flow[v2] += f;
        if (gen_all_nz_mod5_flows(graph, cur_edge_idx + 1)) {
            return true;
        }
        nz_mod5_vertex_flow[v1] -= mod - f;
        nz_mod5_vertex_flow[v2] -= f;

    }
    nz_mod5_edge_flow[cur_edge] = 0;
    return false;
}

bool gen_all_fancy_4colourings(Graph& graph, int v1) {
  if (v1 == graph.number_of_vertices) {
    if ((fancy_colour_count[1] != fancy_colour_count[4]) ||
        (fancy_colour_count[2] != fancy_colour_count[3])) {
      return false;
    }
    vector<int> colour_vec;
    for (int v = 0; v < graph.number_of_vertices; ++v) {
      colour_vec.push_back(fancy_4colouring[v]);
    }
    graph.all_fancy_4colourings.push_back(colour_vec);
    return false;
  }

  int v1_neib_counts[5];
  for (int c = 0; c < 5; ++c) {
    v1_neib_counts[c] = 0;
  }
  for (int j = 0; j < MAX_DEG; ++j) {
    int v2 = graph.v2v[v1][j];
    if (fancy_4colouring[v2] != 0) {
      int v2_col = fancy_4colouring[v2];
      v1_neib_counts[v2_col] += 1;
    }
  }

  for (int v1_col = 1; v1_col <= 4; ++v1_col) {
    fancy_4colouring[v1] = v1_col;
    ;
    if (v1_neib_counts[fancy_single_neibs[v1_col]] > 1 ||
        v1_neib_counts[5 - fancy_single_neibs[v1_col]] == 3) {
      continue;
    }

    bool good_colour = true;
    for (int j = 0; j < MAX_DEG; ++j) {
      int v2 = graph.v2v[v1][j];
      if (fancy_4colouring[v2] != 0) {
        int v2_col = fancy_4colouring[v2];
        if ((v2_col == v1_col) ||
            (fancy_single_neibs[v2_col] == v1_col && fancy_neib_counts[v2][v1_col] == 1) ||
            (fancy_single_neibs[v2_col] == 5 - v1_col && fancy_neib_counts[v2][v1_col] == 2)) {
          good_colour = false;
          break;
        }
      }
    }
    if (good_colour) {
      fancy_colour_count[v1_col] += 1;
      for (int j = 0; j < MAX_DEG; ++j) {
        int v2 = graph.v2v[v1][j];
        fancy_neib_counts[v2][v1_col] += 1;
      }
      if (gen_all_fancy_4colourings(graph, v1 + 1)) {
        return true;
      }
      fancy_colour_count[v1_col] -= 1;
      for (int j = 0; j < MAX_DEG; ++j) {
        int v2 = graph.v2v[v1][j];
        fancy_neib_counts[v2][v1_col] -= 1;
      }
    }
  }
  fancy_4colouring[v1] = 0;
  return false;
}

void compare_nz_mod5_flows_and_fancy_4colourings(Graph& graph) {
  bool printed = false;
  set<string> nz_mod5_colourings_str;
  set<string> fancy_4colourings_str;
  size_t min_e = graph.number_of_edges;
  size_t max_e = 0;
  for (const auto& f : graph.all_nz_mod5_flows) {
    string s;
    set<int> edges;
    for (int v = 0; v < graph.number_of_vertices; ++v) {
      int counts[5];
      for (int i = 0; i < 5; ++i) {
        counts[i] = 0;
      }
      for (int j = 0; j < MAX_DEG; ++j) {
        int flow = f[graph.v2e[v][j]];
        if (graph.v2v[v][j] < v) {
          flow = 5 - flow;
        }
        counts[flow] += 1;
      }
      int type = -1;
      for (int i = 1; i < 5; ++i) {
        if (counts[i] == 1) {
          type = i;
          for (int j = 0; j < MAX_DEG; ++j) {
            int flow = f[graph.v2e[v][j]];
            if (graph.v2v[v][j] < v) {
              flow = 5 - flow;
            }
            if (flow == type) {
              edges.insert(graph.v2e[v][j]);
            }
          }
          break;
        }
      }
      assert(type != -1);
      s += to_string(type);
    }
    min_e = min(min_e, edges.size());
    max_e = max(max_e, edges.size());
    nz_mod5_colourings_str.insert(s);
  }
  cerr << "min: " << min_e << "; max_e: " << max_e << endl;

  for (const auto& c : graph.all_fancy_4colourings) {
    string s;
    for (int v = 0; v < graph.number_of_vertices; ++v) {
      s += to_string(c[v]);
    }
    fancy_4colourings_str.insert(s);
    if (nz_mod5_colourings_str.find(s) != nz_mod5_colourings_str.end()) {
      continue;
    }

    if (!printed) {
      printed = true;
      graph.print();
      for (int v = 0; v < graph.number_of_vertices; ++v) {
        cerr << v << ": " << c[v] << endl;
      }
    }
  }

  for (const auto& s : nz_mod5_colourings_str) {
    if (fancy_4colourings_str.find(s) == fancy_4colourings_str.end()) {
      cerr << "wtf: " << s << endl;
    }
    assert(fancy_4colourings_str.find(s) != fancy_4colourings_str.end());
  }
}

void find_all_nz_mod5_flows(Graph& graph) {
    for (int v = 0; v < graph.number_of_vertices; ++v) {
        nz_mod5_vertex_flow[v] = 0;
    }
    for (int e = 0; e < graph.number_of_edges; ++e) {
        nz_mod5_edge_flow[e] = 0;
    }
    for (int v = 0; v < graph.number_of_vertices; ++v) {
      fancy_4colouring[v] = 0;
      for (int c = 1; c <= 4; ++c) {
        fancy_neib_counts[v][c] = 0;
      }
    }

    gen_all_nz_mod5_flows(graph, 0);
    //gen_all_fancy_4colourings(graph, 0);
    cerr << "number of nz-mod5 flows: " << graph.all_nz_mod5_flows.size() << " vs " << (float)(graph.all_nz_mod5_flows.size()) / 240 << endl;

    map<int, int> type_counts;
    for (int i = 0; i < graph.number_of_vertices / 2; ++i) {
      type_counts[i] = 0;
    }
    for (const auto& f : graph.all_nz_mod5_flows) {
      int type1_count = 0;
      for (int v = 0; v < graph.number_of_vertices; ++v) {
        int counts[5];
        for (int i = 0; i < 5; ++i) {
          counts[i] = 0;
        }
        for (int j = 0; j < MAX_DEG; ++j) {
          int flow = f[graph.v2e[v][j]];
          if (graph.v2v[v][j] < v) {
            flow = 5 - flow;
          }
          counts[flow] += 1;
        }
        if (counts[1] == 1) {
            ++type1_count;
        }
      }
      type_counts[type1_count] += 1;
    }
    cerr << "counts by type 1 vertices:" << endl;
    for (int i = 0; i < graph.number_of_vertices / 2; ++i) {
      if (type_counts[i] != 0) {
        cerr << "#type1 = " << i << "; counts = " << type_counts[i] << endl;
      }
    }
    //cerr << "number of fancy 4-vertex-colourings: " << graph.all_fancy_4colourings.size() << endl;
    //compare_nz_mod5_flows_and_fancy_4colourings(graph);
    // verdict:
    // not fancy enough,
    // works for Petersen graph,
    // not enough for other snarks
}

bool gen_all_nz_complex4_flows(Graph& graph, int cur_edge_idx) {
  if (cur_edge_idx == graph.number_of_edges) {
    vector<vector<int>> flow_vec;
    for (int e = 0; e < graph.number_of_edges; ++e) {
      flow_vec.push_back(nz_complex4_edge_flow[e]);
    }
    graph.all_nz_complex4_flows.push_back(flow_vec);
    return false;
  }

  int cur_edge = graph.faster_edge_order[cur_edge_idx];

  nz_complex4_edge_flow[cur_edge] = {0, 0};
  int v1 = graph.e2v[cur_edge][0];
  int v2 = graph.e2v[cur_edge][1];
  int not_flowed1 = 0;
  for (int j = 0; j < MAX_DEG; ++j) {
    if (nz_complex4_edge_flow[graph.v2e[v1][j]][0] == 0 &&
        nz_complex4_edge_flow[graph.v2e[v1][j]][1] == 0) {
      ++not_flowed1;
    }
  }
  int not_flowed2 = 0;
  for (int j = 0; j < MAX_DEG; ++j) {
    if (nz_complex4_edge_flow[graph.v2e[v2][j]][0] == 0 &&
        nz_complex4_edge_flow[graph.v2e[v2][j]][1] == 0) {
      ++not_flowed2;
    }
  }

  for (int f_idx = 0; f_idx < 12; ++f_idx) {
    vector<int> f = nz_complex4_flow_vals[f_idx];

    if ((not_flowed1 == 1) &&
        !(nz_complex4_vertex_flow[v1][0] == f[0] && nz_complex4_vertex_flow[v1][1] == f[1])) {
      continue;
    }

    if ((not_flowed2 == 1) &&
        !(nz_complex4_vertex_flow[v2][0] == -f[0] && nz_complex4_vertex_flow[v2][1] == -f[1])) {
      continue;
    }

    nz_complex4_edge_flow[cur_edge] = f;

    nz_complex4_vertex_flow[v1][0] -= f[0];
    nz_complex4_vertex_flow[v1][1] -= f[1];
    nz_complex4_vertex_flow[v2][0] += f[0];
    nz_complex4_vertex_flow[v2][1] += f[1];
    if (gen_all_nz_complex4_flows(graph, cur_edge_idx + 1)) {
      return true;
    }
    nz_complex4_vertex_flow[v1][0] += f[0];
    nz_complex4_vertex_flow[v1][1] += f[1];
    nz_complex4_vertex_flow[v2][0] -= f[0];
    nz_complex4_vertex_flow[v2][1] -= f[1];
  }
  nz_complex4_edge_flow[cur_edge] = {0, 0};
  return false;
}

void find_all_nz_complex4_flows(Graph& graph) {
  for (int v = 0; v < graph.number_of_vertices; ++v) {
    nz_complex4_vertex_flow[v] = {0, 0};
  }
  for (int e = 0; e < graph.number_of_edges; ++e) {
    nz_complex4_edge_flow[e] = {0, 0};
  }

  gen_all_nz_complex4_flows(graph, 0);
  cerr << "number of nz complex4 flows: " << graph.all_nz_complex4_flows.size() << endl;
}

} // ExpNZ5
