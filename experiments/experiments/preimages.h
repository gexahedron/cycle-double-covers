/*
 * File:   preimages.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
 *
 */

#pragma once

#include "graph.h"


namespace ExpPreimages {

void preimage_full_cycles(const Graph& petersen_graph, Graph& graph);
void preimage_all_cycles(const Graph& petersen_graph, Graph& graph);
void preimage_dominating_circuits(const Graph& petersen_graph, Graph& graph);
void preimage_6c4c_5cdc_cycles(const Graph& petersen_graph, Graph& graph);

} // ExpPreimages
