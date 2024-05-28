/*
 * File:   petersen_colouring.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 */

#pragma once

#include "graph.h"


namespace ExpPetersenColouring {

void find_all_petersen_colourings(Graph& graph);

void create_cc_mapping(const Graph& petersen_graph);

} // ExpPetersenColouring
