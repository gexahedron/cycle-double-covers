/*
 * File:   flow_parity_pairs.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * TODO:
 *
 *
 */

#pragma once

#include "graph.h"
#include "util/flows.h"

#include <cassert>
#include <set>
#include <iostream>
#include <cstdlib>
#include <vector>

using namespace std;


namespace ExpFlowParityPairs {
void find_all_33pp(Graph& graph);
} // ExpFlowParityPairs
