/*
 * File:   main.h
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * Created on 21 August 2019
 *
 */

#pragma once

#include "args.h"
#include "graph.h"


void init();
void process_petersen_graph(Graph& petersen_graph);
bool read_graph(const Args& args, Graph& graph);
void run_experiments(const Graph& petersen_graph, Graph& graph);
