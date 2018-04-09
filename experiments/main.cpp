/*
 * File:   experiments.cpp
 * Author: Nikolay Ulyanov (ulyanick@gmail.com)
 *
 * Created on 14 August 2017
 * Study of properties and relations
 * between various constructions
 * relevant to snarks:
 * normal (Petersen) colouring, o6c4c, o5cdc, nz5 (nz-mod5, 33-pp)
 *
 */

// TODO: add some kind of documentation, e. g., doxygen-style
// TODO: maybe add in every experiment structure, which will hold all additional data

#include "graph.h"

// include all experiments
#include "experiments/petersen_colouring.h"
#include "experiments/tree_cycle_matching.h"
#include "experiments/cycles.h"
#include "experiments/preimages.h"
#include "experiments/mnk_flows.h" // TODO: this is a subexperiment of 6c4c; needs rethinking a bit
#include "experiments/o6c4c.h"
#include "experiments/o5cdc.h"
#include "experiments/nz5.h"
#include "experiments/flow_parity_pairs.h"
#include "experiments/matching_5cdc_and_6c4c.h"
#include "experiments/tc3_joining.h"

#include <iostream>
#include <cstdlib>

using namespace std;

//const bool PRINT_SOLUTIONS = false;

int main(int argc, char** argv) {
    srand(time(NULL)); // TODO: add function init

    // TODO: wrap up in TArgs and read_args function
    int to_skip = 0;
    if (argc < 2 || argc > 3) {
        cerr << "Error: invalid number of arguments" << endl;
        cerr << "Usage: " << argv[0] << " <path_to_petersen_graph> <number_of_graphs_to_skip>" << endl;
        exit(1);
    } else {
        if (argc >= 3) {
            to_skip = atoi(argv[2]);
        }
    }

    // initialize petersen graph
    TGraph petersen_graph;
    FILE* petersen_file;
    petersen_file = fopen(argv[1], "rb");
    decode_multicode(petersen_file, petersen_graph);
    petersen_graph.print();

    NExpPetersenColouring::find_all_petersen_colourings(petersen_graph);
    NExpPetersenColouring::create_cc_mapping(petersen_graph);
    NExpCycles::prepare_build_cycle(petersen_graph);
    cerr << "petersen cycles: " << petersen_graph.all_cycles.size() << endl;
    cerr << "petersen perfect matchings: " << petersen_graph.all_full_cycles.size() << endl;
    //cerr << "petersen tree-cycle-matching solutions: " << petersen_graph.tree_cycle_matchings.size() << endl;

    // run some useful experiments for Petersen graph
    NExp5cdc::find_all_o5cdc(petersen_graph, true);
    NExp6c4c::find_all_o6c4c(petersen_graph, true);
    cerr << "petersen 5cdc: " << petersen_graph.all_5cdc.size() << endl;
    cerr << "petersen 6c4c: " << petersen_graph.all_6c4c.size() << endl;

    NExp6c4c::NExpMNKFlows::gen_333flows_combinations();

    size_t number_of_graphs_without_solution = 0;
    size_t number_of_graphs_read = 0;
    while (true) {
        // initialize reference graph
        TGraph ref_graph;
        if (!decode_multicode(stdin, ref_graph)) {
            break;
        }
        ++number_of_graphs_read;
        if (to_skip >= number_of_graphs_read) {
            continue;
        }
        cerr << "g" << number_of_graphs_read << "\t" << endl << flush;
        cout << "g" << number_of_graphs_read << "\t" << endl << flush;

        // building cycles
        NExpCycles::prepare_build_cycle(ref_graph);

        // TODO: delete code
        for (int v = 0; v < ref_graph.number_of_vertices; ++v) {
            TMask m = BIT(v);
            for (int i = 0; i < REG; ++i) {
                m += BIT(ref_graph.v2v[v][i]);
            }
            ref_graph.all_vertex_neib_masks.insert(m);
        }

        // experimenting
        NExpTC3Joining::check_tc3_joinability(ref_graph);
    }
    cerr << "fin" << endl;
    return(EXIT_SUCCESS);

    // #include <functional>

    // void ReadGraphs(const std::function<bool()>& experiment, ...
}
