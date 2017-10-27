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
    print_graph(petersen_graph);

    NExpPetersenColouring::find_all_petersen_colourings(petersen_graph);
    NExpPetersenColouring::create_cc_mapping(petersen_graph);
    NExpCycles::prepare_build_cycle(petersen_graph);
    cerr << "petersen cycles: " << petersen_graph.all_cycles.size() << endl;
    cerr << "petersen perfect matchings: " << petersen_graph.all_full_cycles.size() << endl;
    //cerr << "petersen tree-cycle-matching solutions: " << petersen_graph.tree_cycle_matchings.size() << endl;

    NExp5cdc::find_all_o5cdc(petersen_graph, true);
    NExp6c4c::find_all_o6c4c(petersen_graph, true);
    cerr << "petersen 5cdc: " << petersen_graph.all_5cdc.size() << endl;
    cerr << "petersen 6c4c: " << petersen_graph.all_6c4c.size() << endl;

    NExp6c4c::NExpMNKFlows::gen_333flows_combinations();

    size_t number_of_graphs_without_solution = 0;
    size_t number_of_graphs_read = 0;
    while (true) {
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

        NExpCycles::prepare_build_cycle(ref_graph);
        //cerr << "graph oddness: " << ref_graph.oddness << endl;
        //cerr << "graph cycles: " << ref_graph.all_cycles.size() << endl;
        //cerr << "graph not full cycles: " << ref_graph.all_cycles.size() - ref_graph.all_full_cycles.size() << endl;
        //cerr << "graph circuits: " << ref_graph.all_circuits.size() << endl;
        //cerr << "graph dominating circuits: " << ref_graph.all_dominating_circuits.size() << endl;
        //cerr << "graph even cycles: " << ref_graph.all_even_cycles.size() << endl;
        //cerr << "graph perfect matchings: " << ref_graph.all_full_cycles.size() << endl;
        //cerr << endl;
        //NPetersenColouring::find_all_petersen_colourings(ref_graph);
        //cerr << "petersen poor: " << petersen_min_poor << "-" << petersen_max_poor << endl;
        /*cerr << "petersen always poor: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << petersen_always_poor[e] << " ";
        }
        cerr << endl;
        cerr << "petersen always rich: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << petersen_always_rich[e] << " ";
        }
        cerr << endl;*/

        //preimage_full_cycles(petersen_graph, ref_graph);
        //preimage_6c4c_5cdc_cycles(petersen_graph, ref_graph);
        //find_o6c4c_compatible_with_preimages(ref_graph);
        NExpNZ5::find_all_nz5_flows(ref_graph);

        cerr << "searching for combination: ";
        NExp6c4c::find_all_o6c4c(ref_graph);
        cerr << endl;

        NExpFlowParityPairs::find_all_33pp(ref_graph);
        
        //cerr << "all 33pp cycles: " << all_33pp_cycles.size() << endl;
        //cerr << "all 33pp not full cycles: " << all_33pp_cycles.size() - all_33pp_full_cycles.size() << endl;
        //cerr << "all 33pp circuits: " << all_33pp_circuits.size() << endl;
        //cerr << "all 33pp full cycles: " << all_33pp_full_cycles.size() << endl;
        //cerr << "all 333pp cycles: " << all_333pp_cycles.size() << endl;
        //cerr << "all 333pp even cycles: " << all_333pp_even_cycles.size() << endl;
        //cerr << "all 333pp full cycles: " << all_333pp_full_cycles.size() << endl;
        //cerr << endl;

        //find_all_o5cdc(ref_graph);
        //cerr << "all 5cdc cycles: " << all_cycles_from_5cdc.size() << endl;
        //cerr << "all 5cdc not full cycles: " << all_cycles_from_5cdc.size() - all_full_cycles_from_5cdc.size() << endl;
        //cerr << "all 5cdc circuits: " << all_circuits_from_5cdc.size() << endl;
        //cerr << "all 5cdc full cycles: " << all_full_cycles_from_5cdc.size() << endl;
        //cerr << "all even cycles from 5cdc: " << all_even_cycles_from_5cdc.size() << endl;

       
        //cerr << "all 6c4c poor: " << u6c4c_min_poor << "-" << u6c4c_max_poor << endl;
        /*cerr << "all 6c4c always poor: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << u6c4c_always_poor[e] << " ";
        }
        cerr << endl;
        cerr << "all 6c4c always rich: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << u6c4c_always_rich[e] << " ";
        }*/

        /*cerr << "   o6c4c always    1: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << o6c4c_always_1[e] << " ";
        }
        cerr << endl;
        cerr << "   o6c4c always    2: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << o6c4c_always_2[e] << " ";
        }
        cerr << endl;
        cerr << "   o6c4c always    3: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << o6c4c_always_3[e] << " ";
        }
        cerr << endl;
        cerr << "   o6c4c always    4: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << o6c4c_always_4[e] << " ";
        }
        cerr << endl;
        cerr << endl;

        cerr << "   o6c4c never     1: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << o6c4c_never_1[e] << " ";
        }
        cerr << endl;
        cerr << "   o6c4c never     2: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << o6c4c_never_2[e] << " ";
        }
        cerr << endl;
        cerr << "   o6c4c never     3: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << o6c4c_never_3[e] << " ";
        }
        cerr << endl;
        cerr << "   o6c4c never     4: ";
        for (int e = 0; e < ref_graph.number_of_edges; ++e) {
            cerr << o6c4c_never_4[e] << " ";
        }
        cerr << endl;*/

        //cerr << "333pp cycles from o5cdc: " << u333pp_cycles_from_o5cdc.size() << endl;
        //cerr << "333pp even cycles from o5cdc: " << u333pp_even_cycles_from_o5cdc.size() << endl;
        //cerr << "comparing" << endl;
        //compare_6c4c_5cdc_pairs(ref_graph);

        cerr << endl;

        //cerr << "333pp even cycles from 6c4c: " << u333pp_cycles_from_o6c4c.size() << endl;
        //cerr << endl;
        //cerr << "33pp pairs: " << u33pp_pairs.size() << endl;
        //cerr << endl;
        //preimage_tree_cycle_matchings(petersen_graph, ref_graph);
        //break;

        //find_construction(ref_graph);
    }
    cerr << "fin" << endl;
    return(EXIT_SUCCESS);

    // #include <functional>

    // void ReadGraphs(const std::function<bool()>& experiment, ...
}

/*

void ReadGraphs(const function<bool()>& experiment, int number_of_graphs_to_skip, int& number_of_vertices, int& number_of_edges, GRAPH graph) {
    unsigned long long int number_of_graphs_read = 0;
    unsigned long long int number_of_graphs_with_failed_experiment = 0;
    while (DecodeNextGraphInMulticode(stdin, number_of_vertices, number_of_edges, graph)) {
        ++number_of_graphs_read;
        if (number_of_graphs_to_skip >= number_of_graphs_read) {
            continue;
        }

        if (!experiment()) {
            ++number_of_graphs_with_failed_experiment;
        }
    }
}

*/
