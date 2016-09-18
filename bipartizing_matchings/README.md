#3 Bipartizing matchings

##How to run the code
To compile the code just print "make" in your console (you need g++ compiler).

To run this code you need to specify the number of vertices in the snarks and provide the file with snarks in multicode format.
(you can read instructions in README.md file in three_four_flows folder)

##Main result here (September 18, 2016)
* The story begins from theorem by Fleischner [1]: by G we denote a cubic graph, and D a dominating circuit (aka a connected cycle in case of a cubic graph). Then: Every (G,D) has a BM (which implies that G has a nowhere zero 6-flow).
* Next he proposed a conjecture [2], which is now disproven by Hoffman-Ostenhoff in [3]: Let (G,D) be a snark. Then (G,D) has two disjoint bipartizing matchings (2dBMs).
* In the same paper [3] Hoffman-Ostenhoff proposed a weakened version of the conjecture: Every snark G has at least one dominating circuit D such that (G,D) has 2dBMs.
* A natural question arises: What about 3 disjoint bipartizing matchings? Which means that any edge which is not in D lies in exactly 1 matching.

TODO: recheck this conclusion
If we still ask for D to be a dominating circuit, then there is a counterexample on 28 vertices: g2151
It seems likely, that:
* if we don't mind having 2 and more circuits, and still having the domination property (every edge lies in the cycle or is connected to it), then there are no counterexamples (checked for all snarks with upto and including 30 vertices).


[1] H. Fleischner, (Some of) the many uses of Eulerian graphs in graph theory (plus some applications), Discrete Math. 230 (2001) 23–43.

[2] H. Fleischner, Bipartizing matchings and Sabidussi’s compatibility conjecture, Discrete Math. 244 (2002) 77–82.

[3] A. Hoffman-Ostenhoff, A counterexample to the bipartizing matching conjecture, Discrete Math. 307 (2007), 2723-2733.
