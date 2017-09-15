#[The three 4-flows conjecture](http://www.openproblemgarden.org/op/three_4_flows_conjecture)

##How to run the code
To compile the code just print "make" in your console (you need g++ compiler).

To run this code you need to specify the number of vertices in the snarks and provide the file with snarks in multicode format.

##How to generate the file with snarks
Download the code of snarkhunter [2] from [here](http://caagt.ugent.be/cubic/) (and for this you would also need to follow the installation instructions, which include downloading [nauty](http://pallini.di.uniroma1.it/) (you can also find older versions here: http://cs.anu.edu.au/~bdm/nauty/) and copying files from nauty [1]).

Then you can generate snarks in multicode format, e. g., if you want cyclically 4-edge-connected snarks with girth 4 and more (I will denote them as 22.04g) with 22 vertices:

`./snarkhunter-64 22 4 s S C4`

or

`./snarkhunter 22 4 s S C4`

If you run on a 32-bit computer. You'll get a file named something like "Generated_graphs.22.04.sn.cyc4".

If you want cyclically 4-edge-connected snarks with girth 5 and more (I will denote them as 22.05g):

`./snarkhunter-64 22 5 s S C4`

You'll get a file named something like "Generated_graphs.22.05.sn.cyc4".

Use this file as an input to the program:

`./three_4_flows < Generated_graphs.22.05.sn.cyc4`

##Main result here (May 20, 2016)
The conjecture about 233-flows, stated in the last section in [openproblemgarden](http://www.openproblemgarden.org/op/three_4_flows_conjecture) text is already false for Petersen graph (234-flows and 235-flow also don't exist for Petersen graph).

Looks like there do exist 333-flows and 244-flows for all snarks (verified for all snarks with 10, 18, 20, 22, 24, 26 vertices).

Also looks like almost all snarks have 233-flows. Here are the exceptions (girth 5 snarks):

* 10 vertices: g1 (no 23k-flows)
* 18 vertices: none
* 20 vertices: g1 (no 23k-flows)
* 22 vertices: g6 (no 233-flows)
* 24 vertices: g7, g12, g26, g29 (no 233-flows)
* 26 vertices: g62, g67, g88, g89, g93, g98, g109, g119, g138, g143, g166, g177, g189, g191, g203, g246, g277 (no 233-flows); g141 (no 23k-flows)

##References
[1] McKay, B.D. and Piperno, A., Practical Graph Isomorphism, II, Journal of Symbolic Computation, 60 (2014), pp. 94-112.

[2] G. Brinkmann, J. Goedgebeur and B.D. McKay, Generation of Cubic graphs, Discrete Mathematics and Theoretical Computer Science, 13(2):69-80, 2011.
