# Cycle double cover conjecture

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.51766.svg)](http://dx.doi.org/10.5281/zenodo.51766)
[![Join the chat at https://gitter.im/gexahedron/cycle-double-covers](https://badges.gitter.im/gexahedron/cycle-double-covers.svg)](https://gitter.im/gexahedron/cycle-double-covers?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://travis-ci.org/gexahedron/cycle-double-covers.svg?branch=master)](https://travis-ci.org/gexahedron/cycle-double-covers)
[![License](http://img.shields.io/:license-mit-blue.svg)](http://doge.mit-license.org)

Here you will find code (in Python and C++) for different unsolved problems related to cycle double cover conjecture.

## Main references
- [Cycle double cover conjecture](http://www.openproblemgarden.org/op/cycle_double_cover_conjecture)
- [5-flow conjecture](http://www.openproblemgarden.org/op/5_flow_conjecture)
- [The Berge-Fulkerson conjecture](http://www.openproblemgarden.org/op/the_berge_fulkerson_conjecture)
- [(m,n)-cycle covers](http://www.openproblemgarden.org/op/m_n_cycle_covers)
- [Bouchet's 6-flow conjecture](http://www.openproblemgarden.org/op/bouchets_6_flow_conjecture)
- [Unit vector flows](http://www.openproblemgarden.org/op/unit_vector_flows)
- [The three 4-flows conjecture](http://www.openproblemgarden.org/op/three_4_flows_conjecture)

# Upcoming code for:
* oriented 5-cycle double cover conjecture (o5cdc)
* strong Petersen colouring conjecture
* (brand new!) oriented 6-cycle 4-cover conjecture (o6c4c) (what's new here is just the 'oriented' part; verified for all snarks upto and including 30 vertices)
* oriented 9-cycle 6-cover conjecture (o9c6c) (but not for Petersen graph, because it doesn't have any 9-cycle 6-cover and it's kind of easy to understand why (TODO: write the proof): any solution for o9c6c consists of 9 layers, which are complementary to perfect matchings; if we go the complementary problem then it says that we need to cover the graph with 9 perfect matchings which cover the graph 3 times; we have 6 different perfect matchings for Petersen graph (which are all used in o6c4c solution), and each edge lies in exactly 2 of them; if we don't use any one of the perfect matchings, then we must use some of them 3 times, which is problematic (because we already cover 5 edges 3 times with it, we are left with 10 edges and no perfect matchings which can cover them); which means that 3 of perfect matchings are used once and 3 of them are used twice which means that there are some edges which are covered 2 or 4 times)
* investigations into common nowhere-zero 5-flows which come from 2BMs, 3BMs, o5cdc, o6c4c
* maybe some code for conjectures related to signed graphs

# Also upcoming:
* oriented 10-cycle 6-cover for Petersen graph (there's plenty of solutions, actually)
* diagram of all related theorems and conjectures, which i could find in the literature

# Actually, one of 10-cycle 6-covers for Petersen graph
Printing graph: 

0 : 4 6 8

1 : 5 6 9

2 : 4 7 9

3 : 5 7 8

4 : 0 2 5

5 : 1 3 4

6 : 0 1 7

7 : 2 3 6

8 : 0 3 9

9 : 1 2 8

cycle 0 (colour 0): 0 4 2 7 3 5 1 6 

cycle 1 (colour 1): 0 4 2 7 3 5 1 9 8 

cycle 2 (colour 2): 0 4 2 7 3 8 9 1 6 

cycle 3 (colour 3): 0 4 5 1 6 7 2 9 8 

cycle 4 (colour 4): 0 4 5 3 7 6 1 9 8 

cycle 5 (colour 5): 0 4 5 3 8 9 2 7 6 

cycle 6 (colour 6): 0 6 7 3 8 

cycle 7 (colour 6): 1 5 4 2 9 

cycle 8 (colour 7): 0 6 1 9 2 4 5 3 8 

cycle 9 (colour 8): 0 6 7 2 9 1 5 3 8 

cycle 10 (colour 9): 1 5 4 2 9 8 3 7 6 

## License

See the [LICENSE](LICENSE.md) file for license rights and limitations (MIT).
