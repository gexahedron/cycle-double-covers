# Computations around Cycle double cover conjecture (and it's friends)

How to run experiments, example:
`./run_experiments mc < ../../multicode/Generated_graphs.10.05.sn.cyc4`

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.51766.svg)](http://dx.doi.org/10.5281/zenodo.51766)
[![Join the chat at https://gitter.im/gexahedron/cycle-double-covers](https://badges.gitter.im/gexahedron/cycle-double-covers.svg)](https://gitter.im/gexahedron/cycle-double-covers?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
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

# TODO:

* oriented 5-cycle double cover conjecture (o5cdc)
* strong Petersen colouring conjecture
* (brand new!) oriented 6-cycle 4-cover conjecture (o6c4c, we can also call it as 'oriented Berge-Fulkerson') (what's new here is just the 'oriented' part; verified for all snarks upto and including 30 vertices)
* oriented 9-cycle 6-cover conjecture (o9c6c) (but not for Petersen graph, because it doesn't have any 9-cycle 6-cover and it's kind of easy to understand why (TODO: write the proof): any solution for 9c6c consists of 9 layers, which are complementary to perfect matchings; if we go the complementary problem then it says that we need to cover the graph with 9 perfect matchings which cover the graph 3 times; we have 6 different perfect matchings for Petersen graph (which are all used in o6c4c solution), and each edge lies in exactly 2 of them; if we don't use any one of the perfect matchings, then we must use some of them 3 times, which is problematic (because we already cover 5 edges 3 times with it, we are left with 10 edges and no perfect matchings which can cover them); which means that 3 of perfect matchings are used once and 3 of them are used twice which means that there are some edges which are covered 2 or 4 times)
* investigations into common nowhere-zero 5-flows which come from 2BMs, 3BMs, o5cdc, o6c4c
* maybe some code for conjectures related to signed graphs
* Updating diagram, and it also needs some helpful notation

## License

See the [LICENSE](LICENSE.md) file for license rights and limitations (MIT).

## Diagram of conjectures (outdated and has some bugs)

![Diagram of conjectures](/images/full_scheme.png)
