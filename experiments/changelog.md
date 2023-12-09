Milestones:
* added Graph struct
* added fast Petersen colouring
* added cc-mapping
* added cycles
* added Hoffman-Ostenhof decomposition
* added check for preimage Hoffman-Ostenhof solutions (FIXME needs recheck)
* added new code for o6c4c
* added new code for o5cdc (FIXME there's bug somewhere)
* added new code for a new construction, which includes 6c4c (244-flows) and 33-pp (which includes nz5 and 5cdc)
* added code for 33-pp and 333-pp flows
* added counts for circuits (in the list of cycles, full cycles, even cycles)
* compared 33-pp and 5cdc
* added code for 333-flows from 6c4c
* added code for o244-flows from 6c4c
* added check for same triple 33-pp and o244-flows
* added code for dominating circuits
* added code for poor-rich edges from 6c4c
* added code for oriented vertices from o6c4c
* added code for poor-rich edges from o6c4c (4 different types of edges)
* compared o6c4c and petersen colouring (FIXME needs recheck)
* added code for (6c4c, 5cdc) pairs from 6c4c-244-flows-33-pp construction
* added code for (6c4c, 5cdc) pairs from petersen colouring (FIXME needs recheck)
* released code to github
* removed input parameter with number of vertices
* added nz5 code
* added code for generating full cycles from nz5
* added flows to 33-pp (FIXME has bug)
* added code for 6c4c from 33-pp
* found out that in nz-mod5 - number of vertices with 1+1+3 equals number of vertices with 4+4+2, and same for pair of 2+2+1 and 3+3+4
* refactoring (split into separate files)

TODO:
* fix bugs in 33-pp
* now that i have (33-pp from 6c4c), (6c4c from 33-pp), (petersen colouring without bug) - i need again to match and compare them all
* search for dominating circuits in nz5 (orientations_from_nz5)
* rewrite oriented-244-flows with masks
* add o9c6c
* add o2233-flows
* rewrite 5cdc so that i won't have duplicates
* fix bug in 5cdc
* refactor code
* combine code for 5cdc and 6c4c and add support for 9c6c
* check code on some graphs from graph6 files
* print various statistics by cycles, by 5cdc, by 6c4c, by nz5, by nz-mod5
* generalize code not only for cubic graphs? won't be easy for some of the constructions
