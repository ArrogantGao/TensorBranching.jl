```@meta
CurrentModule = TensorBranching
```

# TensorBranching.jl

Documentation for [TensorBranching.jl](https://github.com/ArrogantGao/TensorBranching.jl).

## Quick start

### Install

Press `]` to enter the package mode in Julia Repl, and type the following command to install the package:

```julia
pkg> add https://github.com/ArrogantGao/TensorBranching.jl
```

### solve the maximum independent set problem

Currently we provide a simple interface to solve the maximum independent set problem. You can use the following code:
```julia
julia> using TensorBranching

julia> branching_strategy = SetCoverBranching()
SetCoverBranching(5, IPSetCoverSolver(), false)

julia> vertex_selector = MinBoundarySelector(2)
MinBoundarySelector(2)

julia> measure = D3Measure()
D3Measure()

julia> table_filter = EnvFilter()
EnvFilter()

julia> cfg = SolverConfig(; branching_strategy, vertex_selector, measure, table_filter)
SolverConfig{TensorNetworkSolver, SetCoverBranching, D3Measure, MinBoundarySelector, EnvFilter}(TensorNetworkSolver(), SetCoverBranching(3), D3Measure(), MinBoundarySelector(2), EnvFilter())

julia> g = smallgraph(:tutte)
{46, 69} undirected simple Int64 graph

julia> solve_mis(g, cfg)
19

julia> count_mis(g, cfg)
(19, 5)
```

One can use the `verbose` option to see the detailed information of the set covering process.
```julia
julia> branching_strategy = SetCoverBranching(max_itr = 1, verbose = true)
SetCoverBranching(5, IPSetCoverSolver(), true)

julia> v = 1; k = 1;

julia> vs,openvertices = TensorBranching.neighbor_cover(g,v,k);

julia> subg, subg_vs = induced_subgraph(g,vs);

julia> openvertices = findall(x -> x != 3, degree(subg));

julia> tbl = reduced_alpha_configs(table_solver, subg, openvertices);

julia> TensorBranching.impl_strategy(g, vs, tbl, branching_strategy, measure);
[ Info: γ0 = 1.1054326208163803
feasible solution found by trivial heuristic after 0.0 seconds, objective value 5.404376e+00
presolving:
   (0.0s) running MILP presolver
   (0.0s) MILP presolver found nothing
(round 1, exhaustive) 0 del vars, 0 del conss, 0 add conss, 0 chg bounds, 0 chg sides, 0 chg coeffs, 5 upgd conss, 0 impls, 0 clqs
   (0.0s) probing cycle finished: starting next cycle
   (0.0s) symmetry computation started: requiring (bin +, int -, cont +), (fixed: bin -, int +, cont -)
   (0.0s) symmetry computation finished: 1 generators found (max: 1500, log10 of symmetry group size: 0.3) (symcode time: 0.00)
(round 2, exhaustive) 0 del vars, 0 del conss, 1 add conss, 0 chg bounds, 0 chg sides, 0 chg coeffs, 5 upgd conss, 0 impls, 0 clqs
presolving (3 rounds: 3 fast, 3 medium, 3 exhaustive):
 0 deleted vars, 0 deleted constraints, 1 added constraints, 0 tightened bounds, 0 added holes, 0 changed sides, 0 changed coefficients
 0 implications, 0 cliques
presolved problem has 15 variables (15 bin, 0 int, 0 impl, 0 cont) and 6 constraints
      1 constraints of type <orbitope>
      5 constraints of type <logicor>
Presolving Time: 0.00
transformed 1/1 original solutions to the transformed problem space

 time | node  | left  |LP iter|LP it/n|mem/heur|mdpt |vars |cons |rows |cuts |sepa|confs|strbr|  dualbound   | primalbound  |  gap   | compl. 
i 0.0s|     1 |     0 |     0 |     - |  oneopt|   0 |  15 |   6 |   5 |   0 |  0 |   0 |   0 | 0.000000e+00 | 1.000000e+00 |    Inf | unknown
* 0.0s|     1 |     0 |     5 |     - |    LP  |   0 |  15 |   6 |   5 |   0 |  0 |   0 |   0 | 8.986148e-01 | 8.986148e-01 |   0.00%| unknown
  0.0s|     1 |     0 |     5 |     - |   785k |   0 |  15 |   6 |   5 |   0 |  0 |   0 |   0 | 8.986148e-01 | 8.986148e-01 |   0.00%| unknown

SCIP Status        : problem is solved [optimal solution found]
Solving Time (sec) : 0.00
Solving Nodes      : 1
Primal Bound       : +8.98614789194314e-01 (3 solutions)
Dual Bound         : +8.98614789194314e-01
Gap                : 0.00 %
[ Info: IP Solver, Iteration 1, complexity = 1.0952137194262848
```

The standard `mis1` and `mis2` strategies are also available.
```Julia
julia> using TensorBranching

julia> counting_mis2(g)
CountingMIS(19, 139)

julia> counting_mis1(g)
CountingMIS(19, 1697552)
```

### Branching Tree

To check the branching process in detail, we store all the branching process as tree.
```Julia
julia> tree, c = branching_tree(g, cfg)
(G{46}
├─ G{41}
│  └─ G{39}
│     └─ G{37}
│        └─ G{35}
│           └─ G{32}
│              ⋮
│              
└─ G{42}
   └─ G{41}
      └─ G{39}
         └─ G{37}
            └─ G{35}
               ⋮
               
, 4)
```
and the graph of current state, the vertices to be removed are stored:
```Julia
julia> tree.graph
{46, 69} undirected simple Int64 graph

julia> tree.removed
3-element Vector{Vector{Int64}}:
 [2, 1, 5, 27]
 [2, 1, 3, 4, 5, 6, 34, 27, 26]
 [2, 1, 5, 6, 34, 27, 26, 3, 11, 12, 4, 19, 20]
```
the vertices to be removed correspond to the children of the current node.
