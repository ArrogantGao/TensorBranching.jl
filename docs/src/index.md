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
SetCoverBranching(3)

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
