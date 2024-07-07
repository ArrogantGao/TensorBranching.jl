```@meta
CurrentModule = TensorBranching
```

# TensorBranching

Documentation for [TensorBranching](https://github.com/ArrogantGao/TensorBranching.jl).

## Quick start

### Install

Press `]` to enter the package mode in Julia Repl, and type the following command to install the package:

```julia
pkg> add https://github.com/ArrogantGao/TensorBranching.jl
```

### solve the maximum independent set problem

Currently we provide a simple interface to solve the maximum independent set problem. You can use the following code:
```julia
julia> g = random_regular_graph(60, 3)
{60, 90} undirected simple Int64 graph

julia> missolve(g, strategy = SetCoverBranching(), kneighbor = 1, show_count = true, usr_rv = true)
(26, 206)
```
where `kneighbor` for layer of neighbor considered, the first result is the size of the maximum independent set, and the second result is the number of branches visited.

The standard `mis1` and `mis2` strategies are also available.
```Julia
julia> using TensorBranching: counting_mis1, counting_mis2

julia> using EliminateGraphs

julia> counting_mis1(EliminateGraph(g))
CountingMIS(26, 58148786)

julia> counting_mis2(EliminateGraph(g))
CountingMIS(26, 188)
```

### Branching Tree

To check the branching process in detail, we store all the branching process as tree.
```Julia
julia> g = random_regular_graph(20, 3)
{20, 30} undirected simple Int64 graph

julia> tree, mis = branching_tree(g, SetCoverBranching(), 1, use_rv = true)
(G{20}
├─ G{16}
│  └─ G{14}
│     ├─ G{11}
│     │  └─ G{9}
│     │     └─ G{7}
│     │        ⋮
│     │        
│     └─ G{7}
│        └─ G{5}
│           └─ G{3}
│              ⋮
│              
├─ G{12}
│  └─ G{10}
│     └─ G{8}
│        ├─ G{5}
│        │  └─ G{3}
│        │     ⋮
│        │     
│        └─ G{3}
│           └─ G{0}
└─ G{16}
   └─ G{14}
      ├─ G{11}
      │  └─ G{9}
      │     └─ G{7}
      │        ⋮
      │        
      └─ G{8}
         └─ G{6}
            └─ G{4}
               ⋮
               
, 6)
```
and the graph of current state, the vertices to be removed are stored:
```Julia
julia> tree.graph
{20, 30} undirected simple Int64 graph

julia> tree.removed
3-element Vector{Vector{Int64}}:
 [1, 5, 13, 20]
 [1, 5, 13, 10, 15, 20, 4, 16]
 [1, 5, 12, 16]
```
the vertices to be removed correspond to the children of the current node.

## Detailed documentation

```@index
```
