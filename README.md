# TensorBranching.jl

[![Build Status](https://github.com/ArrogantGao/TensorBranching.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ArrogantGao/TensorBranching.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://app.travis-ci.com/ArrogantGao/TensorBranching.jl.svg?branch=main)](https://app.travis-ci.com/ArrogantGao/TensorBranching.jl)
[![Coverage](https://codecov.io/gh/ArrogantGao/TensorBranching.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ArrogantGao/TensorBranching.jl)


## How to use

### Basic Data Structures

- `BitStr`: A bit string, for more details see [BitBasis.jl](https://yaoquantum.org/BitBasis.jl/dev/man.html).
```julia
julia> a = bit"10101"
10101 ₍₂₎
```

- `Clause`: store both `mask` and `val` as bit strings, the $i$-th bit of `mask` is `1` if and only if the $i$-th bit of `val` is valid.
```julia
julia> Clause(bit"1110", bit"1010")
Clause{4, 3, Int64}: mask: 1110 ₍₂₎, val: 1010 ₍₂₎
```

- BranchingTable{N, C}:

A table of branching configurations. The table is a vector of vectors of `StaticBitVector{N, C}`. Type parameters are:
- `N`: The number of bits in the bit vector.
- `C`: The number of integers as the storage.


```jldoctest
julia> graph_sat = graph_from_tuples(5, [(1, 2), (2,3), (2,4), (1,3), (3, 4), (4, 5), (2,5)])
{5, 7} undirected simple Int64 graph

julia> tbl = reduced_alpha_configs(graph_sat, [1, 4, 5])
BranchingTable{N}
00100, 01000
10010
00101
```


- `SubCover`: store the clause, and record the truth table cover by the clause, and `n_rm` is a number record the vertices removed by the clause.
```julia
julia> graph_sat = graph_from_tuples(5, [(1, 2), (2,3), (2,4), (1,3), (3, 4), (4, 5), (2,5)])
{5, 7} undirected simple Int64 graph

julia> vertices = [1, 2, 3, 4, 5];

julia> tbl = reduced_alpha_configs(graph_sat, [1, 4, 5])
BranchingTable{N}
00100, 01000
10010
00101

julia> bss = Tbl2BitStrs(tbl);

julia> sc = subcovers(bss, vertices, graph_sat)
9-element Vector{SubCover{5, Int64}}:
 SubCover{5, Int64}: ids: Set([2, 1]), mask: 10100 ₍₂₎, val: 00000 ₍₂₎, n_rm: 2
 SubCover{5, Int64}: ids: Set([2, 1]), mask: 10010 ₍₂₎, val: 00000 ₍₂₎, n_rm: 2
 ......
 SubCover{5, Int64}: ids: Set([3]), mask: 11111 ₍₂₎, val: 10100 ₍₂₎, n_rm: 5
```

### load graph from artifacts

The graphs provided by PACE2019 have been preprocessed and saved as artifacts. You can load them using the following code:
```julia
julia> using TensorBranching

julia> graph_from_artifact(1)
{6160, 40207} undirected simple Int64 graph
```
where `1` is the graph id. The graph id can be found in the [PACE2019](https://pacechallenge.org/2019/vc/index) website.


### solve maximum independent set problem

Currently we provide a simple interface to solve the maximum independent set problem. You can use the following code:
```julia
julia> g = random_regular_graph(60, 3)
{60, 90} undirected simple Int64 graph

julia> missolve(g, strategy = SetCoverBranching(), kneighbor = 1, show_count = true)
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

julia> tree, mis = branching_tree(g, SetCoverBranching(), 1)
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