# TensorBranching.jl

[![Coverage](https://codecov.io/gh/ArrogantGao/TensorBranching.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ArrogantGao/TensorBranching.jl)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ArrogantGao.github.io/TensorBranching.jl/dev/)

This package is for combine the branching algorithm and tensor network methods to solve the maximum independent set problem.

## Usage

```julia
using TensorBranching
```

For fast contraction on CPU, use
```julia
using TropicalGEMM
```

For fast contraction on GPU, use
```julia
using CUDA, CuTropicalGEMM
```

## Example

```julia
julia> using TensorBranching, Graphs, TropicalGEMM

julia> g = random_regular_graph(100, 3)
{100, 150} undirected simple Int64 graph

julia> slicer = ContractionTreeSlicer(sc_target = 10)
ContractionTreeSlicer(10, :sc_score, 100.0:1.0:100.0, 1, 100, MaxIntersectRS(20, :mincut, :num_uniques), OptimalBranchingMIS.TensorNetworkSolver(true), OptimalBranchingCore.IPSolver(HiGHS.Optimizer, 20, 2.0, false))

julia> dynamic_ob_mis(g, slicer = slicer)
[ Info: kernelized graph: {84, 133} simple graph
[ Info: branches: 2
44

julia> using EliminateGraphs

julia> mis2(EliminateGraph(g))
44

```
