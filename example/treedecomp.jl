using TensorBranching
using TreeWidthSolver
using TensorBranching.OMEinsum, TensorBranching.GenericTensorNetworks, TensorBranching.AbstractTrees
using TensorBranching.GenericTensorNetworks.Graphs

using Plots

g0 = random_regular_graph(150, 3)

prob = GenericTensorNetwork(IndependentSet(g0); optimizer=TreeSA(; ntrials=1, niters=10))
code = prob.code.eins

cc = contraction_complexity(code, uniformsize(code, 2))
tree = decompose(code)

treebags = [node.bag for node in PostOrderDFS(tree)]

length(treebags)
bagsizes = sort(length.(treebags))

histogram(bagsizes)

