using TensorBranching
using Graphs
using OptimalBranching.OptimalBranchingMIS
using OptimalBranching.OptimalBranchingMIS.EliminateGraphs
using GenericTensorNetworks
using Random

n_sqrt = 60
filling = 0.8
for seed in 1:5
    Random.seed!(seed)
    reducer = TensorNetworkReducer()
    brancher = FixedPointBrancher()
    refiner = TreeSARefiner()
    sc_target = 25
    search_order = :bfs
    slicer_maxintersect = ContractionTreeSlicer(sc_target = sc_target, brancher = brancher, refiner = refiner, search_order = search_order, region_selector = ScoreRS(loss = :num_uniques))
    slicer_scscore = ContractionTreeSlicer(sc_target = sc_target, brancher = brancher, refiner = refiner, search_order = search_order, region_selector = ScoreRS(loss = :sc_score))

    graph = GenericTensorNetworks.random_diagonal_coupled_graph(n_sqrt, n_sqrt, filling)
    graph = SimpleGraph(graph)
    result_maxintersect = dynamic_ob_mis(graph, slicer = slicer_maxintersect, reducer = reducer, verbose = 1) 
    result_scscore = dynamic_ob_mis(graph, slicer = slicer_scscore, reducer = reducer, verbose = 1) 
    @assert result_maxintersect == result_scscore
end