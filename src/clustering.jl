function bithclust(bitstrings::Vector{BitStr{N, T}}) where {N, T}
    dm = [bdistance(b1, b2) for b1 in bitstrings, b2 in bitstrings]
    return hclust(dm, linkage = :average)
end

function clustering(bitstrings::Vector{BitStr{N, T}}, k::Int) where {N, T}
    tree = cutree(bithclust(bitstrings), k = k)
    clustered_bitstrings = Vector{Vector{BitStr{N, T}}}()
    for i in 1:k
        mask = tree .== i
        push!(clustered_bitstrings, bitstrings[mask])
    end
    return clustered_bitstrings
end