module CUDAExt

using CUDA
import TensorBranching: item, gpu_tensors

item(t::CuArray{T, 0}) where {T} = Array(t)[]

function gpu_tensors(tensors::Vector{T}) where {T}
    return CuArray.(tensors)
end

end