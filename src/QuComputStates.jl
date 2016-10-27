module QuComputStates

import QuBase: AbstractQuVector,AbstractFiniteBasis,FiniteBasis,Orthonormal,rawcoeffs,coefftype,rawbases,show,similar_type,bases
import Base: copy
export ComputState, CHPState, rawcoeffs, rawbases, AbstractQC, AbstractQuCircuit, ComputBasis, comput_basis, bases, BellState

abstract AbstractQC{N}
abstract AbstractQuState{T,N}
abstract AbstractQuCircuit{N} <: AbstractQC{N}

include("computbasis.jl")
include("computstate.jl")
include("stlzstate.jl")
include("graphstate.jl")

# This file is for types in chp.c

type CHPState{N}<:AbstractQuVector{N}
    X::AbstractMatrix # (2n+1)*n matrix for stabilizer/destabilizer x bits
    Z::AbstractMatrix # (2n+1)*n matrix for z bits
    R::AbstractVector # phase bits: 0 for +1, 1 for i, 2 for -1, 3 for -i. Normally either 0 or 2
end

CHPState(X::BitMatrix,Z::BitMatrix,R::BitVector) = CHPState{length(R)/2}(X,Z,R)

function CHPState(n::Int)
    X = spzeros(Bool,2n,n)
    Z = spzeros(Bool,2n,n)
    R = spzeros(Bool,2n)

    X[1:n,:] = speye(n)
    Z[n+1:2n,:] = speye(n)

    return CHPState{n}(X,Z,R)
end

function copy!(A::CHPState,B::CHPState)
    copy!(A.X,B.X)
    copy!(A.Z,B.Z)
    copy!(A.R,B.R)

    return A
end

end # module
