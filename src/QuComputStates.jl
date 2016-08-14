module QuComputStates

import QuBase: AbstractQuVector,AbstractFiniteBasis,FiniteBasis,Orthonormal,rawcoeffs,coefftype,rawbases,show
import Base.copy
export ComputState, stlzState, rawcoeffs, rawbases

# general computing basis
typealias ComputBasis FiniteBasis{Orthonormal}
comput_basis(n::Int) = FiniteBasis(ntuple(x->2,n))

type ComputState{A<:AbstractVector,B<:AbstractFiniteBasis,T,N} <: AbstractQuVector{B,T}
    coeffs::A
    bases::B
end

# constructors
ComputState{A<:AbstractVector}(state_vec::A) = ComputState{A,ComputBasis,eltype(A),length(state_vec)|>log2|>Int}(state_vec,state_vec|>length|>comput_basis)

# QuBase functions overload
coefftype{A,B,T,N}(state::ComputState{A,B,T,N}) = A
rawcoeffs(state::ComputState) = state.coeffs
rawbases(state::ComputState) = state.bases

# Base functions overload
copy(state::ComputState) = ComputState(copy(state.coeffs),copy(state.bases))
similar_type{Q<:ComputState}(::Type{Q}) = ComputBasis


# This file is for types in chp.c

type stlzState{N}<:AbstractQuVector{N}
    X::AbstractMatrix # (2n+1)*n matrix for stabilizer/destabilizer x bits
    Z::AbstractMatrix # (2n+1)*n matrix for z bits
    R::AbstractVector # phase bits: 0 for +1, 1 for i, 2 for -1, 3 for -i. Normally either 0 or 2
end

stlzState(X::BitMatrix,Z::BitMatrix,R::BitVector) = stlzState{length(R)/2}(X,Z,R)

function stlzState(n::Int)
    X = spzeros(Bool,2n,n)
    Z = spzeros(Bool,2n,n)
    R = spzeros(Bool,2n)

    X[1:n,:] = speye(n)
    Z[n+1:2n,:] = speye(n)

    return stlzState{n}(X,Z,R)
end

function copy!(A::stlzState,B::stlzState)
    copy!(A.X,B.X)
    copy!(A.Z,B.Z)
    copy!(A.R,B.R)

    return A
end

end # module
