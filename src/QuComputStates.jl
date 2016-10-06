module QuComputStates

import QuBase: AbstractQuVector,AbstractFiniteBasis,FiniteBasis,Orthonormal,rawcoeffs,coefftype,rawbases,show,similar_type,bases
import Base: copy
export ComputState, CHPState, rawcoeffs, rawbases, AbstractQC, AbstractQuCircuit, ComputBasis, comput_basis, bases, BellState

abstract AbstractQC{N}
abstract AbstractQuCircuit{N} <: AbstractQC{N}

# general computing basis
typealias ComputBasis FiniteBasis{Orthonormal}
comput_basis(n::Int) = FiniteBasis(ntuple(x->2,n))

type ComputState{A<:AbstractVector,B<:AbstractFiniteBasis,T,N} <: AbstractQuVector{B,T}
    coeffs::A
    bases::B
end

# constructors
ComputState{A<:AbstractVector,B<:AbstractFiniteBasis}(state_vec::A,bases::B) = ComputState{A,B,eltype(A),length(state_vec)|>log2|>Int}(state_vec,bases)
ComputState{A<:AbstractVector}(state_vec::A) = ComputState{A,ComputBasis,eltype(A),length(state_vec)|>log2|>Int}(state_vec,state_vec|>length|>comput_basis)

# Hints
BellState(n::Int) = ComputState(Complex128[1/sqrt(2^n) for i=1:2^n],comput_basis(n))

# QuBase functions overload
coefftype{A,B,T,N}(state::ComputState{A,B,T,N}) = A
rawcoeffs(state::ComputState) = state.coeffs
rawbases(state::ComputState) = state.bases
rawbases(state::ComputState, i) = state.bases

# Base functions overload
copy(state::ComputState) = ComputState(copy(state.coeffs),copy(state.bases))
similar_type{Q<:ComputState}(::Type{Q}) = ComputState


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

# Graph States

end # module
