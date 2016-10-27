# general computing basis
immutable ComputBasis{N}<:AbstractFiniteBasis{Orthonormal}
end

ComputBasis(n::Int) = ComputBasis{n}()
ComputBasis(n::Int...) = ComputBasis(sum(n))

bases{T,N}(::AbstractQuState{T,N}) = ComputBasis(N)
