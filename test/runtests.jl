using QuComputStates
using Base.Test

t = ComputState([1,2,3,4])
@test rawcoeffs(t) == [1,2,3,4]

# write your own tests here
@test 1 == 1
