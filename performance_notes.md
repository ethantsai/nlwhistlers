## Notes on Performance in Julia

## Using Static Arrays



## Types with values-as-parameters

Say you have a function parameter with unknown dimension AND you want to run a function inside that function that uses that ambiguous parameter. In this case, let's use a 1d array with unknown length. The function inside will be slow so you can pass in the length or type as a parameter, thereby fully eliminating type instability. For example, the syntax looks like:

```
function filter3(A::AbstractArray{T,N}) where {T,N}
    kernel = array3(1, Val(N))
    filter(A, kernel)
end
```

which can take it any type or vector size and still be very happy.

## Function barriers

If you can't avoid type instability, you can hide it with function barriers, which can still avoid allocations. Julia's compiler specializes code for argument types at function boundaries, so 



## fastmath

Use `@fastmath` macro to allow floating point optimizations that are usually pretty correct for real numbers.
https://docs.julialang.org/en/v1/base/math/#Base.FastMath.@fastmath


## loop vectorization

Use special vectorized functions like LoopVectorization when broadcasting trig or exp or sqrt functions. Utilizes `@avx` macro to do this wicked fast.
https://github.com/JuliaSIMD/LoopVectorization.jl