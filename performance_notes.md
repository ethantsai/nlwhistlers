## Notes on Performance in Julia

## Using Static Arrays

Allocations are only expensive if they are heap allocations. A good working definition is that heap allocations are variable-sized slabs of memory which have to be pointed to, and this pointer indirection costs time. Additionally, the heap has to be managed and the garbage controllers has to actively keep track of what's on the heap.

Alternatively, we want to use stack allocations. The stack is statically-sized (known at compile time) and thus its accesses are quick. Additionally, the exact block of memory is known in advance by the compiler, and thus re-using the memory is cheap. This means that allocating on the stack has essentially no cost!

Arrays have to be heap allocated because their size (and thus the amount of memory they take up) is determined at runtime. But there are structures in Julia which are stack-allocated. structs for example are stack-allocated "value-type"s. Tuples are a stack-allocated collection. The most useful data structure for DiffEq though is the StaticArray from the package StaticArrays.jl.


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

## Managing Allocations with Broadcast Fusion

When the system is sufficiently large, or you have to make use of a non-native Julia algorithm, you have to make use of Arrays. In order to use arrays in the most efficient manner, you need to be careful about temporary allocations. Vectorized calculations naturally have plenty of temporary array allocations. We can fix this by pre-allocating the output, and mutating the braoddcast with `.=`. For example, here's a function that calculates `b`:

```
julia> lambda = [rand() for i in 1:100000];

julia> calcb(lambda::Float64) = sqrt(1+3*sin(lambda)^2)/(cos(lambda)^6);

julia> @benchmark calcb.(lambda)
BenchmarkTools.Trial: 
  memory estimate:  781.39 KiB
  allocs estimate:  5
  --------------
  minimum time:     8.309 ms (0.00% GC)
  median time:      8.546 ms (0.00% GC)
  mean time:        8.583 ms (0.05% GC)
  maximum time:     10.246 ms (0.00% GC)
  --------------
  samples:          582
  evals/sample:     1
```

And here's the one where stuff is preallocated, and the output is passed in.

```
julia> b = zeros(length(lambda));

julia> calcb2!(b::Vector{Float64},lambda::Vector{Float64}) = @. b =  sqrt(1+3*sin(lambda)^2)/(cos(lambda)^6);

julia> @benchmark calcb2!(b,lambda)
BenchmarkTools.Trial: 
  memory estimate:  0 bytes
  allocs estimate:  0
  --------------
  minimum time:     8.107 ms (0.00% GC)
  median time:      8.488 ms (0.00% GC)
  mean time:        8.554 ms (0.00% GC)
  maximum time:     10.066 ms (0.00% GC)
  --------------
  samples:          584
  evals/sample:     1
```


# Last Resorts

## fastmath

Use `@fastmath` macro to allow floating point optimizations that are usually pretty correct for real numbers.
https://docs.julialang.org/en/v1/base/math/#Base.FastMath.@fastmath


## loop vectorization

Use special vectorized functions like LoopVectorization when broadcasting trig or exp or sqrt functions. Utilizes `@avx` macro to do this wicked fast.
https://github.com/JuliaSIMD/LoopVectorization.jl