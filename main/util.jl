

logrange(x1, x2, n::Int64) = [10^y for y in range(log10(x1), log10(x2), length=n)]

obtain_normalizer(f::Function) = maximum(f.(0:0.01:90))^-1
obtain_normalizer(f::Function,t) = maximum(f.(0:0.01:90).-f(t))^-1

"""
    last_index_before_zero(x::Vector{Float64})

Takes in a vector of floats, returns the index of the last non-zero value

# Example
```julia-repl
julia> last_index_before_zero([1,0,3,0,5,0,0,0.])
5
```
"""
last_index_before_zero(x::Vector{Float64}) = findlast(x->!iszero(x),x)

"""
    smooth(signal, filterwidth::Int, polydegree::Int)

Smooths a function using the savitzky golay filter.
Specify filderwidth and polynomial order.

# Examples
```julia-repl
x = 0:0.001:2π
r = rand()
y = @. (sin(10*r*x) + 0.01 * randn()) * (4*r*sin(2.5*r*x) + 0.01 * rand()) * (3*sin(x) + 0.01 * randn())
plot(x,[y, smooth(y, 100, 5)], label=["original" "smoothed"])
```
"""
# https://github.com/vincent-picaud/DirectConvolution.jl
# https://pixorblog.wordpress.com/2016/07/13/savitzky-golay-filters-julia/
function smooth(signal, filterwidth::Int, polydegree::Int)
    s = Float64[i for i in signal]
    sg = SG_Filter(halfWidth=filterwidth,degree=polydegree)
    ss = apply_SG_filter(s, sg)
    return ss # 1d savitzky-golay smoothed signal
end