

# https://github.com/vincent-picaud/DirectConvolution.jl
# https://pixorblog.wordpress.com/2016/07/13/savitzky-golay-filters-julia/
function smooth(signal, filterwidth::Int, polydegree::Int)
    s = Float64[i for i in signal]
    sg = SG_Filter(halfWidth=filterwidth,degree=polydegree)
    ss = apply_SG_filter(s, sg)
    return ss # 1d savitzky-golay smoothed signal
end