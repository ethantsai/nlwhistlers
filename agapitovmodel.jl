using StaticArrays, Plots

const α_ij_MLT_234_L15 = SMatrix{5,5}([0 0           1             2             3          ;
                                       1 4.1877e−1   −9.0608e−3    5.2279e−2     −6.5120e−3 ;
                                       1 1.5863      1.4375e−1     −9.3020e−2    9.3822e−3  ;
                                       1 7.4082e−2   −6.8803e−3    6.3612e−3     −6.8278e−4 ;
                                       1 −2.4300     −2.1121       3.6348e−1     −2.0494e−2 ])

const α_ij_MLT_234_L56 = SMatrix{5,5}([0 0           1             2             3          ;
                                       1 3.2203e−1   8.2359e−3     4.6659e−2     −5.4544e−3 ;
                                       2 1.7732      6.0086e−2     −8.5062e−2    9.6633e−3  ;
                                       3 4.5626e−2   9.2193e−3     3.4485e−3     −4.2794e−4 ;
                                       4 −5.7878     −5.2970       1.7597        −1.5317e−1 ])

const α_ij_MLT_412_L15 = SMatrix{5,5}([0 0           1             2             3          ;
                                       1 5.5884e−1   −1.3011e−1    7.9472e-002   −9.2103e−3 ;
                                       2 1.4989      4.5649e−2     −2.1041e−2    3.7914e−4  ;
                                       3 9.9792e−2   −3.3988e−2    1.0392e−2     −9.0381e−4 ;
                                       4 −1.9441     −3.3455       8.5083e−1     −7.2328e−2 ])

const α_ij_MLT_412_L56 = SMatrix{5,5}([0 0           1             2             3          ;
                                       1 2.4510e−1   9.6416e−2     −8.3788e−3    3.5775e−4  ;
                                       2 1.7767      −6.0993e−2    2.4017e−3     −5.5076e−5 ;
                                       3 2.7606e−2   5.9185e−3     −1.3959e−4    −4.0953e−5 ;
                                       4 −1.2235e1   −1.2452       2.8956e−1     −1.4832e−2 ])

const α_ij_MLT_1223_L15 = SMatrix{5,5}([0 0           1            2            3          ;
                                        1 4.8179e−1   2.2206e−2    −2.1949e−2   4.0690e−3  ;
                                        3 1.5997      −9.9165e−2   7.1753e−2    −1.0798e−2 ;
                                        4 6.1269e−2   −1.2257e−3   −3.7738e−3   7.7571e−4  ;
                                        5 −2.7227     −1.2788      −4.9845e−1   1.2404e−1  ])

const α_ij_MLT_1223_L56 = SMatrix{5,5}([0 0          1             2            3          ;
                                        1 2.1277e−1  2.0086e−1     −8.4907e−2   1.1070e−2  ;
                                        2 1.7252     −9.2292e−2    5.7363e−2    −9.2312e−3 ;
                                        3 2.4257e−2  1.0763e−2     −5.4592e−3   9.8882e−4  ;
                                        4 −9.2179    −5.3855e−1    1.7238e−1    2.0337e−2  ])


function α_ij_matrix(L, MLT)::SMatrix{5, 5, Float64, 25}
    if L >= 5
        if MLT >= 18 || MLT < 4 # since 12-23 is dominated by dusk flank, do not use those coeffs for even 18-23
            return α_ij_MLT_234_L56
        elseif 4 <= MLT <= 15 # use dawn-day for even up to 12-15
            return α_ij_MLT_412_L56
        else
            return α_ij_MLT_1223_L56
        end
    elseif L < 5
       if MLT > 18 || MLT < 4
        return α_ij_MLT_234_L15
       elseif 4 <= MLT <= 15
        return α_ij_MLT_412_L15
       else
        return α_ij_MLT_1223_L15
       end
    end
end
b_i(i, Kp, α_ij::SMatrix{5, 5, Float64, 25})::Float64 = sum( α_ij[i+2,2:end] .* ( Kp .^ α_ij[1,2:end] ))
B_w(lambda, Kp, α_ij::SMatrix{5, 5, Float64, 25}) = 10 ^ abs( b_i(0,Kp,α_ij) * (abs(lambda) - b_i(3,Kp,α_ij)) * exp(-abs(lambda) * b_i(2,Kp,α_ij) - b_i(1,Kp,α_ij)))

function agapitov_coeffs(Kp, α_ij::SMatrix{5, 5, Float64, 25})
    return @SArray [b_i(0,Kp,α_ij), b_i(1,Kp,α_ij), b_i(2,Kp,α_ij), b_i(3,Kp,α_ij)]
end

#Ethan's Mods
# forced B_w to be zero at lambda = 0
# forced it to be an odd function as well
# should I normalize this funciton?


# Usage, define an anonymous function using B_w:
# Kp = 5
# L = 6
# MLT = 23
# agapitov_model(lambda) = B_w(lambda, Kp, α_ij_matrix(L, MLT)) * tanh(lambda) * 0.1216439



# plot(agapitov_model, 0, 90, yscale=:log10, ylim = (1e-3,10), xlim = (0,45), title = "L = $L, MLT = $MLT, Kp = $Kp")

# # function test_agapitov_model(L, Kp, MLT)
# L = 6.0;
# p = Vector{Plots.Plot{Plots.GRBackend}}()
# Kp = 5.5; MLT = 23;
# push!(p, plot(agapitov_model, 0, 90, yscale=:log10, ylim = (1e-3,10), xlim = (0,45), title = "L = $L, MLT = $MLT, Kp = $Kp", legend=false))
# Kp = 5.5; MLT = 5;
# push!(p, plot(agapitov_model, 0, 90, yscale=:log10, ylim = (1e-3,10), xlim = (0,45), title = "L = $L, MLT = $MLT, Kp = $Kp", legend=false))
# Kp = 5.5; MLT = 15;
# push!(p, plot(agapitov_model, 0, 90, yscale=:log10, ylim = (1e-3,10), xlim = (0,45), title = "L = $L, MLT = $MLT, Kp = $Kp", legend=false))
# Kp = 3; MLT = 23;
# push!(p, plot(agapitov_model, 0, 90, yscale=:log10, ylim = (1e-3,10), xlim = (0,45), title = "L = $L, MLT = $MLT, Kp = $Kp", legend=false))
# Kp = 3; MLT = 5;
# push!(p, plot(agapitov_model, 0, 90, yscale=:log10, ylim = (1e-3,10), xlim = (0,45), title = "L = $L, MLT = $MLT, Kp = $Kp", legend=false))
# Kp = 3; MLT = 15;
# push!(p, plot(agapitov_model, 0, 90, yscale=:log10, ylim = (1e-3,10), xlim = (0,45), title = "L = $L, MLT = $MLT, Kp = $Kp", legend=false))
# Kp = 1; MLT = 23;
# push!(p, plot(agapitov_model, 0, 90, yscale=:log10, ylim = (1e-3,10), xlim = (0,45), title = "L = $L, MLT = $MLT, Kp = $Kp", legend=false))
# Kp = 1; MLT = 5;
# push!(p, plot(agapitov_model, 0, 90, yscale=:log10, ylim = (1e-3,10), xlim = (0,45), title = "L = $L, MLT = $MLT, Kp = $Kp", legend=false))
# Kp = 1; MLT = 15;
# push!(p, plot(agapitov_model, 0, 90, yscale=:log10, ylim = (1e-3,10), xlim = (0,45), title = "L = $L, MLT = $MLT, Kp = $Kp", legend=false))
# plot(p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9], layout=(3,3), dpi = 96, size=(1500,1000))




# L=5.1
# Kp = 3; MLT = 21.7;
# plot(agapitov_model, 0, 90, yscale=:log10, ylim = (1e-1,10), xlim = (0,45), label = "L = $L, MLT = $MLT, Kp = $Kp")


# L=5.1
# Kp = 3; MLT = 23;
# plot!(agapitov_model, 0, 90, yscale=:log10, ylim = (1e-1,10), xlim = (0,45), label = "L = $L, MLT = $MLT, Kp = $Kp")

