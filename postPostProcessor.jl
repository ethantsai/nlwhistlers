include("plotHelpers.jl")
Egrid, PAgrid = logrange(10,1000,21), 6:4:90

# Plots for 4/29/21
@time @load "results/210429_lat_storage.jld2" equatorial_fluxes_042921 themis_hilat_042921 themis_midlat_042921 themis_lolat_042921
elfin_measurements_042921, elfin_error_042921 = extract_idl_csv("042921_time.csv", "042921_prec.csv",
                                                "042921_precerror.csv", "042921_precerror_time.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2021,4,29,3,14,45),DateTime(2021,4,29,3,15,0)); # time to sample from ELFIN measurements                                                                

x1 = themis_lolat_042921.precipitating_fluxes_mean.-themis_lolat_042921.precipitating_fluxes_minus .+ 1.
x2 = themis_lolat_042921.precipitating_fluxes_mean.+themis_lolat_042921.precipitating_fluxes_plus .+ 1.
y1 = themis_midlat_042921.precipitating_fluxes_mean.-themis_midlat_042921.precipitating_fluxes_minus .+ 1.
y2 = themis_midlat_042921.precipitating_fluxes_mean.+themis_midlat_042921.precipitating_fluxes_plus .+ 1.
z1 = themis_hilat_042921.precipitating_fluxes_mean.-themis_hilat_042921.precipitating_fluxes_minus .+ 1.
z2 = themis_hilat_042921.precipitating_fluxes_mean.+themis_hilat_042921.precipitating_fluxes_plus .+ 1.
# elfin_error_042921[8] += 2000. 
e1 = elfin_measurements_042921[2].-elfin_error_042921
e2 = elfin_measurements_042921[2].+elfin_error_042921
sum_of_means = @. themis_hilat_042921.precipitating_fluxes_mean + themis_midlat_042921.precipitating_fluxes_mean + themis_lolat_042921.precipitating_fluxes_mean
error_sum_minus = @. sqrt(themis_hilat_042921.precipitating_fluxes_minus^2 + themis_midlat_042921.precipitating_fluxes_minus^2 + themis_lolat_042921.precipitating_fluxes_minus^2)
error_sum_plus = @. sqrt(themis_hilat_042921.precipitating_fluxes_plus^2 + themis_midlat_042921.precipitating_fluxes_plus^2 + themis_lolat_042921.precipitating_fluxes_plus^2)
s1 = sum_of_means.-error_sum_minus .+ 1.
s2 = sum_of_means.+error_sum_plus .+ 1. 

plot(Egrid, equatorial_fluxes_042921,label = L"\mathrm{equatorial\ flux\ (THEMIS\ A\ measured)}", color = bipride_pink, linewidth=2, markershape=:circle);
plot!(Egrid,themis_lolat_042921.precipitating_fluxes_mean, fillrange=x1, fillalpha = 0.2, color = c4, label=false)
plot!(Egrid,themis_lolat_042921.precipitating_fluxes_mean, fillrange=x2, fillalpha = 0.2, color = c4, label=false)
plot!(Egrid, themis_lolat_042921.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\lambda_2 = 20\degree\mathrm{,\ occurence\ rate} = 30\%)",
        # yerror=(themis_lolat_042921.precipitating_fluxes_minus, themis_lolat_042921.precipitating_fluxes_plus),
        color = c4, alpha = .7, markerstrokecolor = c4, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(Egrid,themis_midlat_042921.precipitating_fluxes_mean, fillrange=y1, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid,themis_midlat_042921.precipitating_fluxes_mean, fillrange=y2, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid, themis_midlat_042921.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\lambda_2 = 30\degree\mathrm{,\ occurence\ rate} = 4\%)",
        # yerror=(themis_midlat_042921.precipitating_fluxes_minus, themis_midlat_042921.precipitating_fluxes_plus),
        color = c3, alpha = .7, markerstrokecolor = c3, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(Egrid,themis_hilat_042921.precipitating_fluxes_mean, fillrange=z1, fillalpha = 0.2, color = bipride_orange, label=false)
plot!(Egrid,themis_hilat_042921.precipitating_fluxes_mean, fillrange=z2, fillalpha = 0.2, color = bipride_orange, label=false)
plot!(Egrid, themis_hilat_042921.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\lambda_2 = 40\degree\mathrm{,\ occurence\ rate} = 2\%)",
        # yerror=(themis_hilat_042921.precipitating_fluxes_minus, themis_hilat_042921.precipitating_fluxes_plus),
        color = bipride_orange, alpha = .7, markerstrokecolor = bipride_orange, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(Egrid,sum_of_means, fillrange=s1, fillalpha = 0.2, color = :black, label=false)
plot!(Egrid,sum_of_means, fillrange=s2, fillalpha = 0.2, color = :black, label=false)
plot!(Egrid, sum_of_means, label=L"\mathrm{precipitating\ flux\ (sum\ of\ simulations)}",
        # yerror=(error_sum_minus, error_sum_plus),
        color = :black, linewidth=3, markershape=:o, markersize = 3);
plot!(elfin_measurements_042921[1][1:end-5], elfin_measurements_042921[2][1:end-5], fillrange = e1[1:end-5], fillalpha = 0.2, color = c5, label=false)        
plot!(elfin_measurements_042921[1][1:end-5], elfin_measurements_042921[2][1:end-5], fillrange = e2[1:end-5], fillalpha = 0.2, color = c5, label=false)        
plot!(elfin_measurements_042921, yerror=elfin_error_042921.+10, label = L"\mathrm{precipitating\ flux\ (ELFIN\ B\ measured)}",
        color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3);
plot!(ylim =(1e3,1e8), xlim=(50,800), yscale=:log10, legend=(0.45,0.7))
plot!(xlabel=L"\mathrm{Energy\ (keV)}", ylabel=L"\mathrm{Flux\ (1/cm^{2}/s/sr/MeV)}", title=L"\mathrm{Precipitating\ Flux\ Comparison\ of\ \delta\lambda_2\ on\ 4/29/21}",
xtickfontsize=12, ytickfontsize=12, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, titlefontsize=16);
plot1 = plot!(dpi = 500,size=(800,450), margin=3mm, bottom_margin=4mm)
savefig(plot1, "images/042921_flux_lat_comparison.pdf")
savefig(plot1, "images/042921_flux_lat_comparison.png")

@time @load "results/210429_pkt_storage.jld2" equatorial_fluxes_042921 themis_1pkt_042921 themis_30pkt_042921 themis_100pkt_042921;

np1 = themis_1pkt_042921.precipitating_fluxes_mean.-themis_1pkt_042921.precipitating_fluxes_minus .+ 1.
np2 = themis_1pkt_042921.precipitating_fluxes_mean.+themis_1pkt_042921.precipitating_fluxes_plus .+ 1.
op1 = themis_30pkt_042921.precipitating_fluxes_mean.-themis_30pkt_042921.precipitating_fluxes_minus .+ 1.
op2 = themis_30pkt_042921.precipitating_fluxes_mean.+themis_30pkt_042921.precipitating_fluxes_plus .+ 1.
pp1 = themis_100pkt_042921.precipitating_fluxes_mean.-themis_100pkt_042921.precipitating_fluxes_minus .+ 1.
pp2 = themis_100pkt_042921.precipitating_fluxes_mean.+themis_100pkt_042921.precipitating_fluxes_plus .+ 1.

plot(Egrid, equatorial_fluxes_042921,label = L"\mathrm{equatorial\ flux\ (THEMIS\ A\ measured)}", color = bipride_pink, linewidth=2, markershape=:circle);
plot!(Egrid,themis_1pkt_042921.precipitating_fluxes_mean, fillrange=np1, fillalpha = 0.2, color = c4, label=false)
plot!(Egrid,themis_1pkt_042921.precipitating_fluxes_mean, fillrange=np2, fillalpha = 0.2, color = c4, label=false)
plot!(Egrid, themis_1pkt_042921.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\phi = 1)",
        # yerror=(themis_nopkt_042921.precipitating_fluxes_minus, themis_nopkt_042921.precipitating_fluxes_plus),
        color = c4, alpha = .7, markerstrokecolor = c4, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(Egrid,themis_30pkt_042921.precipitating_fluxes_mean, fillrange=op1, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid,themis_30pkt_042921.precipitating_fluxes_mean, fillrange=op2, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid, themis_30pkt_042921.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\phi = 30)",
        # yerror=(themis_100pkt_042921.precipitating_fluxes_minus, themis_100pkt_042921.precipitating_fluxes_plus),
        color = c3, alpha = .7, markerstrokecolor = c3, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(Egrid,themis_100pkt_042921.precipitating_fluxes_mean, fillrange=pp1, fillalpha = 0.2, color = bipride_orange, label=false)
plot!(Egrid,themis_100pkt_042921.precipitating_fluxes_mean, fillrange=pp2, fillalpha = 0.2, color = bipride_orange, label=false)
plot!(Egrid, themis_100pkt_042921.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\phi = 100)",
        # yerror=(themis_100pkt_042921.precipitating_fluxes_minus, themis_100pkt_042921.precipitating_fluxes_plus),
        color = bipride_orange, alpha = .7, markerstrokecolor = bipride_orange, markerstrokewidth = 2, markershape=:o, markersize = 2);
# plot!(elfin_measurements_042921[1][1:end-5], elfin_measurements_042921[2][1:end-5], fillrange = e1[1:end-5], fillalpha = 0.2, color = c5, label=false)        
# plot!(elfin_measurements_042921[1][1:end-5], elfin_measurements_042921[2][1:end-5], fillrange = e2[1:end-5], fillalpha = 0.2, color = c5, label=false)        
# plot!(elfin_measurements_042921, yerror=elfin_error_042921.+10, label = L"\mathrm{precipitating\ flux\ (ELFIN\ B\ measured)}",
#         color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3);
plot!(ylim =(1e3,1e8), xlim=(50,800), yscale=:log10, legend=(0.45,0.7));
plot!(xlabel=L"\mathrm{Energy\ (keV)}", ylabel=L"\mathrm{Flux\ (1/cm^{2}/s/sr/MeV)}", title=L"\mathrm{Precipitating\ Flux\ Comparison\ of\ \delta\phi\ on\ 4/29/21}",
xtickfontsize=12, ytickfontsize=12, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, titlefontsize=16);
plot2 = plot!(dpi = 500,size=(800,450), margin=3mm, bottom_margin=4mm)
# savefig(plot2, "images/042921_flux_pkt_comparison.pdf")
# savefig(plot2, "images/042921_flux_pkt_comparison.png")



include("plotHelpers.jl")
Egrid, PAgrid = logrange(10,1000,21), 6:4:90

# Plots for 9/22/20
@load "results/200922_lat_storage.jld2" equatorial_fluxes_092220 mms_lowerlat_092220 mms_lolat_092220 mms_midlat_092220 mms_hilat_092220;
@time elfin_measurements_092220, elfin_error_092220 = extract_idl_csv("092220_time.csv", "092220_prec.csv",
                                                "092220_precerror.csv", "092220_precerror_time.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2020,9,22,9,16,38), DateTime(2020,9,22,9,16,45)) # time to sample from ELFIN measurements
# @time elfin_measurements_092220 = extract_idl_csv("092220_time.csv", "092220_prec.csv", "ebins.csv", # csvs containing ELFIN measurements
#                                                 DateTime(2020,9,22,9,16,25), DateTime(2020,9,22,9,16,34)) # time to sample from ELFIN measurements
# @time elfin_measurements_092220 = extract_idl_csv("092220_time.csv", "092220_prec.csv", "ebins.csv", # csvs containing ELFIN measurements
#                                                 DateTime(2020,9,22,9,15,48), DateTime(2020,9,22,9,16,3)) # time to sample from ELFIN measurements
# @time elfin_measurements_all_092220 = extract_idl_csv("092220_time.csv", "092220_prec.csv", "ebins.csv", # csvs containing ELFIN measurements
#                                                 DateTime(2020,9,22,9,15,48), DateTime(2020,9,22,9,16,45)) # time to sample from ELFIN measurements                                                


x1 = mms_lowerlat_092220.precipitating_fluxes_mean.-mms_lowerlat_092220.precipitating_fluxes_minus .+ 1.
x2 = mms_lowerlat_092220.precipitating_fluxes_mean.+mms_lowerlat_092220.precipitating_fluxes_plus .+ 1.
y1 = mms_lolat_092220.precipitating_fluxes_mean.-mms_lolat_092220.precipitating_fluxes_minus .+ 1.
y2 = mms_lolat_092220.precipitating_fluxes_mean.+mms_lolat_092220.precipitating_fluxes_plus .+ 1.
z1 = mms_hilat_092220.precipitating_fluxes_mean.-mms_hilat_092220.precipitating_fluxes_minus .+ 1.
z2 = mms_hilat_092220.precipitating_fluxes_mean.+mms_hilat_092220.precipitating_fluxes_plus .+ 1.
e1 = elfin_measurements_092220[2].-elfin_error_092220
e2 = elfin_measurements_092220[2].+elfin_error_092220

plot(Egrid, equatorial_fluxes_092220,label = L"\mathrm{equatorial\ flux\ (MMS\ 1\ measured)}", color = bipride_pink, linewidth=2, markershape=:circle);
plot!(Egrid,mms_lowerlat_092220.precipitating_fluxes_mean, fillrange=x1, fillalpha = 0.2, color = c4, label=false)
plot!(Egrid,mms_lowerlat_092220.precipitating_fluxes_mean, fillrange=x2, fillalpha = 0.2, color = c4, label=false)
plot!(Egrid, mms_lowerlat_092220.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\lambda_2 = 15\degree\mathrm{,\ occurence\ rate} = 6\%)",
        # yerror=(mms_lowerlat_092220.precipitating_fluxes_minus, mms_lowerlat_092220.precipitating_fluxes_plus),
        color = c4, alpha = .7, markerstrokecolor = c4, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(Egrid,mms_lolat_092220.precipitating_fluxes_mean, fillrange=y1, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid,mms_lolat_092220.precipitating_fluxes_mean, fillrange=y2, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid, mms_lolat_092220.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\lambda_2 = 20\degree\mathrm{,\ occurence\ rate} = 6\%)",
        # yerror=(mms_lolat_092220.precipitating_fluxes_minus, mms_lolat_092220.precipitating_fluxes_plus),
        color = c3, alpha = .7, markerstrokecolor = c3, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(Egrid,mms_hilat_092220.precipitating_fluxes_mean, fillrange=z1, fillalpha = 0.2, color = bipride_orange, label=false)
plot!(Egrid,mms_hilat_092220.precipitating_fluxes_mean, fillrange=z2, fillalpha = 0.2, color = bipride_orange, label=false)
plot!(Egrid, mms_hilat_092220.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\lambda_2 = 30\degree\mathrm{,\ occurence\ rate} = 2\%)",
        # yerror=(mms_hilat_092220.precipitating_fluxes_minus, mms_hilat_092220.precipitating_fluxes_plus),
        color = bipride_orange, alpha = .7, markerstrokecolor = bipride_orange, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(elfin_measurements_092220[1][1:6], elfin_measurements_092220[2][1:6], fillrange = e1[1:end-5], fillalpha = 0.2, color = c5, label=false)        
plot!(elfin_measurements_092220[1][1:6], elfin_measurements_092220[2][1:6], fillrange = e2[1:end-5], fillalpha = 0.2, color = c5, label=false)        
plot!(elfin_measurements_092220[1][1:6], elfin_measurements_092220[2][1:6],
        yerror=elfin_error_092220[1:6].+10,
        label = L"\mathrm{precipitating\ flux\ (ELFIN\ A\ measured)}",
        color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3);
plot!(ylim =(1e1,1e10), xlim=(50,800), yscale=:log10, legend=:topright)
plot!(xlabel=L"\mathrm{Energy\ (keV)}", ylabel=L"\mathrm{Flux\ (1/cm^{2}/s/sr/MeV)}", title=L"\mathrm{Precipitating\ Flux\ Comparison\ of\ \delta\lambda_2\ on\ 9/22/20}",
xtickfontsize=12, ytickfontsize=12, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, titlefontsize=16);
plot3 = plot!(dpi = 500,size=(800,450), margin=3mm, bottom_margin=4mm)
savefig(plot3, "images/092220_flux_comparison.pdf")
savefig(plot3, "images/092220_flux_comparison.png")


# plot PA and Energy trajectories for different dLambda
@time @load "results/demo_themis_dlambda_20.jld2" demo_E_6_20
@time @load "results/demo_themis_dlambda_30.jld2" demo_E_6_30
@time @load "results/demo_themis_dlambda_40.jld2" demo_E_6_40
L=5.2
twenty = trajectoryChecking(demo_E_6_20, [4,19,23,29,34,48], L"\delta \lambda_2 = "*"20 deg")
thirty = trajectoryChecking(demo_E_6_30, [2,18,21,29,34,50], L"\delta \lambda_2 = "*"30 deg")
forty = trajectoryChecking(demo_E_6_40, [1,19,21,27,34,49], L"\delta \lambda_2 = "*"40 deg")
bigplot = plot(twenty,thirty,forty, dpi = 96, layout = (1,3), left_margin = 8mm, size=(1500,500))


savefig(bigplot, "images/trajectories.pdf")
savefig(bigplot, "images/trajectories.png")



# plot PA and Energy for different dPhi
@time @load "results/demo_mms_dPhi_3.jld2" demo_mms_dPhi_3
@time @load "results/demo_mms_dPhi_30.jld2" demo_mms_dPhi_30
@time @load "results/demo_mms_dPhi_300.jld2" demo_mms_dPhi_300
L=6
three = trajectoryChecking(demo_mms_dPhi_3, [1,4,5,6,7,9], L"\delta \phi = "*"3")
thirty = trajectoryChecking(demo_mms_dPhi_30, [1,4,5,6,7,9], L"\delta \phi = "*"30")
threehundred = trajectoryChecking(demo_mms_dPhi_300, [1,4,5,6,7,9], L"\delta \phi = "*"300")
bigphiplot = plot(three,thirty,threehundred, dpi = 96, layout = (1,3), left_margin = 8mm, size=(1500,500))
savefig(bigphiplot, "images/dphitrajectories.pdf")
savefig(bigphiplot, "images/dphitrajectories.png")