include("plotHelpers.jl")
Egrid, PAgrid = logrange(10,1000,21), 6:4:90

# Plots for 4/29/21
@time @load "results/210429_lat_storage.jld2" equatorial_fluxes_042921 themis_hilat_042921 themis_midlat_042921 themis_lolat_042921
elfin_measurements_042921, elfin_error_042921 = extract_idl_csv("042921_time.csv", "042921_prec.csv",
                                                "042921_precerror.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2021,4,29,3,14,45),DateTime(2021,4,29,3,15,0)); # time to sample from ELFIN measurements                                                                

x1 = themis_lolat_042921.precipitating_fluxes_mean.-themis_lolat_042921.precipitating_fluxes_minus .+ 1.
x2 = themis_lolat_042921.precipitating_fluxes_mean.+themis_lolat_042921.precipitating_fluxes_plus .+ 1.
y1 = themis_midlat_042921.precipitating_fluxes_mean.-themis_midlat_042921.precipitating_fluxes_minus .+ 1.
y2 = themis_midlat_042921.precipitating_fluxes_mean.+themis_midlat_042921.precipitating_fluxes_plus .+ 1.
z1 = themis_hilat_042921.precipitating_fluxes_mean.-themis_hilat_042921.precipitating_fluxes_minus .+ 1.
z2 = themis_hilat_042921.precipitating_fluxes_mean.+themis_hilat_042921.precipitating_fluxes_plus .+ 1.
elfin_error_042921[8] += 2000. 
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
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\lambda_2 = 20\degree\mathrm{,\ occurence\ rate} = 40\%)",
        # yerror=(themis_lolat_042921.precipitating_fluxes_minus, themis_lolat_042921.precipitating_fluxes_plus),
        color = c4, alpha = .7, markerstrokecolor = c4, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(Egrid,themis_midlat_042921.precipitating_fluxes_mean, fillrange=y1, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid,themis_midlat_042921.precipitating_fluxes_mean, fillrange=y2, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid, themis_midlat_042921.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\lambda_2 = 30\degree\mathrm{,\ occurence\ rate} = 5\%)",
        # yerror=(themis_midlat_042921.precipitating_fluxes_minus, themis_midlat_042921.precipitating_fluxes_plus),
        color = c3, alpha = .7, markerstrokecolor = c3, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(Egrid,themis_hilat_042921.precipitating_fluxes_mean, fillrange=z1, fillalpha = 0.2, color = bipride_orange, label=false)
plot!(Egrid,themis_hilat_042921.precipitating_fluxes_mean, fillrange=z2, fillalpha = 0.2, color = bipride_orange, label=false)
plot!(Egrid, themis_hilat_042921.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\lambda_2 = 40\degree\mathrm{,\ occurence\ rate} = 3\%)",
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
plot1 = plot!(dpi = 300,size=(800,450), margin=3mm, bottom_margin=4mm)
savefig(plot1, "images/042921_flux_lat_comparison.pdf")
savefig(plot1, "images/042921_flux_lat_comparison.png")

@time @load "results/210429_pkt_storage.jld2" themis_hilat_042921 themis_nopkt_042921 themis_100pkt_042921;

np1 = themis_nopkt_042921.precipitating_fluxes_mean.-themis_nopkt_042921.precipitating_fluxes_minus .+ 1.
np2 = themis_nopkt_042921.precipitating_fluxes_mean.+themis_nopkt_042921.precipitating_fluxes_plus .+ 1.
pp1 = themis_100pkt_042921.precipitating_fluxes_mean.-themis_100pkt_042921.precipitating_fluxes_minus .+ 1.
pp2 = themis_100pkt_042921.precipitating_fluxes_mean.+themis_100pkt_042921.precipitating_fluxes_plus .+ 1.

plot(Egrid, equatorial_fluxes_042921,label = L"\mathrm{equatorial\ flux\ (THEMIS\ A\ measured)}", color = bipride_pink, linewidth=2, markershape=:circle);
plot!(Egrid,themis_nopkt_042921.precipitating_fluxes_mean, fillrange=np1, fillalpha = 0.2, color = c4, label=false)
plot!(Egrid,themis_nopkt_042921.precipitating_fluxes_mean, fillrange=np2, fillalpha = 0.2, color = c4, label=false)
plot!(Egrid, themis_nopkt_042921.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\phi = 1)",
        # yerror=(themis_nopkt_042921.precipitating_fluxes_minus, themis_nopkt_042921.precipitating_fluxes_plus),
        color = c4, alpha = .7, markerstrokecolor = c4, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(Egrid,themis_100pkt_042921.precipitating_fluxes_mean, fillrange=pp1, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid,themis_100pkt_042921.precipitating_fluxes_mean, fillrange=pp2, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid, themis_100pkt_042921.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\phi_2 = 100)",
        # yerror=(themis_100pkt_042921.precipitating_fluxes_minus, themis_100pkt_042921.precipitating_fluxes_plus),
        color = c3, alpha = .7, markerstrokecolor = c3, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(Egrid,themis_hilat_042921.precipitating_fluxes_mean, fillrange=z1, fillalpha = 0.2, color = bipride_orange, label=false)
plot!(Egrid,themis_hilat_042921.precipitating_fluxes_mean, fillrange=z2, fillalpha = 0.2, color = bipride_orange, label=false)
plot!(Egrid, themis_hilat_042921.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\phi = 30)",
        # yerror=(themis_hilat_042921.precipitating_fluxes_minus, themis_hilat_042921.precipitating_fluxes_plus),
        color = bipride_orange, alpha = .7, markerstrokecolor = bipride_orange, markerstrokewidth = 2, markershape=:o, markersize = 2);
# plot!(elfin_measurements_042921[1][1:end-5], elfin_measurements_042921[2][1:end-5], fillrange = e1[1:end-5], fillalpha = 0.2, color = c5, label=false)        
# plot!(elfin_measurements_042921[1][1:end-5], elfin_measurements_042921[2][1:end-5], fillrange = e2[1:end-5], fillalpha = 0.2, color = c5, label=false)        
# plot!(elfin_measurements_042921, yerror=elfin_error_042921.+10, label = L"\mathrm{precipitating\ flux\ (ELFIN\ B\ measured)}",
#         color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3);
plot!(ylim =(1e3,1e8), xlim=(50,800), yscale=:log10, legend=(0.45,0.7));
plot!(xlabel=L"\mathrm{Energy\ (keV)}", ylabel=L"\mathrm{Flux\ (1/cm^{2}/s/sr/MeV)}", title=L"\mathrm{Precipitating\ Flux\ Comparison\ of\ \delta\lambda_2\ on\ 4/29/21}",
xtickfontsize=12, ytickfontsize=12, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, titlefontsize=16);
plot2 = plot!(dpi = 300,size=(800,450), margin=3mm, bottom_margin=4mm)
savefig(plot2, "images/042921_flux_pkt_comparison.pdf")
savefig(plot2, "images/042921_flux_pkt_comparison.png")


equatorial_fluxes_092220
mms_lolat_092220
elfin_measurements_092220

include("plotHelpers.jl")
Egrid, PAgrid = logrange(10,1000,21), 6:4:90

# Plots for 9/22/20
@load "results/200922_lat_storage.jld2" equatorial_fluxes_092220 mms_lolat_092220 mms_midlat_092220 mms_hilat_092220;
@time elfin_measurements_3_092220 = extract_idl_csv("092220_time.csv", "092220_prec.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2020,9,22,9,16,38), DateTime(2020,9,22,9,16,45)) # time to sample from ELFIN measurements
@time elfin_measurements_2_092220 = extract_idl_csv("092220_time.csv", "092220_prec.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2020,9,22,9,16,25), DateTime(2020,9,22,9,16,34)) # time to sample from ELFIN measurements
@time elfin_measurements_1_092220 = extract_idl_csv("092220_time.csv", "092220_prec.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2020,9,22,9,15,48), DateTime(2020,9,22,9,16,3)) # time to sample from ELFIN measurements
@time elfin_measurements_all_092220 = extract_idl_csv("092220_time.csv", "092220_prec.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2020,9,22,9,15,48), DateTime(2020,9,22,9,16,45)) # time to sample from ELFIN measurements                                                


plot(Egrid, equatorial_fluxes_092220, label = "Equatorial Flux", color = bipride_blue, linewidth=2, markershape=:circle)
plot!(Egrid, mms_lolat_092220.precipitating_fluxes_mean,
        yerror=(mms_lolat_092220.precipitating_fluxes_minus, mms_lolat_092220.precipitating_fluxes_plus),
        color = bipride_orange, alpha = .5, label="dLambda2 = 20 deg", markershape=:x)
plot!(Egrid, mms_midlat_092220.precipitating_fluxes_mean,
        yerror=(mms_midlat_092220.precipitating_fluxes_minus, mms_midlat_092220.precipitating_fluxes_plus),
        color = bipride_orange, alpha = .5, label="dLambda2 = 30 deg", markershape=:x)
plot!(Egrid, mms_hilat_092220.precipitating_fluxes_mean,
        yerror=(mms_hilat_092220.precipitating_fluxes_minus, mms_hilat_092220.precipitating_fluxes_plus),
        color = bipride_orange, alpha = .5, label="dLambda2 = 40 deg", markershape=:x)
plot!(elfin_measurements_3_092220, label = "ELFIN Measured Precipitating Flux", color = bipride_pink, linewidth=3, markershape=:dtriangle);
plot!(ylim =(1e1,1e10), xlim=(50,800), yscale=:log10);
plot!(xlabel=L"\mathrm{Energy\ (keV)}", ylabel=L"\mathrm{Flux\ (1/cm^{2}/s/sr/MeV)}", title=L"\mathrm{Flux\ Comparison\ of\ Precipitating\ Particles\ on\ 09/22/20}",
xtickfontsize=12, ytickfontsize=12, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, titlefontsize=16);
plot2 = plot!(dpi = 300,size=(800,450), margin=3mm, bottom_margin=4mm)
# savefig(plot2, "images/092220_flux_comparison.pdf")
# savefig(plot2, "images/092220_flux_comparison.png")