include("plotHelpers.jl")
Egrid, PAgrid = logrange(10,1000,21), 6:4:90


# Plots for 4/29/21
@time @load "results/210429_lat_storage_cleaned.jld2" equatorial_fluxes_042921 themis_hilat_combined_042921 themis_midlat_042921 themis_lolat_042921
elfin_measurements_042921, elfin_error_042921 = extract_idl_csv("042921_time.csv", "042921_prec.csv",
                                                "042921_precerror.csv", "042921_precerror_time.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2021,4,29,3,14,45),DateTime(2021,4,29,3,15,0)); # time to sample from ELFIN measurements                                                                

x1 = themis_lolat_042921.precipitating_fluxes_mean.-themis_lolat_042921.precipitating_fluxes_minus .+ 1.
x2 = themis_lolat_042921.precipitating_fluxes_mean.+themis_lolat_042921.precipitating_fluxes_plus .+ 1.
y1 = themis_midlat_042921.precipitating_fluxes_mean.-themis_midlat_042921.precipitating_fluxes_minus .+ 1.
y2 = themis_midlat_042921.precipitating_fluxes_mean.+themis_midlat_042921.precipitating_fluxes_plus .+ 1.
z1 = themis_hilat_combined_042921.precipitating_fluxes_mean.-themis_hilat_combined_042921.precipitating_fluxes_minus .+ 1.
z2 = themis_hilat_combined_042921.precipitating_fluxes_mean.+themis_hilat_combined_042921.precipitating_fluxes_plus .+ 1.
# elfin_error_042921[8] += 2000. 
e1 = elfin_measurements_042921[2].-elfin_error_042921
e2 = elfin_measurements_042921[2].+elfin_error_042921
sum_of_means = @. themis_hilat_combined_042921.precipitating_fluxes_mean + themis_midlat_042921.precipitating_fluxes_mean + themis_lolat_042921.precipitating_fluxes_mean
error_sum_minus = @. sqrt(themis_hilat_combined_042921.precipitating_fluxes_minus^2 + themis_midlat_042921.precipitating_fluxes_minus^2 + themis_lolat_042921.precipitating_fluxes_minus^2)
error_sum_plus = @. sqrt(themis_hilat_combined_042921.precipitating_fluxes_plus^2 + themis_midlat_042921.precipitating_fluxes_plus^2 + themis_lolat_042921.precipitating_fluxes_plus^2)
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
plot!(Egrid,themis_hilat_combined_042921.precipitating_fluxes_mean, fillrange=z1, fillalpha = 0.2, color = bipride_orange, label=false)
plot!(Egrid,themis_hilat_combined_042921.precipitating_fluxes_mean, fillrange=z2, fillalpha = 0.2, color = bipride_orange, label=false)
plot!(Egrid, themis_hilat_combined_042921.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\lambda_2 = 40\degree\mathrm{,\ occurence\ rate} = 2\%)",
        # yerror=(themis_hilat_combined_042921.precipitating_fluxes_minus, themis_hilat_combined_042921.precipitating_fluxes_plus),
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





include("plotHelpers.jl")
Egrid, PAgrid = logrange(10,1000,21), 6:4:90

# Plots for 9/22/20
@time @load "results/200922_pkt_a3_storage.jld2" mms_med_a3_092220;
@time @load "results/200922_pkt_storage.jld2" mms_med_a5_092220;
@time @load "results/200922_pkt_a7_storage.jld2" mms_med_a7_092220;
@load "results/200922_lat_storage.jld2" equatorial_fluxes_092220;
@time elfin_measurements_092220, elfin_error_092220 = extract_idl_csv("092220_time.csv", "092220_prec.csv",
                                                "092220_precerror.csv", "092220_precerror_time.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2020,9,22,9,16,38), DateTime(2020,9,22,9,16,45)) # time to sample from ELFIN measurements
# @time elfin_measurements_092220 = extract_idl_csv("092220_time.csv", "092220_prec.csv", "ebins.csv", # csvs containing ELFIN measurements
#                                                 DateTime(2020,9,22,9,16,25), DateTime(2020,9,22,9,16,34)) # time to sample from ELFIN measurements
# @time elfin_measurements_092220 = extract_idl_csv("092220_time.csv", "092220_prec.csv", "ebins.csv", # csvs containing ELFIN measurements
#                                                 DateTime(2020,9,22,9,15,48), DateTime(2020,9,22,9,16,3)) # time to sample from ELFIN measurements
# @time elfin_measurements_all_092220 = extract_idl_csv("092220_time.csv", "092220_prec.csv", "ebins.csv", # csvs containing ELFIN measurements
#                                                 DateTime(2020,9,22,9,15,48), DateTime(2020,9,22,9,16,45)) # time to sample from ELFIN measurements                                                

x1 = mms_med_a3_092220.precipitating_fluxes_mean.-mms_med_a3_092220.precipitating_fluxes_minus .+ 1.
x2 = mms_med_a3_092220.precipitating_fluxes_mean.+mms_med_a3_092220.precipitating_fluxes_plus .+ 1.
y1 = mms_med_a5_092220.precipitating_fluxes_mean.-mms_med_a5_092220.precipitating_fluxes_minus .+ 1.
y2 = mms_med_a5_092220.precipitating_fluxes_mean.+mms_med_a5_092220.precipitating_fluxes_plus .+ 1.
z1 = mms_med_a7_092220.precipitating_fluxes_mean.-mms_med_a7_092220.precipitating_fluxes_minus .+ 1.
z2 = mms_med_a7_092220.precipitating_fluxes_mean.+mms_med_a7_092220.precipitating_fluxes_plus .+ 1.
e1 = elfin_measurements_092220[2].-elfin_error_092220
e2 = elfin_measurements_092220[2].+elfin_error_092220

plot(Egrid, equatorial_fluxes_092220,label = L"\mathrm{equatorial\ flux\ (MMS\ 1\ measured)}", color = bipride_pink, linewidth=2, markershape=:circle);
plot!(Egrid,mms_med_a3_092220.precipitating_fluxes_mean, fillrange=x1, fillalpha = 0.2, color = c4, label=false)
plot!(Egrid,mms_med_a3_092220.precipitating_fluxes_mean, fillrange=x2, fillalpha = 0.2, color = c4, label=false)
plot!(Egrid, mms_med_a3_092220.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }a=3\mathrm{,\delta\lambda_2} = 20\degree)",
        # yerror=(mms_med_a3_092220.precipitating_fluxes_minus, mms_med_a3_092220.precipitating_fluxes_plus),
        color = c4, alpha = .7, markerstrokecolor = c4, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(Egrid,mms_med_a5_092220.precipitating_fluxes_mean, fillrange=y1, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid,mms_med_a5_092220.precipitating_fluxes_mean, fillrange=y2, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid, mms_med_a5_092220.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }a=5\mathrm{,\delta\lambda_2} = 20\degree)",
        # yerror=(mms_med_a5_092220.precipitating_fluxes_minus, mms_med_a5_092220.precipitating_fluxes_plus),
        color = c3, alpha = .7, markerstrokecolor = c3, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(Egrid,mms_med_a7_092220.precipitating_fluxes_mean, fillrange=z1, fillalpha = 0.2, color = bipride_orange, label=false)
plot!(Egrid,mms_med_a7_092220.precipitating_fluxes_mean, fillrange=z2, fillalpha = 0.2, color = bipride_orange, label=false)
plot!(Egrid, mms_med_a7_092220.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }a=7\mathrm{,\delta\lambda_2} = 20\degree)",
        # yerror=(mms_med_a7_092220.precipitating_fluxes_minus, mms_med_a7_092220.precipitating_fluxes_plus),
        color = bipride_orange, alpha = .7, markerstrokecolor = bipride_orange, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(elfin_measurements_092220[1][1:6], elfin_measurements_092220[2][1:6], fillrange = e1[1:end-5], fillalpha = 0.2, color = c5, label=false)        
plot!(elfin_measurements_092220[1][1:6], elfin_measurements_092220[2][1:6], fillrange = e2[1:end-5], fillalpha = 0.2, color = c5, label=false)        
plot!(elfin_measurements_092220[1][1:6], elfin_measurements_092220[2][1:6],
        yerror=elfin_error_092220[1:6].+10,
        label = L"\mathrm{precipitating\ flux\ (ELFIN\ A\ measured)}",
        color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3);
plot!(ylim =(1e1,1e10), xlim=(50,800), yscale=:log10, legend=:topright)
plot!(xlabel=L"\mathrm{Energy\ (keV)}", ylabel=L"\mathrm{Flux\ (1/cm^{2}/s/sr/MeV)}", title=L"\mathrm{Precipitating\ Flux\ Comparison\ of\ Wave\ Modulation\ Depth\ on\ 9/22/20}",
xtickfontsize=12, ytickfontsize=12, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, titlefontsize=16);
plot2 = plot!(dpi = 500,size=(800,450), margin=3mm, bottom_margin=4mm, right_margin=8mm)
savefig(plot2, "images/092220_a_flux_comparison.pdf")
savefig(plot2, "images/092220_a_flux_comparison.png")


include("plotHelpers.jl")
Egrid, PAgrid = logrange(10,1000,21), 6:4:90

# Plots for 9/22/20
@time @load "results/200922_pkt_a3_storage.jld2" mms_med_a3_092220;
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
y1 = mms_med_a3_092220.precipitating_fluxes_mean.-mms_med_a3_092220.precipitating_fluxes_minus .+ 1.
y2 = mms_med_a3_092220.precipitating_fluxes_mean.+mms_med_a3_092220.precipitating_fluxes_plus .+ 1.
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
plot!(Egrid,mms_med_a3_092220.precipitating_fluxes_mean, fillrange=y1, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid,mms_med_a3_092220.precipitating_fluxes_mean, fillrange=y2, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid, mms_med_a3_092220.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\lambda_2 = 20\degree\mathrm{,\ occurence\ rate} = 6\%)",
        # yerror=(mms_med_a3_092220.precipitating_fluxes_minus, mms_med_a3_092220.precipitating_fluxes_plus),
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
savefig(plot3, "images/092220_dl2_flux_comparison.pdf")
savefig(plot3, "images/092220_dl2_flux_comparison.png")







include("plotHelpers.jl")
Egrid, PAgrid = logrange(10,1000,21), 6:4:90

# Plots packet comparison on 9/22/20 w/ a=5
@time @load "results/200922_pkt_storage.jld2" equatorial_fluxes_092220 mms_short_a5_092220 mms_med_a5_092220 mms_long_a5_092220;
@time elfin_measurements_092220, elfin_error_092220 = extract_idl_csv("092220_time.csv", "092220_prec.csv",
                                                "092220_precerror.csv", "092220_precerror_time.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2020,9,22,9,16,38), DateTime(2020,9,22,9,16,45)) # time to sample from ELFIN measurements

x1 = mms_short_a5_092220.precipitating_fluxes_mean.-mms_short_a5_092220.precipitating_fluxes_minus .+ 1.
x2 = mms_short_a5_092220.precipitating_fluxes_mean.+mms_short_a5_092220.precipitating_fluxes_plus .+ 1.
y1 = mms_med_a5_092220.precipitating_fluxes_mean.-mms_med_a5_092220.precipitating_fluxes_minus .+ 1.
y2 = mms_med_a5_092220.precipitating_fluxes_mean.+mms_med_a5_092220.precipitating_fluxes_plus .+ 1.
z1 = mms_long_a5_092220.precipitating_fluxes_mean.-mms_long_a5_092220.precipitating_fluxes_minus .+ 1.
z2 = mms_long_a5_092220.precipitating_fluxes_mean.+mms_long_a5_092220.precipitating_fluxes_plus .+ 1.
e1 = elfin_measurements_092220[2].-elfin_error_092220
e2 = elfin_measurements_092220[2].+elfin_error_092220

plot(Egrid, equatorial_fluxes_092220,label = L"\mathrm{equatorial\ flux\ (MMS\ 1\ measured)}", color = bipride_pink, linewidth=2, markershape=:circle);
plot!(Egrid,mms_short_a5_092220.precipitating_fluxes_mean, fillrange=x1, fillalpha = 0.2, color = c4, label=false)
plot!(Egrid,mms_short_a5_092220.precipitating_fluxes_mean, fillrange=x2, fillalpha = 0.2, color = c4, label=false)
plot!(Egrid, mms_short_a5_092220.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\phi = 3, a=5)",
        # yerror=(mms_short_a5_092220.precipitating_fluxes_minus, mms_short_a5_092220.precipitating_fluxes_plus),
        color = c4, alpha = .7, markerstrokecolor = c4, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(Egrid,mms_med_a5_092220.precipitating_fluxes_mean, fillrange=y1, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid,mms_med_a5_092220.precipitating_fluxes_mean, fillrange=y2, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid, mms_med_a5_092220.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\phi = 30, a=5)",
        # yerror=(mms_med_a5_092220.precipitating_fluxes_minus, mms_med_a5_092220.precipitating_fluxes_plus),
        color = c3, alpha = .7, markerstrokecolor = c3, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(Egrid,mms_long_a5_092220.precipitating_fluxes_mean, fillrange=z1, fillalpha = 0.2, color = bipride_orange, label=false)
plot!(Egrid,mms_long_a5_092220.precipitating_fluxes_mean, fillrange=z2, fillalpha = 0.2, color = bipride_orange, label=false)
plot!(Egrid, mms_long_a5_092220.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\phi = 300, a=5)",
        # yerror=(mms_long_a5_092220.precipitating_fluxes_minus, mms_long_a5_092220.precipitating_fluxes_plus),
        color = bipride_orange, alpha = .7, markerstrokecolor = bipride_orange, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(elfin_measurements_092220[1][1:6], elfin_measurements_092220[2][1:6], fillrange = e1[1:end-5], fillalpha = 0.2, color = c5, label=false)        
plot!(elfin_measurements_092220[1][1:6], elfin_measurements_092220[2][1:6], fillrange = e2[1:end-5], fillalpha = 0.2, color = c5, label=false)        
plot!(elfin_measurements_092220[1][1:6], elfin_measurements_092220[2][1:6],
        yerror=elfin_error_092220[1:6].+10,
        label = L"\mathrm{precipitating\ flux\ (ELFIN\ A\ measured)}",
        color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3);
plot!(ylim =(1e1,1e10), xlim=(50,800), yscale=:log10, legend=:topright)
plot!(xlabel=L"\mathrm{Energy\ (keV)}", ylabel=L"\mathrm{Flux\ (1/cm^{2}/s/sr/MeV)}", title=L"\mathrm{Precipitating\ Flux\ Comparison\ of\ \delta\phi\ on\ 9/22/20}",
xtickfontsize=12, ytickfontsize=12, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, titlefontsize=16);
plot4 = plot!(dpi = 500,size=(800,450), margin=3mm, bottom_margin=4mm)
savefig(plot4, "images/092220_flux_pkt_a5_comparison.pdf")
savefig(plot4, "images/092220_flux_pkt_a5_comparison.png")



include("plotHelpers.jl")
Egrid, PAgrid = logrange(10,1000,21), 6:4:90

# Plots packet comparison on 9/22/20 w/ a=3
@time @load "results/200922_pkt_a3_storage.jld2" equatorial_fluxes_092220 mms_short_a3_092220 mms_med_a3_092220 mms_long_a3_092220;
@time elfin_measurements_092220, elfin_error_092220 = extract_idl_csv("092220_time.csv", "092220_prec.csv",
                                                "092220_precerror.csv", "092220_precerror_time.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2020,9,22,9,16,38), DateTime(2020,9,22,9,16,45)) # time to sample from ELFIN measurements

x1 = mms_short_a3_092220.precipitating_fluxes_mean.-mms_short_a3_092220.precipitating_fluxes_minus .+ 1.
x2 = mms_short_a3_092220.precipitating_fluxes_mean.+mms_short_a3_092220.precipitating_fluxes_plus .+ 1.
y1 = mms_med_a3_092220.precipitating_fluxes_mean.-mms_med_a3_092220.precipitating_fluxes_minus .+ 1.
y2 = mms_med_a3_092220.precipitating_fluxes_mean.+mms_med_a3_092220.precipitating_fluxes_plus .+ 1.
z1 = mms_long_a3_092220.precipitating_fluxes_mean.-mms_long_a3_092220.precipitating_fluxes_minus .+ 1.
z2 = mms_long_a3_092220.precipitating_fluxes_mean.+mms_long_a3_092220.precipitating_fluxes_plus .+ 1.
e1 = elfin_measurements_092220[2].-elfin_error_092220
e2 = elfin_measurements_092220[2].+elfin_error_092220

plot(Egrid, equatorial_fluxes_092220,label = L"\mathrm{equatorial\ flux\ (MMS\ 1\ measured)}", color = bipride_pink, linewidth=2, markershape=:circle);
plot!(Egrid,mms_short_a3_092220.precipitating_fluxes_mean, fillrange=x1, fillalpha = 0.2, color = c4, label=false)
plot!(Egrid,mms_short_a3_092220.precipitating_fluxes_mean, fillrange=x2, fillalpha = 0.2, color = c4, label=false)
plot!(Egrid, mms_short_a3_092220.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\phi = 3, a=3)",
        # yerror=(mms_short_a3_092220.precipitating_fluxes_minus, mms_short_a3_092220.precipitating_fluxes_plus),
        color = c4, alpha = .7, markerstrokecolor = c4, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(Egrid,mms_med_a3_092220.precipitating_fluxes_mean, fillrange=y1, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid,mms_med_a3_092220.precipitating_fluxes_mean, fillrange=y2, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid, mms_med_a3_092220.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\phi = 30, a=3)",
        # yerror=(mms_med_a3_092220.precipitating_fluxes_minus, mms_med_a3_092220.precipitating_fluxes_plus),
        color = c3, alpha = .7, markerstrokecolor = c3, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(Egrid,mms_long_a3_092220.precipitating_fluxes_mean, fillrange=z1, fillalpha = 0.2, color = bipride_orange, label=false)
plot!(Egrid,mms_long_a3_092220.precipitating_fluxes_mean, fillrange=z2, fillalpha = 0.2, color = bipride_orange, label=false)
plot!(Egrid, mms_long_a3_092220.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\phi = 300, a=3)",
        # yerror=(mms_long_a3_092220.precipitating_fluxes_minus, mms_long_a3_092220.precipitating_fluxes_plus),
        color = bipride_orange, alpha = .7, markerstrokecolor = bipride_orange, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(elfin_measurements_092220[1][1:6], elfin_measurements_092220[2][1:6], fillrange = e1[1:end-5], fillalpha = 0.2, color = c5, label=false)        
plot!(elfin_measurements_092220[1][1:6], elfin_measurements_092220[2][1:6], fillrange = e2[1:end-5], fillalpha = 0.2, color = c5, label=false)        
plot!(elfin_measurements_092220[1][1:6], elfin_measurements_092220[2][1:6],
        yerror=elfin_error_092220[1:6].+10,
        label = L"\mathrm{precipitating\ flux\ (ELFIN\ A\ measured)}",
        color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3);
plot!(ylim =(1e1,1e10), xlim=(50,800), yscale=:log10, legend=:topright)
plot!(xlabel=L"\mathrm{Energy\ (keV)}", ylabel=L"\mathrm{Flux\ (1/cm^{2}/s/sr/MeV)}", title=L"\mathrm{Precipitating\ Flux\ Comparison\ of\ \delta\phi\ on\ 9/22/20}",
xtickfontsize=12, ytickfontsize=12, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, titlefontsize=16);
plot5 = plot!(dpi = 500,size=(800,450), margin=3mm, bottom_margin=4mm)
savefig(plot5, "images/092220_flux_pkt_a3_comparison.pdf")
savefig(plot5, "images/092220_flux_pkt_a3_comparison.png")

include("plotHelpers.jl")
Egrid, PAgrid = logrange(10,1000,21), 6:4:90

# Plots packet comparison on 9/22/20 w/ a=7
@time @load "results/200922_pkt_a7_storage.jld2" equatorial_fluxes_092220 mms_short_a7_092220 mms_med_a7_092220 mms_long_a7_092220;
@time elfin_measurements_092220, elfin_error_092220 = extract_idl_csv("092220_time.csv", "092220_prec.csv",
                                                "092220_precerror.csv", "092220_precerror_time.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2020,9,22,9,16,38), DateTime(2020,9,22,9,16,45)) # time to sample from ELFIN measurements

x1 = mms_short_a7_092220.precipitating_fluxes_mean.-mms_short_a7_092220.precipitating_fluxes_minus .+ 1.
x2 = mms_short_a7_092220.precipitating_fluxes_mean.+mms_short_a7_092220.precipitating_fluxes_plus .+ 1.
y1 = mms_med_a7_092220.precipitating_fluxes_mean.-mms_med_a7_092220.precipitating_fluxes_minus .+ 1.
y2 = mms_med_a7_092220.precipitating_fluxes_mean.+mms_med_a7_092220.precipitating_fluxes_plus .+ 1.
z1 = mms_long_a7_092220.precipitating_fluxes_mean.-mms_long_a7_092220.precipitating_fluxes_minus .+ 1.
z2 = mms_long_a7_092220.precipitating_fluxes_mean.+mms_long_a7_092220.precipitating_fluxes_plus .+ 1.
e1 = elfin_measurements_092220[2].-elfin_error_092220
e2 = elfin_measurements_092220[2].+elfin_error_092220

plot(Egrid, equatorial_fluxes_092220,label = L"\mathrm{equatorial\ flux\ (MMS\ 1\ measured)}", color = bipride_pink, linewidth=2, markershape=:circle);
plot!(Egrid,mms_short_a7_092220.precipitating_fluxes_mean, fillrange=x1, fillalpha = 0.2, color = c4, label=false)
plot!(Egrid,mms_short_a7_092220.precipitating_fluxes_mean, fillrange=x2, fillalpha = 0.2, color = c4, label=false)
plot!(Egrid, mms_short_a7_092220.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\phi = 3, a=7)",
        # yerror=(mms_short_a7_092220.precipitating_fluxes_minus, mms_short_a7_092220.precipitating_fluxes_plus),
        color = c4, alpha = .7, markerstrokecolor = c4, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(Egrid,mms_med_a7_092220.precipitating_fluxes_mean, fillrange=y1, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid,mms_med_a7_092220.precipitating_fluxes_mean, fillrange=y2, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid, mms_med_a7_092220.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\phi = 30, a=7)",
        # yerror=(mms_med_a7_092220.precipitating_fluxes_minus, mms_med_a7_092220.precipitating_fluxes_plus),
        color = c3, alpha = .7, markerstrokecolor = c3, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(Egrid,mms_long_a7_092220.precipitating_fluxes_mean, fillrange=z1, fillalpha = 0.2, color = bipride_orange, label=false)
plot!(Egrid,mms_long_a7_092220.precipitating_fluxes_mean, fillrange=z2, fillalpha = 0.2, color = bipride_orange, label=false)
plot!(Egrid, mms_long_a7_092220.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\phi = 300, a=7)",
        # yerror=(mms_long_a7_092220.precipitating_fluxes_minus, mms_long_a7_092220.precipitating_fluxes_plus),
        color = bipride_orange, alpha = .7, markerstrokecolor = bipride_orange, markerstrokewidth = 2, markershape=:o, markersize = 2);
plot!(elfin_measurements_092220[1][1:6], elfin_measurements_092220[2][1:6], fillrange = e1[1:end-5], fillalpha = 0.2, color = c5, label=false)        
plot!(elfin_measurements_092220[1][1:6], elfin_measurements_092220[2][1:6], fillrange = e2[1:end-5], fillalpha = 0.2, color = c5, label=false)        
plot!(elfin_measurements_092220[1][1:6], elfin_measurements_092220[2][1:6],
        yerror=elfin_error_092220[1:6].+10,
        label = L"\mathrm{precipitating\ flux\ (ELFIN\ A\ measured)}",
        color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3);
plot!(ylim =(1e1,1e10), xlim=(50,800), yscale=:log10, legend=:topright)
plot!(xlabel=L"\mathrm{Energy\ (keV)}", ylabel=L"\mathrm{Flux\ (1/cm^{2}/s/sr/MeV)}", title=L"\mathrm{Precipitating\ Flux\ Comparison\ of\ \delta\phi\ on\ 9/22/20}",
xtickfontsize=12, ytickfontsize=12, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, titlefontsize=16);
plot6 = plot!(dpi = 500,size=(800,450), margin=3mm, bottom_margin=4mm)
savefig(plot6, "images/092220_flux_pkt_a7_comparison.pdf")
savefig(plot6, "images/092220_flux_pkt_a7_comparison.png")



include("plotHelpers.jl")
# Plots packet comparison on 9/22/20 w/ a=3 and short packets dphi=3 
@time @load "results/200922_short_a3_hires_compare_storage.jld2" equatorial_fluxes_092220 mms_short_a3_hires_092220 mms_short_a3_092220;
@time @load "results/200922_long_a7_hires_compare_storage.jld2" equatorial_fluxes_092220 mms_long_a7_hires_092220 mms_long_a7_092220;
@time elfin_measurements_092220, elfin_error_092220 = extract_idl_csv("092220_time.csv", "092220_prec.csv",
                                                "092220_precerror.csv", "092220_precerror_time.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2020,9,22,9,16,38), DateTime(2020,9,22,9,16,45)) # time to sample from ELFIN measurements

x1 = mms_short_a3_hires_092220.precipitating_fluxes_mean.-mms_short_a3_hires_092220.precipitating_fluxes_minus .+ 1.
x2 = mms_short_a3_hires_092220.precipitating_fluxes_mean.+mms_short_a3_hires_092220.precipitating_fluxes_plus .+ 1.
y1 = mms_short_a3_092220.precipitating_fluxes_mean.-mms_short_a3_092220.precipitating_fluxes_minus .+ 1.
y2 = mms_short_a3_092220.precipitating_fluxes_mean.+mms_short_a3_092220.precipitating_fluxes_plus .+ 1.
xx1 = mms_long_a7_hires_092220.precipitating_fluxes_mean.-mms_long_a7_hires_092220.precipitating_fluxes_minus .+ 1.
xx2 = mms_long_a7_hires_092220.precipitating_fluxes_mean.+mms_long_a7_hires_092220.precipitating_fluxes_plus .+ 1.
yy1 = mms_long_a7_092220.precipitating_fluxes_mean.-mms_long_a7_092220.precipitating_fluxes_minus .+ 1.
yy2 = mms_long_a7_092220.precipitating_fluxes_mean.+mms_long_a7_092220.precipitating_fluxes_plus .+ 1.
e1 = elfin_measurements_092220[2].-elfin_error_092220
e2 = elfin_measurements_092220[2].+elfin_error_092220

Egrid, PAgrid = logrange(10,1000,41), 6:4:90
plot(Egrid, equatorial_fluxes_092220,label = L"\mathrm{equatorial\ flux\ (MMS\ 1\ measured)}", color = bipride_pink, linewidth=2, markershape=:circle);
# plot!(Egrid,mms_short_a3_hires_092220.precipitating_fluxes_mean, fillrange=x1, fillalpha = 0.2, color = c4, label=false)
# plot!(Egrid,mms_short_a3_hires_092220.precipitating_fluxes_mean, fillrange=x2, fillalpha = 0.2, color = c4, label=false)
# plot!(Egrid, mms_short_a3_hires_092220.precipitating_fluxes_mean,
#         label=L"\mathrm{precipitating\ flux\ hires\ }(\mathrm{simulated\ with\ }\delta\phi = 3, a=3)",
#         # yerror=(mms_short_a3_hires_092220.precipitating_fluxes_minus, mms_short_a3_hires_092220.precipitating_fluxes_plus),
#         color = c4, alpha = .7, markerstrokecolor = c4, markerstrokewidth = 2, markershape=:o, markersize = 2);
# Egrid, PAgrid = logrange(10,1000,21), 6:4:90
# plot!(Egrid,mms_short_a3_092220.precipitating_fluxes_mean, fillrange=y1, fillalpha = 0.2, color = c3, label=false)
# plot!(Egrid,mms_short_a3_092220.precipitating_fluxes_mean, fillrange=y2, fillalpha = 0.2, color = c3, label=false)
# plot!(Egrid, mms_short_a3_092220.precipitating_fluxes_mean,
#         label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\phi = 3, a=3)",
#         # yerror=(mms_short_a3_092220.precipitating_fluxes_minus, mms_short_a3_092220.precipitating_fluxes_plus),
#         color = c3, alpha = .7, markerstrokecolor = c3, markerstrokewidth = 2, markershape=:o, markersize = 2);
Egrid, PAgrid = logrange(10,1000,41), 6:4:90
plot!(Egrid,mms_long_a7_hires_092220.precipitating_fluxes_mean, fillrange=xx1, fillalpha = 0.2, color = c4, label=false)
plot!(Egrid,mms_long_a7_hires_092220.precipitating_fluxes_mean, fillrange=xx2, fillalpha = 0.2, color = c4, label=false)
plot!(Egrid, mms_long_a7_hires_092220.precipitating_fluxes_mean,
                label=L"\mathrm{precipitating\ flux\ hires\ }(\mathrm{simulated\ with\ }\delta\phi = 300, a=7)",
                # yerror=(mms_long_a7_hires_092220.precipitating_fluxes_minus, mms_long_a7_hires_092220.precipitating_fluxes_plus),
                color = c4, alpha = .7, markerstrokecolor = c4, markerstrokewidth = 2, markershape=:o, markersize = 2);
Egrid, PAgrid = logrange(10,1000,21), 6:4:90
plot!(Egrid,mms_long_a7_092220.precipitating_fluxes_mean, fillrange=yy1, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid,mms_long_a7_092220.precipitating_fluxes_mean, fillrange=yy2, fillalpha = 0.2, color = c3, label=false)
plot!(Egrid, mms_long_a7_092220.precipitating_fluxes_mean,
        label=L"\mathrm{precipitating\ flux\ }(\mathrm{simulated\ with\ }\delta\phi = 300, a=7)",
        # yerror=(mms_long_a7_092220.precipitating_fluxes_minus, mms_long_a7_092220.precipitating_fluxes_plus),
        color = c3, alpha = .7, markerstrokecolor = c3, markerstrokewidth = 2, markershape=:o, markersize = 2);        
plot!(elfin_measurements_092220[1][1:6], elfin_measurements_092220[2][1:6], fillrange = e1[1:end-5], fillalpha = 0.2, color = c5, label=false)        
plot!(elfin_measurements_092220[1][1:6], elfin_measurements_092220[2][1:6], fillrange = e2[1:end-5], fillalpha = 0.2, color = c5, label=false)        
plot!(elfin_measurements_092220[1][1:6], elfin_measurements_092220[2][1:6],
        yerror=elfin_error_092220[1:6].+10,
        label = L"\mathrm{precipitating\ flux\ (ELFIN\ A\ measured)}",
        color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3);
plot!(ylim =(1e1,1e10), xlim=(50,800), yscale=:log10, legend=:topright)
plot!(xlabel=L"\mathrm{Energy\ (keV)}", ylabel=L"\mathrm{Flux\ (1/cm^{2}/s/sr/MeV)}", title=L"\mathrm{Precipitating\ Flux\ Comparison\ of\ \delta\phi\ on\ 9/22/20}",
xtickfontsize=12, ytickfontsize=12, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, titlefontsize=16);
plot7 = plot!(dpi = 500,size=(800,450), margin=3mm, bottom_margin=4mm)
savefig(plot7, "images/092220_hires_comparison.pdf")
savefig(plot7, "images/092220_hires_comparison.png")









include("plotHelpers.jl")
Egrid, PAgrid = logrange(10,1000,21), 6:4:90

# plot PA and Energy trajectories for different dLambda
@time @load "results/demo_themis_dlambda_20.jld2" demo_E_6_20
@time @load "results/demo_themis_dlambda_30.jld2" demo_E_6_30
@time @load "results/demo_themis_dlambda_40.jld2" demo_E_6_40
L=5.2
twenty = trajectoryChecking(demo_E_6_20, [4,19,23,29,34,48], L"\delta \lambda_2 = "*"20"*L"^\circ")
thirty = trajectoryChecking(demo_E_6_30, [2,18,21,29,34,50], L"\delta \lambda_2 = "*"30"*L"^\circ")
forty = trajectoryChecking(demo_E_6_40, [1,19,21,27,34,49], L"\delta \lambda_2 = "*"40"*L"^\circ")
biglambdaplot = plot(twenty,thirty,forty, dpi = 96, layout = (1,3),
        left_margin = 10mm, top_margin = 2mm, bottom_margin=10mm,size=(1500,500))
savefig(biglambdaplot, "images/dlambda_trajectories.pdf")
savefig(biglambdaplot, "images/dlambda_trajectories.png")



# plot PA and Energy for different dPhi
@time @load "results/demo_mms_dPhi_3.jld2" demo_mms_dPhi_3
@time @load "results/demo_mms_dPhi_30.jld2" demo_mms_dPhi_30
@time @load "results/demo_mms_dPhi_300.jld2" demo_mms_dPhi_300
L=6
short = trajectoryChecking(demo_mms_dPhi_3, [1,4,5,6,7,9], L"\delta \phi = "*"3")
med = trajectoryChecking(demo_mms_dPhi_30, [1,4,5,6,7,9], L"\delta \phi = "*"30")
long = trajectoryChecking(demo_mms_dPhi_300, [1,4,5,6,7,9], L"\delta \phi = "*"300")
bigphiplot = plot(short,med,long, dpi = 96, layout = (1,3),
        left_margin = 10mm, top_margin = 2mm, bottom_margin=10mm,size=(1500,500))
savefig(bigphiplot, "images/dphi_trajectories.pdf")
savefig(bigphiplot, "images/dphi_trajectories.png")


hugeplot = plot(twenty,short,thirty,med,forty,long, dpi = 96, layout = (3,2), size=(1000,2000))
savefig(hugeplot, "images/trajectories.pdf")
savefig(hugeplot, "images/trajectories.png")
