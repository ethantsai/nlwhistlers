include("plotHelpers.jl")
Egrid, PAgrid = logrange(10,1000,21), 6:4:90

# Plots for 4/29/21
@load "210429_data_storage.jld2" themis_hilat_100wav themis_hilat_nowav themis_hilat themis_midlat themis_lolat
elfin_measurements_042921 = extract_idl_csv("042921_time.csv", "042921_prec.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2021,4,29,3,14,45),DateTime(2021,4,29,3,15,0)) # time to sample from ELFIN measurements                                                                

plot(Egrid, themis_hilat.equatorial_fluxes, label = "Equatorial Flux", color = bipride_blue, linewidth=2, markershape=:circle);
plot!(Egrid, themis_hilat.precipitating_fluxes_mean,
        yerror=(themis_hilat.precipitating_fluxes_minus, themis_hilat.precipitating_fluxes_plus),
        color = bipride_orange, alpha = .5, label="dLambda2 = 40 deg", markershape=:x)
plot!(Egrid, themis_midlat.precipitating_fluxes_mean,
        yerror=(themis_lolat.precipitating_fluxes_minus, themis_lolat.precipitating_fluxes_plus),
        color = bipride_lavender, alpha = .5, label="dLambda2 = 30 deg", markershape=:x)
plot!(Egrid, themis_lolat.precipitating_fluxes_mean,
        yerror=(themis_lolat.precipitating_fluxes_minus, themis_lolat.precipitating_fluxes_plus),
        color = bipride_lavender, alpha = .5, label="dLambda2 = 20 deg", markershape=:x)
plot!(elfin_measurements_042921, label = "ELFIN Measured Precipitating Flux", color = bipride_pink, linewidth=3, markershape=:dtriangle);
plot!(ylim =(1e2,1e9), xlim=(50,800), yscale=:log10);
plot!(xlabel=L"\mathrm{Energy\ (keV)}", ylabel=L"\mathrm{Flux\ (1/cm^{2}/s/sr/MeV)}", title=L"\mathrm{Flux\ Comparison\ of\ Precipitating\ Particles\ on\ 4/29/21}",
xtickfontsize=12, ytickfontsize=12, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, titlefontsize=16);
plot1 = plot!(dpi = 300,size=(800,450), margin=3mm, bottom_margin=4mm)
savefig(plot1, "042921_flux_lat_comparison.pdf")
savefig(plot1, "042921_flux_lat_comparison.png")


plot(Egrid, themis_hilat.equatorial_fluxes, label = "Equatorial Flux", color = bipride_blue, linewidth=2, markershape=:circle);
plot!(Egrid, themis_hilat_100wav.precipitating_fluxes_mean,
        yerror=(themis_hilat_100wav.precipitating_fluxes_minus, themis_hilat_100wav.precipitating_fluxes_plus),
        color = bipride_orange, alpha = .5, label="dPhi = 100", markershape=:x)
plot!(Egrid, themis_hilat.precipitating_fluxes_mean,
        yerror=(themis_hilat.precipitating_fluxes_minus, themis_hilat.precipitating_fluxes_plus),
        color = bipride_lavender, alpha = .5, label="dPhi = 30", markershape=:x)
plot!(Egrid, themis_hilat_nowav.precipitating_fluxes_mean,
        yerror=(themis_hilat_nowav.precipitating_fluxes_minus, themis_hilat_nowav.precipitating_fluxes_plus),
        color = bipride_lavender, alpha = .5, label="dPhi = 1", markershape=:x)
plot!(elfin_measurements_042921, label = "ELFIN Measured Precipitating Flux", color = bipride_pink, linewidth=3, markershape=:dtriangle);
plot!(ylim =(1e2,1e9), xlim=(50,800), yscale=:log10);
plot!(xlabel=L"\mathrm{Energy\ (keV)}", ylabel=L"\mathrm{Flux\ (1/cm^{2}/s/sr/MeV)}", title=L"\mathrm{Flux\ Comparison\ of\ Precipitating\ Particles\ on\ 4/29/21}",
xtickfontsize=12, ytickfontsize=12, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, titlefontsize=16);
plot2= plot!(dpi = 300,size=(800,450), margin=3mm, bottom_margin=4mm)
savefig(plot2, "042921_flux_pkt_comparison.pdf")
savefig(plot2, "042921_flux_pkt_comparison.png")





# Plots for 9/22/20
@load "200922_data_storage.jld2" mms_midlat mms_hilat
elfin_measurements_092220 = extract_idl_csv("092220_time.csv", "092220_prec.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2020,9,22,9,16,15), DateTime(2020,9,22,9,16,50)) # time to sample from ELFIN measurements

plot(Egrid, mms_midlat.equatorial_fluxes, label = "Equatorial Flux", color = bipride_blue, linewidth=2, markershape=:circle);
plot!(Egrid, 2*mms_hilat.precipitating_fluxes_mean,
        yerror=(mms_hilat.precipitating_fluxes_minus, mms_hilat.precipitating_fluxes_plus),
        color = bipride_orange, alpha = .5, label="dLambda2 = 40 deg", markershape=:x)
plot!(Egrid, mms_midlat.precipitating_fluxes_mean,
        yerror=(mms_midlat.precipitating_fluxes_minus, mms_midlat.precipitating_fluxes_plus),
        color = bipride_lavender, alpha = .5, label="dLambda2 = 30 deg", markershape=:x)
plot!(elfin_measurements_092220, label = "ELFIN Measured Precipitating Flux", color = bipride_pink, linewidth=3, markershape=:dtriangle);
plot!(ylim =(1e1,1e9), xlim=(50,800), yscale=:log10);
plot!(xlabel=L"\mathrm{Energy\ (keV)}", ylabel=L"\mathrm{Flux\ (1/cm^{2}/s/sr/MeV)}", title=L"\mathrm{Flux\ Comparison\ of\ Precipitating\ Particles\ on\ 09/22/20}",
xtickfontsize=12, ytickfontsize=12, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, titlefontsize=16);
plot2 = plot!(dpi = 300,size=(800,450), margin=3mm, bottom_margin=4mm)
savefig(plot2, "092220_flux_comparison.pdf")
savefig(plot2, "092220_flux_comparison.png")