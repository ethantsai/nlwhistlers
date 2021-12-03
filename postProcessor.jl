include("plotHelpers.jl")
Egrid, PAgrid = logrange(10,1000,21), 6:4:90

# themis_lolat = load_resultant_matrix("210429_themis_lolat", "results/themis_lolat_m1pro", "210429_themis_lolat_96600", "setupasrun.conf", 161);
# themis_midlat = load_resultant_matrix("210429_themis_midlat", "results/themis_midlat_m1pro", "210429_themis_midlat_96600", "setupasrun.conf", 161);
# themis_hilat = load_resultant_matrix("210429_themis_hilat", "results/themis_hilat_m1pro", "210429_themis_hilat_96600", "setupasrun.conf", 161);

# themis_1pkt = load_resultant_matrix("210429_themis_1pkt", "results/210429_themis_1pkt2", "210429_themis_1pkt2_96600", "setupasrun.conf", 161);
# themis_30pkt = load_resultant_matrix("210429_themis_30pkt", "results/210429_themis_30pkt2", "210429_themis_30pkt2_96600", "setupasrun.conf", 161);
# themis_100pkt = load_resultant_matrix("210429_themis_100pkt", "results/210429_themis_100pkt2", "210429_themis_100pkt2_96600", "setupasrun.conf", 161);

# mms_lowerlat = load_resultant_matrix("200922_mms_lowerlat", "results/mms_lowerlat_m1pro", "200922_mms_lowerlat_96600", "setupasrun.conf", 161);
# mms_lolat = load_resultant_matrix("200922_mms_lolat", "results/mms_lolat_m1pro", "200922_mms_lolat_96600", "setupasrun.conf", 161);
# mms_midlat = load_resultant_matrix("200922_mms_midlat", "results/mms_midlat_m1pro", "200922_mms_midlat_96600", "setupasrun.conf", 161);
# mms_hilat = load_resultant_matrix("200922_mms_hilat", "results/mms_hilat_m1pro", "200922_mms_hilat_96600", "setupasrun.conf", 161);

# @time @save "results/210429_themis_lolat.jld2" themis_lolat 
# @time @save "results/210429_themis_midlat.jld2" themis_midlat 
# @time @save "results/210429_themis_hilat.jld2" themis_hilat 

# @time @save "results/210429_themis_1pkt.jld2" themis_1pkt 
# @time @save "results/210429_themis_30pkt.jld2" themis_30pkt 
# @time @save "results/210429_themis_100pkt.jld2" themis_100pkt

# @time @save "results/200922_mms_lowerlat.jld2" mms_lowerlat;
# @time @save "results/200922_mms_lolat.jld2" mms_lolat
# @time @save "results/200922_mms_midlat.jld2" mms_midlat;
# @time @save "results/200922_mms_hilat.jld2" mms_hilat;


@time @load "results/210429_themis_hilat.jld2" themis_hilat;
@time @load "results/210429_themis_midlat.jld2" themis_midlat; 
@time @load "results/210429_themis_lolat.jld2" themis_lolat;

@time @load "results/210429_themis_1pkt.jld2" themis_1pkt 
@time @load "results/210429_themis_30pkt.jld2" themis_30pkt 
@time @load "results/210429_themis_100pkt.jld2" themis_100pkt

@time @load "results/200922_mms_lowerlat.jld2" mms_lowerlat;
@time @load "results/200922_mms_lolat.jld2" mms_lolat;
@time @load "results/200922_mms_midlat.jld2" mms_midlat;
@time @load "results/200922_mms_hilat.jld2" mms_hilat;
@info "Loaded result matrices."

@time equatorial_fluxes_042921 = calc_equatorial_fluxes(themis_hilat, f0_042921);
@time themis_hilat_042921  = export_results("042921_themis_hilat",  calc_prec_flux(themis_hilat, 10,f0_042921,.02));
@time themis_midlat_042921 = export_results("042921_themis_midlat", calc_prec_flux(themis_midlat,10,f0_042921,.04));
@time themis_lolat_042921  = export_results("042921_themis_lolat",  calc_prec_flux(themis_lolat, 10,f0_042921,5*.06));
@time themis_1pkt_042921   = export_results("042921_themis_1pkt",   calc_prec_flux(themis_1pkt, 10,f0_042921,.02));
@time themis_30pkt_042921  = export_results("042921_themis_30pkt",  calc_prec_flux(themis_30pkt,10,f0_042921,.02));
@time themis_100pkt_042921 = export_results("042921_themis_100pkt", calc_prec_flux(themis_100pkt,10,f0_042921,.02));
@time elfin_measurements_042921, elfin_error_042921 = extract_idl_csv("042921_time.csv", "042921_prec.csv",
                                                "042921_precerror.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2021,4,29,3,14,45),DateTime(2021,4,29,3,15,0)) # time to sample from ELFIN measurements
@time @save "results/210429_lat_storage.jld2" equatorial_fluxes_042921 themis_hilat_042921 themis_midlat_042921 themis_lolat_042921;
@time @save "results/210429_pkt_storage.jld2" equatorial_fluxes_042921 themis_1pkt_042921 themis_30pkt_042921 themis_100pkt_042921;


@time equatorial_fluxes_092220 = calc_equatorial_fluxes(mms_lolat, f0_092220);
@time mms_lowerlat_092220  = export_results("092220_mms_lowerlat",  calc_prec_flux(mms_lowerlat, 10,f0_092220,0.06));
@time mms_lolat_092220  = export_results("092220_mms_lolat",  calc_prec_flux(mms_lolat, 10,f0_092220,0.06));
@time mms_midlat_092220 = export_results("092220_mms_midlat", calc_prec_flux(mms_midlat,10,f0_092220,0.02));
@time mms_hilat_092220  = export_results("092220_mms_hilat",  calc_prec_flux(mms_hilat, 10,f0_092220,0.02));
@time @save "results/200922_lat_storage.jld2" equatorial_fluxes_092220 mms_lowerlat_092220 mms_lolat_092220 mms_midlat_092220 mms_hilat_092220;




equatorial_fluxes_092220, elfin_measurements_092220, prec_flux_timeseries_092220 = generate_flux_comparison(mms_midlat,
                                                                10, f0_092220, 1.2,     # timebins, dist_func, whistler occurence rate
                                                                Egrid, PAgrid, # Ebins and PA bins to use
                                                                "092220_time.csv", "092220_prec.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                                DateTime(2020,9,22,9,16,15), DateTime(2020,9,22,9,16,50)) # time to sample from ELFIN measurements

# Plots for 9/22/20
plot(Egrid, equatorial_fluxes_092220, label = "Equatorial Flux", color = bipride_blue, linewidth=2, markershape=:circle);
plot!(Egrid, prec_flux_timeseries_092220[2:end-1], color = bipride_lavender, alpha = .5, label=false, markershape=:x);
plot!(Egrid, prec_flux_timeseries_092220[1], color = bipride_lavender, alpha = .5, label="Modelled Precipitating Flux", markershape=:x);
plot!(elfin_measurements_092220, label = "ELFIN Measured Precipitating Flux", color = bipride_pink, linewidth=3, markershape=:dtriangle);
plot!(ylim =(1e1,1e9), xlim=(50,800), yscale=:log10);
plot!(xlabel=L"\mathrm{Energy\ (keV)}", ylabel=L"\mathrm{Flux\ (1/cm^{2}/s/sr/MeV)}", title=L"\mathrm{Flux\ Comparison\ of\ Precipitating\ Particles\ on\ 09/22/20}",
xtickfontsize=12, ytickfontsize=12, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, titlefontsize=16);
plot1 = plot!(dpi = 300,size=(800,450), margin=3mm, bottom_margin=4mm)
#savefig(plot1, "092220_flux_comparison.pdf")
equatorial_fluxes_102720, elfin_measurements_102720, prec_flux_timeseries_102720 = generate_flux_comparison(themis_hilat,
                                                                10, f0_102720, .05,     # timebins, dist_func, whistler occurence rate
                                                                Egrid, PAgrid, # Ebins and PA bins to use
                                                                "102720_time.csv", "102720_prec.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                                DateTime(2020,10,27,10,34,7),DateTime(2020,10,27,10,34,40)) # time to sample from ELFIN measurements

# Plots for 10/27/20
plot(Egrid, equatorial_fluxes_102720, label = "Equatorial Flux", color = bipride_blue, linewidth=2, markershape=:circle);
plot!(Egrid, prec_flux_timeseries_102720[2:end-1], color = bipride_lavender, alpha = .5, label=false, markershape=:x);
plot!(Egrid, prec_flux_timeseries_102720[1], color = bipride_lavender, alpha = .5, label="Modelled Precipitating Flux", markershape=:x);
plot!(elfin_measurements_102720, label = "ELFIN Measured Precipitating Flux", color = bipride_pink, linewidth=3, markershape=:dtriangle);
plot!(ylim =(1e4,1e11), xlim=(50,800), yscale=:log10);
plot!(xlabel=L"\mathrm{Energy\ (keV)}", ylabel=L"\mathrm{Flux\ (1/cm^{2}/s/sr/MeV)}", title=L"\mathrm{Flux\ Comparison\ of\ Precipitating\ Particles\ on\ 10/27/20}",
xtickfontsize=12, ytickfontsize=12, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, titlefontsize=16);
plot2 = plot!(dpi = 300,size=(800,450), margin=3mm, bottom_margin=4mm)
#savefig(plot2, "102720_flux_comparison.pdf")

equatorial_fluxes_042921, elfin_measurements_042921, prec_flux_timeseries_042921 = generate_flux_comparison(themis_midlat,
                                                                10, f0_042921, .05,     # timebins, dist_func, whistler occurence rate
                                                                Egrid, PAgrid, # Ebins and PA bins to use
                                                                "042921_time.csv", "042921_prec.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                                DateTime(2021,4,29,3,14,45),DateTime(2021,4,29,3,15,0)) # time to sample from ELFIN measurements                                                                

# Plots for 04/29/21
plot(Egrid, equatorial_fluxes_042921, label = "Equatorial Flux", color = bipride_blue, linewidth=2, markershape=:circle);
plot!(Egrid, prec_flux_timeseries_042921[2:end-1], color = bipride_lavender, alpha = .5, label=false, markershape=:x);
plot!(Egrid, prec_flux_timeseries_042921[1], color = bipride_lavender, alpha = .5, label="Modelled Precipitating Flux", markershape=:x);
plot!(elfin_measurements_042921, label = "ELFIN Measured Precipitating Flux", color = bipride_pink, linewidth=3, markershape=:dtriangle);
plot!(ylim =(1e2,1e9), xlim=(50,800), yscale=:log10);
plot!(xlabel=L"\mathrm{Energy\ (keV)}", ylabel=L"\mathrm{Flux\ (1/cm^{2}/s/sr/MeV)}", title=L"\mathrm{Flux\ Comparison\ of\ Precipitating\ Particles\ on\ 4/29/21}",
xtickfontsize=12, ytickfontsize=12, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, titlefontsize=16);
plot3 = plot!(dpi = 300,size=(800,450), margin=3mm, bottom_margin=4mm)
# savefig(plot3, "042921_flux_comparison.pdf")




mms_hilat = export_results("200922_mms_hilat", equatorial_fluxes_092220, prec_flux_timeseries_092220)
plot(Egrid, mms_hilat.precipitating_fluxes_mean, yerror=(mms_hilat.precipitating_fluxes_minus, mms_hilat.precipitating_fluxes_plus), ylim =(1e2,1e9), xlim=(50,800), yscale=:log10)
plot!(Egrid, mms_hilat.equatorial_fluxes)
@save "200922_data_storage.jld2" mms_midlat

themis_midlat = export_results("210429_themis_midlat", equatorial_fluxes_042921, prec_flux_timeseries_042921)
plot(Egrid, themis_midlat.precipitating_fluxes_mean, yerror=(themis_midlat.precipitating_fluxes_minus, themis_midlat.precipitating_fluxes_plus), ylim =(1e2,1e9), xlim=(50,800), yscale=:log10)
plot!(Egrid, themis_midlat.equatorial_fluxes)
@save "210429_data_storage.jld2" themis_hilat_100wav themis_hilat_nowav themis_hilat themis_midlat themis_lolat








trapped_fluxes = generate_trapped_psd(f0_092220, 1.0)
simulated_ratio = [prec_flux_timeseries[i]./((trapped_fluxes[i]).-prec_flux_timeseries[i]) for i in 1:length(trapped_fluxes)]
plot(Egrid, simulated_ratio, xlim=(50,1000),legend=false, xscale=:log10)
elfin_prec = extract_idl_csv("092220_time.csv", "092220_prec.csv", "ebins.csv", DateTime(2020,9,22,9,16,15), DateTime(2020,9,22,9,16,50))
elfin_trap = extract_idl_csv("092220_time.csv", "092220_perp.csv", "ebins.csv", DateTime(2020,9,22,9,16,15), DateTime(2020,9,22,9,16,50))
plot!(elfin_prec[1], elfin_prec[2] ./ elfin_trap[2], xlim=(50,1000), yscale=:log10)

trapped_fluxes = generate_trapped_psd(f0_102720, 0.1)
simulated_ratio = [prec_flux_timeseries[i]./((trapped_fluxes[i]./45.).-prec_flux_timeseries[i]) for i in 1:length(trapped_fluxes)]
plot(Egrid, simulated_ratio, xlim=(50,1000),legend=false, xscale=:log10)
elfin_prec = extract_idl_csv("102720_time.csv", "102720_prec.csv", "ebins.csv", DateTime(2020,10,27,10,34,7),DateTime(2020,10,27,10,34,40))
elfin_trap = extract_idl_csv("102720_time.csv", "102720_perp.csv", "ebins.csv", DateTime(2020,10,27,10,34,7),DateTime(2020,10,27,10,34,40))
plot!(elfin_prec[1], elfin_prec[2] ./ elfin_trap[2], xlim=(50,1000), yscale=:log10)


plot(Egrid,binned_psd_prec_timeseries[14].*psd_0i.*Egrid./1000,xlim=(Egrid[1], Egrid[end]), xscale=:log10)
extract_idl_csv("time.csv", "prec.csv", "ebins.csv", DateTime(2020,9,22,9,16,30),DateTime(2020,9,22,9,16,35))
plot(ebins,flux, xlim = (ebins[1], 1e3), xscale =:log10, ylim = (1, 1e8), yscale = :log10)



plot!(xlim=(10,1000), ylim=(1e2,1e11), yscale=:log10)




animateNewPSD("21000_animation.gif", Egrid, PAgrid)

plot_numLostParticles = plot(lostParticles[:,1]*Re/(c), lostParticles[:,2], legend=false,
    xlim = [0,maximum(maximum(allT))*Re/(c)], ylim = [0,numParticles],
    title = "Number of Lost Particles", xlabel="Time (s)", ylabel="Number of Particles")
savefig(plot_numLostParticles, string("particleLosses.png"))

plot3 = plot(
    (1000*allT*Re)./(c), allPA, title = "Pitch Angles",
    legend=false,xlabel="Time (ms)", ylabel="Pitch Angle (deg)",
    xlim = [0,50], ylim = [0,90]);
plot4 = plot(
    (1000*allT*Re)./(c), allE, title = "Energies",
    legend=false,xlabel="Time (ms)", ylabel="Energy (keV)",
    xlim = [0,50], ylim = [0,3000]);
bigplot = plot(plot3,plot4, dpi = 300, layout = (2,1))




animatePAD("63000PADevolution.gif", 2000, 5)
animateESD("4200ESDevolution.gif", 4200, 50)
animatePSD("63000PSDevolution2.gif", 5, 50, 450)

# makes 
allPrecip, indexArray = precipitatingParticles(themis_100pkt, 10);
animatePrecipitatingParticles("test_precipanimation.gif", themis_100pkt, 10, 0.015)


initial, final = 1, 5
Egrid = logrange(10,1000,20)
PAgrid = 2:4:90
f, psd_init, psd_final = recalcDistFunc(Ematrix,PAmatrix,initial,final,f0_102720, Egrid, PAgrid, 1.0);
checkDistFunc(f, psd_init, psd_final, initial, final, Egrid, PAgrid)
savefig(plot_numLostParticles, string("particleLosses.png"))


# function f0(E, alpha)
#     Estar = 30; # kev
#     return (E/Estar * (sin(deg2rad(alpha)) - sin(deg2rad(lossConeAngle))))
# end

# check ICs
# plot(1:numParticles,sort(Ematrix[1,:]))
# plot(1:numParticles,sort(PAmatrix[1,:]))


# energyDistributionofPrecipiatingParticles
# energy and time phase space dist at bin of time
# can use large bin of time 
# 1 simulation time = 0.1 s
# each bin = 10 simulation time

# trajectory of conjunctions for mms and elfin
# plot whistler wave spectrum
# elfin data


# plot trajectory
trajectoryChecking(45)
trajectoryTracing([7300,7309,6900,6905,6910,6918,6927], 10)





animDec = 1; # make a png for animation every 10 points
animScale = 10; # i.e. animscale = 10 means every 10 seconds in animation is 1 second of simulation time (increase for longer animation)
anim = @animate for i in eachindex(tVec)
    p1 = plot!(tVec[1:i], PAmatrix[:,j][1:i], label = "")
    scatter!((tVec[i], PAmatrix[:,j][i]), color = 1, label = "")
    
    # annotate!(20, 0.1*maxFraction, text("t = $(round(tVec[indexArray[i]]*Re*L/(c),digits=4)) s"), :left)
end every animDec
savename = string("results/","trajectories")
gif(anim, savename, fps = (length(tVec)/animDec)/(animScale*endTime*Re*L/(c)))
