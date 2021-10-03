directoryname = "results/jld2_210823_21"
conffile = "setupasrun.conf"
basename = "ghostlog_96600"
num_batches = 161

include("plotHelpers.jl")

# Load in data
@time allZ, allPZ, allT, allPA, allE = loadData(directoryname, basename, num_batches);
# Count lost particles
@time lostParticles = countLostParticles(allT);
# Convert into a workable matrix
@time tVec, Zmatrix, PZmatrix, Ematrix, PAmatrix = postProcessor(allT, allZ, allPZ, allE, allPA);
# Makes an array of lost particles 
@time allPrecip, indexArray, allPrecipInitial = precipitatingParticles(tVec, Ematrix, 20);

# define dist function here
f0 = ((E::Float64, PA::Float64) -> ((1. * (E/70)^-3)*(sin(deg2rad(PA)))))
f0_092220 = function (E::Float64, PA)
    A = 3
    B0 = 1.0
    B1 = 2.5
    C = 1e4
    E0 = 300
    E1 = 50
    if E < E1
        B = B0
    elseif E <= E0
        B = B0 + (B1-B0)*(E-E0)/(E1-E0)
    elseif E > E0
        B = B1
    end
    return ((C*(E/70)^-A)*(sin(deg2rad(PA))-sin(deg2rad(8)))^B)
end
f0_102720 = function (E::Float64, PA)
    A = 3
    B0 = .2
    B1 = .3
    C = 1.5e5
    E0 = 300
    E1 = 50
    if E < E1
        B = B0
    elseif E <= E0
        B = B0 + (B1-B0)*(E-E0)/(E1-E0)
    elseif E > E0
        B = B1
    end
    return ((C*(E/70)^-A)*(sin(deg2rad(PA)))^B)
end




Egrid, PAgrid = logrange(10,1000,21), 2:4:90

equatorial_fluxes_092220, elfin_measurements_092220, prec_flux_timeseries_092220 = generate_flux_comparison(
                                                                20, f0_092220, 1.4,     # timebins, dist_func, whistler occurence rate
                                                                Egrid, PAgrid, # Ebins and PA bins to use
                                                                "092220_time.csv", "092220_prec.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                                DateTime(2020,9,22,9,16,15), DateTime(2020,9,22,9,16,50)) # time to sample from ELFIN measurements

equatorial_fluxes_102720, elfin_measurements_102720, prec_flux_timeseries_102720 = generate_flux_comparison(
                                                                20, f0_102720, .05,     # timebins, dist_func, whistler occurence rate
                                                                Egrid, PAgrid, # Ebins and PA bins to use
                                                                "102720_time.csv", "102720_prec.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                                DateTime(2020,10,27,10,34,7),DateTime(2020,10,27,10,34,40)) # time to sample from ELFIN measurements

using Plots.PlotMeasures
bipride_pink = RGB(234/255, 2/255, 112/255);
bipride_lavender = RGB(155/255, 79/255, 150/255);
bipride_blue = RGB(0/255, 56/255, 168/255);

# Plots for 9/22/20
plot(Egrid, equatorial_fluxes_092220, label = "Equatorial Flux", color = bipride_blue, linewidth=2, markershape=:circle);
plot!(Egrid, prec_flux_timeseries_092220[2:end-1], color = bipride_lavender, alpha = .5, label=false, markershape=:x);
plot!(Egrid, prec_flux_timeseries_092220[1], color = bipride_lavender, alpha = .5, label="Modelled Precipitating Flux", markershape=:x);
plot!(elfin_measurements_092220, label = "ELFIN Measured Precipitating Flux", color = bipride_pink, linewidth=3, markershape=:dtriangle);
plot!(ylim =(1e1,1e9), xlim=(50,800), yscale=:log10);
plot!(xlabel=L"\mathrm{Energy\ (keV)}", ylabel=L"\mathrm{Flux\ (keV/cm^{2}/s/sr/MeV)}", title=L"\mathrm{Flux\ Comparison\ of\ Precipitating\ Particles\ on\ 09/22/20}",
xtickfontsize=12, ytickfontsize=12, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, titlefontsize=16);
plot1 = plot!(dpi = 300,size=(800,450), margin=3mm, bottom_margin=4mm)
savefig(plot1, "092220_flux_comparison.png")

# Plots for 10/27/20
plot(Egrid, equatorial_fluxes_102720, label = "Equatorial Flux", color = bipride_blue, linewidth=2, markershape=:circle);
plot!(Egrid, prec_flux_timeseries_102720[2:end-1], color = bipride_lavender, alpha = .5, label=false, markershape=:x);
plot!(Egrid, prec_flux_timeseries_102720[1], color = bipride_lavender, alpha = .5, label="Modelled Precipitating Flux", markershape=:x);
plot!(elfin_measurements_102720, label = "ELFIN Measured Precipitating Flux", color = bipride_pink, linewidth=3, markershape=:dtriangle);
plot!(ylim =(1e4,1e11), xlim=(50,800), yscale=:log10);
plot!(xlabel=L"\mathrm{Energy\ (keV)}", ylabel=L"\mathrm{Flux\ (keV/cm^{2}/s/sr/MeV)}", title=L"\mathrm{Flux\ Comparison\ of\ Precipitating\ Particles\ on\ 10/27/20}",
xtickfontsize=12, ytickfontsize=12, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, titlefontsize=16);
plot2 = plot!(dpi = 300,size=(800,450), margin=3mm, bottom_margin=4mm)
savefig(plot2, "102720_flux_comparison.png")

# Plots for 04/29/21
elfin_measurements_042921 = extract_idl_csv("042921_time.csv", "042921_prec.csv", "ebins.csv", DateTime(2021,4,29,3,14,44),DateTime(2021,4,29,3,15,00))
plot(elfin_measurements, label = "ELFIN Measured Precipitating Flux", linewidth=3, markershape=:dtriangle)
plot!(ylim =(1e4,1e11), xlim=(50,800), yscale=:log10);


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

f_timeseries, psd_timeseries, psd_prec_timeseries =  make_psd_timeseries(Ematrix,PAmatrix,tVec, f0_102720, Egrid, PAgrid, 1.)
animatePSD("96600_animation.gif", psd_timeseries)



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
allPrecip, indexArray = precipitatingParticles(tVec, Ematrix, 10);
animatePrecipitatingParticles("test_precipanimation.gif", allPrecip, indexArray)


initial, final = 1, 1
Egrid = logrange(10,1000,16)
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
