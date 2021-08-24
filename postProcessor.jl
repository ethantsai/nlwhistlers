directoryname = "results/jld2_210821_19"
conffile = "setupasrun.conf"
basename = "ghostlog_63000"
num_batches = 126

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
f0 = ((C::Float64, E::Float64, PA::Float64) -> ((C*(E/70)^-3)*(sin(deg2rad(PA)))))
f0_092220 = function (E::Float64, PA)
    A = 3
    B0 = 1
    B1 = 2.5
    C = 1e4
    E0 = 300
    E1 = 70
    if E < 70
        B = B0
    elseif E <= 300
        B = B0 + (B1-B0)*(E-E0)/(E1-E0)
    elseif E > 300
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
    E1 = 70
    if E < 70
        B = B0
    elseif E <= 300
        B = B0 + (B1-B0)*(E-E0)/(E1-E0)
    elseif E > 300
        B = B1
    end
    return ((C*(E/70)^-A)*(sin(deg2rad(PA)))^B)
end


Egrid, PAgrid = logrange(10,1000,21), 2:4:90

equatorial_fluxes, elfin_measurements, prec_flux_timeseries = generate_flux_comparison(
                                                                20, f0_092220, 1.0,     # timebins, dist_func, whistler occurence rate
                                                                Egrid, PAgrid, # Ebins and PA bins to use
                                                                "092220_time.csv", "092220_prec.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                                DateTime(2020,9,22,9,16,15), DateTime(2020,9,22,9,16,50)) # time to sample from ELFIN measurements
animate_flux_comparison("recalc_prec_flux_092220_3.gif", equatorial_fluxes, elfin_measurements, prec_flux_timeseries,1e2,1e11)

equatorial_fluxes, elfin_measurements, prec_flux_timeseries = equatorial_fluxes, elfin_measurements, prec_flux_timeseries = generate_flux_comparison(
                                                                20, f0_102720, 0.1,     # timebins, dist_func, whistler occurence rate
                                                                Egrid, PAgrid, # Ebins and PA bins to use
                                                                "102720_time.csv", "102720_prec.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                                DateTime(2020,10,27,10,34,7),DateTime(2020,10,27,10,34,40)) # time to sample from ELFIN measurements
animate_flux_comparison("recalc_prec_flux_102720.gif", equatorial_fluxes, elfin_measurements, prec_flux_timeseries,1e2,1e11)



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
allPrecip, indexArray = precipitatingParticles(tVec, Ematrix, 10);
animatePrecipitatingParticles("test_precipanimation.gif", allPrecip, indexArray)


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
