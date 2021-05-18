dirname = "results/4200testrun"
conffile = "setupasrun.conf"
basename = "ghost_4200"
num_batches = 20

include("plotHelpers.jl")

# Load in data
@time allZ, allPZ, allT, allPA, allE = loadData(dirname, basename, num_batches);


# Count lost particles
lostParticles = countLostParticles(allT);

tVec, Zmatrix, PZmatrix, Ematrix, PAmatrix = postProcessor(allT, allZ, allPZ, allE, allPA);

#define dist function here
f0 = ((C::Float64, E::Float64, PA::Float64) -> ((C*(E/10)^-3)*sin(deg2rad(PA))))

initial, final = 1, 50
f, psd_init, psd_final = recalcDistFunc(Ematrix,PAmatrix,initial,final,f0,0:50:1000,0:5:90);
checkDistFunc(f, psd_init, psd_final, initial, final)

animateNewPSD("recalculatedAnimation.gif")

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
    xlim = [0,50], ylim = [0,1000]);
bigplot = plot(plot3,plot4, dpi = 300, layout = (2,1))

# animatePAD("10000PADevolution.gif", 2000, 5)
# animateESD("10000ESDevolution.gif", 1500, 50)
# animatePSD("10000PSDevolution.gif", 5, 50, 150)


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
