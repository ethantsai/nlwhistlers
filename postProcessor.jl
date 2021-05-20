directoryname = "results/63000run2"
conffile = "setupasrun.conf"
basename = "63000_2"
num_batches = 126

include("plotHelpers.jl")

# Load in data
@time allZ, allPZ, allT, allPA, allE = loadData(directoryname, basename, num_batches);
# jldsave("example.jld2"; allZ, allPZ, allT, allPA, allE)
jldsave("63000run2.jld2"; allZ, allPZ, allT, allPA, allE)
# JLD2.@load "results/63000run/63000run.jld2" allT
# JLD2.@load "results/63000run/63000run.jld2" allZ
# JLD2.@load "results/63000run/63000run.jld2" allPZ

JLD2.@load "results/63000run/63000run.jld2" allT
JLD2.@load "results/63000run/63000run.jld2" allPA
JLD2.@load "results/63000run/63000run.jld2" allE



@time tVec,  Ematrix, PAmatrix = postProcessor2(allT, allE, allPA);

# Count lost particles
@time lostParticles = countLostParticles(allT);

@time tVec, Zmatrix, PZmatrix, Ematrix, PAmatrix = postProcessor(allT, allZ, allPZ, allE, allPA);

#define dist function here
f0 = ((C::Float64, E::Float64, PA::Float64) -> ((C*(E/10)^-3)*sin(deg2rad(PA))))

initial, final = 1, 1
Egrid, PAgrid = 0:50:2000,2:4:90
f, psd_init, psd_final = recalcDistFunc(Ematrix,PAmatrix,initial,final,f0, Egrid, PAgrid);
checkDistFunc(f, psd_init, psd_final, initial, final, Egrid, PAgrid)
# savefig(plot_numLostParticles, string("particleLosses.png"))

animateNewPSD("630002recalculatedAnimation.gif", Egrid, PAgrid)

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

animatePAD("63000PADevolution.gif", 2000, 5)
animateESD("4200ESDevolution.gif", 4200, 50)
animatePSD("63000PSDevolution2.gif", 5, 50, 450)


allPrecip, indexArray = precipitatingParticles(tVec, Ematrix, 10);
animatePrecipitatingParticles("NewwwwPRECIPanimation.gif", allPrecip, indexArray)


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
trajectoryTracing(tVec, 45:2:57) # a nice 7 trajectories with varying properties to model





animDec = 1; # make a png for animation every 10 points
animScale = 10; # i.e. animscale = 10 means every 10 seconds in animation is 1 second of simulation time (increase for longer animation)
anim = @animate for i in eachindex(tVec)
    p1 = plot!(tVec[1:i], PAmatrix[:,j][1:i], label = "")
    scatter!((tVec[i], PAmatrix[:,j][i]), color = 1, label = "")
    
    # annotate!(20, 0.1*maxFraction, text("t = $(round(tVec[indexArray[i]]*Re*L/(c),digits=4)) s"), :left)
end every animDec
savename = string("results/","trajectories")
gif(anim, savename, fps = (length(tVec)/animDec)/(animScale*endTime*Re*L/(c)))
