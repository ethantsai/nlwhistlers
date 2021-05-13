@info "Hello! It's sim time!"
conffile = "setup.conf";
include("helperFunctions.jl")
@info "Imported helper functions!"

outFileBaseName = basename*"_$numParticles";
outFileDir = string("jld2_", Dates.format(now(), DateFormat("yymmdd_HH")))
mkpath(outFileDir)

conf = ConfParse(conffile)
parse_conf!(conf)
numThreads = parse(Int64, retrieve(conf, "numberOfThreads"))
Threads.nthreads() = numThreads #parse(Int64, retrieve(ConfParse("setup.conf"), "numberOfThreads"))
@info "using $(Threads.nthreads()) CPU threads"
@info "All prepped, running model"

h0, f0, eta, epsilon, resolution = generateFlatParticleDistribution(numParticles, ICrange, z0, lambda0);
numParticles = length(h0[:,1]); # new total number of particles may have changed.

## Actual model
prob = ODEProblem(eom!, h0[1,:], tspan, [eta, epsilon]);
palostcondition(H,t,integrator) = (rad2deg(asin(sqrt((2*H[4])/((calcGamma(H[2], H[4], calcb(H[5])))^2 -1)))))<(lossConeAngle); # condition: if particle enters loss cone
ixlostcondition(H,t,integrator) = 2*H[4]*calcb(H[5])<(1/saveDecimation) # condition: if I_x approaches 0
affect!(integrator) = terminate!(integrator); # terminate if condition reached
cb1 = DiscreteCallback(palostcondition,affect!);
cb2 = DiscreteCallback(ixlostcondition,affect!);

probGeneratorList, nPerBatch, percentage = generateModifiableFunction(batches); # consider making all of these constant
tick()
for i in 1:batches
    ensemble_prob = EnsembleProblem(prob,prob_func=probGeneratorList[i])
    @time sol = solve(
        ensemble_prob, Tsit5(), EnsembleThreads(),
        callback = CallbackSet(cb1, cb2), trajectories=nPerBatch,
        dtmax = resolution, dense = false, maxiters = 1e8, saveat = saveDecimation*resolution)
    @save string(outFileDir,"/",outFileBaseName,"$i.jld2") sol
    @info "$(i*percentage)% complete..."
end
@info "this finished in..."
tock()


# after that, see if we can use parallel GPU

