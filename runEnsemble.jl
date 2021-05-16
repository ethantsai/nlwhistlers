@info "Hello! It's sim time!"
conffile = "setup.conf";
include("helperFunctions.jl")
@info "Imported helper functions!"

outFileBaseName = basename*"_$numParticles";
outFileDir = string("jld2_", Dates.format(now(), DateFormat("yymmdd_HH")))
mkpath(outFileDir)

Threads.nthreads() = numThreads
@info "using $(Threads.nthreads()) CPU threads"
@info "All prepped, running model"

h0, f0, eta, epsilon, resolution = generateFlatParticleDistribution(numParticles, ICrange, z0, lambda0);
numParticles = length(h0[:,1]);      # new total number of particles may have changed.
# hf0 = [h0  f0][shuffle(1:end), :];   # shuffle so batches all take same-ish time

## Actual model
prob = ODEProblem(eom!, ~, tspan, @SVector [eta, epsilon, Omegape, omegam]);

@everywhere probGeneratorList, nPerBatch, percentage = generateModifiableFunction(batches); # consider making all of these constant

tick()
@time ensemble()
@info "this finished in..."
tock()


# after that, see if we can use parallel GPU

