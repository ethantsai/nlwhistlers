@info "Hello! It's sim time!"
include("helperFunctions.jl")
@info "Imported helperFunctions.jl"
flush(io)


# Setup multi threading
threadsAvailable = length(Sys.cpu_info())
if numThreads > threadsAvailable
    Threads.nthreads() = threadsAvailable
    @warn "Assigned more threads than available."
elseif numThreads < threadsAvailable
    Threads.nthreads() = numThreads
end
@info "using $(Threads.nthreads()) CPU threads of $threadsAvailable available threads."
flush(io)

## Initialize particles
h0, f0, η, ε, resolution = generateFlatParticleDistribution(numParticles, ICrange, z0, λ0);
numParticles = length(h0[:,1]);      # new total number of particles may have changed.
@info "All prepped, running model on $numParticles particles."
flush(io)
# hf0 = [h0  f0][shuffle(1:end), :];   # shuffle so batches all take same-ish time

## Setup model
prob = ODEProblem(eom!, ~, tspan, @SVector [η, ε, Omegape, omegam]);
# _prob, _alg = auto_optimize(prob) TODO: get auto-optimize to work


@everywhere probGeneratorList, nPerBatch, percentage = generateModifiableFunction(batches); # consider making all of these constant
flush(io)



## Run simulation
tick()
ensemble()
@info "Total sim time:"
tock()
flush(io)
