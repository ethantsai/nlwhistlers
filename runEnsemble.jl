@info "Hello! It's sim time!"
include("helperFunctions.jl")
@info "Imported helperFunctions.jl"
flush(io)

threadsAvailable = length(Sys.cpu_info())
if numThreads > threadsAvailable
    Threads.nthreads() = threadsAvailable
    @warn "Assigned more threads than available."
elseif numThreads < threadsAvailable
    Threads.nthreads() = numThreads
end
@info "using $(Threads.nthreads()) CPU threads of $threadsAvailable available threads."
flush(io)

h0, f0, eta, epsilon, resolution = generateFlatParticleDistribution(numParticles, ICrange, z0, lambda0);
numParticles = length(h0[:,1]);      # new total number of particles may have changed.
@info "All prepped, running model on $numParticles particles."
flush(io)
# hf0 = [h0  f0][shuffle(1:end), :];   # shuffle so batches all take same-ish time

## Actual model
prob = ODEProblem(eom!, ~, tspan, @SVector [eta, epsilon, Omegape, omegam]);

@everywhere probGeneratorList, nPerBatch, percentage = generateModifiableFunction(batches); # consider making all of these constant
flush(io)


tick()
ensemble()
@info "Total sim time:"
tock()
flush(io)



# after that, see if we can use parallel GPU

