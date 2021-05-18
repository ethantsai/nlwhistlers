@info "Compiling packages for running model..."
#meta
using TickTock
using ConfParser
using Logging
using Profile
#sim
using Dates
using Random
using StaticArrays
using Distributed
using OrdinaryDiffEq
using JLD2
using Plots
@info "Packages compiled, running model"

#######################
## Constants n stuff ##
#######################

const Re   = 6370e3;        # Earth radius, f64
const c    = 3e8;           # speedo lite, f64
const Beq  = 3.e-5;         # B field at equator (T), f64

# Parsing the conf file
conffile = "setup.conf";
conf = ConfParse(conffile)
parse_conf!(conf)
const basename = retrieve(conf, "basename");
const directoryname = retrieve(conf, "directoryname");
numParticles = parse(Int64, retrieve(conf, "numberOfParticles"));
startTime = parse(Float64, retrieve(conf, "startTime"));
endTime = parse(Float64, retrieve(conf, "endTime"));
tspan = (startTime, endTime); # integration time
const lossConeAngle = parse(Float64, retrieve(conf, "lossConeAngle"));
const saveDecimation = parse(Float64, retrieve(conf, "saveDecimation"));
const L = parse(Float64, retrieve(conf, "L"));
const omegam = parse(Float64, retrieve(conf, "omegam"));
const Omegape = parse(Float64, retrieve(conf, "Omegape"));
const z0 = parse(Float64, retrieve(conf, "z0"));
const lambda0 = parse(Float64, retrieve(conf, "lambda0"));
const waveAmplitudeModifier = parse(Float64, retrieve(conf, "waveAmplitudeModifier"));
const ELo = parse(Float64, retrieve(conf, "ELo"));
const EHi = parse(Float64, retrieve(conf, "EHi"));
const Esteps = parse(Float64, retrieve(conf, "Esteps"));
const PALo = parse(Float64, retrieve(conf, "PALo"));
const PAHi = parse(Float64, retrieve(conf, "PAHi"));
const PAsteps = parse(Float64, retrieve(conf, "PAsteps"));
ICrange = [ELo, EHi, Esteps, PALo, PAHi, PAsteps];
const batches = parse(Int64, retrieve(conf, "batches"));
const numThreads = parse(Int64, retrieve(conf, "numberOfThreads"))
@info "Parsed Config file: $conffile:"

###################
## Setup Logging ##
###################

function setupDirectories(directoryname::String)
    #=
    directory name
     └ jld2_yymmdd_HH
      | setup_asrun.conf    
      | yymmdd_HHMMSS.log
      └ basename_numparticles_numbatch.jld2
    =#
    mkpath(directoryname)
    outputFileDirectory = string(directoryname, "/jld2_", Dates.format(now(), DateFormat("yymmdd_HH")))
    mkpath(outputFileDirectory)
    outputFileBaseName = outputFileDirectory*"/"*basename*"_$numParticles";
    loggingFileName = string(outputFileDirectory,"/",Dates.format(now(), DateFormat("yymmdd_HHMMSS")),".log")
    asrunConfFileName = outputFileDirectory*"/setup asrun.conf"
    run(`cp setup.conf $asrunConfFileName`) # copy setup over to as run file
    run(`cp setup.conf $loggingFileName`) # copy setup over to log file

    @info "View logs here: "*loggingFileName
    io = open(loggingFileName, "a")
    global_logger(ConsoleLogger(io))
    @info "^^Used this setup configuration file^^^"
    flush(io)

    return outputFileBaseName, io
end
outputFileBaseName, io = setupDirectories(directoryname)


####################
## Initialization ##
####################

function generateFlatParticleDistribution(numParticles::Int64, ICrange, z0=0::Float64, lambda0=0::Float64)
    ELo, EHi, Esteps, PALo, PAHi, PAsteps = ICrange
    @info "Generating a flat particle distribution with"
    @info "Energy from $ELo KeV to $EHi KeV in $Esteps KeV increments"
    @info "PA from $PALo deg to $PAHi deg in $PAsteps deg increments"
    flush(io)

    nBins = ((EHi-ELo)÷Esteps + 1) * ((PAHi-PALo)÷PAsteps + 1)
    N = numParticles ÷ nBins # num of particles per bin
    @info "Flat distribution with $N particles/bin in $nBins bins"
    flush(io)

    if numParticles%nBins != 0
        @warn "Truncating $(numParticles%nBins) particles for an even distribution"
        flush(io)
    end
    if iszero(N)
        N = 1;
        @warn "Use higher number of particles next time. Simulating 1 trajectory/bin."
        @warn "Minimum number of particles to simulate is 1 particles/bin."
        flush(io)
    end

    @views f0 = [[(E+511.)/511. deg2rad(PA)] for PA in PALo:PAsteps:PAHi for E in ELo:Esteps:EHi for i in 1:N] # creates a 2xN array with initial PA and Energy
    #####       [[z0 pz0                          zeta0       mu0                          lambda0]             ]
    @views h0 = [[z0 sqrt(IC[1]^2 - 1)*cos(IC[2]) rand()*2*pi .5*(IC[1]^2-1)*sin(IC[2])^2  lambda0] for IC in f0] # creates a 5xN array with inital h0 terms
    f0 = vcat(f0...) # convert Array{Array{Float64,2},1} to Array{Float64,2}
    h0 = vcat(h0...) # since i used list comprehension it is now a nested list

    # Other ICs that are important
    # Define basic ICs and parameters
    B0          = Beq*sqrt(1. +3. *sin(lambda0)^2.)/L^3.;     # starting B field at eq
    Omegace0    = (1.6e-19*B0)/(9.11e-31);                    # electron gyrofreq @ the equator
    # todo, make lambda IC a funtion of z and L, get rid of dep on Beq
    eta         = Omegace0*L*Re/c;              # should be like 10^3
    epsilon     = waveAmplitudeModifier/eta;    # normalized wave large amplitude, .1 for small, 15 for large
    resolution  = .1/eta;                       # determines max step size of the integrator
    @info "Created Initial Conditions for $(length(h0[:,1])) particles"
    flush(io)
    
    return h0, f0, eta, epsilon, resolution;
end

@everywhere function generateModifiableFunction(batches)
    #=
    Takes in the initial condition and splits them into batches.
    This way, we can feed each one in and get a percent completeness during sim.
    =#
    probGeneratorList = [];
    nPerBatch = numParticles÷batches;
    for j in 0:batches-1
        truncatedIC =  h0[nPerBatch*j+1:(nPerBatch*j+nPerBatch),:]
        push!(probGeneratorList, ((prob,i,repeat) -> remake(prob, u0 = truncatedIC[i,:], p = @SVector [eta, epsilon, Omegape, omegam])))
    end
    percentage = (round(100/batches))
    @info "Each batch will simulate $nPerBatch particles for $(endTime-startTime) dt and correspond with $percentage%"
    flush(io)
    return probGeneratorList, nPerBatch, percentage
end

##################
## Math n stuff ##
##################

function eom!(dH,H,p::SVector{4, Float64},t::Float64)
    # These equations define the motion.

    # z, pz, zeta, mu, lambda = H
    # eta, epsilon, Omegape, omegam = p
    sinlambda = sin(H[5]);
    coslambda = cos(H[5]);
    sinzeta = sin(H[3]);
    coszeta = cos(H[3]);

    u = .5*(tanh(H[5]/deg2rad(1))+1);
    b = sqrt(1+3*sinlambda^2)/(coslambda^6);
    db = (3*(27*sinlambda-5*sin(3*H[5])))/(coslambda^8*(4+12*sinlambda^2));
    gamma = sqrt(1 + H[2]^2 + 2*H[4]*b);
    K = (p[3] * (coslambda^(-5/2)))/sqrt(b/p[4] - 1);
    psi = p[1]*p[2]*u*sqrt(2*H[4]*b)/gamma;

    dH1 = H[2]/gamma;
    dH2 = -(H[4]*db)/gamma - (psi*coszeta);
    dH3 = p[1]*(K*dH1 - p[4] + b/gamma) + (psi*sinzeta)/(2*H[4]*K);
    dH4 = -(psi*coszeta)/K;
    dH5 = H[2]/(gamma*coslambda*sqrt(1+3*sinlambda^2)); 

    dH .= SizedVector{5}([ dH1, dH2, dH3, dH4, dH5 ]);
end


function palostcondition(H,t,integrator)
    # condition: if particle enters loss cone
    b = sqrt(1+3*sin(H[5])^2)/(cos(H[5])^6);
    gamma = sqrt(1 + H[2]^2 + 2*H[4]*b);
    return  (rad2deg(asin(sqrt((2*H[4])/(gamma^2 -1))))) < (lossConeAngle)
end

function ixlostcondition(H,t,integrator)
    # condition: if I_x approaches 0
    b = sqrt(1+3*sin(H[5])^2)/(cos(H[5])^6);
    return 2*H[4]*b < (1/saveDecimation)
end

affect!(integrator) = terminate!(integrator); # terminate if condition reached
cb1 = DiscreteCallback(palostcondition,affect!);
cb2 = DiscreteCallback(ixlostcondition,affect!);

function ensemble()
    for i in 1:batches
        ensemble_prob = EnsembleProblem(prob::ODEProblem,prob_func=probGeneratorList[i])
        tick()
        sol = solve(ensemble_prob, Tsit5(), EnsembleThreads(), save_everystep=false;
                            callback=CallbackSet(cb1, cb2), trajectories=nPerBatch,
                            dtmax=resolution, linear_solver=:LapackDense, maxiters=1e8, 
                            saveat = saveDecimation*resolution)
        @save outputFileBaseName*"_$(i).jld2" sol
        @info "$nPerBatch particles simulated in:"
        tock()
        @info "$(i*percentage)% complete..."
        flush(io)
    end
end


# Simple calcs
calcb(b::Vector{Float64},lambda::Vector{Float64}) = @. b = sqrt(1+3*sin(lambda)^2)/(cos(lambda)^6)
calcdb(db::Vector{Float64},lambda::Vector{Float64}) = @. db = (3*(27*sin(lambda)-5*sin(3*lambda)))/(cos(lambda)^8*(4+12*sin(lambda)^2))
calcGamma(gamma::Vector{Float64},pz::Vector{Float64},mu::Vector{Float64},b::Vector{Float64}) = @. gamma = sqrt(1 + pz^2 + 2*mu*b)
calcK(K::Vector{Float64},b::Vector{Float64},lambda::Vector{Float64}) = @. K = (Omegape * (cos(lambda)^-(5/2)))/sqrt(b/omegam - 1)
calcAlpha(alpha::Vector{Float64},mu::Vector{Float64}, gamma::Vector{Float64}) = @. alpha = rad2deg(asin(sqrt((2*mu)/(gamma^2 - 1))))



# ```
# Useful links
# Community stuff: https://discourse.julialang.org/c/domain/models/21
# Docs for ensemble stuff: https://diffeq.sciml.ai/dev/features/ensemble/
# Marginal Histograms: http://docs.juliaplots.org/latest/recipes/#Marginal-Histograms
# Dumb colorbar ticks issue: https://github.com/JuliaPlots/Plots.jl/issues/2308
# ```