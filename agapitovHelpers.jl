@info "Compiling packages for running model..."
#meta
using TickTock
using ConfParser
using Logging
using Profile
#sim
using Dates
using Random
using JLD2
using Plots
# parallelize
using Distributed
addprocs();
@everywhere using OrdinaryDiffEq
@everywhere using StaticArrays
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
# const basename = retrieve(conf, "basename");
# const directoryname = retrieve(conf, "directoryname");
const numParticles = 105*2;
const startTime = 0;
const endTime = 5;
tspan = (startTime, endTime); # integration time
const lossConeAngle = parse(Float64, retrieve(conf, "lossConeAngle"));
const saveDecimation = parse(Float64, retrieve(conf, "saveDecimation"));
const L = parse(Float64, retrieve(conf, "L"));
const omegam = parse(Float64, retrieve(conf, "omegam"));
const Omegape = parse(Float64, retrieve(conf, "Omegape"));
const z0 = parse(Float64, retrieve(conf, "z0"));
const λ0 = parse(Float64, retrieve(conf, "lambda0"));
const dλ1 = parse(Float64, retrieve(conf, "dLambda1"));
const dλ2 = parse(Float64, retrieve(conf, "dLambda2"));
const a = parse(Float64, retrieve(conf, "a"));
const dPhi = parse(Float64, retrieve(conf, "dPhi"));
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
const B_eq_measured = parse(Float64, retrieve(conf, "B_eq_measured")); # in nT
@info "Parsed Config file: $conffile:"




function generateFlatParticleDistribution(numParticles::Int64, ICrange, z0=0::Float64, λ0=0::Float64)
    ELo, EHi, Esteps, PALo, PAHi, PAsteps = ICrange
    @info "Generating a flat particle distribution with"
    @info "$Esteps steps of energy from $ELo KeV to $EHi KeV"
    @info "$PAsteps steps of pitch angles from $PALo deg to $PAHi deg"

    nBins = PAsteps*Esteps
    N = numParticles ÷ nBins # num of particles per bin
    E_bins = logrange(ELo,EHi, Int64(Esteps))
    PA_bins = range(PALo, PAHi, length = Int64(PAsteps))
    @info "Flat distribution with $N particles/bin in $nBins bins"

    if numParticles%nBins != 0
        @warn "Truncating $(numParticles%nBins) particles for an even distribution"
    end
    if iszero(N)
        N = 1;
        @warn "Use higher number of particles next time. Simulating 1 trajectory/bin."
        @warn "Minimum number of particles to simulate is 1 particles/bin."
    end
    
    @views f0 = [[(E+511.)/511. deg2rad(PA)] for PA in PA_bins for E in E_bins for i in 1:N] # creates a 2xN array with initial PA and Energy

    #####       [[z0 pz0                          ζ0               mu0                          λ0 Φ0              ]             ]
    @views h0 = [[z0 sqrt(IC[1]^2 - 1)*cos(IC[2]) rand()*2*pi*dPhi .5*(IC[1]^2-1)*sin(IC[2])^2  λ0 rand()*2*pi*dPhi] for IC in f0] # creates a 5xN array with inital h0 terms
    f0 = vcat(f0...) # convert Array{Array{Float64,2},1} to Array{Float64,2}
    h0 = vcat(h0...) # since i used list comprehension it is now a nested list

    # Other ICs that are important
    # Define basic ICs and parameters
    # B0          = Beq*sqrt(1. +3. *sin(λ0)^2.)/L^3.;     # starting B field at eq
    B0          = B_eq_measured*1e-9                        # measured equatorial field strength
    Omegace0    = (1.6e-19*B0)/(9.11e-31);                    # electron gyrofreq @ the equator
    η           = Omegace0*L*Re/c;              # should be like 10^3
    ε           = waveAmplitudeModifier/η;    # normalized wave large amplitude, .1 for small, 15 for large
    resolution  = .1/η;                       # determines max step size of the integrator
    @info "Min integration step of $resolution"
    @info "Created Initial Conditions for $(length(h0[:,1])) particles"
    
    return h0, f0, η, ε, resolution;
end

function eom!(dH,H,p::SVector{8},t::Float64)
    # These equations define the motion.

    # z, pz, zeta, mu, lambda, phi = H
    # p[1] p[2]     p[3]     p[4]    p[5] p[6]  p[7] p[8]
    # eta, epsilon, Omegape, omegam, a,   dPhi, B_w, B_w_normalizer = p

    sinλ = sin(H[5]);
    cosλ = cos(H[5]);
    g = exp(-p[5] * (cos(H[6]/(2*π*p[6]))^2)) +  exp(-p[5] * (sin(H[6]/(2*π*p[6]))^2))  
    sinζ = g*sin(H[3]);
    cosζ = g*cos(H[3]);
    
    # helper variables
    b = sqrt(1+3*sinλ^2)/(cosλ^6);
    db = (3*(27*sinλ-5*sin(3*H[5])))/(cosλ^8*(4+12*sinλ^2));
    γ = sqrt(1 + H[2]^2 + 2*H[4]*b);
    K = copysign((p[3] * (cosλ^(-5/2)))/sqrt(b/p[4] - 1), H[5]);

    #     eta * epsilon * u * sqrt (2 mu b) / gamma
    # u = 10 ^ abs( p[7][1] * (abs(H[5]) - p[7][4]) * exp(-abs(H[5]) * p[7][3] - p[7][2])) * tanh(H[5])
    psi = p[1]*p[2]* (10 ^ abs( p[7][1] * (abs(H[5]) - p[7][4]) * exp(-abs(H[5]) * p[7][3] - p[7][2])) * tanh(H[5])) *sqrt(2*H[4]*b)/γ;
    
    # actual integration vars
    dH1 = H[2]/γ;
    dH2 = -(H[4]*db)/γ - (psi*cosζ);
    dH3 = p[1]*(K*dH1 - p[4] + b/γ) + (psi*sinζ)/(2*H[4]*K);
    dH4 = -(psi*cosζ)/K;
    dH5 = H[2]/(γ*cosλ*sqrt(1+3*sinλ^2));
    dH6 = p[1]*(K*dH1 - p[4]);

    dH .= SizedVector{6}([ dH1, dH2, dH3, dH4, dH5, dH6 ]);
end

function palostcondition(H,t,integrator)
    # condition: if particle enters loss cone
    b = sqrt(1+3*sin(H[5])^2)/(cos(H[5])^6);
    γ = sqrt(1 + H[2]^2 + 2*H[4]*b);
    return (rad2deg(asin(sqrt((2*H[4])/(γ^2 -1))))) < (lossConeAngle)
end

function ixlostcondition(H,t,integrator)
    # condition: if I_x approaches 0
    b = sqrt(1+3*sin(H[5])^2)/(cos(H[5])^6);
    return 2*H[4]*b < 10*resolution
end

affect!(integrator) = terminate!(integrator); # terminate if condition reached
cb1 = DiscreteCallback(palostcondition,affect!);
cb2 = DiscreteCallback(ixlostcondition,affect!);

function ensemble()
    ensemble_prob = EnsembleProblem(prob::ODEProblem,prob_func=probGeneratorList[i])
    sol = solve(ensemble_prob, Tsit5(), EnsembleThreads(), save_everystep=false;
                        callback=CallbackSet(cb1, cb2), trajectories=nPerBatch,
                        dtmax=resolution, linear_solver=:LapackDense, maxiters=1e8, 
                        saveat = saveDecimation*resolution)
    @save outputFileBaseName*"_$(i).jld2" sol
end

# Simple calcs
calcb(b::Vector{Float64},λ::Vector{Float64}) = @. b = sqrt(1+3*sin(λ)^2)/(cos(λ)^6)
calcdb(db::Vector{Float64},λ::Vector{Float64}) = @. db = (3*(27*sin(λ)-5*sin(3*λ)))/(cos(λ)^8*(4+12*sin(λ)^2))
calcGamma(γ::Vector{Float64},pz::Vector{Float64},μ::Vector{Float64},b::Vector{Float64}) = @. γ = sqrt(1 + pz^2 + 2*μ*b)
calcK(K::Vector{Float64},b::Vector{Float64},λ::Vector{Float64}) = @. K = (Omegape * (cos(λ)^-(5/2)))/sqrt(b/omegam - 1)
calcAlpha(α::Vector{Float64},μ::Vector{Float64}, γ::Vector{Float64}) = @. α = rad2deg(asin(sqrt((2*μ)/(γ^2 - 1))))

# Useful helpers
logrange(x1, x2, n::Int64) = [10^y for y in range(log10(x1), log10(x2), length=n)]

obtain_normalizer(f::Function) = maximum(f.(0:0.01:90))^-1

calcb!(b::Vector{Float64}, lambda::SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}) = @. b = sqrt(1+3*sin(lambda)^2)/(cos(lambda)^6)
calcGamma!(gamma::Vector{Float64}, pz::SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, mu::SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, b::Vector{Float64}) = @. gamma = sqrt(1 + pz^2 + 2*mu*b)
calcAlpha!(alpha::Vector{Float64}, mu::SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, gamma::Vector{Float64}) = @. alpha = rad2deg(asin(sqrt((2*mu)/(gamma^2 - 1))))





# plot helpers


function extract(sol::EnsembleSolution)
    allZ = Vector{Vector{Float64}}();
    allPZ = Vector{Vector{Float64}}();
    allE = Vector{Vector{Float64}}();
    allPA = Vector{Vector{Float64}}();
    allT = Vector{Vector{Float64}}();
    for traj in sol
        vars = Array(traj');
        timesteps = length(traj.t);
        b = zeros(timesteps);
        gamma = zeros(timesteps);
        Alpha = zeros(timesteps);

        @views calcb!(b,vars[:,5]);
        @views calcGamma!(gamma,vars[:,2],vars[:,4],b);
        @views calcAlpha!(Alpha,vars[:,4],gamma);
        @views push!(allT, traj.t);
        @views push!(allZ, vars[:,1]);
        @views push!(allPZ, vars[:,2]);
        @views push!(allPA, Alpha);
        @views push!(allE, @. (511*(gamma - 1)));
    end
    @info "$(length(sol)) particles loaded in..."
    return allT, allZ, allPZ, allE, allPA;
end

function postProcessor(allT::Vector{Vector{Float64}}, allZ::Vector{Vector{Float64}}, allPZ::Vector{Vector{Float64}}, allE::Vector{Vector{Float64}}, allPA::Vector{Vector{Float64}})
    #=
    This function will take the output of the model and convert them into usable m x n matrices
    where m is max number of timesteps for the longest trajectory and N is number of particles
    Arrays that are not m long are filled with NaNs
    =#
    N = length(allT); # num particles from data
    tVec = allT[findall(i->i==maximum(length.(allT)),length.(allT))[1]]; # turns all time vectors into a single time vector spanning over the longest trajectory
    timeseriesLength = length(tVec); # all vectors must be this tall to ride
    Zmatrix = fill(NaN,timeseriesLength,N); 
    PZmatrix = fill(NaN,timeseriesLength,N); 
    Ematrix = fill(NaN,timeseriesLength,N); 
    PAmatrix = fill(NaN,timeseriesLength,N); 
    # iterate over each matrix column and fill with each vector
    for i = 1:N
        @views Zmatrix[1:length(allT[i]),i] = allZ[i]
        @views PZmatrix[1:length(allT[i]),i] = allPZ[i]
        @views Ematrix[1:length(allT[i]),i] = allE[i]
        @views PAmatrix[1:length(allT[i]),i] = allPA[i]
    end
    @info "Matrices generated..."
    return tVec, Zmatrix, PZmatrix, Ematrix, PAmatrix
end