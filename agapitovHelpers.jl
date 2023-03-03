@info "Compiling packages for simulation..."
#meta
using TickTock
using ConfParser
using Logging
using Profile
#sim
using Dates
using Random
using OrdinaryDiffEq
using StaticArrays
#process
using JLD2
using Plots
using StatsBase
@info "Packages compiled."


#######################
## Constants n stuff ##
#######################
@info "Loading constants..."
save_dir = "results_ducting/"
folder = "run21/"
mkpath(save_dir*folder)

# case specific
             #L   MLT  Kp name
test_cases = [4.5 23   3  "LO_NITE_MODEL";
              4.5 16.5 3  "LO_DUSK_MODEL";
              4.5 8.0  3  "LO_DAWN_MODEL";
              6.5 23   3  "HI_NITE_MODEL";
              6.5 16.5 3  "HI_DUSK_MODEL";
              6.5 8.0  3  "HI_DAWN_MODEL";
              ]

omega_m_cases = [0.3] # these are the different frequencies to test
L_array = test_cases[:,1]

const numParticles = 250000;
const startTime = 0;
const endTime = 15;
tspan = (startTime, endTime); # integration time

const ELo = 52;
const EHi = 1000;
const Esteps = 32; # double ELFIN E bins
const PALo = 3;
const PAHi = 15;
const PAsteps = 1300; # only used for flat particle distribution
const factor = 40; #only used for skewed particle distribution
# num particles in highest energy bin = factor * num particles in lowest energy bin
ICrange = [ELo, EHi, Esteps, PALo, PAHi, PAsteps];

const z0 = 0; # start at eq
const λ0 = 0; # start at eq

const lossConeAngle = 3;

const Bw = 300;  # pT
const a = 7;     # exp(-a * (cos(Φ/dΦ)^2))
const dPhi = 3; # exp(-a * (cos(Φ/dΦ)^2)) number of waves in each packet

const Re   = 6370e3;        # Earth radius, f64
const c    = 3e8;           # speedo lite, f64
const Beq  = 3.e-5;         # B field at equator (T), f64

const saveDecimation = 40000; # really only need first and last point
@info "Done."



function generateFlatParticleDistribution(numParticles::Int64, ICrange, L)
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

    # mu0 = .5*(IC[1]^2-1)*sin(IC[2])^2
    # q0 = sqrt(2*mu0) = sqrt((IC[1]^2-1)*sin(IC[2])^2)
    #####       [[z0 pz0                          ζ0               q0                              λ0 Φ0                 B_w0 ]             ]
    @views h0 = [[z0 sqrt(IC[1]^2 - 1)*cos(IC[2]) rand()*2*pi*dPhi sqrt((IC[1]^2-1)*sin(IC[2])^2)  λ0 rand()*2*pi*2*dPhi 0    ] for IC in f0] # creates a 5xN array with inital h0 terms
    f0 = vcat(f0...) # convert Array{Array{Float64,2},1} to Array{Float64,2}
    h0 = vcat(h0...) # since i used list comprehension it is now a nested list

    # Other ICs that are important
    # Define basic ICs and parameters
    B0          = Beq*sqrt(1. +3. *sin(λ0)^2.)/L^3.;     # starting B field at eq
    # B0          = B_eq_measured*1e-9                        # measured equatorial field strength
    Omegace0    = (1.6e-19*B0)/(9.11e-31);                    # electron gyrofreq @ the equator
    Omegape     = L;
    @info "Omegape = $Omegape"
    ε           = (Bw*1e-12)/B0;                      # Bw = 300 pT
    η           = Omegace0*L*Re/c;              # should be like 10^3
    @info "L shell of $L w/ wave amplitude of $Bw pT"
    @info "Yields ε = $ε and η = $η"

    resolution  = (1/η) / 15;  # determines max step size of the integrator
    # due to issues accuracy issues around sqrt(mu)~0, experimentally 
    # found that 15x smaller is sufficient to yield stable results
                                                    
    @info "Min integration step of $resolution"
    @info "Created Initial Conditions for $(length(h0[:,1])) particles"
    
    return h0, f0, η, ε, Omegape, resolution;
end

function generateSkewedParticleDistribution(numParticles::Int64, ICrange, L, factor)
    #=
    generates the initial conditions that biases more particles to higher energy
    for better statistics.

    Uses factor to determine how many more particles are wanted at the highest
    energy compared to the lowest energy. 
    =#
    ELo, EHi, Esteps, PALo, PAHi, PAsteps = ICrange

    @info "Generating a exponential particle distribution such that"
    @info "$Esteps steps of energy from $ELo KeV to $EHi KeV"
    @info "with $(factor)x more particles at $EHi KeV than at $ELo KeV"
    @info "$PAsteps steps of pitch angles from $PALo deg to $PAHi deg"
    
    E_bins = logrange(ELo,EHi, Int64(Esteps))

    N = numParticles
    exponent = log(factor)/log(Esteps)
    num_particles_per_E_bin = floor.((N / sum((1:Esteps).^exponent)) * (1:Esteps).^exponent)
    num_particles_per_E_bin[1:Int64(N-sum(num_particles_per_E_bin))] .+= 1
    PA_bins = [range(PALo, PAHi, length = Int64(steps)) for steps in num_particles_per_E_bin]
    
    @info "Skewed distribution with:"
    @info "$(num_particles_per_E_bin[1]) particles in first energy bin"
    @info "$(num_particles_per_E_bin[end]) particles in last energy bin"
    @info "Total particles to simulate: $(sum(num_particles_per_E_bin))."
    
    f0 = [[(E_bins[i]+511.)/511. deg2rad(PA)] for i in 1:Esteps for PA in PA_bins[i] ]

    #####       [[z0 pz0                          ζ0               q0                              λ0 Φ0               B_w0 ]             ]
    @views h0 = [[z0 sqrt(IC[1]^2 - 1)*cos(IC[2]) rand()*2*pi*dPhi sqrt((IC[1]^2-1)*sin(IC[2])^2)  λ0 rand()*2*pi*dPhi 0    ] for IC in f0] # creates a 5xN array with inital h0 terms
    f0 = vcat(f0...) # convert Array{Array{Float64,2},1} to Array{Float64,2}
    h0 = vcat(h0...) # since i used list comprehension it is now a nested list

    # Other ICs that are important
    # Define basic ICs and parameters
    B0          = Beq*sqrt(1. +3. *sin(λ0)^2.)/L^3.;     # starting B field at eq
    # B0          = B_eq_measured*1e-9                        # measured equatorial field strength
    Omegace0    = (1.6e-19*B0)/(9.11e-31);                    # electron gyrofreq @ the equator
    Omegape     = L;
    @info "Omegape = $Omegape"
    ε           = (Bw*1e-12)/B0;                      # Bw = 300 pT
    η           = Omegace0*L*Re/c;              # should be like 10^3
    @info "L shell of $L w/ wave amplitude of $Bw pT"
    @info "Yields ε = $ε and η = $η"

    resolution  = (1/η) / 15;  # determines max step size of the integrator
    # due to issues accuracy issues around sqrt(mu)~0, experimentally 
    # found that 15x smaller is sufficient to yield stable results
                                                    
    @info "Min integration step of $resolution"
    @info "Created Initial Conditions for $(length(h0[:,1])) particles"
    
    return h0, f0, η, ε, Omegape, resolution;
end


function setup_wave_model(test_cases)
    # take L, MLT, and Kp from test cases
    # return array of functions, normalizers, and coefficients
    wave_model_array = Vector{Function}()
    wave_model_coeff_array = Vector{SVector{4, Float64}}()
    wave_model_normalizer_array = Vector{Float64}()
    wave_model_shifter_array = Vector{Float64}()
    threshold = 40; #degrees, when wave model should be invalidated
    for case in eachrow(test_cases)
        wave_model(lambda) = B_w(lambda, case[3], α_ij_matrix(case[1], case[2]))
        push!(wave_model_array, wave_model)
        push!(wave_model_normalizer_array, obtain_normalizer(wave_model))
        # push!(wave_model_normalizer_array, obtain_normalizer(wave_model, threshold))
        push!(wave_model_coeff_array, agapitov_coeffs(case[3], α_ij_matrix(case[1], case[2]))  )
        push!(wave_model_shifter_array)
        # push!(wave_model_shifter_array, wave_model(threshold))
    end 
    wave_normalizer = minimum(wave_model_normalizer_array)
    return wave_model_array, wave_model_coeff_array, wave_normalizer, wave_model_shifter_array
end

function eom!(dH,H,p::SVector{9},t::Float64)
    # These equations define the motion.
    #                  lambda in radians                     
    # z, pz, zeta, mu, lambda, phi = H
    # p[1] p[2]     p[3]     p[4]    p[5] p[6]  p[7] p[8]           p[9]
    # eta, epsilon, Omegape, omegam, a,   dPhi, B_w, B_w_normalizer B_w_shifter = p

    sinλ = sin(H[5]);
    cosλ = cos(H[5]);
    g = exp(-p[5] * (cos(H[6]/(2*π*p[6]))^2))
    dg = (p[5] / (2 * π * p[6])) * sin(H[6]/(π * p[6]))
    sinζ = sin(H[3])*g;
    cosζ = cos(H[3])*g;
    
    # helper variables
    b = sqrt(1+3*sinλ^2)/(cosλ^6);
    db = (3*(27*sinλ-5*sin(3*H[5])))/(cosλ^8*(4+12*sinλ^2));
    γ = sqrt(1 + H[2]^2 + H[4]^2*b);
    K = copysign(1, H[5]) * (p[3] * (cosλ^(-5/2)))/sqrt(b/p[4] - 1);

    # B_w
    H[7] = p[8] * ((10 ^ abs( p[7][1] * (abs(rad2deg(H[5])) - p[7][4]) * exp(-abs(rad2deg(H[5])) * p[7][3] - p[7][2]))) - p[9]) * tanh(rad2deg(H[5]/deg2rad(1)))
    
    # if H[5]>deg2rad(20)
        # H[7] = 0
    # else
        # H[7] = tanh((H[5]/deg2rad(1)))
    # end
                           # 1 when >20 deg, 0 when <20 deg
                           
    # old psi 
    #     eta  * epsilon * B_w  * sqrt(2 mu   b)/gamma
    # psi = p[1] * p[2]    * H[7] * sqrt(2*H[4]*b)/γ;
    
    # new_psi = sqrt(2mu) * old_psi = q * old_psi 
    #     eta  * epsilon * B_w  * sqrt(b)/gamma
    psi = p[1] * p[2]    * H[7] * sqrt(b)/γ;
    

    # actual integration vars
    dH1 = H[2]/γ;
    dH2 = -(0.5 * H[4]^2 * db)/γ - (H[4] * psi * cosζ) - (H[4] * psi * cosζ * dg);
    dH3 = p[1]*(K*dH1 - p[4] + b/γ) + (psi*sinζ)/(H[4]*K);
    dH4 = -(psi*cosζ)/K;
    dH5 = H[2]/(γ*cosλ*sqrt(1+3*sinλ^2));
    dH6 = p[1]*(K*dH1 - p[4]);

    dH .= SizedVector{7}([ dH1, dH2, dH3, dH4, dH5, dH6, 0 ]);
    
end

function palostcondition(H,t,integrator)
    # condition: if particle enters loss cone
    b = sqrt(1+3*sin(H[5])^2)/(cos(H[5])^6);
    γ = sqrt(1 + H[2]^2 + 2*H[4]*b);
    return (rad2deg(asin(sqrt( (2*H[4])/(γ^2 -1) )))) < (lossConeAngle/2)
end

function ixlostcondition(H,t,integrator)
    # condition: if I_x approaches 0
    #      2*mu  * b 
    return 2*H[4]*sqrt(1+3*sin(H[5])^2)/(cos(H[5])^6) < 3e-5
end

function eqlostcondition(H,t,integrator)
    # condition: if particle crosses eq in negative direction
    return H[1]<=0 && H[2]<0
end

affect!(integrator) = terminate!(integrator); # terminate if condition reached
cb1 = DiscreteCallback(palostcondition,affect!);
cb2 = DiscreteCallback(ixlostcondition,affect!);
cb3 = DiscreteCallback(eqlostcondition,affect!);


# Simple calcs
calcb(b::Vector{Float64},λ::Vector{Float64}) = @. b = sqrt(1+3*sin(λ)^2)/(cos(λ)^6)
calcdb(db::Vector{Float64},λ::Vector{Float64}) = @. db = (3*(27*sin(λ)-5*sin(3*λ)))/(cos(λ)^8*(4+12*sin(λ)^2))
calcGamma(γ::Vector{Float64},pz::Vector{Float64},μ::Vector{Float64},b::Vector{Float64}) = @. γ = sqrt(1 + pz^2 + 2*μ*b)
calcK(K::Vector{Float64},b::Vector{Float64},λ::Vector{Float64}) = @. K = (Omegape * (cos(λ)^-(5/2)))/sqrt(b/omegam - 1)
calcAlpha(α::Vector{Float64},μ::Vector{Float64}, γ::Vector{Float64}) = @. α = rad2deg(asin(sqrt((2*μ)/(γ^2 - 1))))

# Useful helpers
logrange(x1, x2, n::Int64) = [10^y for y in range(log10(x1), log10(x2), length=n)]
const E_bins = logrange(ELo,EHi, Int64(Esteps))

obtain_normalizer(f::Function) = maximum(f.(0:0.01:90))^-1
obtain_normalizer(f::Function,t) = maximum(f.(0:0.01:90).-f(t))^-1

calcb!(b::Vector{Float64}, lambda) = @. b = sqrt(1+3*sin(lambda)^2)/(cos(lambda)^6)
calcGamma!(gamma::Vector{Float64}, pz, mu, b::Vector{Float64}) = @. gamma = sqrt(1 + pz^2 + 2*mu*b)
calcAlpha!(alpha::Vector{Float64}, mu, gamma::Vector{Float64}) = @. alpha = rad2deg(asin(sqrt((2*mu)/(gamma^2 - 1))))



# plot helpers
function extract(sol::EnsembleSolution)
    allZ = Vector{Vector{Float64}}();
    allPZ = Vector{Vector{Float64}}();
    allQ = Vector{Vector{Float64}}();
    allE = Vector{Vector{Float64}}();
    allZeta = Vector{Vector{Float64}}();
    allPhi = Vector{Vector{Float64}}();
    allPA = Vector{Vector{Float64}}();
    allT = Vector{Vector{Float64}}();
    allLambda = Vector{Vector{Float64}}();
    allBw = Vector{Vector{Float64}}();
    for traj in sol
    # for i in eachindex(sol)
    #     traj = sol[i];
    #     @info i
        
        vars = Array(traj');
        # try
            timesteps = length(traj.t);
            b = zeros(timesteps);
            gamma = zeros(timesteps);
            Alpha = zeros(timesteps);

            @views calcb!(b,vars[:,5]);
            @views calcGamma!(gamma,vars[:,2],0.5.*vars[:,4].^2,b);
            @views calcAlpha!(Alpha,0.5.*vars[:,4].^2,gamma);
            @views push!(allT, traj.t);
            @views push!(allZ, vars[:,1]);
            @views push!(allPZ, vars[:,2]);
            @views push!(allQ, vars[:,4]);
            @views push!(allZeta, vars[:,3]);
            @views push!(allPhi, vars[:,6]);
            @views push!(allPA, Alpha);
            @views push!(allE, @. (511*(gamma - 1)));
            @views push!(allLambda, vars[:,5]);
            @views push!(allBw, vars[:,7]);
        # catch
        #     last_positive_index = minimum(findall(x->x<=0,vars[:,4])) -1 
        #     @info "Caught negative mu"
        #     @info "index $(length(vars[:,4])-last_positive_index) from end"

        #     timesteps = length(traj.t[1:last_positive_index]);
        #     b = zeros(timesteps);
        #     gamma = zeros(timesteps);
        #     Alpha = zeros(timesteps);

        #     @views calcb!(b,vars[1:last_positive_index,5]);
        #     @views calcGamma!(gamma,vars[1:last_positive_index,2],0.5*vars[1:last_positive_index,4]^2,b);
        #     @views calcAlpha!(Alpha,0.5*vars[1:last_positive_index,4]^2,gamma);
        #     @views push!(allT, traj.t[1:last_positive_index]);
        #     @views push!(allZ, vars[1:last_positive_index,1]);
        #     @views push!(allPZ, vars[1:last_positive_index,2]);
        #     @views push!(allPA, Alpha);
        #     @views push!(allE, @. (511*(gamma - 1)));
        #     @views push!(allLambda, vars[:,5]);
        #     @views push!(allBw, vars[:,7]);
        # end

    end
    @info "$(length(sol)) particles loaded in."
    return allT, allZ, allPZ, allQ, allZeta, allPhi, allE, allPA, allLambda, allBw;
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
    @info "Matrices generated."
    return tVec, Zmatrix, PZmatrix, Ematrix, PAmatrix
end

function countLostParticles(allT::Vector{Vector{Float64}}, endTime::Float64)
    #=
    Based on time vectors, counts which ones were lost
    and at what time. Returns a Nx2 array where the first
    column is the time at which the particle was lost, and
    the 2nd column denotes the number lost at that point.
    =#
    lossCounter = []; # initialize vector

    for ntime in allT # loop thru each particle
        if maximum(ntime) != endTime # particle lost if ended early
            push!(lossCounter, maximum(ntime)); # the final entry in the time vector is when the particle got lost
        end
    end
    if isempty(lossCounter) # if particle wasn't lost, then throw in a NaN
        push!(lossCounter, NaN) 
        @warn "All particles trapped!" # since all particles trapped
    end
    lostParticles = [sort(lossCounter) collect(1:length(lossCounter))] # sort it by time, so you can see when particles got lost

    if maximum(lostParticles[:,1]) != endTime # this adds in a final hline from last particle lost to end of simulation
        lostParticles = vcat(lostParticles, [endTime (maximum(lostParticles[:,2]))]);
    end 

    @info "Total of $(lostParticles[end,end]) particles lost during sim"
    return lostParticles
end

struct Resultant_Matrix
    label::String
    numParticles::Int64
    endTime::Float64
    allZ::Vector{Vector{Float64}}
    allPZ::Vector{Vector{Float64}}
    allQ::Vector{Vector{Float64}}
    allZeta::Vector{Vector{Float64}}
    allPhi::Vector{Vector{Float64}}
    allT::Vector{Vector{Float64}}
    allPA::Vector{Vector{Float64}}
    allE::Vector{Vector{Float64}}
    allLambda::Vector{Vector{Float64}}
    allBw::Vector{Vector{Float64}}
    lostParticles::Matrix{Float64}
    tVec::Vector{Float64}
    Zmatrix::Matrix{Float64}
    PZmatrix::Matrix{Float64}
    Ematrix::Matrix{Float64}
    PAmatrix::Matrix{Float64}
end

function sol2rm(sol, label)
    @info "Extracting data..."
    @time allT, allZ, allPZ, allQ, allZeta, allPhi, allE, allPA, allLambda, allBw = extract(sol);
    @info "Processing data..."
    @time tVec, Zmatrix, PZmatrix, PAmatrix, Ematrix = postProcessor(allT, allZ, allPZ, allPA, allE);
    return Resultant_Matrix(label, length(sol), tVec[end], allZ, allPZ, allQ, allZeta, allPhi, allT, allPA, allE, allLambda, allBw, countLostParticles(allT, tVec[end]), tVec, Zmatrix, PZmatrix, Ematrix, PAmatrix)
end

function last_index_before_zero(x::Vector{Float64})
    k = 1
    @inbounds for i in eachindex(x)
        if iszero(x[i])
            k = i-1
            return k
        end
    end
    return
end