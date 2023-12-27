# functions for setting up the simulations
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
    #####       [[z0 pz0                          ζ0               q0                              λ0 Φ0                 B_w0 K0]             ]
    @views h0 = [[z0 sqrt(IC[1]^2 - 1)*cos(IC[2]) rand()*2*pi*dPhi sqrt((IC[1]^2-1)*sin(IC[2])^2)  λ0 rand()*2*pi*2*dPhi 0    0] for IC in f0] # creates a 5xN array with inital h0 terms
    f0 = vcat(f0...) # convert Array{Array{Float64,2},1} to Array{Float64,2}
    h0 = vcat(h0...) # since i used list comprehension it is now a nested list

    # Other ICs that are important
    # Define basic ICs and parameters
    B0          = Beq*sqrt(1. +3. *sin(λ0)^2.)/L^3.;     # starting B field at eq
    # B0          = B_eq_measured*1e-9                        # measured equatorial field strength
    Omegace0    = (1.6e-19*B0)/(9.11e-31);                    # electron gyrofreq @ the equator
    Omegape     = L;
    @info "Ω_pe = $Omegape"
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

    #####       [[z0 pz0                          ζ0               q0                              λ0 Φ0               B_w0 K0]             ]
    @views h0 = [[z0 sqrt(IC[1]^2 - 1)*cos(IC[2]) rand()*2*pi*dPhi sqrt((IC[1]^2-1)*sin(IC[2])^2)  λ0 rand()*2*pi*dPhi 0    0] for IC in f0] # creates a 5xN array with inital h0 terms
    f0 = vcat(f0...) # convert Array{Array{Float64,2},1} to Array{Float64,2}
    h0 = vcat(h0...) # since i used list comprehension it is now a nested list

    # Other ICs that are important
    # Define basic ICs and parameters
    B0          = Beq*sqrt(1. +3. *sin(λ0)^2.)/L^3.;     # starting B field at eq
    # B0          = B_eq_measured*1e-9                        # measured equatorial field strength
    Omegace0    = (1.6e-19*B0)/(9.11e-31);                    # electron gyrofreq @ the equator
    Omegape     = L;
    @info "Ω_pe = $Omegape"
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

function setup_single_wave_model(L, MLT, Kp)
    # take L, MLT, and Kp from test cases
    # return array of functions, normalizers, and coefficients
    wave_model(lambda) = B_w(lambda, Kp, α_ij_matrix(L, MLT))
    wave_model_normalizer = obtain_normalizer(wave_model)
    wave_model_coeff = agapitov_coeffs(Kp, α_ij_matrix(L, MLT))
    wave_model_shifter = 0
    # threshold = 40; #degrees, when wave model should be invalidated
    # wave_model_normalizer = obtain_normalizer(wave_model, threshold)
    # wave_model_shifter = wave_model(threshold)
    return wave_model, wave_model_coeff, wave_model_normalizer, wave_model_shifter
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
        push!(wave_model_shifter_array, 0)
        # push!(wave_model_shifter_array, wave_model(threshold))
    end 
    wave_normalizer = minimum(wave_model_normalizer_array)
    return wave_model_array, wave_model_coeff_array, wave_normalizer, wave_model_shifter_array
end

#=
Here are the different particle tracing models:
eom_og! -> original equations of motion with only dphi and dlambda mods
eom_Bw! -> equations of motion with modified wave amplitude as a function of latitude
eom_Bwωm! -> eom with modified wave amplitude and wave frequency as a function of latitude
eom_wna1! -> eom w/ modded Bw with a slightly oblique wave normal angle model
eom_wna1! -> eom w/ modded Bw with a moderately oblique wave normal angle model
eom_wna3! -> eom w/ modded Bw with a very oblique wave normal angle model

The run_model arguments are described below:
numParticles = number of particles to simulate in total
ICrange = [ELo, EHi, Esteps, PALo, PAHi, PAsteps]
L = what L shell to use
MLT = what MLT to use for wave model
Kp = what Kp to use for wave model
omega_m = what wave frequency to use
save_decimation = decimation factor for data saved, typically 
use ~40000 when you care about only final vs initial values...
use ~40 when you care about details of the trajectories
particle_distribution = X --> generate skewed particle distribution with X
times more particles at highest energy bin as compared to lowest energy bin,
typically X>30 or so to get better statistics at higher energy bins, but
if X=1, then will default to flat particle distribution

Optional args:
model = either freq, wna1, wna2, or wna3; if unused, but previous 2 are specified, then will default to Bw
=#
function run_model(numParticles::Int64,
    ICrange,
    L::Float64,
    MLT::Float64,
    Kp::Float64,
    omega_m::Float64,
    save_decimation::Int64,
    particle_distribution::Int64,
    model::Any
    )

    if particle_distribution==1
        h0, f0, η, ε, Omegape, resolution = generateFlatParticleDistribution(numParticles, ICrange, L);
    elseif particle_distribution>1
        h0, f0, η, ε, Omegape, resolution = generateSkewedParticleDistribution(numParticles, ICrange, L, particle_distribution);
    else
        @error "Particle_distribution must be greater than 1"
        return
    end

    @info "ω_m = $omega_m"

    wave_model, wave_model_coeff, wave_model_normalizer, wave_model_shifter = setup_single_wave_model(L, MLT, Kp)
    @info "Generating wave model at L=$L, MLT=$MLT, and Kp=$Kp"

    if model[1] == "og"
        params = (η, ε, Omegape, omega_m, model[2], model[3], model[4], model[5], model[6]);
        prob = ODEProblem(eom_og!, ~, tspan, params);
    elseif model == "bw"
        params = (η, ε, Omegape, omega_m, a, dPhi, wave_model_coeff, wave_model_normalizer, wave_model_shifter);
        prob = ODEProblem(eom_Bw!, ~, tspan, params);
    elseif model == "freq"
        params = (η, ε, Omegape, omega_m, a, dPhi, wave_model_coeff, wave_model_normalizer, wave_model_shifter);
        prob = ODEProblem(eom_Bwωm!, ~, tspan, params);
    elseif model == "wna1"
        params = (η, ε, Omegape, omega_m, a, dPhi, wave_model_coeff, wave_model_normalizer, wave_model_shifter);
        prob = ODEProblem(eom_wna1!, ~, tspan, params);
    elseif model == "wna2"
        params = (η, ε, Omegape, omega_m, a, dPhi, wave_model_coeff, wave_model_normalizer, wave_model_shifter);
        prob = ODEProblem(eom_wna1!, ~, tspan, params);
    elseif model == "wna3"
        params = (η, ε, Omegape, omega_m, a, dPhi, wave_model_coeff, wave_model_normalizer, wave_model_shifter);
        prob = ODEProblem(eom_wna3!, ~, tspan, params);
    else
        @error "model needs to be either og, bw, freq, wna1, wna2, or wna3"
        return
    end
    prob_func = ((prob,i,repeat) -> remake(prob, u0 = h0[i,:], p = params))
    ensemble_prob = EnsembleProblem(prob::ODEProblem,prob_func=prob_func)
    
    sol = solve(ensemble_prob, Tsit5(), EnsembleThreads(), save_everystep=false;
                            callback=CallbackSet(cb3), trajectories=numParticles,
                            dtmax=resolution, maxiters=1e8,
                            saveat = save_decimation*resolution, kwargshandle=KeywordArgSilent)
    
    @info "Done!"

    println()
    return sol
end
