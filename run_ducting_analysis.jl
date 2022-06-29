# ingest multiple conjunction events
# need L, MLT, Kp, produce a wave function for each model
# need B_mag, B_freq
# run model for 5 bounce periods for 300000 particles
# produce flux as a function of pitch angle
# find j_precip particles (< PA_LC)
# find j_trapped particles (from PA_LC to PA_2LC)
# plot j_precip/j_trapped as a function of energy from sim and ELFIN
# ???
# profit

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

include("agapitovHelpers.jl")

       #L   MLT  Kp B_eq
test = [6   23.5 2 ;
        6.5 1    1  ]

include("agapitovmodel.jl")

wave_model_array = Vector{Function}()
wave_model_normalizer_array = Vector{Float64}()

for case in eachrow(test)
    wave_model(lambda) = B_w(lambda, case[3], α_ij_matrix(case[1], case[2]))
    push!(wave_model_normalizer_array, obtain_normalizer(wave_model))
    push!(wave_model_array, wave_model)
end

@everywhere would_be_nice(lambda) = wave_model_array[1](lambda)


h0, f0, η, ε, resolution = generateFlatParticleDistribution(numParticles, ICrange, z0, λ0);
numParticles = length(h0[:,1]);    
  # new total number of particles may have changed.
params = @SVector [η, ε, Omegape, omegam, a, dPhi, wave_model_array[1], obtain_normalizer(wave_model_array[1])];
prob = ODEProblem(eom!, ~, tspan, params);
prob_func = ((prob,i,repeat) -> remake(prob, u0 = h0[i,:], p = @SVector [η, ε, Omegape, omegam, a, dPhi, wave_model_array[1], obtain_normalizer(wave_model_array[1])]))
ensemble_prob = EnsembleProblem(prob::ODEProblem,prob_func=prob_func)

@info "Solving..."
@time sol = solve(ensemble_prob, Tsit5(), EnsembleThreads(), save_everystep=false;
                        callback=CallbackSet(cb1, cb2), trajectories=numParticles,
                        dtmax=resolution, linear_solver=:LapackDense, maxiters=1e8, 
                        saveat = saveDecimation*resolution, kwargshandle=KeywordArgSilent)
@info "Sim complete! Plotting..."
@time quicklook(sol)

@info "Saving..."


@time allT, allZ, allPZ, allE, allPA = extract(sol);
@time tVec, Zmatrix, PZmatrix, Ematrix, PAmatrix = postProcessor(allT, allZ, allPZ, allPA, allE);


# using copysign
# 29.2 sec

# using passed in function as parameter
# 52.75

# using an externally called function
# 54.922

# using an externally called function defined @everywhere
# 51.520930

