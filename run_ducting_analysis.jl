#=
run by typing julia run_ducting_analysis.jl into terminal
after specifying cases and simulation parameters in agapitovHelpers.jl
=#

include("agapitovHelpers.jl")
include("agapitovmodel.jl")
tick()

wave_model_array, wave_model_coeff_array, wave_normalizer = setup_wave_model(test_cases)
total_sims = length(L_array)*length(omega_m_cases)
@info "Performing simulations on $(length(L_array)) scenarios with $(length(omega_m_cases)) cases each."
@info "For a total of $total_sims sims."

println()
# loop for simulation
for case_index in eachindex(L_array)
    h0, f0, η, ε, Omegape, resolution = generateFlatParticleDistribution(numParticles, ICrange, L_array[case_index]);

    wave_model_coeffs = wave_model_coeff_array[case_index]

    for n in eachindex(omega_m_cases)
        omega_m = omega_m_cases[n]
        @info "Starting new sim: $(test_cases[:,end][case_index]) with omega_m = $omega_m"
        savename = save_dir*folder*test_cases[:,end][case_index]*"_"*string(omega_m)[end]*"_$numParticles.jld2"

        params = @SVector [η, ε, Omegape, omega_m, a, dPhi, wave_model_coeffs, wave_normalizer];
        prob = ODEProblem(eom!, ~, tspan, params);
        prob_func = ((prob,i,repeat) -> remake(prob, u0 = h0[i,:], p = @SVector [η, ε, Omegape, omega_m, a, dPhi, wave_model_coeffs, wave_normalizer]))
        ensemble_prob = EnsembleProblem(prob::ODEProblem,prob_func=prob_func)
        
        tick()
        sol = solve(ensemble_prob, Tsit5(), EnsembleThreads(), save_everystep=false;
                                callback=CallbackSet(cb1, cb2, cb3), trajectories=numParticles,
                                dtmax=resolution, maxiters=1e8,
                                saveat = saveDecimation*resolution, kwargshandle=KeywordArgSilent)
        tock()
        @info "Saving solution to $savename"
        @save savename sol
        
        
        @info "Done! ($(((case_index-1)*length(omega_m_cases))+n)/$total_sims)"

        println()
    end
end
tock()