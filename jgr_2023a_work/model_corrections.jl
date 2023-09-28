include("agapitovHelpers.jl")
include("agapitovmodel.jl")
include("agapitovPlotHelpers.jl")
include("modelCorrectionHelpers.jl")

############
### Demo ###
############

############
scenario = "ELA_SA_211102T2218"
start = DateTime(2021,11,2,22,18,21)
stop = DateTime(2021,11,2,22,18,35)
############
@time model_bins, model_ratio, ELFIN_bins, ELFIN_ratio, elfin_p2t_error, sim_ratio_sm, normalizer, L, MLT, Kp, ELFIN_label, index = extract_elfin_comparison(scenario,start,stop);
@time E_new, spl_model, spl_elfin, k = obtain_K_E(ELFIN_bins, ELFIN_ratio, model_bins, model_ratio);
@time k_mult, degrees, E2λ, E_of_λ, E_of_λ_model_bins, E_of_λ_ELFIN_bins = obtain_E2λ(k, E_new, model_bins, ELFIN_bins, 0.0001, omega_m_cases[1], L);
@time B_w_og, B_w_mod, B_w_mod_ip, B_w_mod_ip_rad = obtain_B_w_mod(test_cases, scenario, k_mult, 200, 10, 20);

###########################
# plots to check B_w mods #
###########################
# check imported ELFIN and model data
plot(E_bins, model_ratio, label=" Model: L=$L, MLT=$MLT, Kp=$Kp", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 6))
plot!(ELFIN_bins, ELFIN_ratio, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 0, label=ELFIN_label)
ELFIN_model_compare_plot = plot!()

# check interpolation
scatter(E_bins, model_ratio)
plot!(E_new, spl_model(E_new))
scatter!(ELFIN_bins, ELFIN_ratio)
plot!(E_new, spl_elfin(E_new))
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 6))
plot!(E_new, k, xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 6), legend=false)
interpolation_plot = plot!()

# check energy to latitude conversion
scatter(E_of_λ_model_bins, model_ratio)
plot!(E_of_λ, spl_model(E_new), label = "model")
scatter!(E_of_λ_ELFIN_bins, ELFIN_ratio)
plot!(E_of_λ, spl_elfin(E_new), label = "ELFIN observation")
plot!(E_of_λ, k, yscale=:log10, xlim=(20,50), ylim=(1e-2, 6), label = "k")
plot!(degrees, k_mult(degrees), label = "k_mult", legend=:bottomright)
conversion_plot = plot!()

# check original B_w, modified B_w, and smoothed/interpolated B_w
plot(B_w_og, 0,90, label = "Original B_w")
plot!(B_w_mod, 0,90, label = "Corrected B_w")
plot!(0:0.01:90, B_w_mod_ip_rad.(deg2rad.(0:0.01:90)), label = "Smoothed + Corrected B_w")

################
### End Demo ###
################



###################
### Sim Attempt ###
###################
############
scenario = "ELA_SA_211102T2218"
start = DateTime(2021,11,2,22,18,21)
stop = DateTime(2021,11,2,22,18,35)
############
@time model_bins, model_ratio, ELFIN_bins, ELFIN_ratio, elfin_p2t_error, sim_ratio_sm, normalizer, L, MLT, Kp, ELFIN_label, index = extract_elfin_comparison(scenario,start,stop);
@time E_new, spl_model, spl_elfin, k = obtain_K_E(ELFIN_bins, ELFIN_ratio, model_bins, model_ratio);
@time k_mult, degrees, E2λ, E_of_λ, E_of_λ_model_bins, E_of_λ_ELFIN_bins = obtain_E2λ(k, E_new, model_bins, ELFIN_bins, 0.0001, omega_m_cases[1], L);
@time B_w_og, B_w_mod, B_w_mod_ip, B_w_mod_ip_rad = obtain_B_w_mod(test_cases, scenario, k_mult, 200, 10, 20);

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
const saveDecimation = 500; # really only need first and last point


h0, f0, η, ε, Omegape, resolution = generateSkewedParticleDistribution(numParticles, ICrange, L, factor);

omega_m = omega_m_cases[1]
@info "Starting new sim: $(test_cases[:,end][index]) with omega_m = $omega_m"
savename = save_dir*folder*test_cases[:,end][index]*"_"*string(omega_m)[end]*"_$numParticles.jld2"

params = @SVector [η, ε, Omegape, omega_m, a, dPhi, B_w_mod_ip_rad];
prob = ODEProblem(eom_mod!, ~, tspan, params);
prob_func = ((prob,i,repeat) -> remake(prob, u0 = h0[i,:], p = @SVector [η, ε, Omegape, omega_m, a, dPhi, B_w_mod_ip_rad]))
ensemble_prob = EnsembleProblem(prob::ODEProblem,prob_func=prob_func)

sol = solve(ensemble_prob, Tsit5(), EnsembleThreads(), save_everystep=false;
                        callback=CallbackSet(cb3), trajectories=numParticles,
                        dtmax=resolution, maxiters=1e8,
                        saveat = saveDecimation*resolution, kwargshandle=KeywordArgSilent)
@info "Saving solution to $savename"
@save savename sol

allT, allZ, allPZ, allQ, allZeta, allPhi, allE, allPA, allLambda, allBw = extract(sol);
plot([rad2deg.(Lambda) for Lambda in allLambda], allBw)