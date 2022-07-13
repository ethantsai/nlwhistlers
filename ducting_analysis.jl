include("agapitovHelpers.jl")
include("agapitovmodel.jl")


# plot all wave models
wave_model_array, wave_model_normalizer_array, wave_model_coeff_array = setup_wave_model(test_cases)
i=1
plot1 = plot(0:.01:90,tanh.(0:.01:90)*wave_model_normalizer_array[i].*wave_model_array[i].(0:.01:90), label = test_cases[:,end][i], linewidth=2)
i=2
plot!(0:.01:90,tanh.(0:.01:90)*wave_model_normalizer_array[i].*wave_model_array[i].(0:.01:90), label = test_cases[:,end][i], linewidth=2)
i=3
plot!(0:.01:90,tanh.(0:.01:90)*wave_model_normalizer_array[i].*wave_model_array[i].(0:.01:90), label = test_cases[:,end][i], linewidth=2, linestyle=:dashdot)
i=4
plot!(0:.01:90,tanh.(0:.01:90)*wave_model_normalizer_array[i].*wave_model_array[i].(0:.01:90), label = test_cases[:,end][i], linewidth=2)
i=5
plot!(0:.01:90,tanh.(0:.01:90)*wave_model_normalizer_array[i].*wave_model_array[i].(0:.01:90), label = test_cases[:,end][i], linewidth=2, linestyle=:dash)

savefig(plot1, "images/wave_models.png")



rm_array = Vector{Resultant_Matrix}()
# loading loop
for case_index in eachindex(L_array)
    for omega_m in omega_m_cases
        loadname = save_dir*folder*test_cases[:,end][case_index]*"_"*string(omega_m)[end]*"_$numParticles.jld2"
        @info "Loading solution from $loadname..."
        @time @load loadname sol
        label = test_cases[:,end][case_index]*"_"*string(omega_m)[end]
        push!(rm_array, sol2rm(sol, label));
        @info "Loaded $label."
    end
end


function prec_to_trap_ratio(rm::Resultant_Matrix)
    initial_E = [rm.allE[i][1] for i = 1:length(rm.allT)]
    final_E = [rm.allE[i][end] for i = 1:length(rm.allT)]
    initial_PA = [rm.allPA[i][1] for i = 1:length(rm.allT)]
    final_PA = [rm.allPA[i][end] for i = 1:length(rm.allT)]
    # E_prec = sort(Ematrix[1,findall(isnan,truncated_matrix[end,:])])
    # E_trap = sort(Ematrix[1,findall(!isnan,truncated_Ematrix[end,:])])
  
    multiplier = 3

    trap_range = findall(x->multiplier*(lossConeAngle+0.002)>x>(lossConeAngle+0.002), final_PA)
    loss_range = findall(x->x<(lossConeAngle+0.002), final_PA)
    E_prec = sort(final_E[loss_range])
    E_trap = sort(final_E[trap_range])
  
    # PA_min = round(minimum(initial_PA))-1
    # PA_max = round(maximum(initial_PA))+1
    # final_PA_dist = fit(Histogram, (final_PA), PA_min:1:PA_max)
    # # initial_PA_dist = fit(Histogram, round.(initial_PA), PA_min:1:PA_max)
    # plot(PA_min:1:(PA_max-1), final_PA_dist.weights, label=false)
  
    j_prec = fit(Histogram, E_prec, logrange(ELo, EHi, Esteps+1)).weights
    j_trap = fit(Histogram, E_trap, logrange(ELo, EHi, Esteps+1)).weights

    f = Vector{Union{Float64,Missing}}(undef, length(j_trap))
    for i in eachindex(j_trap)
        if !iszero(j_prec)
            f[i] = (multiplier-1)*j_prec[i] ./ j_trap[i]
        end
    end

    return f, j_prec, j_trap
end


############
scenario = 1
###########^
num_omega_m = length(omega_m_cases)

index_array = 1:3 #( 3*(scenario-1) + 1 ) : ( 3*(scenario-1) + num_omega_m )
if maximum(index_array) > length(L_array)*length(index_array)
    @error "There's only $(length(L_array)) scenarios; pick a valid ID"
end
results_array = [prec_to_trap_ratio(x) for x in rm_array]


normalizer = 1/results_array[2][1][7] # normalize to dayside precip @ 100 keV to 100%
i = 1
plot2 = plot(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(100,2000), ylim=(1e-2, 1))
i = 2
plot!(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(100,2000), ylim=(1e-2, 1))
i = 3
plot!(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(100,2000), ylim=(1e-2, 1))
i = 4
plot!(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(100,2000), ylim=(1e-2, 1))

savefig(plot2, "images/ratios.png")


rm = rm_array[3]

initial_E = [rm.allE[i][1] for i = 1:length(rm.allT)]
final_E = [rm.allE[i][end] for i = 1:length(rm.allT)]
initial_PA = [rm.allPA[i][1] for i = 1:length(rm.allT)]
final_PA = [rm.allPA[i][end] for i = 1:length(rm.allT)]

trap_range = findall(x->2*(lossConeAngle+0.002)>x>(lossConeAngle+0.002), final_PA)
PA_inclusion_range = findall(x->(5<=x<=10), initial_PA)
E_inclusion_range = findall(x->200<x<210, initial_E)
loss_range = findall(x->x<(lossConeAngle+0.002), final_PA)

x_range = intersect(PA_inclusion_range, E_inclusion_range)


E_prec = sort(final_E[loss_range])
E_trap = sort(final_E[trap_range])

j_prec = fit(Histogram, E_prec, logrange(ELo, EHi, Esteps+1))
j_trap = fit(Histogram, E_trap, logrange(ELo, EHi, Esteps+1))
f = j_prec.weights ./ j_trap.weights





plot(E_bins, 1e-10.+results_array[1][1], label="omega_m = 0.2", xscale=:log10, yscale=:log10)
plot!(E_bins, 1e-10.+results_array[2][1], label="omega_m = 0.3")
plot!(E_bins, 1e-10.+results_array[3][1], label="omega_m = 0.4", ylim=(1e-3, 1))









@load "results_ducting/run2/ELA_ND_210105T1454_4_38400.jld2" sol

allT = Vector{Vector{Float64}}();
allBw = Vector{Vector{Float64}}();
allLambda = Vector{Vector{Float64}}();
for traj in sol

    vars = Array(traj');
    timesteps = length(traj.t);
    
    @views push!(allT, traj.t);
    @views push!(allLambda, rad2deg.(vars[:,5]));
    @views push!(allBw, wave_model_normalizer_array[1]*wave_model_array[1].(rad2deg.(vars[:,5])))
end


