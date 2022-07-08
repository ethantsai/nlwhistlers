include("agapitovHelpers.jl")
include("agapitovmodel.jl")





# plot all wave models
wave_model_array, wave_model_normalizer_array, wave_model_coeff_array = setup_wave_model(test_cases)
for i in eachindex(wave_model_array)
    if isone(i)
        plot(0:.01:90, wave_model_normalizer_array[i].*wave_model_array[i].(0:.01:90), label = test_cases[:,end][i])
    else
        plot!(0:.01:90, wave_model_normalizer_array[i].*wave_model_array[i].(0:.01:90), label = test_cases[:,end][i])
    end
end




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
  
    trap_range = findall(x->2*(lossConeAngle+0.002)>x>(lossConeAngle+0.002), final_PA)
    loss_range = findall(x->x<(lossConeAngle+0.002), final_PA)
    E_prec = sort(final_E[loss_range])
    E_trap = sort(initial_E[trap_range])
  
    # PA_min = round(minimum(initial_PA))-1
    # PA_max = round(maximum(initial_PA))+1
    # final_PA_dist = fit(Histogram, (final_PA), PA_min:1:PA_max)
    # # initial_PA_dist = fit(Histogram, round.(initial_PA), PA_min:1:PA_max)
    # plot(PA_min:1:(PA_max-1), final_PA_dist.weights, label=false)
  
    j_prec = fit(Histogram, E_prec, logrange(ELo, EHi, Esteps+1))
    j_trap = fit(Histogram, E_trap, logrange(ELo, EHi, Esteps+1))
    f = j_prec.weights ./ j_trap.weights
    return f, j_prec.weights, j_trap.weights
  end


############
scenario = 1
###########^
num_omega_m = length(omega_m_cases)

index_array = ( 3*(scenario-1) + 1 ) : ( 3*(scenario-1) + num_omega_m )
if maximum(index_array) > length(L_array)*num_omega_m
    @error "There's only $(length(L_array)) scenarios; pick a valid ID"
end
results_array = [prec_to_trap_ratio(x) for x in rm_array]



rm = rm_array[3]

initial_E = [rm.allE[i][1] for i = 1:length(rm.allT)]
final_E = [rm.allE[i][end] for i = 1:length(rm.allT)]
initial_PA = [rm.allPA[i][1] for i = 1:length(rm.allT)]
final_PA = [rm.allPA[i][end] for i = 1:length(rm.allT)]

trap_range = findall(x->2*(lossConeAngle+0.002)>x>(lossConeAngle+0.002), final_PA)
inclusion_range = findall(x->4.002<x<5.002, initial_PA)
loss_range = findall(x->x<(lossConeAngle+0.002), final_PA)


E_prec = sort(final_E[loss_range])
E_trap = sort(final_E[trap_range])

j_prec = fit(Histogram, E_prec, logrange(ELo, EHi, Esteps+1))
j_trap = fit(Histogram, E_trap, logrange(ELo, EHi, Esteps+1))
f = j_prec.weights ./ j_trap.weights





