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


############
scenario = 1
###########^
num_omega_m = length(omega_m_cases)

index_array = ( 3*(scenario-1) + 1 ) : ( 3*(scenario-1) + num_omega_m )
if maximum(index_array) > length(L_array)*num_omega_m
    @error "There's only $(length(L_array)) scenarios; pick a valid ID"
end
[prec_to_trap_ratio(x) for x in rm_array




#     @info "Extracting data..."
#     @time allT, allZ, allPZ, allE, allPA = extract(sol);
#     @time tVec, Zmatrix, PZmatrix, PAmatrix, Ematrix = postProcessor(allT, allZ, allPZ, allPA, allE);

#     @info "Saving to $savename"
#     @time event220112T0226_3 = Resultant_Matrix(test_cases[:,end][case_index], length(sol), tVec[end], allZ, allPZ, allT, allPA, allE,countLostParticles(allT, tVec[end]), tVec, Zmatrix, PZmatrix, Ematrix, PAmatrix);    
    

    




# @time event220112T0226_3 = Resultant_Matrix("ducttest", length(sol), tVec[end], allZ, allPZ, allT, allPA, allE,countLostParticles(allT, tVec[end]), tVec, Zmatrix, PZmatrix, Ematrix, PAmatrix);
# @save "220112T0226_3.jld2" event220112T0226_3


# # at each energy range, fine j_precip/j_trapped
# # obtain all lost particles and energy at which they're lost
# # then obtain all trapped particles and their respective energies (where LC < PA < 2LC)

# @time test = Resultant_Matrix("ducttest", length(sol), tVec[end], allZ, allPZ, allT, allPA, allE,countLostParticles(allT, tVec[end]), tVec, Zmatrix, PZmatrix, Ematrix, PAmatrix);


# @time oneenergy = Resultant_Matrix("ducttest", length(sol), tVec[end], allZ, allPZ, allT, allPA, allE,countLostParticles(allT, tVec[end]), tVec, Zmatrix, PZmatrix, Ematrix, PAmatrix);
# @save "oneenergy.jld2" oneenergy

# @time oneenergyhr = Resultant_Matrix("ducttest", length(sol), tVec[end], allZ, allPZ, allT, allPA, allE,countLostParticles(allT, tVec[end]), tVec, Zmatrix, PZmatrix, Ematrix, PAmatrix);
# @save "oneenergyhr.jld2" oneenergyhr

# @load "220112T0226_2.jld2" event220112T0226_2
# @load "220112T0226_3.jld2" event220112T0226_3
# @load "220112T0226_4.jld2" event220112T0226_4




# fobs2, j_prec2, j_trap2 = prec_to_trap_ratio(event220112T0226_2)
# fobs3, j_prec3, j_trap3 = prec_to_trap_ratio(event220112T0226_3)
# fobs4, j_prec4, j_trap4 = prec_to_trap_ratio(event220112T0226_4)

# plot(E_bins, fobs2, xscale=:log10, label="omega_m=0.2")
# plot!(E_bins, fobs3, xscale=:log10, label="omega_m=0.3")
# plot!(E_bins, fobs4, xscale=:log10, label="omega_m=0.4")
