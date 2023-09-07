include("agapitovHelpers.jl")
include("agapitovmodel.jl")

# case specific
             #L   MLT  Kp name
test_cases = [7.1 8.4  2  "ELB_SA_210106T1154"; # ELB SA 01/06 11:54
              6.5 19.8 0  "ELB_ND_210108T0646"; # ELB ND 01/08 06:46
              4.8 19.0 3  "ELA_SD_210111T1750"; # ELA SD 01/11 17:50
              6   8.4  3  "ELA_NA_210112T0226"; # ELA NA 01/12 02:26
              6.5 3.8  3  "ELA_ND_200904T0112"; # ELA ND 09/04 01:12
              4.8 2.6  4  "ELB_ND_200926T0101"; # ELB ND 09/26 01:01
              5.1 20.2 3  "ELA_SD_210203T0902"; # ELA SD 02/03 09:02
              6.6 19.3 3  "ELA_SD_210203T1342"; # ELA SD 02/03 13:42
              6.2 5.8  3  "ELA_NA_210203T2046"; # ELA NA 02/03 20:46
              7.1 10.2 4  "ELA_SD_211001T0501"; # ELA SD 10/01 05:01
              6.1 9.0  3  "ELA_SD_211001T0810"; # ELA SD 10/01 08:10
              6.1 13.1 3  "ELA_SA_211101T0424"; # ELA SA 11/01 04:24
              4.5 20.8 3  "ELA_SA_211102T2218"] # ELA SA 11/02 22:18


# plot all wave models
wave_model_array, wave_model_normalizer_array, wave_model_coeff_array = setup_wave_model(test_cases)
# i=1
# plot1 = plot(0:.01:90,tanh.(0:.01:90)*wave_model_normalizer_array[i].*wave_model_array[i].(0:.01:90), label = test_cases[:,end][i], linewidth=2)
i=1
plot(0:.01:90,tanh.(0:.01:90)*wave_model_normalizer_array[i].*wave_model_array[i].(0:.01:90), label = test_cases[:,end][i], linewidth=2)
i=2
plot!(0:.01:90,tanh.(0:.01:90)*wave_model_normalizer_array[i].*wave_model_array[i].(0:.01:90), label = test_cases[:,end][i], linewidth=2)
i=3
plot!(0:.01:90,tanh.(0:.01:90)*wave_model_normalizer_array[i].*wave_model_array[i].(0:.01:90), label = test_cases[:,end][i], linewidth=2, linestyle=:dash)
i=4
plot!(0:.01:90,tanh.(0:.01:90)*wave_model_normalizer_array[i].*wave_model_array[i].(0:.01:90), label = test_cases[:,end][i], linewidth=2, linestyle=:dash)
plot1 = plot!(ylabel = "B_w Multiplier", xlabel = "Latitude (degs)", title = "Synthetic Lower Band Chorus Wave Model", dpi = 500,size =(800,450), margin=3mm, bottom_margin=4mm)
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
  
    multiplier = 2

    trap_range = findall(x->multiplier*(lossConeAngle)>x>(lossConeAngle), initial_PA)
    loss_range = findall(x->x<(lossConeAngle), final_PA)
    E_prec = sort(final_E[loss_range])
    E_trap = sort(initial_E[trap_range])
  
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
scenario = 2
###########^
num_omega_m = length(omega_m_cases)

index_array = ( 3*(scenario-1) + 1 ) : ( 3*(scenario-1) + num_omega_m )
if maximum(index_array) > length(L_array)*length(index_array)
    @error "There's only $(length(L_array)) scenarios; pick a valid ID"
end
results_array = [prec_to_trap_ratio(x) for x in rm_array]


i = 1
plot2 = plot(E_bins, (1/results_array[i][1][7])*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(100,2000), ylim=(1e-2, 2))
i = 2
plot!(E_bins, (1/results_array[i][1][7])*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(100,2000), ylim=(1e-2, 2))
i = 3
plot!(E_bins, (1/results_array[i][1][7])*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(100,2000), ylim=(1e-2, 2))
i = 4
plot!(E_bins, (1/results_array[i][1][7])*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(100,2000), ylim=(1e-2, 2))
i = 5
plot!(E_bins, (1/results_array[i][1][7])*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(100,2000), ylim=(1e-2, 2))

normalizer = 1/results_array[2][1][7] # normalize to dayside precip @ 100 keV to 100%
i = 1
plot2 = plot(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(100,2000), ylim=(1e-2, 1))
i = 2
plot!(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(100,2000), ylim=(1e-2, 1))
i = 3
plot!(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(100,2000), ylim=(1e-2, 1))
i = 4
plot!(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(100,2000), ylim=(1e-2, 1))
# i = 5
# plot!(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(100,2000), ylim=(1e-2, 1))


i = 1
normalizer = 1/results_array[i][1][7] # normalize to dayside precip @ 100 keV to 100%
plot2 = plot(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(100,2000), ylim=(1e-2, 1))
i = 2
normalizer = 1/results_array[i][1][7] # normalize to dayside precip @ 100 keV to 100%
plot!(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(100,2000), ylim=(1e-2, 1))
i = 3
normalizer = 1/results_array[i][1][7] # normalize to dayside precip @ 100 keV to 100%
plot!(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(100,2000), ylim=(1e-2, 1))
i = 4
normalizer = 1/results_array[i][1][7] # normalize to dayside precip @ 100 keV to 100%
plot!(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(100,2000), ylim=(1e-2, 1))
# i = 5
# plot!(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(100,2000), ylim=(1e-2, 1))



i = 1
normalizer = 1/results_array[2][1][7] # normalize to dayside precip @ 100 keV to 100%
@info normalizer
plot2 = plot(E_bins, results_array[i][2], label=rm_array[i].label, xscale=:log10, xlim=(100,2000))
i = 2
normalizer = 1/results_array[i][1][7] # normalize to dayside precip @ 100 keV to 100%
@info normalizer
plot!(E_bins, results_array[i][2], label=rm_array[i].label, xscale=:log10, xlim=(100,2000))
i = 3
normalizer = 1/results_array[i][1][7] # normalize to dayside precip @ 100 keV to 100%
@info normalizer
plot!(E_bins, results_array[i][2], label=rm_array[i].label, xscale=:log10, xlim=(100,2000))
i = 4
normalizer = 1/results_array[i][1][7] # normalize to dayside precip @ 100 keV to 100%
plot!(E_bins, results_array[i][2], label=rm_array[i].label, xscale=:log10, xlim=(100,2000))
i = 5
normalizer = 1/results_array[i][1][7] # normalize to dayside precip @ 100 keV to 100%
plot!(E_bins, results_array[i][2], label=rm_array[i].label, xscale=:log10, xlim=(100,2000))



i = 1
normalizer = 1/results_array[i][1][7] # normalize to dayside precip @ 100 keV to 100%
@info normalizer
plot2 = plot(E_bins, results_array[i][3], label=rm_array[i].label, xscale=:log10, xlim=(100,2000))
i = 2
normalizer = 1/results_array[i][1][7] # normalize to dayside precip @ 100 keV to 100%
@info normalizer
plot!(E_bins, results_array[i][3], label=rm_array[i].label, xscale=:log10, xlim=(100,2000))
i = 3
normalizer = 1/results_array[i][1][7] # normalize to dayside precip @ 100 keV to 100%
@info normalizer
plot!(E_bins, results_array[i][3], label=rm_array[i].label, xscale=:log10, xlim=(100,2000))
i = 4
normalizer = 1/results_array[i][1][7] # normalize to dayside precip @ 100 keV to 100%
plot!(E_bins, results_array[i][3], label=rm_array[i].label, xscale=:log10, xlim=(100,2000))
i = 5
normalizer = 1/results_array[i][1][7] # normalize to dayside precip @ 100 keV to 100%
plot!(E_bins, results_array[i][3], label=rm_array[i].label, xscale=:log10, xlim=(100,2000))


savefig(plot2, "images/ratios.png")


rm = rm_array[3]

initial_E = [rm2.allE[i][1] for i = 1:length(rm2.allT)]
final_E = [rm2.allE[i][end] for i = 1:length(rm2.allT)]
initial_PA = [rm2.allPA[i][1] for i = 1:length(rm2.allT)]
final_PA = [rm2.allPA[i][end] for i = 1:length(rm2.allT)]

trap_range = findall(x->2*(lossConeAngle+0.002)>x>(lossConeAngle+0.002), final_PA)
findall(x->(x>1999), initial_E[trap_range])


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





plot(E_bins, results_array[1][2])
plot!(E_bins, results_array[2][2])
plot!(E_bins, results_array[3][2])
plot!(E_bins, results_array[4][2])
plot!(E_bins, results_array[5][2])

plot(E_bins, results_array[1][3])
plot!(E_bins, results_array[2][3])
plot!(E_bins, results_array[3][3])
plot!(E_bins, results_array[4][3])
plot!(E_bins, results_array[5][3])





# for case #1, omega_m = 0.3, L = 5.1, mlt = 21.7
case = 1
rm_case = rm_array[case]
# first, obtain what B_w looks like:
full_length = length(rm_case.tVec)
time = rm_case.tVec
i=1
while full_length != length(rm_case.allBw[i])
    i+=1
end
b_w = rm_case.allBw[i]
lambda = rad2deg.(rm_case.allLambda[i])
plot_bw = plot(time, b_w, xlabel="time", ylabel="B_w", title="Wave amplitude as function of time", legend=false)
plot_bw = plot(lambda, b_w, xlabel="Latitude (deg)", ylabel="B_w", title="Wave amplitude as function of Latitude", legend=false)

# next, isolate just those with starting energy of 63 and 138 keV
initial_E = [rm_case.allE[i][1] for i = 1:length(rm_case.allT)]
final_E = [rm_case.allE[i][end] for i = 1:length(rm_case.allT)]
initial_PA = [rm_case.allPA[i][1] for i = 1:length(rm_case.allT)]
final_PA = [rm_case.allPA[i][end] for i = 1:length(rm_case.allT)]
prec_range = findall(x->4<=x<=4.1, final_PA) 
E63_range = findall(x->62<x<64, initial_E)
E145_range = findall(x->144<x<147, initial_E)

plot1 = scatter(xlabel="Initial Pitch Angle (deg)", ylabel="Final Pitch Angle (deg)", title="L = 5.1, mlt = 21.7", legend=:topleft)
plot1 = scatter!([initial_PA[E63_range], initial_PA[E145_range]], [final_PA[E63_range], final_PA[E145_range]], label = ["63 keV" "145 keV"], xlim=[0,30], ylim=[0,90], ms=1, msw=0)
plot1 = plot!(0:30, 0:30, label=false, ls=:dash, lc=:red)


# for case #1, omega_m = 0.3, L = 7.1, mlt = 8.4
case = 2
rm_case = rm_array[case]
# first, obtain what B_w looks like:
full_length = length(rm_case.tVec)
time = rm_case.tVec
i=1
while full_length != length(rm_case.allBw[i])
    i+=1
end
b_w = rm_case.allBw[i]
lambda = rad2deg.(rm_case.allLambda[i])
plot_bw = plot(time, b_w, xlabel="time", ylabel="B_w", title="Wave amplitude as function of time", legend=false)
plot_bw = plot(lambda, b_w, xlabel="Latitude (deg)", ylabel="B_w", title="Wave amplitude as function of Latitude", legend=false)

# next, isolate just those with starting energy of 63 and 138 keV
initial_E = [rm_case.allE[i][1] for i = 1:length(rm_case.allT)]
final_E = [rm_case.allE[i][end] for i = 1:length(rm_case.allT)]
initial_PA = [rm_case.allPA[i][1] for i = 1:length(rm_case.allT)]
final_PA = [rm_case.allPA[i][end] for i = 1:length(rm_case.allT)]
prec_range = findall(x->4<=x<=4.1, final_PA) 
E63_range = findall(x->62<x<64, initial_E)
E145_range = findall(x->144<x<147, initial_E)

plot2 = scatter(xlabel="Initial Pitch Angle (deg)", ylabel="Final Pitch Angle (deg)", title="L = 7.1, mlt = 8.4", legend=:topleft)
plot2 = scatter!([initial_PA[E63_range], initial_PA[E145_range]], [final_PA[E63_range], final_PA[E145_range]], label = ["63 keV" "145 keV"], xlim=[0,30], ylim=[0,90], ms=1, msw=0)
plot2 = plot!(0:30, 0:30, label=false, ls=:dash, lc=:red)

savefig(plot1, "images/resonance_analysis_1.png")
savefig(plot2, "images/resonance_analysis_2.png")
