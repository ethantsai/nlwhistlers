include("util.jl")
include("constants.jl")
include("sim_setup.jl")
include("data_processor.jl")

# most simulation results were run with these bins
const ELo = 52;
const EHi = 1000;
const Esteps = 32; # double ELFIN E bins
const PALo = 3;
const PAHi = 15;
const PAsteps = 1300;
sim_plot_E_bins = logrange(52,1000,32) 

save_density_plot = false
save_frequency_plot = false
save_trajectory_comparison = false
visualize_frequency_models = false
save_frequency_model_comparison = false

## Import data
# Import ELFIN statistics CSVs
stats_1_csv_name = "stats_1_norm+.csv"
stats_2_csv_name = "stats_2_norm+.csv"
stats_3_csv_name = "stats_3_norm+.csv"
stats_1 =  CSV.File("external_data/stats_csvs/$stats_1_csv_name"; header=true, delim=' ', types=Float64) |> DataFrame
stats_2 =  CSV.File("external_data/stats_csvs/$stats_2_csv_name"; header=true, delim=' ', types=Float64) |> DataFrame
stats_3 =  CSV.File("external_data/stats_csvs/$stats_3_csv_name"; header=true, delim=' ', types=Float64) |> DataFrame
energy = stats_1.en
ll_dawn_lo = stats_1."L<5_day"
ll_dawn_md = stats_2."L<5_day"
ll_dawn_hi = stats_3."L<5_day"
hl_dawn_lo = stats_1."L>5_day"
hl_dawn_md = stats_2."L>5_day"
hl_dawn_hi = stats_3."L>5_day"
ll_dusk_lo = stats_1."L<5_dusk"
ll_dusk_md = stats_2."L<5_dusk"
ll_dusk_hi = stats_3."L<5_dusk"
hl_dusk_lo = stats_1."L>5_dusk"
hl_dusk_md = stats_2."L>5_dusk"
hl_dusk_hi = stats_3."L>5_dusk"
ll_nite_lo = stats_1."L<5_night"
ll_nite_md = stats_2."L<5_night"
ll_nite_hi = stats_3."L<5_night"
hl_nite_lo = stats_1."L>5_night"
hl_nite_md = stats_2."L>5_night"
hl_nite_hi = stats_3."L>5_night"

###########################
# Density Comparison Plot #
###########################
test_cases = [6.5 23.0 3 "main/results/original/" "HI_NITE_MODEL_1m" c4 "Sheeley";
              6.5 23.0 3 "main/results/plasma_density/" "HI_NITE_MODEL_omegam3" c5 "Low Density Bound";
             ]

L, MLT, Kp, dir, scenario, colour, label = test_cases[1,:]
@time @load "$dir$scenario.jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_1 = smooth(sim_ratio[1], 8, 5)
norm_1 = normalize_to_elfin(hl_nite_md, sim_ratio_sm_1)

L, MLT, Kp, dir, scenario, colour, label = test_cases[2,:]
@time @load "$dir$scenario.jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_2 = smooth(sim_ratio[1], 8, 5)
norm_2 = normalize_to_elfin(hl_nite_md, sim_ratio_sm_2)

density_comparison_plot = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)

L, MLT, Kp, dir, scenario, colour, label = test_cases[1,:]
density_comparison_plot = plot!(sim_plot_E_bins, norm_1*sim_ratio_sm_1, label="$label Model: L=$L, MLT=$MLT, Ω_pe=6.5", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
L, MLT, Kp, dir, scenario, colour, label = test_cases[2,:]
density_comparison_plot = plot!(sim_plot_E_bins, norm_2*sim_ratio_sm_2, label="$label: L=$L, MLT=$MLT, Ω_pe=3.0", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)

density_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = c1, label=false)
density_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = c1, label=false)
density_comparison_plot = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = c1, linewidth=2, markershape=:circle);
density_comparison_plot = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
density_comparison_plot = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)

if save_density_plot
    savefig(density_comparison_plot, "main/images/density_comparison.png")
    savefig(density_comparison_plot, "main/images/density_comparison.pdf")
end


#############################
# Frequency Comparison Plot #
#############################
test_cases = [6.5 23.0 3  "main/results/frequency/" "HI_NITE_MODEL_15" c2 "0.15";
              6.5 23.0 3  "main/results/frequency/" "HI_NITE_MODEL_20" c3 "0.20";
              6.5 23.0 3  "main/results/frequency/" "HI_NITE_MODEL_25" c4 "0.25";
              6.5 23.0 3  "main/results/frequency/" "HI_NITE_MODEL_30" c5 "0.30";
              6.5 23.0 3  "main/results/frequency/" "HI_NITE_MODEL_35" c6 "0.35";
              6.5 23.0 3  "main/results/frequency/" "HI_NITE_MODEL_40" c7 "0.40";
              6.5 23.0 3  "main/results/frequency/" "HI_NITE_MODEL_45" c8 "0.45";
              ]

frequency_comparison_plot = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)

# I want to normalize to omega_m = 0.35 which is what I have been using in previous plots
L, MLT, Kp, dir, scenario, colour, label = test_cases[5,:]
@time @load "$dir$scenario.jld2" rm
@info "loaded $scenario.jld2 from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 8, 5)
norm_35 = normalize_to_elfin(hl_nite_md, sim_ratio_sm)

for i in 1:length(test_cases[:,1])
    L, MLT, Kp, dir, scenario, colour, label = test_cases[i,:]
    @time @load "$dir$scenario.jld2" rm
    @info "loaded $scenario.jld2 from $dir"
    sim_ratio = prec_to_trap_ratio(rm)
    sim_ratio_sm = smooth(sim_ratio[1], 8, 5)
    norm = normalize_to_elfin(hl_nite_md, sim_ratio_sm)
    frequency_comparison_plot = plot!(sim_plot_E_bins, norm_35*sim_ratio_sm, label="ω_m = $label, L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
end

frequency_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = c1, label=false)
frequency_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = c1, label=false)
frequency_comparison_plot = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = c1, linewidth=2, markershape=:circle);
frequency_comparison_plot = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
frequency_comparison_plot = plot!(dpi = 500,size=(1000,600), margin=20px, bottom_margin=12px)

if save_frequency_plot
    savefig(frequency_comparison_plot, "main/images/freq_comparison.png")
    savefig(frequency_comparison_plot, "main/images/freq_comparison.pdf")
end


####################################################
# Trajectory comparison with different frequencies #
# This will actually run the simulation for N=500  #
# This takes ~5 min on my laptop per run. It then  #
# plots the trajectories for ω_m = 0.15 and 0.45   #
####################################################
include("agapitov_setup.jl")
include("particle_tracing.jl")
include("sim_setup.jl")
# run a sim w/ 500 particles with E=500 keV and PA=3 deg
ICrange = [500, 500, 1, 3, 3, 1]; 
wave_model_array, wave_model_coeff_array, wave_normalizer, wave_shifter_array = setup_wave_model(test_cases)
sol_15 = @time run_model(500, ICrange, 6.5, 0.15, 40, 1, wave_model_coeff_array[1], wave_shifter_array[1]);
sol_45 = @time run_model(500, ICrange, 6.5, 0.45, 40, 1, wave_model_coeff_array[1], wave_shifter_array[1]);

rm_15 = @time sol2rm(sol_15, test_cases[1,4]);
rm_45 = @time sol2rm(sol_45, test_cases[2,4]);

p1 = plot(rm_15.allT, rm_15.allE, legend=false, title = "Energies ω_m = 0.15", ylim = (492,508))
p2 = plot(rm_15.allT, rm_15.allPA, legend=false, title = "PA ω_m = 0.15", ylim = (2.85,3.15))
p3 = plot(rm_45.allT, rm_45.allE, legend=false, title = "Energies ω_m = 0.45", ylim = (492,508))
p4 = plot(rm_45.allT, rm_45.allPA, legend=false, title = "PA ω_m = 0.45", ylim = (2.85,3.15))

bigplot = plot(p1,p3,p2,p4,
        dpi = 96, layout = (2,2), size=(1000,1000),
        xtickfontsize=14, ytickfontsize=14, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, titlefontsize=16)

if save_trajectory_comparison
    savefig(bigplot, "main/images/freq_traj_comparison.png")
    savefig(bigplot, "main/images/freq_traj_comparison.pdf")
end

if visualize_frequency_models
    # let's add frequency dependence on Latitude
    function omega_m(lambda)
        if lambda < 20
            return 0.41 - 0.0125*lambda
        else
            return 0.2
        end
    end
    function omega_m_2(lambda)
        if lambda < 16.8
            return 0.41 - 0.0125*lambda
        else
            return 0.2
        end
    end
    function omega_m_3(lambda)
        if lambda < 20
            return 0.41 - 0.0125*lambda
        else
            return 0.16
        end
    end
    function omega_m_4(lambda)
        if lambda < 10
            return 0.41 - 0.025*lambda
        else
            return 0.16
        end
    end
    plot(omega_m, 0, 90, ylim = (0.1,0.5), label=["agapitov+2018 model"], ylabel = "omega_m = f_m/f_ce", xlabel = "latitude (deg)")
    plot([omega_m_2,omega_m_3,omega_m_4], 0, 90, ylim = (0.1,0.5), label=["omega_m_1" "omega_m_2" "omega_m_3"], ylabel = "omega_m = f_m/f_ce", xlabel = "latitude (deg)")
end


##############################
# Frequency model comparison #
##############################
test_cases = [6.5 23.0 3 "main/results/frequency/" "HI_NITE_MODEL_omega_mod_1" c2 "ω_mod_1";
              6.5 23.0 3 "main/results/frequency/" "HI_NITE_MODEL_omega_mod_2" c5 "ω_mod_2";
              6.5 23.0 3 "main/results/frequency/" "HI_NITE_MODEL_omega_mod_3" c7 "ω_mod_3";
              6.5 23.0 3 "main/results/frequency/" "HI_NITE_MODEL_15" c3 "ω_m = 0.15";
              6.5 23.0 3 "main/results/frequency/" "HI_NITE_MODEL_45" c4 "ω_m = 0.45";
              ]

frequency_model_comparison_plot = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)

@time @load "main/results/frequency/HI_NITE_MODEL_35.jld2" rm
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 8, 5)
norm_35 = normalize_to_elfin(hl_nite_md, sim_ratio_sm)

for i in 1:length(test_cases[:,1])
    L, MLT, Kp, dir, scenario, colour, label = test_cases[i,:]
    @time @load "$dir$scenario.jld2" rm
    @info "loaded $scenario from $dir"
    sim_ratio = prec_to_trap_ratio(rm)
    sim_ratio_sm = smooth(sim_ratio[1], 8, 5)
    norm = normalize_to_elfin(hl_nite_md, sim_ratio_sm)
    frequency_model_comparison_plot = plot!(sim_plot_E_bins, norm_35*sim_ratio_sm, label="$label: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
end

frequency_model_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = c1, label=false)
frequency_model_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = c1, label=false)
frequency_model_comparison_plot = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = c1, linewidth=2, markershape=:circle);
frequency_model_comparison_plot = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
frequency_model_comparison_plot = plot!(dpi = 500,size=(1000,600), margin=20px, bottom_margin=12px)

if save_frequency_model_comparison
    savefig(frequency_model_comparison_plot, "main/images/freq_model_comparison.png")
    savefig(frequency_model_comparison_plot, "main/images/freq_model_comparison.pdf")
end

###########################
# Oblique Comparison Plot #
###########################

test_cases = [6.5 23.0 3  "HI_NITE_MODEL_1m" c4 "Night";
            6.5 23.0 3  "HI_NITE_WNA1_5m" c5 "WNA1";
            6.5 23.0 3  "HI_NITE_WNA2_5m" c6 "WNA2";
            6.5 23.0 3  "HI_NITE_WNA3_5m" c3 "WNA3";
            ]

L, MLT, Kp, scenario, colour, label = test_cases[1,:]
@time @load "result_matrix_stats/"*scenario*".jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_1 = smooth(sim_ratio[1], 8, 5)
norm_1 = normalize_to_elfin(hl_nite_md, sim_ratio_sm_1)

L, MLT, Kp, scenario, colour, label = test_cases[2,:]
@time @load "result_matrix_oblique_4/"*scenario*".jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_2 = smooth(sim_ratio[1], 12, 6)
norm_2 = normalize_to_elfin(hl_nite_md, sim_ratio_sm_2)

L, MLT, Kp, scenario, colour, label = test_cases[3,:]
@time @load "result_matrix_oblique_4/"*scenario*".jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_3 = smooth(sim_ratio[1], 12, 6)
norm_3 = normalize_to_elfin(hl_nite_md, sim_ratio_sm_3)

L, MLT, Kp, scenario, colour, label = test_cases[4,:]
@time @load "result_matrix_oblique_4/"*scenario*".jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_4 = smooth(sim_ratio[1], 12, 6)
norm_4 = normalize_to_elfin(hl_nite_md, sim_ratio_sm_4)

nite_plot_hil = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)

L, MLT, Kp, scenario, colour, label = test_cases[1,:]
nite_plot_hil = plot!(sim_plot_E_bins, norm_1*sim_ratio_sm_1, label="label Model: L=L, MLT=MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
L, MLT, Kp, scenario, colour, label = test_cases[2,:]
nite_plot_hil = plot!(sim_plot_E_bins, norm_2*sim_ratio_sm_2, label="label Model: L=L, MLT=MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
L, MLT, Kp, scenario, colour, label = test_cases[3,:]
nite_plot_hil = plot!(sim_plot_E_bins, norm_3*sim_ratio_sm_3, label="label Model: L=L, MLT=MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
L, MLT, Kp, scenario, colour, label = test_cases[4,:]
nite_plot_hil = plot!(sim_plot_E_bins, norm_4*sim_ratio_sm_4, label="label Model: L=L, MLT=MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)

nite_plot_hil = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = c4, label=false)
nite_plot_hil = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = c4, label=false)
nite_plot_hil = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = bipride_pink, linewidth=2, markershape=:circle);
nite_plot_hil = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
nite_plot_hil = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)

if multi_oblique plot
    # try to get results to match at higher Energies
    # f3 = f1*a + f2*(1-a), w/ a<1; f1 is wna3, f2 is wna2
    f3(a) = @. norm_4*sim_ratio_sm_4*a + norm_3*sim_ratio_sm_3*(1-a)
    nite_plot_hil = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
                xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)
    L, MLT, Kp, scenario, colour, label = test_cases[3,:]
    nite_plot_hil = plot!(sim_plot_E_bins, norm_3*sim_ratio_sm_3, label="label Model: L=L, MLT=MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
    L, MLT, Kp, scenario, colour, label = test_cases[4,:]
    for i in 0.1:0.1:0.9
        nite_plot_hil = plot!(sim_plot_E_bins, f3(i), label="")
    end
    nite_plot_hil = plot!(sim_plot_E_bins, norm_4*sim_ratio_sm_4, label="label Model: L=L, MLT=MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
    nite_plot_hil = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = c4, label=false)
    nite_plot_hil = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = c4, label=false)
    nite_plot_hil = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = bipride_pink, linewidth=2, markershape=:circle);
    nite_plot_hil = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
    nite_plot_hil = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
end

if save_oblique_plot
    savefig(nite_plot_hil, "images/oblique_comparison_5m.png")
    savefig(nite_plot_hil, "images/oblique_comparison_5m.pdf")
end





save_dir = "results_ducting/"
folder = "run25/"
test_cases = [6.5 23.0 3  "HI_NITE_MODEL_1m" c4 "Night";
              6.5 23.0 3  "HI_NITE_WNA1" c5 "WNA1";
              6.5 23.0 3  "HI_NITE_WNA2" c6 "WNA2";
              6.5 23.0 3  "HI_NITE_WNA3" c3 "WNA3";
              ]


case = "HI_NITE_WNA3_3_1000000"
@time @load save_dir*folder*case*".jld2" sol;
rm = sol2rm2(sol, case);
@time @save "result_matrix_oblique_3/"*case[1:12]*".jld2" rm

i = 2
L, MLT, Kp, scenario, colour, label = test_cases[i,:]
@time @load "result_matrix_oblique_2/"*scenario*".jld2" rm
@info "loaded $scenario from $dir"
rm.allK[1]

sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_1 = smooth(sim_ratio[1], 8, 5)
norm_1 = normalize_to_elfin(hl_nite_md, sim_ratio_sm_1)

L, MLT, Kp, scenario, colour, label = test_cases[2,:]
@time @load "result_matrix_oblique/"*scenario*".jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_2 = smooth(sim_ratio[1], 8, 5)
norm_2 = normalize_to_elfin(hl_nite_md, sim_ratio_sm_2)

L, MLT, Kp, scenario, colour, label = test_cases[3,:]
@time @load "result_matrix_oblique/"*scenario*".jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_3 = smooth(sim_ratio[1], 8, 5)
norm_3 = normalize_to_elfin(hl_nite_md, sim_ratio_sm_3)

L, MLT, Kp, scenario, colour, label = test_cases[4,:]
@time @load "result_matrix_oblique/"*scenario*".jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_4 = smooth(sim_ratio[1], 8, 5)
norm_4 = normalize_to_elfin(hl_nite_md, sim_ratio_sm_4)








### MLT comparison, FIGURE 4 in JGR2023 except now with diffusion code 

# Dawn plot and generate normalizers
test_cases = [4.5 8.0  3  "hr_LO_DAWN_MODEL" c5 "Dawn-Noon";
              6.5 8.0  3  "hr_HI_DAWN_MODEL" c4 "Dawn-Noon";
              ]

dawn_plot_lol = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)

L, MLT, Kp, scenario, colour, label = test_cases[1,:]
@time @load "result_matrix_stats/"*scenario*".jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 8, 5)
norm = normalize_to_elfin(ll_dawn_md, sim_ratio_sm)
b1 = 1/(norm * maximum(sim_ratio_sm))
normalizer = b1*norm

E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("LO", "DAWN", 1)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("LO", "DAWN", 2)
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("LO", "DAWN", 3)
dc_norm = normalizer * normalize_to(80, E_bins, E, sim_ratio_sm, prec_ratio_wna1)

dawn_plot_lol = plot!(sim_plot_E_bins, normalizer .* sim_ratio_sm, label="label Model: L=L, MLT=MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
dawn_plot_lol = plot!(energy, b1*ll_dawn_md, fillrange=b1*ll_dawn_hi, fillalpha = 0.2, color = c5, label=false)
dawn_plot_lol = plot!(energy, b1*ll_dawn_md, fillrange=b1*ll_dawn_lo, fillalpha = 0.2, color = c5, label=false)
dawn_plot_lol = plot!(energy, b1*ll_dawn_md, label = "ELFIN Dawn-Noon: L<5, 4<MLT<13 ", color = c6, linewidth=2, markershape=:circle);
dawn_plot_lol = plot!(E, dc_norm * prec_ratio_wna1, label = "WNA1: Parallel Waves")
dawn_plot_lol = plot!(E, dc_norm * prec_ratio_wna2, label = "WNA2: initially parallel -> resonance cone")
dawn_plot_lol = plot!(E, dc_norm * prec_ratio_wna3, label = "WNA3: gendrin angle -> resonance cone")
dawn_plot_lol = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
dawn_plot_lol = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(dawn_plot_lol, "images/dawn_plot_lol.png")
savefig(dawn_plot_lol, "images/dawn_plot_lol.pdf")

dawn_plot_hil = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)
L, MLT, Kp, scenario, colour, label = test_cases[2,:]
@time @load "result_matrix_stats/"*scenario*".jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 8, 5)
norm = normalize_to_elfin(hl_dawn_md, sim_ratio_sm)
b1 = 1/(norm * maximum(sim_ratio_sm))
normalizer = b1*norm

E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("HI", "DAWN", 1)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("HI", "DAWN", 2)
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("HI", "DAWN", 3)
dc_norm = normalizer * normalize_to(80, E_bins, E, sim_ratio_sm, prec_ratio_wna1)

dawn_plot_hil = plot!(sim_plot_E_bins, normalizer .* sim_ratio_sm, label="label Model: L=L, MLT=MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
dawn_plot_hil = plot!(energy, b1*hl_dawn_md, fillrange=b1*hl_dawn_hi, fillalpha = 0.2, color = c4, label=false)
dawn_plot_hil = plot!(energy, b1*hl_dawn_md, fillrange=b1*hl_dawn_lo, fillalpha = 0.2, color = c4, label=false)
dawn_plot_hil = plot!(energy, b1*hl_dawn_md, label = "ELFIN Dawn-Noon: L>5, 4<MLT<13 ", color = bipride_pink, linewidth=2, markershape=:circle);
dawn_plot_hil = plot!(E, dc_norm * prec_ratio_wna1, label = "WNA1: Parallel Waves")
dawn_plot_hil = plot!(E, dc_norm * prec_ratio_wna2, label = "WNA2: initially parallel -> resonance cone")
dawn_plot_hil = plot!(E, dc_norm * prec_ratio_wna3, label = "WNA3: gendrin angle -> resonance cone")
dawn_plot_hil = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
dawn_plot_hil = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(dawn_plot_hil, "images/dawn_plot_hil.png")
savefig(dawn_plot_hil, "images/dawn_plot_hil.pdf")

# Dusk plot
test_cases = [4.5 16.5 3  "hr_LO_DUSK_MODEL" c5 "Noon-Dusk";
              6.5 16.5 3  "hr_HI_DUSK_MODEL" c4 "Noon-Dusk";
              ]


dusk_plot_lol = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)
L, MLT, Kp, scenario, colour, label = test_cases[1,:]
@time @load "result_matrix_stats/"*scenario*".jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 8, 5)
norm = normalize_to_elfin(ll_dusk_md, sim_ratio_sm)
b1 = 1/(norm * maximum(sim_ratio_sm))
normalizer = b1*norm

E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("LO", "DUSK", 1)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("LO", "DUSK", 2)
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("LO", "DUSK", 3)
dc_norm = normalizer * normalize_to(80, E_bins, E, sim_ratio_sm, prec_ratio_wna1)

dusk_plot_lol = plot!(sim_plot_E_bins, normalizer .* sim_ratio_sm, label="label Model: L=L, MLT=MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
dusk_plot_lol = plot!(energy, b1*ll_dusk_md, fillrange=b1*ll_dusk_hi, fillalpha = 0.2, color = c5, label=false)
dusk_plot_lol = plot!(energy, b1*ll_dusk_md, fillrange=b1*ll_dusk_lo, fillalpha = 0.2, color = c5, label=false)
dusk_plot_lol = plot!(energy, b1*ll_dusk_md, label = "ELFIN Noon-Dusk: L<5, 13<MLT<18", color = c6, linewidth=2, markershape=:circle);
dusk_plot_lol = plot!(E, dc_norm * prec_ratio_wna1, label = "WNA1: Parallel Waves")
dusk_plot_lol = plot!(E, dc_norm * prec_ratio_wna2, label = "WNA2: initially parallel -> resonance cone")
dusk_plot_lol = plot!(E, dc_norm * prec_ratio_wna3, label = "WNA3: gendrin angle -> resonance cone")
dusk_plot_lol = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
dusk_plot_lol = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(dusk_plot_lol, "images/dusk_compare_lol.png")
savefig(dusk_plot_lol, "images/dusk_compare_lol.pdf")

dusk_plot_hil = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)
L, MLT, Kp, scenario, colour, label = test_cases[2,:]
# @time @load "result_matrix_stats/"*scenario*".jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 5, 5)
norm = normalize_to_elfin(hl_dusk_md, sim_ratio_sm)
b1 = 1/(norm * maximum(sim_ratio_sm))
normalizer = b1*norm

E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("HI", "DUSK", 1)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("HI", "DUSK", 2)
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("HI", "DUSK", 3)
dc_norm = normalizer * normalize_to(80, E_bins, E, sim_ratio_sm, prec_ratio_wna1)

dusk_plot_hil = plot!(sim_plot_E_bins, normalizer .* sim_ratio_sm, label="label Model: L=L, MLT=MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
dusk_plot_hil = plot!(energy, b1*hl_dusk_md, fillrange=b1*hl_dusk_hi, fillalpha = 0.2, color = c4, label=false)
dusk_plot_hil = plot!(energy, b1*hl_dusk_md, fillrange=b1*hl_dusk_lo, fillalpha = 0.2, color = c4, label=false)
dusk_plot_hil = plot!(energy, b1*hl_dusk_md, label = "ELFIN Noon-Dusk: L>5, 13<MLT<18", color = bipride_pink, linewidth=2, markershape=:circle);
dusk_plot_hil = plot!(E, dc_norm * prec_ratio_wna1, label = "WNA1: Parallel Waves")
dusk_plot_hil = plot!(E, dc_norm * prec_ratio_wna2, label = "WNA2: initially parallel -> resonance cone")
dusk_plot_hil = plot!(E, dc_norm * prec_ratio_wna3, label = "WNA3: gendrin angle -> resonance cone")
dusk_plot_hil = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
dusk_plot_hil = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(dusk_plot_hil, "images/dusk_compare_hil.png")
savefig(dusk_plot_hil, "images/dusk_compare_hil.pdf")

# Night plot
test_cases = [4.5 23.0 3  "hr_LO_NITE_MODEL" c5 "Night";
              6.5 23.0 3  "HI_NITE_MODEL_1m" c4 "Night";
              ]

nite_plot_lol = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)
L, MLT, Kp, scenario, colour, label = test_cases[1,:]
@time @load "result_matrix_stats/"*scenario*".jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
norm = normalize_to_elfin(ll_nite_md, sim_ratio_sm)
b1 = 1/(norm * maximum(sim_ratio_sm))
normalizer = b1*norm

E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("LO", "NITE", 1)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("LO", "NITE", 2)
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("LO", "NITE", 3)
dc_norm = normalizer * normalize_to(80, E_bins, E, sim_ratio_sm, prec_ratio_wna1)

nite_plot_lol = plot!(sim_plot_E_bins, normalizer .* sim_ratio_sm, label="label Model: L=L, MLT=MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
nite_plot_lol = plot!(energy, b1*ll_nite_md, fillrange=b1*ll_nite_hi, fillalpha = 0.2, color = c5, label=false)
nite_plot_lol = plot!(energy, b1*ll_nite_md, fillrange=b1*ll_nite_lo, fillalpha = 0.2, color = c5, label=false)
nite_plot_lol = plot!(energy, b1*ll_nite_md, label = "ELFIN Night: L<5, 18<MLT<4 ", color = c6, linewidth=2, markershape=:circle);
nite_plot_lol = plot!(E, dc_norm * prec_ratio_wna1, label = "WNA1: Parallel Waves")
nite_plot_lol = plot!(E, dc_norm * prec_ratio_wna2, label = "WNA2: initially parallel -> resonance cone")
nite_plot_lol = plot!(E, dc_norm * prec_ratio_wna3, label = "WNA3: gendrin angle -> resonance cone")
nite_plot_lol = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
nite_plot_lol = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(nite_plot_lol, "images/nite_plot_lol.png")
savefig(nite_plot_lol, "images/nite_plot_lol.pdf")

nite_plot_hil = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)
L, MLT, Kp, scenario, colour, label = test_cases[2,:]
@time @load "result_matrix_stats/"*scenario*".jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 8, 5)
norm = normalize_to_elfin(hl_nite_md, sim_ratio_sm)
b1 = 1/(norm * maximum(sim_ratio_sm))
normalizer = b1*norm

E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("HI", "NITE", 1, true)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("HI", "NITE", 2, true)
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("HI", "NITE", 3, true)
dc_norm = normalizer * normalize_to(80, E_bins, E, sim_ratio_sm, prec_ratio_wna1)

nite_plot_hil = plot!(sim_plot_E_bins, normalizer .* sim_ratio_sm, label="label Model: L=L, MLT=MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
nite_plot_hil = plot!(energy, b1*hl_nite_md, fillrange=b1*hl_nite_hi, fillalpha = 0.2, color = c4, label=false)
nite_plot_hil = plot!(energy, b1*hl_nite_md, fillrange=b1*hl_nite_lo, fillalpha = 0.2, color = c4, label=false)
nite_plot_hil = plot!(energy, b1*hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = bipride_pink, linewidth=2, markershape=:circle);
nite_plot_hil = plot!(E, dc_norm * prec_ratio_wna1, label = "WNA1: Parallel Waves")
nite_plot_hil = plot!(E, dc_norm * prec_ratio_wna2, label = "WNA2: initially parallel -> resonance cone")
nite_plot_hil = plot!(E, dc_norm * prec_ratio_wna3, label = "WNA3: gendrin angle -> resonance cone")
nite_plot_hil = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
nite_plot_hil = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(nite_plot_hil, "images/nite_plot_hil.png")
savefig(nite_plot_hil, "images/nite_plot_hil.pdf")

statistical_comparison_plot = plot(dawn_plot_lol, dawn_plot_hil, dusk_plot_lol, dusk_plot_hil, nite_plot_lol, nite_plot_hil, layout=(3,2),size=(1200,1600),dpi=500)
savefig(statistical_comparison_plot, "images/statistical_comparison.png")
savefig(statistical_comparison_plot, "images/statistical_comparison.pdf")



#= 
Figure 3
Precipitating-to-trapped electron flux ratio j_prec/j_trap
measured by ELFIN at L>5 on the night side (18-4 MLT) as a
function of electron energy E (black curve). The corresponding
ratio j_prec/j_trap obtained from test particle simulations
is displayed for parallel (WNA1 model) lower-band chorus waves,
using frequency Model 1 of constant frequency chorus waves
and either a typical Ω_pe/Ω_ce=6.5 at L=6.5 and 23 MLT
(solid red curve). The ratios j_prec/j_trap obtained from
the quasi-linear diffusion FDC code for the same parameters,
and from test particle simulations in the case of a reduced
density Ω_pe/Ω_ce=3, are also shown (dashed red and solid
magenta curves, respectively). All simulation results are
normalized to observations at 90 keV.
=#



#=
Figure 4
Precipitating-to-trapped electron flux ratio j_prec/j_trap
measured by ELFIN at L>5 on the night side (18-4 MLT) as a
function of electron energy E (black). The corresponding
ratios j_prec/j_trap obtained from test particle simulations
(solid curves) and from the quasi-linear diffusion code
(dashed curves) are displayed for lower-band chorus waves,
using frequency Model 1 of constant frequency, and four
wave normal angle models: WNA1 (red), WNA2 (blue), WNA3
(green), and WNA4 (purple), with a normalization to
observations at 90 keV, adopting a typical Ω_pe/Ω_ce=6.5
at L=6.5 and 23 MLT.
=#



#=
Figure 5
Precipitating-to-trapped electron flux ratio j_prec/j_trap
measured by ELFIN at L>5 on the night side (18-4 MLT) as a
function of electron energy E (black). The corresponding
ratio j_prec/j_trap obtained from test particle simulations
is displayed for parallel lower-band chorus waves (WNA1
model), using frequency Models 1 and 2 of constant frequency
chorus waves (red and orange, respectively) and frequency
Model 3 of realistic chorus waves of decreasing frequency
toward higher latitudes (blue), with a normalization to
observations at 90 keV, adopting a typical Ω_pe/Ω_ce=6.5
at L=6.5 and 23 MLT.
=#