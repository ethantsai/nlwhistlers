include("agapitovHelpers.jl")
include("agapitovmodel.jl")
include("agapitovPlotHelpers.jl")


###############
# Stats Plots #
###############
# Import CSV
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

# Oleksiy plot
test_cases = [6.5 23.0 3  "hr_HI_NITE_MODEL" c4 "Night";
              6.5 23.0 3  "oleksiy_HI_NITE_MODEL" c5 "Oleksiy";
              ]

energy = stats_1.en
hl_nite_lo = stats_1."L>5_night"
hl_nite_md = stats_2."L>5_night"
hl_nite_hi = stats_3."L>5_night"

L, MLT, Kp, scenario, colour, label = test_cases[1,:]
@time @load "result_matrix_stats/"*scenario*".jld2" rm
@info "loaded $scenario.jld2"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_1 = smooth(sim_ratio[1], 8, 5)
norm_1 = normalize_to_elfin(hl_nite_md, sim_ratio_sm_1)


L, MLT, Kp, scenario, colour, label = test_cases[2,:]
@time @load "result_matrix_oblique/"*scenario*".jld2" rm
@info "loaded $scenario.jld2"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_2 = smooth(sim_ratio[1], 8, 5)
norm_2 = normalize_to_elfin(hl_nite_md, sim_ratio_sm_2)


nite_plot_hil = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)


L, MLT, Kp, scenario, colour, label = test_cases[1,:]
nite_plot_hil = plot!(E_bins, norm_1*sim_ratio_sm_1, label="$label Model: L=$L, MLT=$MLT, omega_pe=6.5", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
L, MLT, Kp, scenario, colour, label = test_cases[2,:]
nite_plot_hil = plot!(E_bins, norm_2*sim_ratio_sm_2, label="$label Model: L=$L, MLT=$MLT, omega_pe=3", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)

nite_plot_hil = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = c4, label=false)
nite_plot_hil = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = c4, label=false)
nite_plot_hil = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = bipride_pink, linewidth=2, markershape=:circle);
nite_plot_hil = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
nite_plot_hil = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)




savefig(statistical_comparison_plot, "images/oleksiy_comparison.png")
savefig(statistical_comparison_plot, "images/oleksiy_comparison.pdf")




# Oblique plot
test_cases = [6.5 23.0 3  "hr_HI_NITE_MODEL" c4 "Night";
              6.5 23.0 3  "HI_NITE_WNA1_5m" c5 "WNA1";
              6.5 23.0 3  "HI_NITE_WNA2_5m" c6 "WNA2";
              6.5 23.0 3  "HI_NITE_WNA3_5m" c3 "WNA3";
              ]

energy = stats_1.en
hl_nite_lo = stats_1."L>5_night"
hl_nite_md = stats_2."L>5_night"
hl_nite_hi = stats_3."L>5_night"

L, MLT, Kp, scenario, colour, label = test_cases[1,:]
@time @load "result_matrix_stats/"*scenario*".jld2" rm
@info "loaded $scenario.jld2"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_1 = smooth(sim_ratio[1], 8, 5)
norm_1 = normalize_to_elfin(hl_nite_md, sim_ratio_sm_1)

L, MLT, Kp, scenario, colour, label = test_cases[2,:]
@time @load "result_matrix_oblique_4/"*scenario*".jld2" rm
@info "loaded $scenario.jld2"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_2 = smooth(sim_ratio[1], 12, 6)
norm_2 = normalize_to_elfin(hl_nite_md, sim_ratio_sm_2)

L, MLT, Kp, scenario, colour, label = test_cases[3,:]
@time @load "result_matrix_oblique_4/"*scenario*".jld2" rm
@info "loaded $scenario.jld2"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_3 = smooth(sim_ratio[1], 12, 6)
norm_3 = normalize_to_elfin(hl_nite_md, sim_ratio_sm_3)

L, MLT, Kp, scenario, colour, label = test_cases[4,:]
@time @load "result_matrix_oblique_4/"*scenario*".jld2" rm
@info "loaded $scenario.jld2"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_4 = smooth(sim_ratio[1], 12, 6)
norm_4 = normalize_to_elfin(hl_nite_md, sim_ratio_sm_4)

nite_plot_hil = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)

L, MLT, Kp, scenario, colour, label = test_cases[1,:]
nite_plot_hil = plot!(E_bins, norm_1*sim_ratio_sm_1, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
L, MLT, Kp, scenario, colour, label = test_cases[2,:]
nite_plot_hil = plot!(E_bins, norm_2*sim_ratio_sm_2, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
L, MLT, Kp, scenario, colour, label = test_cases[3,:]
nite_plot_hil = plot!(E_bins, norm_3*sim_ratio_sm_3, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
L, MLT, Kp, scenario, colour, label = test_cases[4,:]
nite_plot_hil = plot!(E_bins, norm_4*sim_ratio_sm_4, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)

nite_plot_hil = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = c4, label=false)
nite_plot_hil = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = c4, label=false)
nite_plot_hil = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = bipride_pink, linewidth=2, markershape=:circle);
nite_plot_hil = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
nite_plot_hil = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)

savefig(nite_plot_hil, "images/oblique_comparison_5m.png")
savefig(nite_plot_hil, "images/oblique_comparison_5m.pdf")

# try to get results to match at higher Energies
# f3 = f1*a + f2*(1-a), w/ a<1; f1 is wna3, f2 is wna2
f3(a) = @. norm_4*sim_ratio_sm_4*a + norm_3*sim_ratio_sm_3*(1-a)
nite_plot_hil = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)
L, MLT, Kp, scenario, colour, label = test_cases[3,:]
nite_plot_hil = plot!(E_bins, norm_3*sim_ratio_sm_3, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
L, MLT, Kp, scenario, colour, label = test_cases[4,:]
for i in 0.1:0.1:0.9
    nite_plot_hil = plot!(E_bins, f3(i), label="")
end
nite_plot_hil = plot!(E_bins, norm_4*sim_ratio_sm_4, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
nite_plot_hil = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = c4, label=false)
nite_plot_hil = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = c4, label=false)
nite_plot_hil = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = bipride_pink, linewidth=2, markershape=:circle);
nite_plot_hil = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
nite_plot_hil = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)


save_dir = "results_ducting/"
folder = "run25/"
test_cases = [6.5 23.0 3  "hr_HI_NITE_MODEL" c4 "Night";
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
@info "loaded $scenario.jld2"
rm.allK[1]

sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_1 = smooth(sim_ratio[1], 8, 5)
norm_1 = normalize_to_elfin(hl_nite_md, sim_ratio_sm_1)

L, MLT, Kp, scenario, colour, label = test_cases[2,:]
@time @load "result_matrix_oblique/"*scenario*".jld2" rm
@info "loaded $scenario.jld2"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_2 = smooth(sim_ratio[1], 8, 5)
norm_2 = normalize_to_elfin(hl_nite_md, sim_ratio_sm_2)

L, MLT, Kp, scenario, colour, label = test_cases[3,:]
@time @load "result_matrix_oblique/"*scenario*".jld2" rm
@info "loaded $scenario.jld2"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_3 = smooth(sim_ratio[1], 8, 5)
norm_3 = normalize_to_elfin(hl_nite_md, sim_ratio_sm_3)

L, MLT, Kp, scenario, colour, label = test_cases[4,:]
@time @load "result_matrix_oblique/"*scenario*".jld2" rm
@info "loaded $scenario.jld2"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_4 = smooth(sim_ratio[1], 8, 5)
norm_4 = normalize_to_elfin(hl_nite_md, sim_ratio_sm_4)







# frequency comparison plot
test_cases = [6.5 23.0 3  "HI_NITE_MODEL_15" c0 "0.15";
              6.5 23.0 3  "HI_NITE_MODEL_20" c1 "0.20";
              6.5 23.0 3  "HI_NITE_MODEL_25" c2 "0.25";
              6.5 23.0 3  "HI_NITE_MODEL_30" c3 "0.30";
              6.5 23.0 3  "HI_NITE_MODEL_35" c5 "0.35";
              6.5 23.0 3  "HI_NITE_MODEL_40" c6 "0.40";
              6.5 23.0 3  "HI_NITE_MODEL_45" c7 "0.45";
              ]

energy = stats_1.en
hl_nite_lo = stats_1."L>5_night"
hl_nite_md = stats_2."L>5_night"
hl_nite_hi = stats_3."L>5_night"

frequency_comparison_plot = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)

L, MLT, Kp, scenario, colour, label = test_cases[5,:]
@time @load "result_matrix_freqs/"*scenario*".jld2" rm
@info "loaded $scenario.jld2"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 8, 5)
norm_35 = normalize_to_elfin(hl_nite_md, sim_ratio_sm)


for i in 1:length(test_cases[:,1])
    L, MLT, Kp, scenario, colour, label = test_cases[i,:]
    @time @load "result_matrix_freqs/"*scenario*".jld2" rm
    @info "loaded $scenario.jld2"
    sim_ratio = prec_to_trap_ratio(rm)
    sim_ratio_sm = smooth(sim_ratio[1], 8, 5)
    norm = normalize_to_elfin(hl_nite_md, sim_ratio_sm)
    frequency_comparison_plot = plot!(E_bins, norm_35*sim_ratio_sm, label="$label Model: L=$L, MLT=$MLT, omega_pe=6.5", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
end

frequency_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = c4, label=false)
frequency_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = c4, label=false)
frequency_comparison_plot = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = bipride_pink, linewidth=2, markershape=:circle);
frequency_comparison_plot = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
frequency_comparison_plot = plot!(dpi = 500,size=(1000,600), margin=20px, bottom_margin=12px)


savefig(frequency_comparison_plot, "images/freq_comparison.png")
savefig(frequency_comparison_plot, "images/freq_comparison.pdf")



# plot a few trajectories for omega_m = 0.15 and 0.45

test_cases = [6.5 23.0 3  "HI_NITE_MODEL_small_15" c0 "0.15";
              6.5 23.0 3  "HI_NITE_MODEL_small_45" c7 "0.45";
              ]

const ELo = 500;
const EHi = 500;
const Esteps = 1;
const PALo = 3;
const PAHi = 3;
const PAsteps = 1; 
const factor = 1; 
# num particles in highest energy bin = factor * num particles in lowest energy bin
ICrange = [ELo, EHi, Esteps, PALo, PAHi, PAsteps];

wave_model_array, wave_model_coeff_array, wave_normalizer, wave_shifter_array = setup_wave_model(test_cases)

sol_15 = run_model(500, ICrange, 6.5, 1, 0.15, wave_model_coeff_array[1], wave_shifter_array[1],40);
sol_45 = run_model(500, ICrange, 6.5, 1, 0.45, wave_model_coeff_array[1], wave_shifter_array[1],40);

rm_15 = sol2rm(sol_15, test_cases[1,4]);
rm_45 = sol2rm(sol_45, test_cases[2,4]);


p1 = plot(rm_15.allT, rm_15.allE, legend=false, title = "Energies omega_m = 0.15", ylim = (499.8,500.2))
p2 = plot(rm_15.allT, rm_15.allPA, legend=false, title = "PA omega_m = 0.15")
p3 = plot(rm_45.allT, rm_45.allE, legend=false, title = "Energies omega_m = 0.45", ylim = (499.8,500.2))
p4 = plot(rm_45.allT, rm_45.allPA, legend=false, title = "PA omega_m = 0.45")

bigplot = plot(p1,p3,p2,p4,
        dpi = 96, layout = (2,2), size=(1000,1000),
        xtickfontsize=14, ytickfontsize=14, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, titlefontsize=16)

savefig(bigplot, "images/freq_traj_comparison.png")




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

# now to import the actual runs

test_cases = [6.5 23.0 3  "HI_NITE_MODEL_omega_mod_1" c5 "omega_mod_1";
              6.5 23.0 3  "HI_NITE_MODEL_omega_mod_2" c6 "omega_mod_2";
              6.5 23.0 3  "HI_NITE_MODEL_omega_mod_3" c6 "omega_mod_3";
              6.5 23.0 3  "HI_NITE_MODEL_15" c0 "omega_m = 0.15";
              6.5 23.0 3  "HI_NITE_MODEL_45" c8 "omega_m = 0.45";
              ]


energy = stats_1.en
hl_nite_lo = stats_1."L>5_night"
hl_nite_md = stats_2."L>5_night"
hl_nite_hi = stats_3."L>5_night"

omega_m_comparison_plot = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)

@time @load "result_matrix_freqs/HI_NITE_MODEL_35.jld2" rm
@info "loaded $scenario.jld2"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 8, 5)
norm_35 = normalize_to_elfin(hl_nite_md, sim_ratio_sm)

for i in 1:length(test_cases[:,1])
    L, MLT, Kp, scenario, colour, label = test_cases[i,:]
    @time @load "result_matrix_freqs/"*scenario*".jld2" rm
    @info "loaded $scenario.jld2"
    sim_ratio = prec_to_trap_ratio(rm)
    sim_ratio_sm = smooth(sim_ratio[1], 8, 5)
    norm = normalize_to_elfin(hl_nite_md, sim_ratio_sm)
    omega_m_comparison_plot = plot!(E_bins, norm_3*sim_ratio_sm, label="$label Model: L=$L, MLT=$MLT, omega_pe=6.5", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
end

omega_m_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = c4, label=false)
omega_m_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = c4, label=false)
omega_m_comparison_plot = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = bipride_pink, linewidth=2, markershape=:circle);
omega_m_comparison_plot = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
omega_m_comparison_plot = plot!(dpi = 500,size=(1000,600), margin=20px, bottom_margin=12px)

savefig(omega_m_comparison_plot, "images/freq_model_comparison.png")
savefig(omega_m_comparison_plot, "images/freq_model_comparison.pdf")


### MLT comparison, FIGURE 4 in JGR2023 except now with diffusion code 

# Dawn plot and generate normalizers
test_cases = [4.5 8.0  3  "hr_LO_DAWN_MODEL" c5 "Dawn-Noon";
              6.5 8.0  3  "hr_HI_DAWN_MODEL" c4 "Dawn-Noon";
              ]
              
energy = stats_1.en
ll_dawn_lo = stats_1."L<5_day"
ll_dawn_md = stats_2."L<5_day"
ll_dawn_hi = stats_3."L<5_day"
hl_dawn_lo = stats_1."L>5_day"
hl_dawn_md = stats_2."L>5_day"
hl_dawn_hi = stats_3."L>5_day"

dawn_plot_lol = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)

L, MLT, Kp, scenario, colour, label = test_cases[1,:]
@time @load "result_matrix_stats/"*scenario*".jld2" rm
@info "loaded $scenario.jld2"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 8, 5)
norm = normalize_to_elfin(ll_dawn_md, sim_ratio_sm)
b1 = 1/(norm * maximum(sim_ratio_sm))
normalizer = b1*norm

E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("LO", "DAWN", 1)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("LO", "DAWN", 2)
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("LO", "DAWN", 3)
dc_norm = normalizer * normalize_to(80, E_bins, E, sim_ratio_sm, prec_ratio_wna1)

dawn_plot_lol = plot!(E_bins, normalizer .* sim_ratio_sm, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
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
@info "loaded $scenario.jld2"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 8, 5)
norm = normalize_to_elfin(hl_dawn_md, sim_ratio_sm)
b1 = 1/(norm * maximum(sim_ratio_sm))
normalizer = b1*norm

E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("HI", "DAWN", 1)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("HI", "DAWN", 2)
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("HI", "DAWN", 3)
dc_norm = normalizer * normalize_to(80, E_bins, E, sim_ratio_sm, prec_ratio_wna1)

dawn_plot_hil = plot!(E_bins, normalizer .* sim_ratio_sm, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
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

energy = stats_1.en
ll_dusk_lo = stats_1."L<5_dusk"
ll_dusk_md = stats_2."L<5_dusk"
ll_dusk_hi = stats_3."L<5_dusk"
hl_dusk_lo = stats_1."L>5_dusk"
hl_dusk_md = stats_2."L>5_dusk"
hl_dusk_hi = stats_3."L>5_dusk"

dusk_plot_lol = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)
L, MLT, Kp, scenario, colour, label = test_cases[1,:]
@time @load "result_matrix_stats/"*scenario*".jld2" rm
@info "loaded $scenario.jld2"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 8, 5)
norm = normalize_to_elfin(ll_dusk_md, sim_ratio_sm)
b1 = 1/(norm * maximum(sim_ratio_sm))
normalizer = b1*norm

E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("LO", "DUSK", 1)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("LO", "DUSK", 2)
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("LO", "DUSK", 3)
dc_norm = normalizer * normalize_to(80, E_bins, E, sim_ratio_sm, prec_ratio_wna1)

dusk_plot_lol = plot!(E_bins, normalizer .* sim_ratio_sm, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
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
@info "loaded $scenario.jld2"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 5, 5)
norm = normalize_to_elfin(hl_dusk_md, sim_ratio_sm)
b1 = 1/(norm * maximum(sim_ratio_sm))
normalizer = b1*norm

E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("HI", "DUSK", 1)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("HI", "DUSK", 2)
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("HI", "DUSK", 3)
dc_norm = normalizer * normalize_to(80, E_bins, E, sim_ratio_sm, prec_ratio_wna1)

dusk_plot_hil = plot!(E_bins, normalizer .* sim_ratio_sm, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
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
              6.5 23.0 3  "hr_HI_NITE_MODEL" c4 "Night";
              ]

energy = stats_1.en
ll_nite_lo = stats_1."L<5_night"
ll_nite_md = stats_2."L<5_night"
ll_nite_hi = stats_3."L<5_night"
hl_nite_lo = stats_1."L>5_night"
hl_nite_md = stats_2."L>5_night"
hl_nite_hi = stats_3."L>5_night"


nite_plot_lol = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)
L, MLT, Kp, scenario, colour, label = test_cases[1,:]
@time @load "result_matrix_stats/"*scenario*".jld2" rm
@info "loaded $scenario.jld2"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
norm = normalize_to_elfin(ll_nite_md, sim_ratio_sm)
b1 = 1/(norm * maximum(sim_ratio_sm))
normalizer = b1*norm

E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("LO", "NITE", 1)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("LO", "NITE", 2)
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("LO", "NITE", 3)
dc_norm = normalizer * normalize_to(80, E_bins, E, sim_ratio_sm, prec_ratio_wna1)

nite_plot_lol = plot!(E_bins, normalizer .* sim_ratio_sm, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
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
@info "loaded $scenario.jld2"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 8, 5)
norm = normalize_to_elfin(hl_nite_md, sim_ratio_sm)
b1 = 1/(norm * maximum(sim_ratio_sm))
normalizer = b1*norm

E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("HI", "NITE", 1)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("HI", "NITE", 2)
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("HI", "NITE", 3)
dc_norm = normalizer * normalize_to(80, E_bins, E, sim_ratio_sm, prec_ratio_wna1)

nite_plot_hil = plot!(E_bins, normalizer .* sim_ratio_sm, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
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