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

save_density_plot = true
save_relative_density_comparison_plot = true
save_oblique_comparison_plot = true
save_frequency_model_compare_plot = true
save_constant_frequency_plot = true
save_trajectory_comparison = false
save_qldc_comparison_plot = true
save_qldc_model_comparison_plot = true
visualize_frequency_models = true
save_frequency_model_comparison = true
save_multi_line_plot_comparison = true
save_tps_QLDC_statistical_comparison_plot = true
save_qldc_scaling_comparison_plot = true
save_all_qldc_comparison_plot = true

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
test_cases = [6.5 23.0 3 "main/results/original/" "HI_NITE_MODEL_1m" red "FAW, Ω_pe = 6.5, TPS";
              6.5 23.0 3 "main/results/plasma_density/" "HI_NITE_MODEL_omegape3" purple "FAW, Ω_pe = 3.0, TPS";
              ]

L, MLT, Kp, dir, scenario, colour, label = test_cases[1,:]
@time @load "$dir$scenario.jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_1 = smooth(sim_ratio[1], 11, 6)
norm_1 = normalize_to(97, energy, sim_plot_E_bins, hl_nite_md, sim_ratio_sm_1)

L, MLT, Kp, dir, scenario, colour, label = test_cases[2,:]
@time @load "$dir$scenario.jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_2 = smooth(sim_ratio[1], 11, 6)
norm_2 = normalize_to(97, energy, sim_plot_E_bins, hl_nite_md, sim_ratio_sm_2)

density_comparison_plot = plot(xscale=:log10, yscale=:log10, xlim=(80,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true)

density_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = black, label=false)
density_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = black, label=false)
density_comparison_plot = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = black, linewidth=5, markershape=:circle);
            
L, MLT, Kp, dir, scenario, colour, label = test_cases[1,:]
density_comparison_plot = plot!(sim_plot_E_bins, norm_1*sim_ratio_sm_1, label="$label: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)

E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("HI", "NITE", 1, "const", "6.5")
dc_norm_1 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna1)
density_comparison_plot = plot!(E, dc_norm_1 * prec_ratio_wna1, label = "FAW (θ≤30), Ω_pe = 6.5, QLDC: L=6.5, Night",  color = colour, linewidth=2, linestyle=:dash)

L, MLT, Kp, dir, scenario, colour, label = test_cases[2,:]
density_comparison_plot = plot!(sim_plot_E_bins, norm_2*sim_ratio_sm_2, label="$label: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)

E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("HI", "NITE", 1, "const", "3")
dc_norm_2 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna2)
density_comparison_plot = plot!(E, dc_norm_2 * prec_ratio_wna2, label = "FAW (θ≤30), Ω_pe = 3.0, QLDC: L=6.5, Night",  color = blue, linewidth=2, linestyle=:dash)

# E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("Chorus_Daa_Precip_FreConst_FpeFce3Eq_Max300pT_L4.50_LO_NITE_WN1_Neq1.txt")
# dc_norm = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna3)
# density_comparison_plot = plot!(E, dc_norm * prec_ratio_wna3, label = "FAW (θ≤5), QLDC, Ω_pe,eq = 3, N=1",  color = green, linewidth=2, linestyle=:dash)

E, Daa_wna4, prec_ratio_wna4 = obtain_diffusion_results("Chorus_Daa_Precip_FreConst_FpeFce3Eq_Max300pT_L4.50_LO_NITE_WN1_Nlt10.txt")
dc_norm_4 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna4)
density_comparison_plot = plot!(E, dc_norm_4 * prec_ratio_wna4, label = "FAW (θ≤5), Ω_pe = 3.0, QLDC: L=6.5, Night",  color = purple, linewidth=2, linestyle=:dash)

density_comparison_plot = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
density_comparison_plot = plot!(dpi = 500,size=(1000,650), margin=20px, bottom_margin=12px)

if save_density_plot
    savefig(density_comparison_plot, "main/images/density_comparison.png")
    savefig(density_comparison_plot, "main/images/density_comparison.pdf")
end


x = dc_norm_1 * prec_ratio_wna1
p = LinearInterpolator(E,x)
relative_density_comparison_plot = plot(xscale=:log10, xlim=(80,1000), ylim=(0.0, 2.1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true, legend=:bottomleft)
relative_density_comparison_plot = plot!(energy, hl_nite_md ./ p.(energy), label = "ELFIN Night: L>5, 18<MLT<4", color = black, linewidth = 3, markershape=:circle)
L, MLT, Kp, dir, scenario, colour, label = test_cases[1,:]
relative_density_comparison_plot = plot!(E, dc_norm_1 * prec_ratio_wna1 ./ x, label = "FAW (θ≤30), Ω_pe = 6.5, QLDC: L=6.5, Night", color = colour, linewidth=2, linestyle=:dash)
relative_density_comparison_plot = plot!(sim_plot_E_bins, norm_1 * sim_ratio_sm_1 ./ p.(sim_plot_E_bins), label="$label: L=$L, MLT=$MLT",  color = colour, linewidth=2)
relative_density_comparison_plot = plot!(sim_plot_E_bins, norm_1 * sim_ratio_sm_1 ./ p.(sim_plot_E_bins), fillrange=p.(sim_plot_E_bins)./p.(sim_plot_E_bins), fillalpha = 0.05, label=false, color = colour, linealpha=0)

relative_density_comparison_plot = plot!(E, dc_norm_2 * prec_ratio_wna2 ./ x, label = "FAW (θ≤30), Ω_pe = 3.0, QLDC: L=6.5, Night", color = blue, linewidth=2, linestyle=:dash)
relative_density_comparison_plot = plot!(E, dc_norm_2 * prec_ratio_wna2 ./ x, fillrange=x./x, fillalpha = 0.05, label=false, color = blue, linealpha=0)

L, MLT, Kp, dir, scenario, colour, label = test_cases[2,:]
relative_density_comparison_plot = plot!(E, dc_norm_4 * prec_ratio_wna4 ./ x, label = "FAW (θ≤5), Ω_pe = 3.0, QLDC: L=6.5, Night", color = colour, linewidth=2, linestyle=:dash)
relative_density_comparison_plot = plot!(E, dc_norm_4 * prec_ratio_wna4 ./ x, fillrange=x./x, fillalpha = 0.05, label=false, color = colour, linealpha=0)
relative_density_comparison_plot = plot!(sim_plot_E_bins, norm_2*sim_ratio_sm_2 ./ p.(sim_plot_E_bins), label="$label: L=$L, MLT=$MLT", color = colour, linewidth=2)
relative_density_comparison_plot = plot!(sim_plot_E_bins, norm_2*sim_ratio_sm_2 ./ p.(sim_plot_E_bins), fillrange=p.(sim_plot_E_bins)./p.(sim_plot_E_bins), fillalpha = 0.05, label=false, color = colour, linealpha=0)

if save_relative_density_comparison_plot
    savefig(relative_density_comparison_plot, "main/images/relative_density_comparison.png")
    savefig(relative_density_comparison_plot, "main/images/relative_density_comparison.pdf")
end


###########################
# Oblique Comparison Plot #
###########################
test_cases = [6.5 23.0 3 "main/results/original/" "HI_NITE_MODEL_1m" red "FAW, TPS";
              6.5 23.0 3  "main/results/oblique/" "HI_NITE_WNA1_5m"  green "WNA1, TPS";
              6.5 23.0 3  "main/results/oblique/" "HI_NITE_WNA2_5m"  blue "WNA2, TPS";
              6.5 23.0 3  "main/results/oblique/" "HI_NITE_WNA3_5m"  purple "WNA3, TPS";
              ]

L, MLT, Kp, dir, scenario, colour, label = test_cases[1,:]
@time @load "$dir$scenario.jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_1 = smooth(sim_ratio[1], 11, 6)
norm_1 = normalize_to(97, energy, sim_plot_E_bins, hl_nite_md, sim_ratio_sm_1)

L, MLT, Kp, dir, scenario, colour, label = test_cases[2,:]
@time @load "$dir$scenario.jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_2 = smooth(sim_ratio[1], 11, 6)
norm_2 = normalize_to(97, energy, sim_plot_E_bins, hl_nite_md, sim_ratio_sm_2)

L, MLT, Kp, dir, scenario, colour, label = test_cases[3,:]
@time @load "$dir$scenario.jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_3 = smooth(sim_ratio[1], 11, 6)
norm_3 = normalize_to(97, energy, sim_plot_E_bins, hl_nite_md, sim_ratio_sm_3)

L, MLT, Kp, dir, scenario, colour, label = test_cases[4,:]
@time @load "$dir$scenario.jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_4 = smooth(sim_ratio[1], 11, 6)
norm_4 = normalize_to(97, energy, sim_plot_E_bins, hl_nite_md, sim_ratio_sm_4)

E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("HI", "NITE", 1, "const", "6.5")
dc_norm_1 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna1)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("HI", "NITE", 6, "const", "6.5")
dc_norm_2 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna2)
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("HI", "NITE", 4, "const", "6.5")
dc_norm_3 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna3)
E, Daa_wna4, prec_ratio_wna4 = obtain_diffusion_results("HI", "NITE", 7, "const", "6.5")
dc_norm_4 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna4)

oblique_comparison_plot = plot(xscale=:log10, yscale=:log10, xlim=(80,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true)

oblique_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = black, label=false)
oblique_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = black, label=false)
oblique_comparison_plot = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = black, linewidth=4, markershape=:circle);

L, MLT, Kp, dir, scenario, colour, label = test_cases[1,:]
oblique_comparison_plot = plot!(sim_plot_E_bins, norm_1*sim_ratio_sm_1, label="$label: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=3, markersize = 3)
oblique_comparison_plot = plot!(E, dc_norm_1 * prec_ratio_wna1, label = "FAW, QLDC: L=6.5, Night",  color = colour, linewidth=2, linestyle=:dash)

L, MLT, Kp, dir, scenario, colour, label = test_cases[2,:]
oblique_comparison_plot = plot!(sim_plot_E_bins, norm_2*sim_ratio_sm_2, label="$label: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=3, markersize = 3)
oblique_comparison_plot = plot!(E, dc_norm_2 * prec_ratio_wna2, label = "WNA1, QLDC: L=6.5, Night",  color = colour, linewidth=2, linestyle=:dash)

L, MLT, Kp, dir, scenario, colour, label = test_cases[3,:]
oblique_comparison_plot = plot!(sim_plot_E_bins, norm_3*sim_ratio_sm_3, label="$label: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=3, markersize = 3)
oblique_comparison_plot = plot!(E, dc_norm_3 * prec_ratio_wna3, label = "WNA2, QLDC: L=6.5, Night",  color = colour, linewidth=2, linestyle=:dashdot)

L, MLT, Kp, dir, scenario, colour, label = test_cases[4,:]
# oblique_comparison_plot = plot!(sim_plot_E_bins, norm_4*sim_ratio_sm_4, label="$label: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=3, markersize = 3)
oblique_comparison_plot = plot!(E, dc_norm_4 * prec_ratio_wna4, label = "WNA3, QLDC: L=6.5, Night",  color = colour, linewidth=2, linestyle=:dash)

oblique_comparison_plot = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
oblique_comparison_plot = plot!(dpi = 500,size=(1000,600), margin=20px, bottom_margin=12px)

if save_oblique_comparison_plot
    savefig(oblique_comparison_plot, "main/images/oblique_comparison.png")
    savefig(oblique_comparison_plot, "main/images/oblique_comparison.pdf")
end


# oblique_comparison_tps_only = plot(xscale=:log10, yscale=:log10, xlim=(80,1000), ylim=(1e-2, 1),
# xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true)
# oblique_comparison_tps_only = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = black, label=false)
# oblique_comparison_tps_only = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = black, label=false)
# oblique_comparison_tps_only = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = black, linewidth=4, markershape=:circle);
# L, MLT, Kp, dir, scenario, colour, label = test_cases[1,:]
# oblique_comparison_tps_only = plot!(sim_plot_E_bins, norm_1*sim_ratio_sm_1, label="$label: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=3, markersize = 3)
# L, MLT, Kp, dir, scenario, colour, label = test_cases[2,:]
# oblique_comparison_tps_only = plot!(sim_plot_E_bins, norm_2*sim_ratio_sm_2, label="$label: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=3, markersize = 3)
# L, MLT, Kp, dir, scenario, colour, label = test_cases[3,:]
# oblique_comparison_tps_only = plot!(sim_plot_E_bins, norm_3*sim_ratio_sm_3, label="$label: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=3, markersize = 3)
# L, MLT, Kp, dir, scenario, colour, label = test_cases[4,:]
# oblique_comparison_tps_only = plot!(sim_plot_E_bins, norm_4*sim_ratio_sm_4, label="$label: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=3, markersize = 3)
# oblique_comparison_tps_only = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
# oblique_comparison_tps_only = plot!(dpi = 500,size=(1000,600), margin=20px, bottom_margin=12px)


E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("HI", "NITE", 1, "const", "3")
dc_norm_1 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna1)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("HI", "NITE", 6, "const", "3")
dc_norm_2 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna2)
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("HI", "NITE", 4, "const", "3")
dc_norm_3 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna3)
E, Daa_wna4, prec_ratio_wna4 = obtain_diffusion_results("HI", "NITE", 7, "const", "3")
dc_norm_4 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna4)

E, Daa_wna5, prec_ratio_wna5 = obtain_diffusion_results("HI", "NITE", 1, "vary", "6.5")
dc_norm_5 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna5)
E, Daa_wna6, prec_ratio_wna6 = obtain_diffusion_results("HI", "NITE", 6, "vary", "6.5")
dc_norm_6 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna6)
E, Daa_wna7, prec_ratio_wna7 = obtain_diffusion_results("HI", "NITE", 4, "vary", "6.5")
dc_norm_7 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna7)
E, Daa_wna8, prec_ratio_wna8 = obtain_diffusion_results("HI", "NITE", 7, "vary", "6.5")
dc_norm_8 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna8)


oblique_comparison_qldc_only = plot(xscale=:log10, yscale=:log10, xlim=(80,1000), ylim=(1e-2, 1),
xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true)
oblique_comparison_qldc_only = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = black, label=false)
oblique_comparison_qldc_only = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = black, label=false)
oblique_comparison_qldc_only = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = black, linewidth=4, markershape=:circle);
L, MLT, Kp, dir, scenario, colour, label = test_cases[1,:]
oblique_comparison_qldc_only = plot!(E, dc_norm_1 * prec_ratio_wna1, label = "FAW, ω_m = const, Ω_pe = 3.0",  color = colour, linewidth=2, linestyle=:dash)
oblique_comparison_qldc_only = plot!(E, dc_norm_5 * prec_ratio_wna5, label = "FAW, ω_m = vary, Ω_pe = 6.5",  color = colour, linewidth=3, linestyle=:dashdot)
L, MLT, Kp, dir, scenario, colour, label = test_cases[2,:]
oblique_comparison_qldc_only = plot!(E, dc_norm_2 * prec_ratio_wna2, label = "WNA1, ω_m = const, Ω_pe = 3.0",  color = colour, linewidth=2, linestyle=:dash)
oblique_comparison_qldc_only = plot!(E, dc_norm_6 * prec_ratio_wna6, label = "WNA1, ω_m = vary, Ω_pe = 6.5",  color = colour, linewidth=3, linestyle=:dashdot)
L, MLT, Kp, dir, scenario, colour, label = test_cases[3,:]
oblique_comparison_qldc_only = plot!(E, dc_norm_3 * prec_ratio_wna3, label = "WNA2, ω_m = const, Ω_pe = 3.0",  color = colour, linewidth=2, linestyle=:dash)
oblique_comparison_qldc_only = plot!(E, dc_norm_7 * prec_ratio_wna7, label = "WNA2, ω_m = vary, Ω_pe = 6.5",  color = colour, linewidth=3, linestyle=:dashdot)
L, MLT, Kp, dir, scenario, colour, label = test_cases[4,:]
oblique_comparison_qldc_only = plot!(E, dc_norm_4 * prec_ratio_wna4, label = "WNA3, ω_m = const, Ω_pe = 3.0",  color = colour, linewidth=2, linestyle=:dash)
oblique_comparison_qldc_only = plot!(E, dc_norm_8 * prec_ratio_wna8, label = "WNA3, ω_m = vary, Ω_pe = 6.5",  color = colour, linewidth=3, linestyle=:dashdot)
oblique_comparison_qldc_only = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
oblique_comparison_qldc_only = plot!(dpi = 500,size=(600,700), margin=20px, bottom_margin=12px)

    
if save_oblique_comparison_plot
    savefig(oblique_comparison_qldc_only, "main/images/oblique_comparison_low_plasma_density.png")
    savefig(oblique_comparison_qldc_only, "main/images/oblique_comparison_low_plasma_density.pdf")
end


##############################
# Frequency model comparison #
##############################
test_cases = [6.5 23.0 3 "main/results/original/" "HI_NITE_MODEL_1m" red "FAW, ω_m = const";
              6.5 23.0 3 "main/results/frequency/" "HI_NITE_MODEL_omega_mod_3" blue "FAW, ω_m(λ)";
            #   6.5 23.0 3 "main/results/combined/" "HI_NITE_WNA1MODEL3" purple "WNA1, ω_m(λ)";
              ]

L, MLT, Kp, dir, scenario, colour, label = test_cases[1,:]
@time @load "$dir$scenario.jld2" rm
@info "loaded $scenario from $dir"
sim_ratio_1 = prec_to_trap_ratio(rm)
sim_ratio_sm_1 = smooth(sim_ratio_1[1], 11, 6)
norm_1 = normalize_to(97, energy, sim_plot_E_bins, hl_nite_md, sim_ratio_sm_1)
L, MLT, Kp, dir, scenario, colour, label = test_cases[2,:]
@time @load "$dir$scenario.jld2" rm
@info "loaded $scenario from $dir"
sim_ratio_2 = prec_to_trap_ratio(rm)
sim_ratio_sm_2 = smooth(sim_ratio_2[1], 11, 6)
norm_2 = normalize_to(97, energy, sim_plot_E_bins, hl_nite_md, sim_ratio_sm_2)
# L, MLT, Kp, dir, scenario, colour, label = test_cases[3,:]
# @time @load "$dir$scenario.jld2" rm
# @info "loaded $scenario from $dir"
# sim_ratio_3 = prec_to_trap_ratio(rm)
# sim_ratio_sm_3 = smooth(sim_ratio_3[1], 11, 6)
# norm_3 = normalize_to(97, energy, sim_plot_E_bins, hl_nite_md, sim_ratio_sm_3)

E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("HI", "NITE", 1, "const", "6.5")
dc_norm_1 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna1)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("HI", "NITE", 1, "vary", "6.5")
dc_norm_2 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna2)
# E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("HI", "NITE", 6, "vary", "6.5")
# dc_norm_3 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna3)

frequency_model_compare_plot = plot(xscale=:log10, yscale=:log10, xlim=(80,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true)

frequency_model_compare_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = black, label=false)
frequency_model_compare_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = black, label=false)
frequency_model_compare_plot = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = black, linewidth=5, markershape=:circle);

L, MLT, Kp, dir, scenario, colour, label = test_cases[1,:]
frequency_model_compare_plot = plot!(sim_plot_E_bins, norm_1*sim_ratio_sm_1, label="$label, TPS: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=3, markersize = 3)
frequency_model_compare_plot = plot!(E, dc_norm_1 * prec_ratio_wna1, label = "FAW, ω_m = const, QLDC: L=$L, Night",  color = colour, linewidth=2, linestyle=:dash)

L, MLT, Kp, dir, scenario, colour, label = test_cases[2,:]
frequency_model_compare_plot = plot!(sim_plot_E_bins, norm_2*sim_ratio_sm_2, label="$label, TPS: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=3, markersize = 3)
frequency_model_compare_plot = plot!(E, dc_norm_2 * prec_ratio_wna2, label = "FAW, ω_m(λ), QLDC: L=$L, Night",  color = blue, linewidth=2, linestyle=:dash)

# L, MLT, Kp, dir, scenario, colour, label = test_cases[3,:]
# frequency_model_compare_plot = plot!(sim_plot_E_bins, norm_3*sim_ratio_sm_3, label="$label, TPS: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=3, markersize = 3)
# frequency_model_compare_plot = plot!(E, dc_norm_3 * prec_ratio_wna3, label = "WNA1, ω_m(λ), QLDC: L=$L, Night",  color = purple, linewidth=2, linestyle=:dash)

frequency_model_compare_plot = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
frequency_model_compare_plot = plot!(dpi = 500,size=(1000,600), margin=20px, bottom_margin=12px)

if save_frequency_model_compare_plot
    savefig(frequency_model_compare_plot, "main/images/freq_model_comparison.png")
    savefig(frequency_model_compare_plot, "main/images/freq_model_comparison.pdf")
end


#############################
# Frequency Comparison Plot #
#############################
test_cases = [6.5 23.0 3  "main/results/frequency/" "HI_NITE_MODEL_15" red "ω_m = 0.15";
              6.5 23.0 3  "main/results/frequency/" "HI_NITE_MODEL_20" orange "ω_m = 0.20";
              6.5 23.0 3  "main/results/frequency/" "HI_NITE_MODEL_25" yellow "ω_m = 0.25";
              6.5 23.0 3  "main/results/frequency/" "HI_NITE_MODEL_30" green "ω_m = 0.30";
              6.5 23.0 3  "main/results/frequency/" "HI_NITE_MODEL_35" blue "ω_m = 0.35";
              6.5 23.0 3  "main/results/frequency/" "HI_NITE_MODEL_40" purplish_pink "ω_m = 0.40";
              6.5 23.0 3  "main/results/frequency/" "HI_NITE_MODEL_45" purple "ω_m = 0.45";
              ]


# I want to normalize to omega_m = 0.35 which is what I have been using in previous plots
L, MLT, Kp, dir, scenario, colour, label = test_cases[5,:]
@time @load "$dir$scenario.jld2" rm
@info "loaded $scenario.jld2 from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 11, 6)
norm_35 = normalize_to(97, energy, sim_plot_E_bins, hl_nite_md, sim_ratio_sm)

constant_frequency_comparison_plot = plot(xscale=:log10, yscale=:log10, xlim=(80,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true)


constant_frequency_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = black, label=false)
constant_frequency_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = black, label=false)
constant_frequency_comparison_plot = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = black, linewidth=2, markershape=:circle);

for i in 1:length(test_cases[:,1])
    L, MLT, Kp, dir, scenario, colour, label = test_cases[i,:]
    @time @load "$dir$scenario.jld2" rm
    @info "loaded $scenario.jld2 from $dir"
    sim_ratio = prec_to_trap_ratio(rm)
    sim_ratio_sm = smooth(sim_ratio[1], 11, 6)
    constant_frequency_comparison_plot = plot!(sim_plot_E_bins, norm_35*sim_ratio_sm, label="$label, L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
end

constant_frequency_comparison_plot = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
constant_frequency_comparison_plot = plot!(dpi = 500,size=(600,700), margin=20px, bottom_margin=12px)

if save_constant_frequency_plot
    savefig(constant_frequency_comparison_plot, "main/images/const_freq_comparison.png")
    savefig(constant_frequency_comparison_plot, "main/images/const_freq_comparison.pdf")
end


###################
# QLDC Comparison #
###################
E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("HI", "NITE", 1, "const", "6.5")
dc_norm_1 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna1)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("HI", "NITE", 1, "vary", "6.5")
dc_norm_2 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna2)
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("HI", "NITE", 7, "const", "6.5")
dc_norm_3 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna3)
E, Daa_wna4, prec_ratio_wna4 = obtain_diffusion_results("HI", "NITE", 7, "vary", "6.5")
dc_norm_4 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna4)
E, Daa_wna5, prec_ratio_wna5 = obtain_diffusion_results("HI", "NITE", 1, "const", "3")
dc_norm_5 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna5)
E, Daa_wna6, prec_ratio_wna6 = obtain_diffusion_results("HI", "NITE", 7, "vary", "3")
dc_norm_6 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna6)

qldc_model_comparison_plot = plot(xscale=:log10, yscale=:log10, xlim=(80,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true)

qldc_model_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = black, label=false)
qldc_model_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = black, label=false)
qldc_model_comparison_plot = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = black, linewidth=5, markershape=:circle);
qldc_model_comparison_plot = plot!(E, dc_norm_1 * prec_ratio_wna1, label = "FAW, ω_m = const, Ω_pe = 6.5",  color = red, linewidth=2, linestyle=:dash)
qldc_model_comparison_plot = plot!(E, dc_norm_2 * prec_ratio_wna2, label = "FAW, ω_m(λ), Ω_pe = 6.5",  color = blue, linewidth=2, linestyle=:dash)
qldc_model_comparison_plot = plot!(E, dc_norm_3 * prec_ratio_wna3, label = "WNA3, ω_m = const, Ω_pe = 6.5",  color = green, linewidth=2, linestyle=:dash)
qldc_model_comparison_plot = plot!(E, dc_norm_4 * prec_ratio_wna4, label = "WNA3, ω_m(λ), Ω_pe = 6.5",  color = purple, linewidth=2, linestyle=:dash)
qldc_model_comparison_plot = plot!(E, dc_norm_5 * prec_ratio_wna5, label = "FAW, ω_m = const, Ω_pe = 3.0",  color = purplish_pink, linewidth=2, linestyle=:dash)
qldc_model_comparison_plot = plot!(E, dc_norm_6 * prec_ratio_wna6, label = "WNA3, ω_m(λ), Ω_pe = 3.0",  color = reddish_orange, linewidth=2, linestyle=:dash)

qldc_model_comparison_plot = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
qldc_model_comparison_plot = plot!(dpi = 500,size=(1000,600), margin=20px, bottom_margin=12px)

if save_qldc_model_comparison_plot
    savefig(qldc_model_comparison_plot, "main/images/qldc_model_comparison.png")
    savefig(qldc_model_comparison_plot, "main/images/qldc_model_comparison.pdf")
end

x = dc_norm_1 * prec_ratio_wna1
p = LinearInterpolator(E,x)
qldc_scaling_comparison_plot = plot(xscale=:log10, xlim=(80,1000), ylim=(0.0, 2.1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true)
qldc_scaling_comparison_plot = plot!(energy, hl_nite_md ./ p.(energy), label = "ELFIN Night: L>5, 18<MLT<4", color = black, linewidth = 3, markershape=:circle)
qldc_scaling_comparison_plot = plot!(E, dc_norm_1 * prec_ratio_wna1 ./ x, label = "FAW, ω_m = const, Ω_pe = 6.5", color = red, linewidth=2, linestyle=:dash)
qldc_scaling_comparison_plot = plot!(E, dc_norm_2 * prec_ratio_wna2 ./ x, label = "FAW, ω_m(λ), Ω_pe = 6.5",  color = blue, linewidth=2, linestyle=:dash)
qldc_scaling_comparison_plot = plot!(E, dc_norm_2 * prec_ratio_wna2 ./ x, fillrange=x./x, fillalpha = 0.05, label=false, color = blue, linealpha=0)
qldc_scaling_comparison_plot = plot!(E, dc_norm_3 * prec_ratio_wna3 ./ x, label = "WNA3, ω_m = const, Ω_pe = 6.5", color = green, linewidth=2, linestyle=:dash)
qldc_scaling_comparison_plot = plot!(E, dc_norm_3 * prec_ratio_wna3 ./ x, fillrange=x./x, fillalpha = 0.05, label=false, color = green, linealpha=0)
qldc_scaling_comparison_plot = plot!(E, dc_norm_4 * prec_ratio_wna4 ./ x, label = "WNA3, ω_m(λ), Ω_pe = 6.5",  color = purple, linewidth=2, linestyle=:dash)
qldc_scaling_comparison_plot = plot!(E, dc_norm_4 * prec_ratio_wna4 ./ x, fillrange=x./x, fillalpha = 0.03, label=false, color = purple, linealpha=0)
qldc_scaling_comparison_plot = plot!(E, dc_norm_5 * prec_ratio_wna5 ./ x, label = "FAW, ω_m = const, Ω_pe = 3.0",  color = purplish_pink, linewidth=2, linestyle=:dash)
qldc_scaling_comparison_plot = plot!(E, dc_norm_5 * prec_ratio_wna5 ./ x, fillrange=x./x, fillalpha = 0.05, label=false, color = purplish_pink, linealpha=0)
qldc_scaling_comparison_plot = plot!(E, dc_norm_6 * prec_ratio_wna6 ./ x, label = "WNA3, ω_m(λ), Ω_pe = 3.0", color = reddish_orange, linewidth=2, linestyle=:dash)
qldc_scaling_comparison_plot = plot!(E, dc_norm_6 * prec_ratio_wna6 ./ x, fillrange=x./x, fillalpha = 0.01, label=false, color = reddish_orange, linealpha=0)
qldc_scaling_comparison_plot = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
qldc_scaling_comparison_plot = plot!(dpi = 500,size=(1000,600), margin=20px, bottom_margin=12px)

if save_qldc_scaling_comparison_plot
    savefig(qldc_scaling_comparison_plot, "main/images/qldc_scaling_comparison.png")
    savefig(qldc_scaling_comparison_plot, "main/images/qldc_scaling_comparison.pdf")
end


#######################
# All QLDC Comparison #
#######################
E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("HI", "NITE", 1, "const", "L") #
dc_norm_1 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna1)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("HI", "NITE", 1, "vary", "L")
dc_norm_2 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna2)
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("HI", "NITE", 6, "const", "6.5")
dc_norm_3 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna3)
E, Daa_wna4, prec_ratio_wna4 = obtain_diffusion_results("HI", "NITE", 6, "vary", "6.5")
dc_norm_4 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna4)
E, Daa_wna5, prec_ratio_wna5 = obtain_diffusion_results("HI", "NITE", 1, "const", "3")
dc_norm_5 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna5)
E, Daa_wna6, prec_ratio_wna6 = obtain_diffusion_results("HI", "NITE", 1, "vary", "3")
dc_norm_6 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna6)
# E, Daa_wna7, prec_ratio_wna7 = obtain_diffusion_results("HI", "NITE", 6, "const", "3")
# dc_norm_7 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna7)
E, Daa_wna8, prec_ratio_wna8 = obtain_diffusion_results("HI", "NITE", 6, "vary", "3")
dc_norm_8 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna8)


all_qldc_comparison_plot = plot(xscale=:log10, yscale=:log10, xlim=(80,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true)

all_qldc_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = black, label=false)
all_qldc_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = black, label=false)
all_qldc_comparison_plot = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = black, linewidth=5, markershape=:circle);

all_qldc_comparison_plot = plot!(E, dc_norm_1 * prec_ratio_wna1, label = "FAW, ω_m = const, Ω_pe = 6.5, Night, QLDC",  color = red, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_2 * prec_ratio_wna2, label = "FAW, ω_m = vary, Ω_pe = 6.5, Night, QLDC",  color = reddish_orange, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_3 * prec_ratio_wna3, label = "WNA1, ω_m = const, Ω_pe = 6.5, Night, QLDC",  color = orange, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_4 * prec_ratio_wna4, label = "WNA1, ω_m = vary, Ω_pe = 6.5, Night, QLDC",  color = yellow, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_5 * prec_ratio_wna5, label = "FAW, ω_m = const, Ω_pe = 3.0, Night, QLDC",  color = green, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_6 * prec_ratio_wna6, label = "FAW, ω_m = vary, Ω_pe = 3.0, Night, QLDC",  color = blue, linewidth=2, linestyle=:dash)
# all_qldc_comparison_plot = plot!(E, dc_norm_7 * prec_ratio_wna7, label = "WNA1, ω_m = const, Ω_pe = 3.0, Night, QLDC",  color = purplish_pink, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_8 * prec_ratio_wna8, label = "WNA1, ω_m = vary, Ω_pe = 3.0, Night, QLDC",  color = purple, linewidth=2, linestyle=:dash)

all_qldc_comparison_plot = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
all_qldc_comparison_plot = plot!(dpi = 500,size=(1000,600), margin=20px, bottom_margin=12px)

if save_all_qldc_comparison_plot
    savefig(all_qldc_comparison_plot, "main/images/all_qldc_comparison_plot.png")
    savefig(all_qldc_comparison_plot, "main/images/all_qldc_comparison_plot.pdf")
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
sol_15 = @time run_model(500, ICrange, 6.5, 23., 3., 0.15, 40, 1, "bw");
sol_45 = @time run_model(500, ICrange, 6.5, 23., 3., 0.45, 40, 1, "bw");

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
test_cases = [6.5 23.0 3 "main/results/frequency/" "HI_NITE_MODEL_omega_mod_1" orange "ω_mod_1";
              6.5 23.0 3 "main/results/frequency/" "HI_NITE_MODEL_omega_mod_2" yellow "ω_mod_2";
              6.5 23.0 3 "main/results/frequency/" "HI_NITE_MODEL_omega_mod_3" red "ω_mod_3";
              6.5 23.0 3 "main/results/frequency/" "HI_NITE_MODEL_15" light_blue "ω_m = 0.15";
              6.5 23.0 3 "main/results/frequency/" "HI_NITE_MODEL_45" green "ω_m = 0.45";
              ]

frequency_model_comparison_plot = plot(xscale=:log10, yscale=:log10, xlim=(80,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true)

@time @load "main/results/frequency/HI_NITE_MODEL_35.jld2" rm
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 8, 5)
norm_35 = normalize_to(97, energy, sim_plot_E_bins, hl_nite_md, sim_ratio_sm)

for i in 1:length(test_cases[:,1])
    L, MLT, Kp, dir, scenario, colour, label = test_cases[i,:]
    @time @load "$dir$scenario.jld2" rm
    @info "loaded $scenario from $dir"
    sim_ratio = prec_to_trap_ratio(rm)
    sim_ratio_sm = smooth(sim_ratio[1], 8, 5)
    norm = normalize_to(97, energy, sim_plot_E_bins, hl_nite_md, sim_ratio_sm)
    frequency_model_comparison_plot = plot!(sim_plot_E_bins, norm_35*sim_ratio_sm, label="$label: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
end

frequency_model_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = black, label=false)
frequency_model_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = black, label=false)
frequency_model_comparison_plot = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = black, linewidth=2, markershape=:circle);
frequency_model_comparison_plot = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
frequency_model_comparison_plot = plot!(dpi = 500,size=(1000,600), margin=20px, bottom_margin=12px)

if save_frequency_model_comparison
    savefig(frequency_model_comparison_plot, "main/images/freq_model_comparison.png")
    savefig(frequency_model_comparison_plot, "main/images/freq_model_comparison.pdf")
end

##############################
# Multi line plot comparison #
##############################
# this blends together model 3 with FAW and with WNA1
test_cases = [6.5 23.0 3 "main/results/frequency/" "HI_NITE_MODEL_omega_mod_2" blue "ω_m(λ), ω_m(λ)";
              ]
L, MLT, Kp, dir, scenario, colour, label = test_cases[1,:]
@time @load "$dir$scenario.jld2" rm
@info "loaded $scenario from $dir"
sim_ratio_1 = prec_to_trap_ratio(rm)
sim_ratio_sm_1 = smooth(sim_ratio_1[1], 11, 6)
norm_1 = normalize_to(97, energy, sim_plot_E_bins, hl_nite_md, sim_ratio_sm_1)

E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("HI", "NITE", 1, "vary", "6.5")
dc_norm_2 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna2)
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("HI", "NITE", 6, "vary", "6.5")
dc_norm_3 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna3)
E, Daa_wna4, prec_ratio_wna4 = obtain_diffusion_results("HI", "NITE", 6, "vary", "6.5")
dc_norm_4 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna4)
# f3 = f1*a + f2*(1-a), w/ a<1; f1 is model3faw, f2 is model3wna1
f3(a) = @. dc_norm_2*prec_ratio_wna2*a + dc_norm_3*prec_ratio_wna3*(1-a)

multi_line_plot = plot(xscale=:log10, yscale=:log10, xlim=(80,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true)
multi_line_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = black, label=false)
multi_line_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = black, label=false)
multi_line_plot = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = black, linewidth=3, markershape=:circle);
for i in 0.1:0.2:0.9
    multi_line_plot = plot!(E, f3(i), label="", linestyle=:dash, linewidth = 1)
end
multi_line_plot = plot!(E, dc_norm_2 * prec_ratio_wna2, label = "ω_m(λ), FAW, Night, Ω_pe = 6.5, QLDC",  color = blue, linewidth=2.5, linestyle=:dash)
multi_line_plot = plot!(E, dc_norm_3 * prec_ratio_wna3, label = "ω_m(λ), WNA1, Night, Ω_pe = 6.5, QLDC",  color = green, linewidth=2.5, linestyle=:dash)
multi_line_plot = plot!(E, dc_norm_4 * prec_ratio_wna4, label = "ω_m(λ), WNA1, Night, Ω_pe = 3, QLDC",  color = purple, linewidth=2.5, linestyle=:dash)

# multi_line_plot = plot!(sim_plot_E_bins, norm_4*sim_ratio_sm_4, label="label Model: L=L, MLT=MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
multi_line_plot = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
multi_line_plot = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)

if save_multi_line_plot_comparison
    savefig(multi_line_plot, "main/images/multi_line_plot.png")
    savefig(multi_line_plot, "main/images/multi_line_plot.pdf")
end


#######################
# All QLDC Comparison #
#######################
############################
# FAW, Freqconst, vary fpe #
############################
E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("HI", "NITE", 1, "const", "6")
dc_norm_1 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna1)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("HI", "NITE", 1, "const", "5.5")
dc_norm_2 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna2)
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("HI", "NITE", 1, "const", "5")
dc_norm_3 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna3)
E, Daa_wna4, prec_ratio_wna4 = obtain_diffusion_results("HI", "NITE", 1, "const", "4.5")
dc_norm_4 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna4)
E, Daa_wna5, prec_ratio_wna5 = obtain_diffusion_results("HI", "NITE", 1, "const", "4")
dc_norm_5 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna5)
E, Daa_wna6, prec_ratio_wna6 = obtain_diffusion_results("HI", "NITE", 1, "const", "3.5")
dc_norm_6 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna6)
E, Daa_wna7, prec_ratio_wna7 = obtain_diffusion_results("HI", "NITE", 1, "const", "3")
dc_norm_7 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna7)
E, Daa_wna8, prec_ratio_wna8 = obtain_diffusion_results("HI", "NITE", 1, "const", "2.5")
dc_norm_8 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna8)

all_qldc_comparison_plot = plot(xscale=:log10, yscale=:log10, xlim=(80,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true)

all_qldc_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = black, label=false)
all_qldc_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = black, label=false)
all_qldc_comparison_plot = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = black, linewidth=5, markershape=:circle);

all_qldc_comparison_plot = plot!(E, dc_norm_1 * prec_ratio_wna1, label = "FAW, ω_m = const, Ω_pe = 6.0, Night, QLDC",  color = red, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_2 * prec_ratio_wna2, label = "FAW, ω_m = const, Ω_pe = 5.5, Night, QLDC",  color = reddish_orange, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_3 * prec_ratio_wna3, label = "FAW, ω_m = const, Ω_pe = 5.0, Night, QLDC",  color = orange, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_4 * prec_ratio_wna4, label = "FAW, ω_m = const, Ω_pe = 4.5, Night, QLDC",  color = yellow, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_5 * prec_ratio_wna5, label = "FAW, ω_m = const, Ω_pe = 4.0, Night, QLDC",  color = green, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_6 * prec_ratio_wna6, label = "FAW, ω_m = const, Ω_pe = 3.5, Night, QLDC",  color = blue, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_7 * prec_ratio_wna7, label = "FAW, ω_m = const, Ω_pe = 3.0, Night, QLDC",  color = purplish_pink, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_8 * prec_ratio_wna8, label = "FAW, ω_m = const, Ω_pe = 2.5, Night, QLDC",  color = purple, linewidth=2, linestyle=:dash)

all_qldc_comparison_plot = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
all_qldc_comparison_plot = plot!(dpi = 500,size=(1200,600), margin=20px, bottom_margin=12px)

savefig(all_qldc_comparison_plot, "main/images/multiplota_faw_freqconst_fpeparam.pdf")

###########################
# FAW, Freqvary, vary fpe #
###########################
E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("HI", "NITE", 1, "vary", "6")
dc_norm_1 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna1)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("HI", "NITE", 1, "vary", "5.5")
dc_norm_2 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna2)
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("HI", "NITE", 1, "vary", "5")
dc_norm_3 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna3)
E, Daa_wna4, prec_ratio_wna4 = obtain_diffusion_results("HI", "NITE", 1, "vary", "4.5")
dc_norm_4 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna4)
E, Daa_wna5, prec_ratio_wna5 = obtain_diffusion_results("HI", "NITE", 1, "vary", "4")
dc_norm_5 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna5)
E, Daa_wna6, prec_ratio_wna6 = obtain_diffusion_results("HI", "NITE", 1, "vary", "3.5")
dc_norm_6 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna6)
E, Daa_wna7, prec_ratio_wna7 = obtain_diffusion_results("HI", "NITE", 1, "vary", "3")
dc_norm_7 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna7)
E, Daa_wna8, prec_ratio_wna8 = obtain_diffusion_results("HI", "NITE", 1, "vary", "2.5")
dc_norm_8 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna8)

all_qldc_comparison_plot = plot(xscale=:log10, yscale=:log10, xlim=(80,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true)

all_qldc_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = black, label=false)
all_qldc_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = black, label=false)
all_qldc_comparison_plot = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = black, linewidth=5, markershape=:circle);

all_qldc_comparison_plot = plot!(E, dc_norm_1 * prec_ratio_wna1, label = "FAW, ω_m = vary, Ω_pe = 6.0, Night, QLDC",  color = red, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_2 * prec_ratio_wna2, label = "FAW, ω_m = vary, Ω_pe = 5.5, Night, QLDC",  color = reddish_orange, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_3 * prec_ratio_wna3, label = "FAW, ω_m = vary, Ω_pe = 5.0, Night, QLDC",  color = orange, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_4 * prec_ratio_wna4, label = "FAW, ω_m = vary, Ω_pe = 4.5, Night, QLDC",  color = yellow, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_5 * prec_ratio_wna5, label = "FAW, ω_m = vary, Ω_pe = 4.0, Night, QLDC",  color = green, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_6 * prec_ratio_wna6, label = "FAW, ω_m = vary, Ω_pe = 3.5, Night, QLDC",  color = blue, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_7 * prec_ratio_wna7, label = "FAW, ω_m = vary, Ω_pe = 3.0, Night, QLDC",  color = purplish_pink, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_8 * prec_ratio_wna8, label = "FAW, ω_m = vary, Ω_pe = 2.5, Night, QLDC",  color = purple, linewidth=2, linestyle=:dash)

all_qldc_comparison_plot = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
all_qldc_comparison_plot = plot!(dpi = 500,size=(1200,600), margin=20px, bottom_margin=12px)

savefig(all_qldc_comparison_plot, "main/images/multiplotb_faw_freqvary_fpeparam.pdf")

############################
# WNA1, Freqvary, vary fpe #
############################
E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("HI", "NITE", 6, "vary", "6")
dc_norm_1 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna1)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("HI", "NITE", 6, "vary", "5.5")
dc_norm_2 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna2)
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("HI", "NITE", 6, "vary", "5")
dc_norm_3 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna3)
E, Daa_wna4, prec_ratio_wna4 = obtain_diffusion_results("HI", "NITE", 6, "vary", "4.5")
dc_norm_4 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna4)
E, Daa_wna5, prec_ratio_wna5 = obtain_diffusion_results("HI", "NITE", 6, "vary", "4")
dc_norm_5 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna5)
E, Daa_wna6, prec_ratio_wna6 = obtain_diffusion_results("HI", "NITE", 6, "vary", "3.5")
dc_norm_6 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna6)
E, Daa_wna7, prec_ratio_wna7 = obtain_diffusion_results("HI", "NITE", 6, "vary", "3")
dc_norm_7 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna7)
E, Daa_wna8, prec_ratio_wna8 = obtain_diffusion_results("HI", "NITE", 6, "vary", "2.5")
dc_norm_8 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna8)

all_qldc_comparison_plot = plot(xscale=:log10, yscale=:log10, xlim=(80,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true)

all_qldc_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = black, label=false)
all_qldc_comparison_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = black, label=false)
all_qldc_comparison_plot = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = black, linewidth=5, markershape=:circle);

all_qldc_comparison_plot = plot!(E, dc_norm_1 * prec_ratio_wna1, label = "WNA1, ω_m = vary, Ω_pe = 6.0, Night, QLDC",  color = red, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_2 * prec_ratio_wna2, label = "WNA1, ω_m = vary, Ω_pe = 5.5, Night, QLDC",  color = reddish_orange, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_3 * prec_ratio_wna3, label = "WNA1, ω_m = vary, Ω_pe = 5.0, Night, QLDC",  color = orange, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_4 * prec_ratio_wna4, label = "WNA1, ω_m = vary, Ω_pe = 4.5, Night, QLDC",  color = yellow, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_5 * prec_ratio_wna5, label = "WNA1, ω_m = vary, Ω_pe = 4.0, Night, QLDC",  color = green, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_6 * prec_ratio_wna6, label = "WNA1, ω_m = vary, Ω_pe = 3.5, Night, QLDC",  color = blue, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_7 * prec_ratio_wna7, label = "WNA1, ω_m = vary, Ω_pe = 3.0, Night, QLDC",  color = purplish_pink, linewidth=2, linestyle=:dash)
all_qldc_comparison_plot = plot!(E, dc_norm_8 * prec_ratio_wna8, label = "WNA1, ω_m = vary, Ω_pe = 2.5, Night, QLDC",  color = purple, linewidth=2, linestyle=:dash)

all_qldc_comparison_plot = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
all_qldc_comparison_plot = plot!(dpi = 500,size=(1200,600), margin=20px, bottom_margin=12px)

savefig(all_qldc_comparison_plot, "main/images/multiplotc_wna1_freqvary_fpeparam.pdf")

#######################
# Final Matching Plot #
#######################
E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("HI", "NITE", 1, "vary", "2.5")
dc_norm_1 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna1)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("HI", "NITE", 6, "vary", "2.5")
dc_norm_2 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna2)
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("HI", "NITE", 4, "vary", "2.5")
dc_norm_3 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna3)
E, Daa_wna4, prec_ratio_wna4 = obtain_diffusion_results("HI", "NITE", 7, "vary", "4")
dc_norm_4 = normalize_to(97, energy, E, hl_nite_md, prec_ratio_wna4)


hil_matching_plot = plot(xscale=:log10, yscale=:log10, xlim=(80,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true);

hil_matching_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = black, label=false)
hil_matching_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = black, label=false)
hil_matching_plot = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = black, linewidth=5, markershape=:circle);

hil_matching_plot = plot!(E, dc_norm_1 * prec_ratio_wna1, label = "FAW, ω_m(λ), Ω_pe = 2.5, Night, QLDC",  color = purple, linewidth=2, linestyle=:dash)
hil_matching_plot = plot!(E, dc_norm_2 * prec_ratio_wna2, label = "WNA1, ω_m(λ), Ω_pe = 2.5, Night, QLDC",  color = green, linewidth=2, linestyle=:dash)
hil_matching_plot = plot!(E, dc_norm_3 * prec_ratio_wna3, label = "WNA2, ω_m(λ), Ω_pe = 2.5, Night, QLDC",  color = blue, linewidth=2, linestyle=:dash)
hil_matching_plot = plot!(E, dc_norm_4 * prec_ratio_wna4, label = "WNA3, ω_m(λ), Ω_pe = 4.0, Night, QLDC",  color = orange, linewidth=2, linestyle=:dash)

hil_matching_plot = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
hil_matching_plot = plot!(dpi = 500,size=(1200,600), margin=20px, bottom_margin=12px)

savefig(hil_matching_plot, "main/images/multiplotd_elfin_matched.pdf")



### MLT comparison, FIGURE 4 in JGR2023 except now with diffusion code 

# Dawn plot and generate normalizers
test_cases = [4.5 8.0  3 "main/results/original/" "LO_DAWN_MODEL_1m" purple "Dawn-Noon";
              6.5 8.0  3 "main/results/original/" "HI_DAWN_MODEL_1m" purple "Dawn-Noon";
              ]

L, MLT, Kp, dir, scenario, colour, label = test_cases[1,:]
@time @load "$dir$scenario.jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 11, 6)
norm = normalize_to(97, energy, sim_plot_E_bins, ll_dawn_md, sim_ratio_sm)
b1 = 1/(norm * maximum(sim_ratio_sm))
normalizer = b1*norm

E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("LO", "DAWN", 1, "const", "L")
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("LO", "DAWN", 2, "const", "L")
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("LO", "DAWN", 4, "const", "L")
dc_norm = normalizer * normalize_to(96, E_bins, E, sim_ratio_sm, prec_ratio_wna1)

dawn_plot_lol = plot(xscale=:log10, yscale=:log10, xlim=(80,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true)
dawn_plot_lol = plot!(energy, b1*ll_dawn_md, fillrange=b1*ll_dawn_hi, fillalpha = 0.2, color = black, label=false)
dawn_plot_lol = plot!(energy, b1*ll_dawn_md, fillrange=b1*ll_dawn_lo, fillalpha = 0.2, color = black, label=false)
dawn_plot_lol = plot!(energy, b1*ll_dawn_md, label = "ELFIN Dawn-Noon: L<5, 4<MLT<13 ", color = black, linewidth=2, markershape=:circle);
dawn_plot_lol = plot!(sim_plot_E_bins, normalizer .* sim_ratio_sm, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
dawn_plot_lol = plot!(E, dc_norm * prec_ratio_wna1, label = "FAW: Parallel Waves", color = blue, linestyle=:dash, linewidth=2)
dawn_plot_lol = plot!(E, dc_norm * prec_ratio_wna2, label = "WNA1: initially parallel -> resonance cone", color = green, linestyle=:dash, linewidth=2)
dawn_plot_lol = plot!(E, dc_norm * prec_ratio_wna3, label = "WNA2: gendrin angle -> resonance cone", color = red, linestyle=:dash, linewidth=2)
dawn_plot_lol = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
dawn_plot_lol = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(dawn_plot_lol, "main/images/dawn_plot_lol.png")
savefig(dawn_plot_lol, "main/images/dawn_plot_lol.pdf")


L, MLT, Kp, dir, scenario, colour, label = test_cases[2,:]
@time @load "$dir$scenario.jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 11, 6)
norm = normalize_to(97, energy, sim_plot_E_bins, hl_dawn_md, sim_ratio_sm)
b1 = 1/(norm * maximum(sim_ratio_sm))
normalizer = b1*norm

E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("HI", "DAWN", 1, "const", "L")
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("HI", "DAWN", 2, "const", "L")
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("HI", "DAWN", 4, "const", "L")
dc_norm = normalizer * normalize_to(96, E_bins, E, sim_ratio_sm, prec_ratio_wna1)

dawn_plot_hil = plot(xscale=:log10, yscale=:log10, xlim=(80,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true)
dawn_plot_hil = plot!(energy, b1*hl_dawn_md, fillrange=b1*hl_dawn_hi, fillalpha = 0.2, color = black, label=false)
dawn_plot_hil = plot!(energy, b1*hl_dawn_md, fillrange=b1*hl_dawn_lo, fillalpha = 0.2, color = black, label=false)
dawn_plot_hil = plot!(energy, b1*hl_dawn_md, label = "ELFIN Dawn-Noon: L>5, 4<MLT<13 ", color = black, linewidth=2, markershape=:circle);
dawn_plot_hil = plot!(sim_plot_E_bins, normalizer .* sim_ratio_sm, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
dawn_plot_hil = plot!(E, dc_norm * prec_ratio_wna1, label = "WNA1: Parallel Waves", color = blue, linestyle=:dash, linewidth=2)
dawn_plot_hil = plot!(E, dc_norm * prec_ratio_wna2, label = "WNA2: initially parallel -> resonance cone", color = green, linestyle=:dash, linewidth=2)
dawn_plot_hil = plot!(E, dc_norm * prec_ratio_wna3, label = "WNA3: gendrin angle -> resonance cone", color = red, linestyle=:dash, linewidth=2)
dawn_plot_hil = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
dawn_plot_hil = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(dawn_plot_hil, "main/images/dawn_plot_hil.png")
savefig(dawn_plot_hil, "main/images/dawn_plot_hil.pdf")

# Dusk plot
test_cases = [4.5 16.5 3 "main/results/original/" "LO_DUSK_MODEL_1m" purple "Noon-Dusk";
              6.5 16.5 3 "main/results/original/" "HI_DUSK_MODEL_1m" purple "Noon-Dusk";
              ]



L, MLT, Kp, dir, scenario, colour, label = test_cases[1,:]
@time @load "$dir$scenario.jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 11, 6)
norm = normalize_to(97, energy, sim_plot_E_bins, ll_dusk_md, sim_ratio_sm)
b1 = 1/(norm * maximum(sim_ratio_sm))
normalizer = b1*norm

E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("LO", "DUSK", 1, "const", "L")
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("LO", "DUSK", 2, "const", "L")
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("LO", "DUSK", 4, "const", "L")
dc_norm = normalizer * normalize_to(96, E_bins, E, sim_ratio_sm, prec_ratio_wna1)

dusk_plot_lol = plot(xscale=:log10, yscale=:log10, xlim=(80,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true)
dusk_plot_lol = plot!(energy, b1*ll_dusk_md, fillrange=b1*ll_dusk_hi, fillalpha = 0.2, color = black, label=false)
dusk_plot_lol = plot!(energy, b1*ll_dusk_md, fillrange=b1*ll_dusk_lo, fillalpha = 0.2, color = black, label=false)
dusk_plot_lol = plot!(energy, b1*ll_dusk_md, label = "ELFIN Noon-Dusk: L<5, 13<MLT<18", color = black, linewidth=2, markershape=:circle);
dusk_plot_lol = plot!(sim_plot_E_bins, normalizer .* sim_ratio_sm, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
dusk_plot_lol = plot!(E, dc_norm * prec_ratio_wna1, label = "WNA1: Parallel Waves", color = blue, linestyle=:dash, linewidth=2)
dusk_plot_lol = plot!(E, dc_norm * prec_ratio_wna2, label = "WNA2: initially parallel -> resonance cone", color = green, linestyle=:dash, linewidth=2)
dusk_plot_lol = plot!(E, dc_norm * prec_ratio_wna3, label = "WNA3: gendrin angle -> resonance cone", color = red, linestyle=:dash, linewidth=2)
dusk_plot_lol = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
dusk_plot_lol = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(dusk_plot_lol, "main/images/dusk_compare_lol.png")
savefig(dusk_plot_lol, "main/images/dusk_compare_lol.pdf")


L, MLT, Kp, dir, scenario, colour, label = test_cases[2,:]
# @time @load "$dir$scenario.jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 5, 5)
norm = normalize_to(97, energy, sim_plot_E_bins, hl_dusk_md, sim_ratio_sm)
b1 = 1/(norm * maximum(sim_ratio_sm))
normalizer = b1*norm

E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("HI", "DUSK", 1, "const", "L")
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("HI", "DUSK", 2, "const", "L")
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("HI", "DUSK", 4, "const", "L")
dc_norm = normalizer * normalize_to(96, E_bins, E, sim_ratio_sm, prec_ratio_wna1)

dusk_plot_hil = plot(xscale=:log10, yscale=:log10, xlim=(80,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true)
dusk_plot_hil = plot!(sim_plot_E_bins, normalizer .* sim_ratio_sm, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
dusk_plot_hil = plot!(energy, b1*hl_dusk_md, fillrange=b1*hl_dusk_hi, fillalpha = 0.2, color = black, label=false)
dusk_plot_hil = plot!(energy, b1*hl_dusk_md, fillrange=b1*hl_dusk_lo, fillalpha = 0.2, color = black, label=false)
dusk_plot_hil = plot!(energy, b1*hl_dusk_md, label = "ELFIN Noon-Dusk: L>5, 13<MLT<18", color = black, linewidth=2, markershape=:circle);
dusk_plot_hil = plot!(E, dc_norm * prec_ratio_wna1, label = "WNA1: Parallel Waves", color = blue, linestyle=:dash, linewidth=2)
dusk_plot_hil = plot!(E, dc_norm * prec_ratio_wna2, label = "WNA2: initially parallel -> resonance cone", color = green, linestyle=:dash, linewidth=2)
dusk_plot_hil = plot!(E, dc_norm * prec_ratio_wna3, label = "WNA3: gendrin angle -> resonance cone", color = red, linestyle=:dash, linewidth=2)
dusk_plot_hil = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
dusk_plot_hil = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(dusk_plot_hil, "main/images/dusk_compare_hil.png")
savefig(dusk_plot_hil, "main/images/dusk_compare_hil.pdf")

# Night plot
test_cases = [4.5 23.0 3 "main/results/original/" "LO_NITE_MODEL_1m" purple "Night";
              6.5 23.0 3 "main/results/original/" "HI_NITE_MODEL_1m" purple "Night";
              ]

L, MLT, Kp, dir, scenario, colour, label = test_cases[1,:]
@time @load "$dir$scenario.jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
norm = normalize_to(97, energy, sim_plot_E_bins, ll_nite_md, sim_ratio_sm)
b1 = 1/(norm * maximum(sim_ratio_sm))
normalizer = b1*norm

E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("LO", "NITE", 1, "const", "L")
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("LO", "NITE", 2, "const", "L")
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("LO", "NITE", 4, "const", "L")
dc_norm = normalizer * normalize_to(96, E_bins, E, sim_ratio_sm, prec_ratio_wna1)

nite_plot_lol = plot(xscale=:log10, yscale=:log10, xlim=(80,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true)
nite_plot_lol = plot!(energy, b1*ll_nite_md, fillrange=b1*ll_nite_hi, fillalpha = 0.2, color = black, label=false)
nite_plot_lol = plot!(energy, b1*ll_nite_md, fillrange=b1*ll_nite_lo, fillalpha = 0.2, color = black, label=false)
nite_plot_lol = plot!(energy, b1*ll_nite_md, label = "ELFIN Night: L<5, 18<MLT<4 ", color = black, linewidth=2, markershape=:circle);
nite_plot_lol = plot!(sim_plot_E_bins, normalizer .* sim_ratio_sm, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
nite_plot_lol = plot!(E, dc_norm * prec_ratio_wna1, label = "WNA1: Parallel Waves", color = blue, linestyle=:dash, linewidth=2)
nite_plot_lol = plot!(E, dc_norm * prec_ratio_wna2, label = "WNA2: initially parallel -> resonance cone", color = green, linestyle=:dash, linewidth=2)
nite_plot_lol = plot!(E, dc_norm * prec_ratio_wna3, label = "WNA3: gendrin angle -> resonance cone", color = red, linestyle=:dash, linewidth=2)
nite_plot_lol = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
nite_plot_lol = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(nite_plot_lol, "main/images/nite_plot_lol.png")
savefig(nite_plot_lol, "main/images/nite_plot_lol.pdf")


L, MLT, Kp, dir, scenario, colour, label = test_cases[2,:]
@time @load "$dir$scenario.jld2" rm
@info "loaded $scenario from $dir"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 11, 6)
norm = normalize_to(97, energy, sim_plot_E_bins, hl_nite_md, sim_ratio_sm)
b1 = 1/(norm * maximum(sim_ratio_sm))
normalizer = b1*norm

E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("HI", "NITE", 1, "const", "L")
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("HI", "NITE", 6, "const", "6.5")
E, Daa_wna3, prec_ratio_wna3 = obtain_diffusion_results("HI", "NITE", 4, "const", "L")
dc_norm = normalizer * normalize_to(96, E_bins, E, sim_ratio_sm, prec_ratio_wna1)

nite_plot_hil = plot(xscale=:log10, yscale=:log10, xlim=(80,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true)
nite_plot_hil = plot!(energy, b1*hl_nite_md, fillrange=b1*hl_nite_hi, fillalpha = 0.2, color = black, label=false)
nite_plot_hil = plot!(energy, b1*hl_nite_md, fillrange=b1*hl_nite_lo, fillalpha = 0.2, color = black, label=false)
nite_plot_hil = plot!(energy, b1*hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = black, linewidth=2, markershape=:circle);
nite_plot_hil = plot!(sim_plot_E_bins, normalizer .* sim_ratio_sm, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
nite_plot_hil = plot!(E, dc_norm * prec_ratio_wna1, label = "WNA1: Parallel Waves", color = blue, linestyle=:dash, linewidth=2)
nite_plot_hil = plot!(E, dc_norm * prec_ratio_wna2, label = "WNA2: initially parallel -> resonance cone", color = green, linestyle=:dash, linewidth=2)
nite_plot_hil = plot!(E, dc_norm * prec_ratio_wna3, label = "WNA3: gendrin angle -> resonance cone", color = red, linestyle=:dash, linewidth=2)
nite_plot_hil = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
nite_plot_hil = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(nite_plot_hil, "main/images/nite_plot_hil.png")
savefig(nite_plot_hil, "main/images/nite_plot_hil.pdf")

if save_tps_QLDC_statistical_comparison_plot
    tps_QLDC_statistical_comparison_plot = plot(dawn_plot_lol, dawn_plot_hil, dusk_plot_lol, dusk_plot_hil, nite_plot_lol, nite_plot_hil, layout=(3,2),size=(1600,1200),dpi=500)
    savefig(tps_QLDC_statistical_comparison_plot, "main/images/tps_QLDC_statistical_comparison.png")
    savefig(tps_QLDC_statistical_comparison_plot, "main/images/tps_QLDC_statistical_comparison.pdf")
end


#check didier's idea
E, Daa_wna1, prec_ratio_wna1 = obtain_diffusion_results("HI", "NITE", 1, "vary", "5") #
dc_norm_1 = normalize_to(97, energy, E, ll_nite_md, prec_ratio_wna1)
E, Daa_wna2, prec_ratio_wna2 = obtain_diffusion_results("Chorus_Daa_Precip_FreVary2_FpeFce5.0_L6.50_HI_NITE_WN1.txt") #
dc_norm_2 = normalize_to(97, energy, E, ll_nite_md, prec_ratio_wna2)

didiers_plot = plot(xscale=:log10, yscale=:log10, xlim=(80,1000), ylim=(1e-2, 2),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10, minorgrid=true);

didiers_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = black, label=false)
didiers_plot = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = black, label=false)
didiers_plot = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = black, linewidth=5, markershape=:circle);

didiers_plot = plot!(E, prec_ratio_wna1, label = "WNA1, Original ω_m(λ), Ω_pe = 5.0, Night, QLDC",  color = blue, linewidth=2, linestyle=:dash)
didiers_plot = plot!(E, prec_ratio_wna2, label = "WNA1, Didier's ω_m(λ), Ω_pe = 5.0, Night, QLDC",  color = reddish_orange, linewidth=2, linestyle=:dash)
didiers_plot = plot!(E, prec_ratio_wna2./prec_ratio_wna1, label = "Ratio of Orange to Blue",  color = green, linewidth=2, linestyle=:dash)

didiers_plot = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
didiers_plot = plot!(dpi = 500,size=(1000,600), margin=20px, bottom_margin=12px)
# savefig(didiers_plot, "main/images/comparison_didier.png")




include("agapitov_setup.jl")
include("particle_tracing.jl")
include("sim_setup.jl")
# run a sim w/ 500 particles with E=500 keV and PA=3 deg
ICrange = [100, 100, 1, 5, 80, 6]; 
sol_og = @time run_model(60, ICrange, 6.5, 23., 3., 0.35, 800, 1, ["og", 0, 300, 1., 40., 5]);
rm_og = @time sol2rm(sol_og, "thesis");
p1 = plot(21 .* rm_og.allT[indexrange], rm_og.allE[indexrange], legend=false, xlabel= "Time (ms)", ylabel = "Energy (keV)", lw = 1, xlim = [0,160], xminorticks = 5, yminorticks = 3, minorgrid=true);
p2 = plot(21 .* rm_og.allT[indexrange], rm_og.allPA[indexrange], legend=false, xlabel= "Time (ms)", ylabel = "Pitch angle α (°)", lw = 1, ylim = [0,90], xlim = [0,160], xminorticks = 5, yminorticks = 3, minorgrid=true);
p = plot(p1, p2, layout = (1,2),dpi = 500, tickfontsize=14, labelfontsize = 16, size=(1500,600), margin=20px, bottom_margin=40px, left_margin=40px)
savefig(p, "main/images/og_traj.pdf")
savefig(p, "main/images/og_traj.png")

sol_short = @time run_model(60, ICrange, 6.5, 23., 3., 0.35, 800, 1, ["og", 5, 30, 1., 40., 5]);
rm_short = @time sol2rm(sol_short, "thesis");
p1 = plot(21 .* rm_short.allT[indexrange], rm_short.allE[indexrange], legend=false, xlabel= "Time (ms)", ylabel = "Energy (keV)", lw = 1, xlim = [0,160], xminorticks = 5, yminorticks = 3, minorgrid=true);
p2 = plot(21 .* rm_short.allT[indexrange], rm_short.allPA[indexrange], legend=false, xlabel= "Time (ms)", ylabel = "Pitch angle α (°)", lw = 1, ylim = [0,90], xlim = [0,160], xminorticks = 5, yminorticks = 3, minorgrid=true);
p = plot(p1, p2, layout = (1,2),dpi = 500, tickfontsize=14, labelfontsize = 16, size=(1500,600), margin=20px, bottom_margin=40px, left_margin=40px)
savefig(p, "main/images/pack_traj.pdf")
savefig(p, "main/images/pack_traj.png")

sol_bw = @time run_model(60, ICrange, 6.5, 23., 3., 0.35, 800, 1, "bw");
rm_bw = @time sol2rm(sol_bw, "thesis");
p1 = plot(21 .* rm_bw.allT[indexrange], rm_bw.allE[indexrange], legend=false, xlabel= "Time (ms)", ylabel = "Energy (keV)", lw = 1, xlim = [0,160], xminorticks = 5, yminorticks = 3, minorgrid=true);
p2 = plot(21 .* rm_bw.allT[indexrange], rm_bw.allPA[indexrange], legend=false, xlabel= "Time (ms)", ylabel = "Pitch angle α (°)", lw = 1, ylim = [0,90], xlim = [0,160], xminorticks = 5, yminorticks = 3, minorgrid=true);
p = plot(p1, p2, layout = (1,2),dpi = 500, tickfontsize=14, labelfontsize = 16, size=(1500,600), margin=20px, bottom_margin=40px, left_margin=40px)
savefig(p, "main/images/bw_traj.pdf")

sol_wna3 = @time run_model(60, ICrange, 6.5, 23., 3., 0.35, 800, 1, "wna3");
rm_wna3 = @time sol2rm(sol_wna3, "thesis");
p1 = plot(21 .* rm_wna3.allT[indexrange], rm_wna3.allE[indexrange], legend=false, xlabel= "Time (ms)", ylabel = "Energy (keV)", lw = 1, xlim = [0,160], xminorticks = 5, yminorticks = 3, minorgrid=true);
p2 = plot(21 .* rm_wna3.allT[indexrange], rm_wna3.allPA[indexrange], legend=false, xlabel= "Time (ms)", ylabel = "Pitch angle α (°)", lw = 1, ylim = [0,90], xlim = [0,160], xminorticks = 5, yminorticks = 3, minorgrid=true);
p = plot(p1, p2, layout = (1,2),dpi = 500, tickfontsize=14, labelfontsize = 16, size=(1500,600), margin=20px, bottom_margin=40px, left_margin=40px)
savefig(p, "main/images/wna3_traj.pdf")


indexrange = 1:1:60
plot(rm_bw.allZ[indexrange], rm_bw.allPZ[indexrange], legend=false, xlim = [-1,1])


p1 = plot(rm_15.allT, rm_15.allE, legend=false, title = "Energies ω_m = 0.15", ylim = (492,508))
p2 = plot(rm_15.allT, rm_15.allPA, legend=false, title = "PA ω_m = 0.15", ylim = (2.85,3.15))
p3 = plot(rm_45.allT, rm_45.allE, legend=false, title = "Energies ω_m = 0.45", ylim = (492,508))
p4 = plot(rm_45.allT, rm_45.allPA, legend=false, title = "PA ω_m = 0.45", ylim = (2.85,3.15))

bigplot = plot(p1,p3,p2,p4,
        dpi = 96, layout = (2,2), size=(1000,1000),
        xtickfontsize=14, ytickfontsize=14, xguidefontsize=16, yguidefontsize=16, legendfontsize=10, titlefontsize=16)
