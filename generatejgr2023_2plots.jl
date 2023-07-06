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
stats_1 =  CSV.File("stats_csvs/$stats_1_csv_name"; header=true, delim=' ', types=Float64) |> DataFrame
stats_2 =  CSV.File("stats_csvs/$stats_2_csv_name"; header=true, delim=' ', types=Float64) |> DataFrame
stats_3 =  CSV.File("stats_csvs/$stats_3_csv_name"; header=true, delim=' ', types=Float64) |> DataFrame


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
              6.5 23.0 3  "HI_NITE_WNA1" c5 "WNA1";
              6.5 23.0 3  "HI_NITE_WNA1" c6 "WNA2";
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

L, MLT, Kp, scenario, colour, label = test_cases[3,:]
@time @load "result_matrix_oblique/"*scenario*".jld2" rm
@info "loaded $scenario.jld2"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm_3 = smooth(sim_ratio[1], 8, 5)
norm_3 = normalize_to_elfin(hl_nite_md, sim_ratio_sm_3)

nite_plot_hil = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)

L, MLT, Kp, scenario, colour, label = test_cases[1,:]
nite_plot_hil = plot!(E_bins, norm_1*sim_ratio_sm_1, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
L, MLT, Kp, scenario, colour, label = test_cases[2,:]
nite_plot_hil = plot!(E_bins, norm_2*sim_ratio_sm_2, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
L, MLT, Kp, scenario, colour, label = test_cases[3,:]
nite_plot_hil = plot!(E_bins, norm_3*sim_ratio_sm_3, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)

nite_plot_hil = plot!(energy, hl_nite_md, fillrange=hl_nite_hi, fillalpha = 0.2, color = c4, label=false)
nite_plot_hil = plot!(energy, hl_nite_md, fillrange=hl_nite_lo, fillalpha = 0.2, color = c4, label=false)
nite_plot_hil = plot!(energy, hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = bipride_pink, linewidth=2, markershape=:circle);
nite_plot_hil = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
nite_plot_hil = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)



