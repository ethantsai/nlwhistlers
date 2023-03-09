include("agapitovHelpers.jl")
include("agapitovmodel.jl")
include("agapitovPlotHelpers.jl")

test_cases = [7.1 8.4  3  "ELB_SA_210106T1154"; # ELB SA 01/06 11:54
              6.5 19.8 3  "ELB_ND_210108T0646"; # ELB ND 01/08 06:46
              4.8 19.0 3  "ELA_SD_210111T1750"; # ELA SD 01/11 17:50
              6   8.4  3  "ELA_NA_210112T0226"; # ELA NA 01/12 02:26
              6.5 3.8  3  "ELA_ND_200904T0112"; # ELA ND 09/04 01:12
              4.8 2.6  4  "ELB_ND_200926T0101"; # ELB ND 09/26 01:01
              5.1 20.2 3  "ELA_SD_210203T0902"; # ELA SD 02/03 09:02
              6.6 19.3 3  "ELA_SD_210203T1342"; # ELA SD 02/03 13:42
              6.2 5.8  3  "ELA_NA_210203T2046"; # ELA NA 02/03 20:46
              7.5 5.7  3  "ELA_NA_210203T2047"; # ELA NA 02/03 20:47
              6.4 10.2 4  "ELA_SD_211001T0459"; # ELA SD 10/01 05:01
              6.6 8.8  3  "ELA_SD_211001T0810"; # ELA SD 10/01 08:10
              6.1 13.1 3  "ELA_SA_211101T0424"; # ELA SA 11/01 04:24
              4.5 20.8 3  "ELA_SA_211102T2218"] # ELA SA 11/02 22:18


####################
# Fig 1 ELFIN demo #
####################

############
scenario = "ELB_SA_210106T1154"
start = DateTime(2021,1,6,11,53,50)
stop = DateTime(2021,1,6,11,54,1)
############
id, yy, mm, dd, HH, MM = extract_idyymmddHHMM(scenario);
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    mm*dd*yy*"_"*HH, start, stop
    )
index = findfirst(x->x==scenario, test_cases[:,4])
L, MLT, Kp = test_cases[index,1:3]
start_string = Dates.format(start, DateFormat("yyyy-mm-dd HH:MM:SS"))
stop_string = Dates.format(stop, DateFormat("HH:MM:SS"))
plot(elfin_p2t, yerror=elfin_p2t_error, color = c2, marker = stroke(3,c2), linewidth=4, markersize = 0, label=" ELFIN-$id $start_string-$stop_string" )
############
scenario = "ELA_ND_200904T0112"
start = DateTime(2020,9,4,1,12,20)
stop = DateTime(2020,9,4,1,12,35)
############
id, yy, mm, dd, HH, MM = extract_idyymmddHHMM(scenario);
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    mm*dd*yy*"_"*HH, start, stop
    )
index = findfirst(x->x==scenario, test_cases[:,4])
L, MLT, Kp = test_cases[index,1:3]
start_string = Dates.format(start, DateFormat("yyyy-mm-dd HH:MM:SS"))
stop_string = Dates.format(stop, DateFormat("HH:MM:SS"))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c4, marker = stroke(3,c4), linewidth=4, markersize = 0, label=" ELFIN-$id $start_string-$stop_string" )
plot!(xscale=:log10, yscale=:log10, xlim=(52,1050), ylim=(1e-1, 4), xticks=([100, 1000], [100, 1000]), xminorticks=10)
elfin_plot = plot!(size=(800,450), dpi=300, title="ELFIN Precipitating Fluxes", titlefontsize=16,
xlabel="Energy (keV)", ylabel = "Precipitating to Trapped Electron Flux Ratio",legendfontsize=12, xtickfontsize=12, ytickfontsize=12,
xlabelfontsize=14, ylabelfontsize=14, margin=20px)

savefig(elfin_plot, "images/elfin_plot.pdf")




# Plot day and nightside comparisons

#############
# Day Cases #
#############

############
scenario = "ELA_NA_210203T2046"
start = DateTime(2021,2,3,20,46,37)
stop = DateTime(2021,2,3,20,46,50)
############
@time @load "result_matrix/"*scenario*".jld2" rm
id, yy, mm, dd, HH, MM = extract_idyymmddHHMM(scenario);
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    mm*dd*yy*"_"*HH, start, stop
    )
@info "Model E_max = $(elfin_p2t[1][last_index_before_zero(elfin_p2t[2])]) keV for scenario $scenario"
normalizer = normalize_to_elfin(elfin_p2t[2], sim_ratio_sm)
index = findfirst(x->x==scenario, test_cases[:,4])
L, MLT, Kp = test_cases[index,1:3]
start_string = Dates.format(start, DateFormat("yyyy-mm-dd HH:MM:SS"))
stop_string = Dates.format(stop, DateFormat("HH:MM:SS"))
plot(E_bins, normalizer*sim_ratio_sm, label=" Model: L=$L, MLT=$MLT, Kp=$Kp", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 8))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 0, label=" ELFIN-$id $start_string-$stop_string" )
p1 = plot!()

############
scenario = "ELB_SA_210106T1154"
start = DateTime(2021,1,6,11,53,50)
stop = DateTime(2021,1,6,11,54,1)
############
@time @load "result_matrix/"*scenario*".jld2" rm
id, yy, mm, dd, HH, MM = extract_idyymmddHHMM(scenario);
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    mm*dd*yy*"_"*HH, start, stop
    )
@info "Model E_max = $(elfin_p2t[1][last_index_before_zero(elfin_p2t[2])]) keV for scenario $scenario"
normalizer = normalize_to_elfin(elfin_p2t[2], sim_ratio_sm)
index = findfirst(x->x==scenario, test_cases[:,4])
L, MLT, Kp = test_cases[index,1:3]
start_string = Dates.format(start, DateFormat("yyyy-mm-dd HH:MM:SS"))
stop_string = Dates.format(stop, DateFormat("HH:MM:SS"))
plot(E_bins, normalizer*sim_ratio_sm, label=" Model: L=$L, MLT=$MLT, Kp=$Kp", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 8))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 0, label=" ELFIN-$id $start_string-$stop_string" )
p2 = plot!()

############
scenario = "ELA_SD_211001T0459"
start = DateTime(2021,10,1,4,59,57)
stop = DateTime(2021,10,1,5,0,18)
############
@time @load "result_matrix/"*scenario*".jld2" rm
id, yy, mm, dd, HH, MM = extract_idyymmddHHMM(scenario)
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    mm*dd*yy*"_"*HH, start, stop
    )
@info "Model E_max = $(elfin_p2t[1][last_index_before_zero(elfin_p2t[2])]) keV for scenario $scenario"
normalizer = normalize_to_elfin(elfin_p2t[2], sim_ratio_sm)
index = findfirst(x->x==scenario, test_cases[:,4])
L, MLT, Kp = test_cases[index,1:3]
start_string = Dates.format(start, DateFormat("yyyy-mm-dd HH:MM:SS"))
stop_string = Dates.format(stop, DateFormat("HH:MM:SS"))
plot(E_bins, normalizer*sim_ratio_sm, label=" Model: L=$L, MLT=$MLT, Kp=$Kp", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 8))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 0, label=" ELFIN-$id $start_string-$stop_string" )
p3 = plot!()

############
scenario = "ELA_NA_210203T2047"
start = DateTime(2021,2,3,20,47,13)
stop = DateTime(2021,2,3,20,47,33)
############
@time @load "result_matrix/"*scenario*".jld2" rm
id, yy, mm, dd, HH, MM = extract_idyymmddHHMM(scenario)
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    mm*dd*yy*"_"*HH, start, stop
    )
@info "Model E_max = $(elfin_p2t[1][last_index_before_zero(elfin_p2t[2])]) keV for scenario $scenario"
normalizer = normalize_to_elfin(elfin_p2t[2], sim_ratio_sm)
index = findfirst(x->x==scenario, test_cases[:,4])
L, MLT, Kp = test_cases[index,1:3]
start_string = Dates.format(start, DateFormat("yyyy-mm-dd HH:MM:SS"))
stop_string = Dates.format(stop, DateFormat("HH:MM:SS"))
plot(E_bins, normalizer*sim_ratio_sm, label=" Model: L=$L, MLT=$MLT, Kp=$Kp", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 8))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 0, label=" ELFIN-$id $start_string-$stop_string" )
p4 = plot!()

############
scenario = "ELA_SD_211001T0810"
start = DateTime(2021,10,1,8,10,2)
stop = DateTime(2021,10,1,8,10,13)
############
@time @load "result_matrix/"*scenario*".jld2" rm
id, yy, mm, dd, HH, MM = extract_idyymmddHHMM(scenario)
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    mm*dd*yy*"_"*HH, start, stop
    )
@info "Model E_max = $(elfin_p2t[1][last_index_before_zero(elfin_p2t[2])]) keV for scenario $scenario"
normalizer = normalize_to_elfin(elfin_p2t[2], sim_ratio_sm)
index = findfirst(x->x==scenario, test_cases[:,4])
L, MLT, Kp = test_cases[index,1:3]
start_string = Dates.format(start, DateFormat("yyyy-mm-dd HH:MM:SS"))
stop_string = Dates.format(stop, DateFormat("HH:MM:SS"))
plot(E_bins, normalizer*sim_ratio_sm, label=" Model: L=$L, MLT=$MLT, Kp=$Kp", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 8))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 0, label=" ELFIN-$id $start_string-$stop_string" )
p5 = plot!()

############
scenario = "ELA_SA_211101T0424"
start = DateTime(2021,11,1,4,23,34)
stop = DateTime(2021,11,1,4,23,44)
############
@time @load "result_matrix/"*scenario*".jld2" rm
id, yy, mm, dd, HH, MM = extract_idyymmddHHMM(scenario)
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    mm*dd*yy*"_"*HH, start, stop
    )
@info "Model E_max = $(elfin_p2t[1][last_index_before_zero(elfin_p2t[2])]) keV for scenario $scenario"
normalizer = normalize_to_elfin(elfin_p2t[2], sim_ratio_sm)
index = findfirst(x->x==scenario, test_cases[:,4])
L, MLT, Kp = test_cases[index,1:3]
start_string = Dates.format(start, DateFormat("yyyy-mm-dd HH:MM:SS"))
stop_string = Dates.format(stop, DateFormat("HH:MM:SS"))
plot(E_bins, normalizer*sim_ratio_sm, label=" Model: L=$L, MLT=$MLT, Kp=$Kp", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 8))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 0, label=" ELFIN-$id $start_string-$stop_string" )
p6 = plot!()

day_plot = plot(p1,p2,p3,p4,p5,p6,layout=(2,3),# plot_title="Dayside Comparisons",
                    legendfontsize=12, xtickfontsize=12, ytickfontsize=12, linewidth=4,
                    size=(1600,1000), dpi=300, margin=8px, xticks=([100, 1000], [100, 1000]), xminorticks=10)
# savefig(day_plot, "images/dayside_comparisons.png")
savefig(day_plot, "images/dayside_comparisons.pdf")


              
###############
# Night Cases #
###############

############
scenario = "ELA_SD_210111T1750"
start = DateTime(2021,1,11,17,50,50)
stop = DateTime(2021,1,11,17,50,58)
############
@time @load "result_matrix/"*scenario*".jld2" rm
id, yy, mm, dd, HH, MM = extract_idyymmddHHMM(scenario);
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    mm*dd*yy*"_"*HH, start, stop
    )
@info "Model E_max = $(elfin_p2t[1][last_index_before_zero(elfin_p2t[2])]) keV for scenario $scenario"
normalizer = normalize_to_elfin(elfin_p2t[2], sim_ratio_sm)
index = findfirst(x->x==scenario, test_cases[:,4])
L, MLT, Kp = test_cases[index,1:3]
start_string = Dates.format(start, DateFormat("yyyy-mm-dd HH:MM:SS"))
stop_string = Dates.format(stop, DateFormat("HH:MM:SS"))
plot(E_bins, normalizer*sim_ratio_sm, label=" Model: L=$L, MLT=$MLT, Kp=$Kp", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 8))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 0, label=" ELFIN-$id $start_string-$stop_string" )
p1 = plot!()

############
scenario = "ELA_SD_210203T0902"
start = DateTime(2021,2,3,9,1,43)
stop = DateTime(2021,2,3,9,1,56)
############
@time @load "result_matrix/"*scenario*".jld2" rm
id, yy, mm, dd, HH, MM = extract_idyymmddHHMM(scenario);
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    mm*dd*yy*"_"*HH, start, stop
    )
@info "Model E_max = $(elfin_p2t[1][last_index_before_zero(elfin_p2t[2])]) keV for scenario $scenario"
normalizer = normalize_to_elfin(elfin_p2t[2], sim_ratio_sm)
index = findfirst(x->x==scenario, test_cases[:,4])
L, MLT, Kp = test_cases[index,1:3]
start_string = Dates.format(start, DateFormat("yyyy-mm-dd HH:MM:SS"))
stop_string = Dates.format(stop, DateFormat("HH:MM:SS"))
plot(E_bins, normalizer*sim_ratio_sm, label=" Model: L=$L, MLT=$MLT, Kp=$Kp", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 8))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 0, label=" ELFIN-$id $start_string-$stop_string" )
p2 = plot!()

############
scenario = "ELB_ND_200926T0101"
start = DateTime(2020,9,26,1,1,12)
stop = DateTime(2020,9,26,1,1,20)
############
@time @load "result_matrix/"*scenario*".jld2" rm
id, yy, mm, dd, HH, MM = extract_idyymmddHHMM(scenario);
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    mm*dd*yy*"_"*HH, start, stop
    )
@info "Model E_max = $(elfin_p2t[1][last_index_before_zero(elfin_p2t[2])]) keV for scenario $scenario"
normalizer = normalize_to_elfin(elfin_p2t[2], sim_ratio_sm)
index = findfirst(x->x==scenario, test_cases[:,4])
L, MLT, Kp = test_cases[index,1:3]
start_string = Dates.format(start, DateFormat("yyyy-mm-dd HH:MM:SS"))
stop_string = Dates.format(stop, DateFormat("HH:MM:SS"))
plot(E_bins, normalizer*sim_ratio_sm, label=" Model: L=$L, MLT=$MLT, Kp=$Kp", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 8))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 0, label=" ELFIN-$id $start_string-$stop_string" )
p3 = plot!()

############
scenario = "ELA_SD_210203T1342"
start = DateTime(2021,2,3,13,42,25)
stop = DateTime(2021,2,3,13,42,32)
############
@time @load "result_matrix/"*scenario*".jld2" rm
id, yy, mm, dd, HH, MM = extract_idyymmddHHMM(scenario);
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    mm*dd*yy*"_"*HH, start, stop
    )
@info "Model E_max = $(elfin_p2t[1][last_index_before_zero(elfin_p2t[2])]) keV for scenario $scenario"
normalizer = normalize_to_elfin(elfin_p2t[2], sim_ratio_sm)
index = findfirst(x->x==scenario, test_cases[:,4])
L, MLT, Kp = test_cases[index,1:3]
start_string = Dates.format(start, DateFormat("yyyy-mm-dd HH:MM:SS"))
stop_string = Dates.format(stop, DateFormat("HH:MM:SS"))
plot(E_bins, normalizer*sim_ratio_sm, label=" Model: L=$L, MLT=$MLT, Kp=$Kp", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 8))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 0, label=" ELFIN-$id $start_string-$stop_string" )
p4 = plot!()

############
scenario = "ELA_SA_211102T2218"
start = DateTime(2021,11,2,22,18,21)
stop = DateTime(2021,11,2,22,18,35)
############
@time @load "result_matrix/"*scenario*".jld2" rm
id, yy, mm, dd, HH, MM = extract_idyymmddHHMM(scenario)
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    mm*dd*yy*"_"*HH, start, stop
    )
@info "Model E_max = $(elfin_p2t[1][last_index_before_zero(elfin_p2t[2])]) keV for scenario $scenario"
normalizer = normalize_to_elfin(elfin_p2t[2], sim_ratio_sm)
index = findfirst(x->x==scenario, test_cases[:,4])
L, MLT, Kp = test_cases[index,1:3]
start_string = Dates.format(start, DateFormat("yyyy-mm-dd HH:MM:SS"))
stop_string = Dates.format(stop, DateFormat("HH:MM:SS"))
plot(E_bins, normalizer*sim_ratio_sm, label=" Model: L=$L, MLT=$MLT, Kp=$Kp", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 8))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 0, label=" ELFIN-$id $start_string-$stop_string" )
p5 = plot!()

############
scenario = "ELA_ND_200904T0112"
start = DateTime(2020,9,4,1,12,20)
stop = DateTime(2020,9,4,1,12,35)
############
@time @load "result_matrix/"*scenario*".jld2" rm
id, yy, mm, dd, HH, MM = extract_idyymmddHHMM(scenario);
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    mm*dd*yy*"_"*HH, start, stop
    )
@info "Model E_max = $(elfin_p2t[1][last_index_before_zero(elfin_p2t[2])]) keV for scenario $scenario"
normalizer = normalize_to_elfin(elfin_p2t[2], sim_ratio_sm)
index = findfirst(x->x==scenario, test_cases[:,4])
L, MLT, Kp = test_cases[index,1:3]
start_string = Dates.format(start, DateFormat("yyyy-mm-dd HH:MM:SS"))
stop_string = Dates.format(stop, DateFormat("HH:MM:SS"))
plot(E_bins, normalizer*sim_ratio_sm, label=" Model: L=$L, MLT=$MLT, Kp=$Kp", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 8))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 0, label=" ELFIN-$id $start_string-$stop_string" )
p6 = plot!()

night_plot = plot(p1,p2,p3,p4,p5,p6,layout=(2,3),# plot_title="Nightside Comparisons",
                    legendfontsize=12, xtickfontsize=12, ytickfontsize=12, linewidth=4,
                    size=(1600,1000), dpi=300, margin=8px, xticks=([100, 1000], [100, 1000]), xminorticks=10)
# savefig(night_plot, "images/nightside_comparisons.png")
savefig(night_plot, "images/nightside_comparisons.pdf")


