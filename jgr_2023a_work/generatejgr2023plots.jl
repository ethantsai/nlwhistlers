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
plot(elfin_p2t, yerror=elfin_p2t_error, color = bipride_pink, marker = stroke(3,bipride_pink), linewidth=4, markersize = 0, label=" ELFIN-$id $start_string-$stop_string" )
############
scenario = "ELA_ND_200904T0112"
start = DateTime(2020,9,4,1,12,24)
stop = DateTime(2020,9,4,1,12,33)
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
plot!(xscale=:log10, yscale=:log10, xlim=(52,1050), ylim=(1e-1, 4), xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)
elfin_plot = plot!(size=(800,450), dpi=300, title="ELFIN Precipitating Fluxes", titlefontsize=16,
xlabel="Energy (keV)", ylabel = "Precipitating to Trapped Electron Flux Ratio",legendfontsize=12, xtickfontsize=12, ytickfontsize=12,
xlabelfontsize=14, ylabelfontsize=14, margin=20px)

## want to remove that particular subspin from statistics and replot
# indices_to_remove = [72] # this one corresponds to microburst
# id, yy, mm, dd, HH, MM = extract_idyymmddHHMM(scenario);
# elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
#     mm*dd*yy*"_"*HH, start, stop, indices_to_remove
#     )
# index = findfirst(x->x==scenario, test_cases[:,4])
# L, MLT, Kp = test_cases[index,1:3]
# start_string = Dates.format(start, DateFormat("yyyy-mm-dd HH:MM:SS"))
# stop_string = Dates.format(stop, DateFormat("HH:MM:SS"))
# plot!(elfin_p2t, yerror=elfin_p2t_error, color = c4, marker = stroke(3,c4), linewidth=4, markersize = 0, label=" ELFIN-$id $start_string-$stop_string" )
# plot!(xscale=:log10, yscale=:log10, xlim=(52,1050), ylim=(1e-1, 4), xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)
# elfin_plot = plot!(size=(800,450), dpi=300, title="ELFIN Precipitating Fluxes", titlefontsize=16,
# xlabel="Energy (keV)", ylabel = "Precipitating to Trapped Electron Flux Ratio",legendfontsize=12, xtickfontsize=12, ytickfontsize=12,
# xlabelfontsize=14, ylabelfontsize=14, margin=20px, legend=false)

savefig(elfin_plot, "images/elfin_plot.pdf")




# Plot day and nightside comparisons
# table for paper
# 2a     & ELA & 21-02-03 20:46:37-20:46:50 & Dawn/Day & 6.2 & 5.8  & 3 & {\green \faCheckSquare} \\
# 2b     & ELA & 21-02-03 20:47:13-20:47:33 & Dawn/Day & 7.5 & 5.7  & 3 & {\green \faCheckSquare} \\
# 1a, 2c & ELB & 21-01-06 11:53:50-11:54:01 & Dawn/Day & 7.1 & 8.4  & 2 & {\green \faCheckSquare} \\
# 2d     & ELA & 21-10-01 08:10:02-08:10:13 & Dawn/Day & 6.6 & 8.8  & 3 & {\red \faTimes} \\
# 2e     & ELA & 21-10-01 04:59:57-05:01:18 & Dawn/Day & 6.4 & 10.2 & 4 & {\green \faCheckSquare} \\
# 2f     & ELA & 21-11-01 04:23:34-04:23:44 & Dawn/Day & 6.1 & 13.1 & 3 & {\green \faCheckSquare} \\
# 3a     & ELA & 21-01-11 17:50:50-17:50:58 & Night    & 4.8 & 19.0 & 3 & {\red \faQuestion} \\
# 3b     & ELA & 21-02-03 13:42:25-13:42:32 & Night    & 6.6 & 19.3 & 3 & {\red \faQuestion} \\
# 3c     & ELA & 21-02-03 09:01:43-09:01:56 & Night    & 5.1 & 20.2 & 2 & {\green \faCheckSquare} \\
# 3d     & ELA & 21-11-02 22:18:21-22:18:35 & Night    & 4.5 & 20.8 & 2 & {\red \faTimes} \\
# 3e     & ELB & 20-09-26 01:01:12-01:01:20 & Night    & 4.8 & 2.6  & 4 & {\red \faTimes} \\
# 1b, 3f & ELA & 20-09-04 01:12:24-01:12:33 & Night    & 6.5 & 3.8  & 3 & {\red \faTimes}

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
plot!(xscale=:log10, yscale=:log10, xlim=(1,10000), ylim=(1e-2, 1))
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
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 6))
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
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 6))
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
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 6))
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
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 6))
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
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 6))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 0, label=" ELFIN-$id $start_string-$stop_string" )
p6 = plot!()

day_plot = plot(p1,p2,p3,p4,p5,p6,layout=(2,3),# plot_title="Dayside Comparisons",
                    legendfontsize=12, xtickfontsize=12, ytickfontsize=12, linewidth=4,
                    size=(1600,1000), dpi=300, margin=8px, xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)
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
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 6))
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
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 6))
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
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 6))
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
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 6))
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
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 6))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 0, label=" ELFIN-$id $start_string-$stop_string" )
p5 = plot!()

############
scenario = "ELA_ND_200904T0112"
start = DateTime(2020,9,4,1,12,24)
stop = DateTime(2020,9,4,1,12,33)
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
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 6))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 0, label=" ELFIN-$id $start_string-$stop_string" )
p6 = plot!()

night_plot = plot(p1,p2,p3,p4,p5,p6,layout=(2,3),# plot_title="Nightside Comparisons",
                    legendfontsize=12, xtickfontsize=12, ytickfontsize=12, linewidth=4,
                    size=(1600,1000), dpi=300, margin=8px, xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)
# savefig(night_plot, "images/nightside_comparisons.png")
savefig(night_plot, "images/nightside_comparisons.pdf")


#####################
# Fig 4 B_w comapre #
#####################

# comparison of B_w with >40 = 0
############
scenario = "ELB_SA_210106T1154"
start = DateTime(2021,1,6,11,53,50)
stop = DateTime(2021,1,6,11,54,1)
############
id, yy, mm, dd, HH, MM = extract_idyymmddHHMM(scenario);
@time @load "result_matrix/"*scenario*".jld2" rm
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    mm*dd*yy*"_"*HH, start, stop
    )
@info "Model E_max = $(elfin_p2t[1][last_index_before_zero(elfin_p2t[2])]) keV for scenario $scenario"
normalizer = normalize_to_elfin(elfin_p2t[2], sim_ratio_sm)
@time @load "result_matrix_bw_compare/"*scenario*"_bwmod.jld2" rm
sim_ratio_new = prec_to_trap_ratio(rm)
sim_ratio_new_sm = smooth(sim_ratio_new[1], 6, 5)
normalizer_new = normalize_to_elfin(elfin_p2t[2], sim_ratio_new_sm)
plot(E_bins, normalizer*sim_ratio_sm, label=scenario[9:end], color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(E_bins, normalizer_new*sim_ratio_new_sm, label=scenario[9:end]*"_bwmod", color = c2, marker = stroke(3,c2), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-"*id*" "*mm*"/"*dd*" "*HH*":"*MM)
plot!(title = mm*"/"*dd*" "*HH*":"*MM*" prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(plot2, "images/"*scenario*"_bw_modcomparison.png")

############
scenario = "ELB_ND_210108T0646"
start = DateTime(2021,1,8,6,46,49)
stop = DateTime(2021,1,8,6,46,56)
############
id, yy, mm, dd, HH, MM = extract_idyymmddHHMM(scenario);
@time @load "result_matrix/"*scenario*".jld2" rm
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    mm*dd*yy*"_"*HH, start, stop
    )
@info "Model E_max = $(elfin_p2t[1][last_index_before_zero(elfin_p2t[2])]) keV for scenario $scenario"
normalizer = normalize_to_elfin(elfin_p2t[2], sim_ratio_sm)
@time @load "result_matrix_bw_compare/"*scenario*"_bwmod.jld2" rm
sim_ratio_new = prec_to_trap_ratio(rm)
sim_ratio_new_sm = smooth(sim_ratio_new[1], 6, 5)
normalizer_new = normalize_to_elfin(elfin_p2t[2], sim_ratio_new_sm)
plot(E_bins, normalizer*sim_ratio_sm, label=scenario[9:end], color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(E_bins, normalizer_new*sim_ratio_new_sm, label=scenario[9:end]*"_bwmod", color = c2, marker = stroke(3,c2), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-"*id*" "*mm*"/"*dd*" "*HH*":"*MM)
plot!(title = mm*"/"*dd*" "*HH*":"*MM*" prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(plot2, "images/"*scenario*"_bw_modcomparison.png")



case = "ELB_SA_210106T1154"
@time @load "hires_result_matrix/$case.jld2" rm;
A = [maximum(rad2deg.(x)) for x in rmx.allLambda]
index_of_interest = findall(A .== maximum(A))
bw_compare = plot([rad2deg.(traj) for traj in rmx.allLambda[index_of_interest]], rmx.allBw[index_of_interest], color = c1, linewidth=4, label = "210106T1154_original")
case = "ELB_ND_210108T0646"
@time @load "hires_result_matrix/$case.jld2" rm;
A = [maximum(rad2deg.(x)) for x in rmx.allLambda]
index_of_interest = findall(A .== maximum(A))
bw_compare = plot!([rad2deg.(traj) for traj in rmx.allLambda[index_of_interest]], rmx.allBw[index_of_interest], color = c4, linewidth=4,  label = "210108T0646_original")
case = "ELB_SA_210106T1154_bwmod"
@time @load "hires_result_matrix/$case.jld2" rm;
A = [maximum(rad2deg.(x)) for x in rmx.allLambda]
index_of_interest = findall(A .== maximum(A))
bw_compare = plot!([rad2deg.(traj) for traj in rmx.allLambda[index_of_interest]], rmx.allBw[index_of_interest], color = c5, linewidth=4,  label = "210106T1154_bwmod")
case = "ELB_ND_210108T0646_bwmod"
@time @load "hires_result_matrix/$case.jld2" rm;
A = [maximum(rad2deg.(x)) for x in rmx.allLambda]
index_of_interest = findall(A .== maximum(A))
bw_compare = plot!([rad2deg.(traj) for traj in rmx.allLambda[index_of_interest]], rmx.allBw[index_of_interest], color = c3, linewidth=4, label = "210108T0646_bwmod")
bw_compare = plot!(title = "B_w profile comparison", xlabel = "Latitude (degrees)", ylabel = "B_w multiplier", dpi = 500, size=(600,300), margin=15px, bottom_margin=12px)
# savefig(bw_compare, "images/bwmod_comparison.png")





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


# low L shell plot
test_cases = [4.5 8.0  3  "LO_DAWN_MODEL" c1 "Dawn";
              4.5 16.5 3  "LO_DUSK_MODEL" c3 "Dusk";
              4.5 23.0 3  "LO_NITE_MODEL" c5 "Night";
              ]

energy = stats_1.en
nite_lo = stats_1."L<5_night"
nite_md = stats_2."L<5_night"
nite_hi = stats_3."L<5_night"
dawn_lo = stats_1."L<5_day"
dawn_md = stats_2."L<5_day"
dawn_hi = stats_3."L<5_day"
dusk_lo = stats_1."L<5_dusk"
dusk_md = stats_2."L<5_dusk"
dusk_hi = stats_3."L<5_dusk"

L, MLT, Kp, scenario, colour, label = test_cases[1,:]
@time @load "result_matrix_stats/"*scenario*".jld2" rm
@info "loaded $scenario.jld2"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
normalizer = dawn_md[1] / sim_ratio_sm[3] 

sp_ll = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)
for case in eachrow(test_cases)
    L, MLT, Kp, scenario, colour, label = case
    @time @load "result_matrix_stats/"*scenario*".jld2" rm
    @info "loaded $scenario.jld2"
    sim_ratio = prec_to_trap_ratio(rm)
    sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
    sp_ll = plot!(E_bins, normalizer .* sim_ratio_sm, label=" $label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
end
sp_ll = plot!(energy, dawn_md, label = "ELFIN L<5 Dawn", color = bipride_orange, linewidth=2, markershape=:circle);
sp_ll = plot!(energy, dawn_md, fillrange=dawn_hi, fillalpha = 0.2, color = c1, label=false)
sp_ll = plot!(energy, dawn_md, fillrange=dawn_lo, fillalpha = 0.2, color = c1, label=false)
sp_ll = plot!(energy, dusk_md, label = "ELFIN L<5 Dusk", color = bipride_pink, linewidth=2, markershape=:circle);
sp_ll = plot!(energy, dusk_md, fillrange=dusk_hi, fillalpha = 0.2, color = c4, label=false)
sp_ll = plot!(energy, dusk_md, fillrange=dusk_lo, fillalpha = 0.2, color = c4, label=false)
sp_ll = plot!(energy, nite_md, label = "ELFIN L<5 Night", color = bipride_blue, linewidth=2, markershape=:circle)
sp_ll = plot!(energy, nite_md, fillrange=nite_hi, fillalpha = 0.2, color = c5, label=false)
sp_ll = plot!(energy, nite_md, fillrange=nite_lo, fillalpha = 0.2, color = c5, label=false)
sp_ll = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px, legend=:bottomleft)

savefig(sp_ll, "images/ll_compare.png")
savefig(sp_ll, "images/ll_compare.pdf")




# high L shell plot
test_cases = [6.5 8.0  3  "HI_DAWN_MODEL" c1 "Dawn";
              6.5 16.5 3  "HI_DUSK_MODEL" c3 "Dusk";
              6.5 23.0 3  "HI_NITE_MODEL" c5 "Night";
              ]

energy = stats_1.en
nite_lo = stats_1."L>5_night"
nite_md = stats_2."L>5_night"
nite_hi = stats_3."L>5_night"
dawn_lo = stats_1."L>5_day"
dawn_md = stats_2."L>5_day"
dawn_hi = stats_3."L>5_day"
dusk_lo = stats_1."L>5_dusk"
dusk_md = stats_2."L>5_dusk"
dusk_hi = stats_3."L>5_dusk"

L, MLT, Kp, scenario, colour, label = test_cases[1,:]
@time @load "result_matrix_stats/"*scenario*".jld2" rm
@info "loaded $scenario.jld2"
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
normalizer = dawn_md[1] / sim_ratio_sm[3] 

sp_hl = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)
for case in eachrow(test_cases)
    L, MLT, Kp, scenario, colour, label = case
    @time @load "result_matrix_stats/"*scenario*".jld2" rm
    @info "loaded $scenario.jld2"
    sim_ratio = prec_to_trap_ratio(rm)
    sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
    sp_hl = plot!(E_bins, normalizer .* sim_ratio_sm, label=" $label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
end
sp_hl = plot!(energy, dawn_md, label = "ELFIN L>5 Dawn", color = bipride_orange, linewidth=2, markershape=:circle);
sp_hl = plot!(energy, dawn_md, fillrange=dawn_hi, fillalpha = 0.2, color = c1, label=false)
sp_hl = plot!(energy, dawn_md, fillrange=dawn_lo, fillalpha = 0.2, color = c1, label=false)
sp_hl = plot!(energy, dusk_md, label = "ELFIN L>5 Dusk", color = bipride_pink, linewidth=2, markershape=:circle);
sp_hl = plot!(energy, dusk_md, fillrange=dusk_hi, fillalpha = 0.2, color = c4, label=false)
sp_hl = plot!(energy, dusk_md, fillrange=dusk_lo, fillalpha = 0.2, color = c4, label=false)
sp_hl = plot!(energy, nite_md, label = "ELFIN L>5 Night", color = bipride_blue, linewidth=2, markershape=:circle)
sp_hl = plot!(energy, nite_md, fillrange=nite_hi, fillalpha = 0.2, color = c5, label=false)
sp_hl = plot!(energy, nite_md, fillrange=nite_lo, fillalpha = 0.2, color = c5, label=false)
sp_hl = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px, legend=:bottomleft)

savefig(sp_hl, "images/hl_compare.png")
savefig(sp_hl, "images/hl_compare.pdf")


### MLT comparison, FIGURE 4 in JGR2023

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
dawn_plot_lol = plot!(E_bins, normalizer .* sim_ratio_sm, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
dawn_plot_lol = plot!(energy, b1*ll_dawn_md, fillrange=b1*ll_dawn_hi, fillalpha = 0.2, color = c5, label=false)
dawn_plot_lol = plot!(energy, b1*ll_dawn_md, fillrange=b1*ll_dawn_lo, fillalpha = 0.2, color = c5, label=false)
dawn_plot_lol = plot!(energy, b1*ll_dawn_md, label = "ELFIN Dawn-Noon: L<5, 4<MLT<13 ", color = c6, linewidth=2, markershape=:circle);
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
dawn_plot_hil = plot!(E_bins, normalizer .* sim_ratio_sm, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
dawn_plot_hil = plot!(energy, b1*hl_dawn_md, fillrange=b1*hl_dawn_hi, fillalpha = 0.2, color = c4, label=false)
dawn_plot_hil = plot!(energy, b1*hl_dawn_md, fillrange=b1*hl_dawn_lo, fillalpha = 0.2, color = c4, label=false)
dawn_plot_hil = plot!(energy, b1*hl_dawn_md, label = "ELFIN Dawn-Noon: L>5, 4<MLT<13 ", color = bipride_pink, linewidth=2, markershape=:circle);
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
dusk_plot_lol = plot!(E_bins, normalizer .* sim_ratio_sm, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
dusk_plot_lol = plot!(energy, b1*ll_dusk_md, fillrange=b1*ll_dusk_hi, fillalpha = 0.2, color = c5, label=false)
dusk_plot_lol = plot!(energy, b1*ll_dusk_md, fillrange=b1*ll_dusk_lo, fillalpha = 0.2, color = c5, label=false)
dusk_plot_lol = plot!(energy, b1*ll_dusk_md, label = "ELFIN Noon-Dusk: L<5, 13<MLT<18", color = c6, linewidth=2, markershape=:circle);
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
dusk_plot_hil = plot!(E_bins, normalizer .* sim_ratio_sm, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
dusk_plot_hil = plot!(energy, b1*hl_dusk_md, fillrange=b1*hl_dusk_hi, fillalpha = 0.2, color = c4, label=false)
dusk_plot_hil = plot!(energy, b1*hl_dusk_md, fillrange=b1*hl_dusk_lo, fillalpha = 0.2, color = c4, label=false)
dusk_plot_hil = plot!(energy, b1*hl_dusk_md, label = "ELFIN Noon-Dusk: L>5, 13<MLT<18", color = bipride_pink, linewidth=2, markershape=:circle);
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
nite_plot_lol = plot!(E_bins, normalizer .* sim_ratio_sm, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
nite_plot_lol = plot!(energy, b1*ll_nite_md, fillrange=b1*ll_nite_hi, fillalpha = 0.2, color = c5, label=false)
nite_plot_lol = plot!(energy, b1*ll_nite_md, fillrange=b1*ll_nite_lo, fillalpha = 0.2, color = c5, label=false)
nite_plot_lol = plot!(energy, b1*ll_nite_md, label = "ELFIN Night: L<5, 18<MLT<4 ", color = c6, linewidth=2, markershape=:circle);
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
nite_plot_hil = plot!(E_bins, normalizer .* sim_ratio_sm, label="$label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
nite_plot_hil = plot!(energy, b1*hl_nite_md, fillrange=b1*hl_nite_hi, fillalpha = 0.2, color = c4, label=false)
nite_plot_hil = plot!(energy, b1*hl_nite_md, fillrange=b1*hl_nite_lo, fillalpha = 0.2, color = c4, label=false)
nite_plot_hil = plot!(energy, b1*hl_nite_md, label = "ELFIN Night: L>5, 18<MLT<4 ", color = bipride_pink, linewidth=2, markershape=:circle);
nite_plot_hil = plot!(legendfontsize=12, tickfontsize=12, legend=:bottomleft)
nite_plot_hil = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(nite_plot_hil, "images/nite_plot_hil.png")
savefig(nite_plot_hil, "images/nite_plot_hil.pdf")

statistical_comparison_plot = plot(dawn_plot_lol, dawn_plot_hil, dusk_plot_lol, dusk_plot_hil, nite_plot_lol, nite_plot_hil, layout=(3,2),size=(1200,1600),dpi=500)
savefig(statistical_comparison_plot, "images/statistical_comparison.png")
savefig(statistical_comparison_plot, "images/statistical_comparison.pdf")




# elfin_only_plot --> EOP

# low L shell plot
test_cases = [4.5 8.0  3  "LO_DAWN_MODEL" c1 "Dawn";
              4.5 16.5 3  "LO_DUSK_MODEL" c3 "Dusk";
              4.5 23.0 3  "LO_NITE_MODEL" c5 "Night";
              ]

energy = stats_1.en
nite_lo = stats_1."L<5_night"
nite_md = stats_2."L<5_night"
nite_hi = stats_3."L<5_night"
dawn_lo = stats_1."L<5_day"
dawn_md = stats_2."L<5_day"
dawn_hi = stats_3."L<5_day"
dusk_lo = stats_1."L<5_dusk"
dusk_md = stats_2."L<5_dusk"
dusk_hi = stats_3."L<5_dusk"


eop = plot(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 1),
            xticks=([100, 1000], [100, 1000]), xminorticks=10, yminorticks=10)
for case in eachrow(test_cases)
    L, MLT, Kp, scenario, colour, label = case
    @time @load "result_matrix_stats/"*scenario*".jld2" rm
    @info "loaded $scenario.jld2"
    sim_ratio = prec_to_trap_ratio(rm)
    sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
    eop = plot!(E_bins, normalizer .* sim_ratio_sm, label=" $label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
end
eop = plot!(energy, dawn_md, label = "ELFIN L<5 Dawn", color = bipride_orange, linewidth=2, markershape=:circle);
eop = plot!(energy, dawn_md, fillrange=dawn_hi, fillalpha = 0.2, color = c1, label=false)
eop = plot!(energy, dawn_md, fillrange=dawn_lo, fillalpha = 0.2, color = c1, label=false)
eop = plot!(energy, dusk_md, label = "ELFIN L<5 Dusk", color = bipride_pink, linewidth=2, markershape=:circle);
eop = plot!(energy, dusk_md, fillrange=dusk_hi, fillalpha = 0.2, color = c4, label=false)
eop = plot!(energy, dusk_md, fillrange=dusk_lo, fillalpha = 0.2, color = c4, label=false)
eop = plot!(energy, nite_md, label = "ELFIN L<5 Night", color = bipride_blue, linewidth=2, markershape=:circle)
eop = plot!(energy, nite_md, fillrange=nite_hi, fillalpha = 0.2, color = c5, label=false)
eop = plot!(energy, nite_md, fillrange=nite_lo, fillalpha = 0.2, color = c5, label=false)
eop = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px, legend=:bottomleft)

# high L shell plot
test_cases = [6.5 8.0  3  "HI_DAWN_MODEL" c1 "Dawn";
              6.5 16.5 3  "HI_DUSK_MODEL" c3 "Dusk";
              6.5 23.0 3  "HI_NITE_MODEL" c5 "Night";
              ]

energy = stats_1.en
nite_lo = stats_1."L>5_night"
nite_md = stats_2."L>5_night"
nite_hi = stats_3."L>5_night"
dawn_lo = stats_1."L>5_day"
dawn_md = stats_2."L>5_day"
dawn_hi = stats_3."L>5_day"
dusk_lo = stats_1."L>5_dusk"
dusk_md = stats_2."L>5_dusk"
dusk_hi = stats_3."L>5_dusk"

for case in eachrow(test_cases)
    L, MLT, Kp, scenario, colour, label = case
    @time @load "result_matrix_stats/"*scenario*".jld2" rm
    @info "loaded $scenario.jld2"
    sim_ratio = prec_to_trap_ratio(rm)
    sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
    eop = plot!(E_bins, normalizer .* sim_ratio_sm, label=" $label Model: L=$L, MLT=$MLT", color = colour, marker = stroke(3,colour), linewidth=4, markersize = 3)
end
eop = plot!(energy, dawn_md, label = "ELFIN L>5 Dawn", color = bipride_orange, linewidth=2, markershape=:circle);
eop = plot!(energy, dawn_md, fillrange=dawn_hi, fillalpha = 0.2, color = c1, label=false)
eop = plot!(energy, dawn_md, fillrange=dawn_lo, fillalpha = 0.2, color = c1, label=false)
eop = plot!(energy, dusk_md, label = "ELFIN L>5 Dusk", color = bipride_pink, linewidth=2, markershape=:circle);
eop = plot!(energy, dusk_md, fillrange=dusk_hi, fillalpha = 0.2, color = c4, label=false)
eop = plot!(energy, dusk_md, fillrange=dusk_lo, fillalpha = 0.2, color = c4, label=false)
eop = plot!(energy, nite_md, label = "ELFIN L>5 Night", color = bipride_blue, linewidth=2, markershape=:circle)
eop = plot!(energy, nite_md, fillrange=nite_hi, fillalpha = 0.2, color = c5, label=false)
eop = plot!(energy, nite_md, fillrange=nite_lo, fillalpha = 0.2, color = c5, label=false)
eop = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px, legend=:bottomleft)



############
start = DateTime(2021,1,12,2,26,27)
stop = DateTime(2021,1,12,2,26,40)
############
id = "A"
yy = "21"
mm = "01"
dd = "12"
HH = "02"
MM = "27"
L, MLT, Kp = 6.0, 8.4, 3
start_string = Dates.format(start, DateFormat("yyyy-mm-dd HH:MM:SS"))
stop_string = Dates.format(stop, DateFormat("HH:MM:SS"))
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    mm*dd*yy*"_"*HH, start, stop
    )
plot(xscale=:log10, yscale=:log10, xlim=(50,1000), ylim=(1e-2, 1),  xminorticks=10, yminorticks=10)
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 0, label=" ELFIN-$id $start_string-$stop_string" )
elfun_210112 = plot!()
# savefig(elfun_210112, "images/elfun_210112.png")
savefig(elfun_210112, "images/elfun_210112.pdf")


# FIGURE 5 in JGR2023

using StatsPlots

# AE<100 L<5 night: Norbits=544.0 Nspins=4210.0
# AE<100 L<5 day: Norbits=276.0 Nspins=2295.0
# AE<100 L<5 dusk: Norbits=190.0 Nspins=2026.0
# AE<100 L>5 night: Norbits=601.0 Nspins=3105.0
# AE<100 L>5 day: Norbits=319.0 Nspins=2587.0
# AE<100 L>5 dusk: Norbits=217.0 Nspins=2149.0
# AE[100,350] L<5 night: Norbits=1182.0 Nspins=8520.0
# AE[100,350] L<5 day: Norbits=552.0 Nspins=4995.0
# AE[100,350] L<5 dusk: Norbits=549.0 Nspins=5638.0
# AE[100,350] L>5 night: Norbits=1232.0 Nspins=6744.0
# AE[100,350] L>5 day: Norbits=607.0 Nspins=6671.0
# AE[100,350] L>5 dusk: Norbits=581.0 Nspins=5752.0
# AE>350 L<5 night: Norbits=448.0 Nspins=2886.0
# AE>350 L<5 day: Norbits=171.0 Nspins=2108.0
# AE>350 L<5 dusk: Norbits=164.0 Nspins=1751.0
# AE>350 L>5 night: Norbits=451.0 Nspins=2485.0
# AE>350 L>5 day: Norbits=181.0 Nspins=2684.0
# AE>350 L>5 dusk: Norbits=163.0 Nspins=1533.0

N_orbits_day = [276,319,552,607,171,181]
N_orbits_dusk = [190,217,549,581,164,163]
N_orbits_night = [544,601,1182,1232,448,451]

N_spins_day = [2295,2587,4995,6671,2108,2684]
N_spins_dusk = [2026,2149,5638,5752,1751,1533]
N_spins_night = [4210,3105,8520,6744,2886,2485]

theme(:bright)
ticklabel = ["L<5", "L>5", "L<5", "L>5", "L<5", "L>5"]
n_orb_plot = groupedbar([N_orbits_day N_orbits_dusk N_orbits_night],
        bar_position = :stack,
        bar_width=0.3,
        xticks=(1:6, ticklabel), label=false, ylim = [0,2500])
savefig(n_orb_plot, "images/n_orb_plot.pdf")

n_spin_plot = groupedbar([N_spins_day/8 N_spins_dusk/8 N_spins_night/8],
        bar_position = :stack,
        bar_width=0.3,
        xticks=(1:6, ticklabel), label=false, ylim = [0,2500])
savefig(n_spin_plot, "images/n_spin_plot.pdf")