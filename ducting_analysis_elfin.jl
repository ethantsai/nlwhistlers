include("agapitovHelpers.jl")
include("agapitovmodel.jl")
include("agapitovPlotHelpers.jl")

#######################
## Constants n stuff ##
#######################
@info "Loading constants..."
save_dir = "results_ducting/"
folder = "run23/"

@info "E_min is $ELo keV"
@info "E_max is $EHi keV"
# E_bins = logrange(ELo,EHi, Int64(Esteps))

# case = "ELB_SA_210106T1154_3_250000", done
# case = "ELB_ND_210108T0646_3_250000", done
# case = "ELA_SD_210111T1750_3_250000", done
# case = "ELA_NA_210112T0226_3_250000", done
# case = "ELA_ND_200904T0112_3_250000", done
# case = "ELB_ND_200926T0101_3_250000", done
# case = "ELA_SD_210203T0902_3_250000", done
# case = "ELA_SD_210203T1342_3_250000", done
# case = "ELA_NA_210203T2046_3_250000", done
# case = "ELA_NA_210203T2047_3_250000", done
# case = "ELA_SD_211001T0501_3_250000", done
# case = "ELA_SD_211001T0810_3_250000", done
# case = "ELA_SA_211101T0424_3_250000", done
# case = "ELA_SA_211102T2218_3_250000"

@time @load save_dir*folder*case*".jld2" sol;
rm = sol2rm(sol, case);
@time @save "result_matrix/"*case[1:18]*".jld2" rm

# case = "HI_DAWN_MODEL_3_1000000"
# @time @load save_dir*folder*case*".jld2" sol;
# rm = sol2rm(sol, case);
# @time @save "result_matrix_stats/hr_"*case[1:13]*".jld2" rm

# case = "HI_DUSK_MODEL_3_1000000"
# @time @load save_dir*folder*case*".jld2" sol;
# rm = sol2rm(sol, case);
# @time @save "result_matrix_stats/hr_"*case[1:13]*".jld2" rm

# case = "HI_NITE_MODEL_3_1000000"
# @time @load save_dir*folder*case*".jld2" sol;
# rm = sol2rm(sol, case);
# @time @save "result_matrix_stats/hr_"*case[1:13]*".jld2" rm

# case = "LO_DAWN_MODEL_3_1000000"
# @time @load save_dir*folder*case*".jld2" sol;
# rm = sol2rm(sol, case);
# @time @save "result_matrix_stats/hr_"*case[1:13]*".jld2" rm

# case = "LO_DUSK_MODEL_3_1000000"
# @time @load save_dir*folder*case*".jld2" sol;
# rm = sol2rm(sol, case);
# @time @save "result_matrix_stats/hr_"*case[1:13]*".jld2" rm

# case = "LO_NITE_MODEL_3_1000000"
# @time @load save_dir*folder*case*".jld2" sol;
# rm = sol2rm(sol, case);
# @time @save "result_matrix_stats/hr_"*case[1:13]*".jld2" rm

case = "HI_NITE_WNA2_3_250000"
@time @load save_dir*folder*case*".jld2" sol;
rm = sol2rm(sol, case);
@time @save "result_matrix_oblique/"*case[1:13]*".jld2" rm


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
plot(E_bins, normalizer*sim_ratio_sm, label=scenario[9:end]*"_model", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-"*id*" "*mm*"/"*dd*" "*HH*":"*MM)
plot!(title = mm*"/"*dd*" "*HH*":"*MM*" prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")


############
scenario = "ELB_ND_210108T0646"
start = DateTime(2021,1,8,6,46,49)
stop = DateTime(2021,1,8,6,46,56)
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
plot(E_bins, normalizer*sim_ratio_sm, label=scenario[9:end]*"_model", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-"*id*" "*mm*"/"*dd*" "*HH*":"*MM)
plot!(title = mm*"/"*dd*" "*HH*":"*MM*" prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")


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
plot(E_bins, normalizer*sim_ratio_sm, label=scenario[9:end]*"_model", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-"*id*" "*mm*"/"*dd*" "*HH*":"*MM)
plot!(title = mm*"/"*dd*" "*HH*":"*MM*" prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")


############
scenario = "ELA_NA_210112T0226"
start = DateTime(2021,1,12,2,26,28)
stop = DateTime(2021,1,12,2,26,38)
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
plot(E_bins, normalizer*sim_ratio_sm, label=scenario[9:end]*"_model", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-"*id*" "*mm*"/"*dd*" "*HH*":"*MM)
plot!(title = mm*"/"*dd*" "*HH*":"*MM*" prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")


############
scenario = "ELA_ND_200904T0112"
start = DateTime(2020,9,4,1,12,28)
stop = DateTime(2020,9,4,1,12,32)
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
plot(E_bins, normalizer*sim_ratio_sm, label=scenario[9:end]*"_model", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-"*id*" "*mm*"/"*dd*" "*HH*":"*MM)
plot!(title = mm*"/"*dd*" "*HH*":"*MM*" prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")


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
plot(E_bins, normalizer*sim_ratio_sm, label=scenario[9:end]*"_model", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-"*id*" "*mm*"/"*dd*" "*HH*":"*MM)
plot!(title = mm*"/"*dd*" "*HH*":"*MM*" prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")


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
plot(E_bins, normalizer*sim_ratio_sm, label=scenario[9:end]*"_model", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-"*id*" "*mm*"/"*dd*" "*HH*":"*MM)
plot!(title = mm*"/"*dd*" "*HH*":"*MM*" prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")


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
plot(E_bins, normalizer*sim_ratio_sm, label=scenario[9:end]*"_model", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-"*id*" "*mm*"/"*dd*" "*HH*":"*MM)
plot!(title = mm*"/"*dd*" "*HH*":"*MM*" prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")


############
scenario = "ELA_NA_210203T2046"
start = DateTime(2021,2,3,20,46,41)
stop = DateTime(2021,2,3,20,46,48)
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
plot(E_bins, normalizer*sim_ratio_sm, label=scenario[9:end]*"_model", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-"*id*" "*mm*"/"*dd*" "*HH*":"*MM)
plot!(title = mm*"/"*dd*" "*HH*":"*MM*" prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")


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
plot(E_bins, normalizer*sim_ratio_sm, label=scenario[9:end]*"_model", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-"*id*" "*mm*"/"*dd*" "*HH*":"*MM)
plot!(title = mm*"/"*dd*" "*HH*":"*MM*" prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")


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
plot(E_bins, normalizer*sim_ratio_sm, label=scenario[9:end]*"_model", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-"*id*" "*mm*"/"*dd*" "*HH*":"*MM)
plot!(title = mm*"/"*dd*" "*HH*":"*MM*" prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")


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
plot(E_bins, normalizer*sim_ratio_sm, label=scenario[9:end]*"_model", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-"*id*" "*mm*"/"*dd*" "*HH*":"*MM)
plot!(title = mm*"/"*dd*" "*HH*":"*MM*" prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")


############
scenario = "ELA_SD_211101T0424"
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
plot(E_bins, normalizer*sim_ratio_sm, label=scenario[9:end]*"_model", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-"*id*" "*mm*"/"*dd*" "*HH*":"*MM)
plot!(title = mm*"/"*dd*" "*HH*":"*MM*" prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")


############
scenario = "ELA_SD_211101T2218"
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
plot(E_bins, normalizer*sim_ratio_sm, label=scenario[9:end]*"_model", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-"*id*" "*mm*"/"*dd*" "*HH*":"*MM)
plot!(title = mm*"/"*dd*" "*HH*":"*MM*" prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=20px, bottom_margin=12px)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")








# dawn comparison
scenario = "ELA_SD_211001T0810"
@time @load "result_matrix/"*scenario*".jld2" rm
id, yy, mm, dd, HH, MM = extract_idyymmddHHMM(scenario)
sim_ratio = prec_to_trap_ratio(rm)
sim_ratio_sm = smooth(sim_ratio[1], 6, 5)

scenario = "HI_DAWN_3_250000"
@time @load "result_matrix/"*scenario*".jld2" rm
sim_ratio_2 = prec_to_trap_ratio(rm)
sim_ratio_sm_2 = smooth(sim_ratio_2[1], 6, 5)

plot(E_bins, sim_ratio_sm, label=scenario[9:end]*"_model", color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(52,1000), ylim=(1e-2, 15))
plot!(E_bins, 1.2*sim_ratio_sm_2, label=scenario[9:end]*"_omegap=10", color = c4, marker = stroke(3,c4), linewidth=2, markersize = 1)





