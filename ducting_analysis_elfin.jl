include("agapitovHelpers.jl")
include("agapitovmodel.jl")
using Plots.PlotMeasures


#######################
## Constants n stuff ##
#######################
@info "Loading constants..."
save_dir = "results_ducting/"
folder = "run18/"



test_cases = ["ELB_SA_210106T1154_3_124800",
              "ELB_ND_210108T0646_3_124800",
              "ELA_SD_210111T1750_3_124800",
              "ELA_NA_210112T0226_3_124800",
              "ELA_ND_200904T0112_3_124800",
              "ELB_ND_200926T0101_3_124800",
              "ELA_SD_210203T0902_3_124800",
              "ELA_SD_210203T1342_3_124800"
              ]

rm_array = Vector{Resultant_Matrix}()
# loading loop™
for case in test_cases
    for omega_m in omega_m_cases
        loadname = save_dir*folder*case*".jld2"
        @info "Loading solution from $loadname..."
        @time @load loadname sol
        push!(rm_array, sol2rm(sol, case));
        @info "Loaded $case."
    end
end

function prec_to_trap_ratio(rm::Resultant_Matrix)
    initial_E = [rm.allE[i][1] for i = 1:length(rm.allT)]
    final_E = [rm.allE[i][end] for i = 1:length(rm.allT)]
    initial_PA = [rm.allPA[i][1] for i = 1:length(rm.allT)]
    final_PA = [rm.allPA[i][end] for i = 1:length(rm.allT)]
    # E_prec = sort(Ematrix[1,findall(isnan,truncated_matrix[end,:])])
    # E_trap = sort(Ematrix[1,findall(!isnan,truncated_Ematrix[end,:])])
  
    multiplier = 3

    trap_range = findall(x->multiplier*(lossConeAngle+0.002)>x>(lossConeAngle+0.002), initial_PA)
    loss_range = findall(x->1<x<(lossConeAngle+0.002), final_PA)
    E_prec = sort(final_E[loss_range])
    E_trap = sort(final_E[trap_range])
  
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

results_array = [prec_to_trap_ratio(x) for x in rm_array]

using CSV, DataFrames

function extract_idl_ratio(
    time_name::String,
    data_name::String,
    error_name::String,
    ebin_name::String,
    start::DateTime, stop::DateTime)

    time_csv_name = "idl_csvs/"*time_name
    data_csv_name = "idl_csvs/"*data_name
    ebins_csv_name = "idl_csvs/"*ebin_name
    error_csv_name = "idl_csvs/"*error_name

    times_df =  CSV.File(time_csv_name; header=false, delim=',', types=Float64) |> DataFrame
    time = unix2datetime.(times_df.Column1)
    indices = findall((time.>start).&(time.<stop)) # these are the indices corresponding to the time range to sum over
    time_of_interest = time[indices]
    @info "Summing over $(time_of_interest[end]-time_of_interest[1])"

    # import ratio data 
    data_df  =  CSV.File(data_csv_name; header=false, delim=',', types=Float64) |> DataFrame
    data = [[data_df[row,col] for col in 1:length(data_df[1,:])] for row in 1:16]
    error_df  =  CSV.File(error_csv_name; header=false, delim=',', types=Float64) |> DataFrame
    error = [[error_df[row,col] for col in 1:length(error_df[1,:])] for row in 1:16]
    # get rid of NaNs and Infs
    data_of_interest = data
    error_of_interest = error
    for i = 1:16
        data_of_interest[i][findall(.!isfinite.(data[i]))] .= 0.0
        error_of_interest[i][findall(.!isfinite.(error[i]))] .= 0.0
    end
    avg_ratio = [mean(data_of_interest[energy][indices]) for energy in 1:16]


    error = [sqrt(mean((data_of_interest[energy][indices].*error_of_interest[energy][indices]).^2)) for energy in 1:16]

    ebins_df = CSV.File(ebins_csv_name; header=false, delim=',', types=Float64) |> DataFrame
    ebins = ebins_df.Column1

    return (ebins, avg_ratio), error
end

function extract_elfin_p2t_ratio(datestring, tstart, tend)
    elfin_p2t, elfin_p2t_error = extract_idl_ratio(datestring*"_time.csv", datestring*"_p2t.csv",
                                                    datestring*"_p2t_err.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                    tstart, tend); # time to sample from ELFIN measurements                                                                
    return elfin_p2t, elfin_p2t_error
end

function obtain_elfin_ratio(elfin_prec, elfin_trap)
    return elfin_prec[2]./elfin_trap[2]
end

function normalize_to_elfin(elfin_ratio, sim_ratio)
    return elfin_ratio[1] / sim_ratio[13]
end

function make_elfin_error_bars(elfin_data, elfin_error)
    lo = elfin_data[2] .- elfin_error
    hi = elfin_data[2] .+ elfin_error
    return lo, hi
end

#################
## Plot Colors ##
#################

using Plots.PlotMeasures
function hexcolor(r::UInt8, g::UInt8, b::UInt8)
    return RGB(Int(r)/255,Int(g)/255,Int(b)/255)
end
bipride_pink = RGB(234/255, 2/255, 112/255);
bipride_orange = RGB(255/255, 79/255, 0/255);
bipride_lavender = RGB(155/255, 79/255, 150/255);
bipride_blue = RGB(0/255, 56/255, 168/255);
c1 = hexcolor(0xff,0x4f,0x00)
c2 = hexcolor(0xf8,0x00,0x4d)
c3 = hexcolor(0xd3,0x00,0x7d)
c4 = hexcolor(0x8f,0x19,0x9f)
c5 = hexcolor(0x00,0x38,0xa8)


############
i = 1 # applies to elfin 1/6 and 1/12
scenario = "010621_11"
############
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    scenario, DateTime(2021,1,6,11,53,50),DateTime(2021,1,6,11,54,1)
    )
normalizer = normalize_to_elfin(elfin_p2t[2], results_array[i][1])
plot(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(10,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-B 1/6")
plot!(title = "1/6 prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=5mm, bottom_margin=3mm)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")


############
i = 2 # applies to 1/8
scenario = "010821_06"
############
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    scenario, DateTime(2021,1,8,6,46,49),DateTime(2021,1,08,6,46,56)
    )
normalizer = normalize_to_elfin(elfin_p2t[2], results_array[i][1])
plot(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(10,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-B 1/8")
plot!(title = "1/8 prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=5mm, bottom_margin=3mm)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")


############
i = 3 # applies to elfin 1/11
scenario = "011121_17"
############
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    scenario, DateTime(2021,1,11,17,50,50),DateTime(2021,1,11,17,50,58)
    )
normalizer = normalize_to_elfin(elfin_p2t[2], results_array[i][1])
plot(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(10,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-A 1/11")
plot!(title = "1/11 prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=5mm, bottom_margin=3mm)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")


############
i = 4 # applies to elfin 1/12
scenario = "011221_02"
############
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    scenario, DateTime(2021,1,12,2,26,28),DateTime(2021,1,12,2,26,38)
    )
normalizer = normalize_to_elfin(elfin_p2t[2], results_array[i][1])
plot(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(10,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-A 1/12")
plot!(title = "1/12 prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=5mm, bottom_margin=3mm)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")


############
i = 5 # applies to elfin 9/4
scenario = "090420_01"
############
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    scenario, DateTime(2020,9,4,1,12,28),DateTime(2020,9,4,1,12,32)
    )
normalizer = normalize_to_elfin(elfin_p2t[2], results_array[i][1])
plot(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(10,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-A 9/4")
plot!(title = "9/4 prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=5mm, bottom_margin=3mm)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")


############
i = 6 # applies to elfin 9/26
scenario = "092620_00"
############
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    scenario, DateTime(2020,9,26,1,1,12),DateTime(2020,9,26,1,1,20)
    )
normalizer = normalize_to_elfin(elfin_p2t[2], results_array[i][1])
plot(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(10,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-B 9/26")
plot!(title = "9/26 prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=5mm, bottom_margin=3mm)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")


############
i = 7 # applies to elfin 2/3 9:xx
scenario = "020321_09"
############
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    scenario, DateTime(2021,2,3,9,1,43),DateTime(2021,2,3,9,1,56)
    )
normalizer = normalize_to_elfin(elfin_p2t[2], results_array[i][1])
plot(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(10,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-A 2/3 09:01")
plot!(title = "2/3 prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=5mm, bottom_margin=3mm)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")



############
i = 8 # applies to elfin 2/3 13:xx
scenario = "020321_13"
############
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    scenario, DateTime(2021,2,3,13,42,25),DateTime(2021,2,3,13,42,32)
    )
normalizer = normalize_to_elfin(elfin_p2t[2], results_array[i][1])
plot(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(10,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-A 2/3 13:42")
plot!(title = "2/3 prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=5mm, bottom_margin=3mm)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")



############
i = 9 # applies to elfin 2/3 13:xx
scenario = "020321_20"
############
elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
    scenario, DateTime(2021,2,3,13,42,25),DateTime(2021,2,3,13,42,32)
    )
normalizer = normalize_to_elfin(elfin_p2t[2], results_array[i][1])
plot(E_bins, normalizer*results_array[i][1], label=rm_array[i].label, color = c1, marker = stroke(3,c1), linewidth=4, markersize = 3)
plot!(xscale=:log10, yscale=:log10, xlim=(10,1000), ylim=(1e-2, 15))
plot!(elfin_p2t, yerror=elfin_p2t_error, color = c5, marker = stroke(3,c5), linewidth=4, markersize = 3, label="ELFIN-A 2/3 13:42")
plot!(title = "2/3 prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(800,450), margin=5mm, bottom_margin=3mm)
savefig(plot2, "images/"*scenario*"_ratiocomparison.png")

