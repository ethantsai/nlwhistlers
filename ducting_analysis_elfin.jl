include("agapitovHelpers.jl")
include("agapitovmodel.jl")
using Plots.PlotMeasures


#######################
## Constants n stuff ##
#######################
@info "Loading constants..."
save_dir = "results_ducting/"
folder = "run18/"



test_cases = ["ELB_SA_210106T1154_3_124800", "ELB_ND_210108T0646_3_124800", "ELA_SD_210111T1750_3_124800"]

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


using CSV, DataFrames

function extract_idl_csv(
    time_name::String,
    data_name::String,
    error_name::String,
    error_time_name::String,
    ebin_name::String,
    start::DateTime, stop::DateTime)

    time_csv_name = "idl_csvs/"*time_name
    data_csv_name = "idl_csvs/"*data_name
    ebins_csv_name = "idl_csvs/"*ebin_name
    error_csv_name = "idl_csvs/"*error_name
    errortime_name = "idl_csvs/"*error_time_name

    times_df =  CSV.File(time_csv_name; header=false, delim=',', types=Float64) |> DataFrame
    time = unix2datetime.(times_df.Column1)
    indices = findall((time.>start).&(time.<stop)) # these are the indices corresponding to the time range to sum over
    time_of_interest = time[indices]
    @info "Summing over $(time_of_interest[end]-time_of_interest[1])"

    errtimes_df = CSV.File(errortime_name; header=false, delim=',', types=Float64) |> DataFrame
    errtime = unix2datetime.(errtimes_df.Column1)
    errindices = findall((errtime.>start).&(errtime.<stop)) # these are the indices corresponding to the time range to sum over
    if length(errindices) != length(indices)
        @error "Length mismatch between time and error time."
        return
    end

    # import particle flux data 
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
    flux = [sum(data_of_interest[energy][indices]) for energy in 1:16]


    error = [sqrt(sum((data_of_interest[energy][indices].*error_of_interest[energy][errindices]).^2)) for energy in 1:16]

    ebins_df = CSV.File(ebins_csv_name; header=false, delim=',', types=Float64) |> DataFrame
    ebins = ebins_df.Column1

    return (ebins, flux), error
end





results_array = [prec_to_trap_ratio(x) for x in rm_array]
############
scenario = 1 # applies to elfin 1/6 and 1/12
############
elfin_prec_010621, elfin_error_010621 = extract_idl_csv("010621_time.csv", "010621_prec.csv",
                                                "010621_precerror.csv", "010621_precerror_time.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2021,1,6,11,53,50),DateTime(2021,1,6,11,54,1)); # time to sample from ELFIN measurements                                                                

elfin_trap_010621, elfin_error_010621 = extract_idl_csv("010621_time.csv", "010621_perp.csv",
                                                "010621_precerror.csv", "010621_precerror_time.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2021,1,6,11,53,50),DateTime(2021,1,6,11,54,1)); # time to sample from ELFIN measurements                                                                

elfin_prec_011221, elfin_error_011221 = extract_idl_csv("011221_time.csv", "011221_prec.csv",
                                                "011221_precerror.csv", "011221_precerror_time.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2021,1,12,2,26,28),DateTime(2021,1,12,2,26,38)); # time to sample from ELFIN measurements                                                                

elfin_trap_011221, elfin_error_011221 = extract_idl_csv("011221_time.csv", "011221_perp.csv",
                                                "011221_precerror.csv", "011221_precerror_time.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2021,1,12,2,26,28),DateTime(2021,1,12,2,26,38)); # time to sample from ELFIN measurements                                                                


#ratio comparison
elfin_ratio_010621 = elfin_prec_010621[2]./elfin_trap_010621[2]
elfin_ratio_011221 = elfin_prec_011221[2]./elfin_trap_011221[2]
i = scenario
plot1 = plot(E_bins, results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(50,2000), ylim=(1e-2, 15))
plot!(elfin_prec_010621[1][1:end-5], elfin_ratio_010621[1:end-5], label="elfin-b 1/6")    
plot!(elfin_prec_011221[1][1:end-5], elfin_ratio_011221[1:end-5], label="elfin-a 1/12")    

#precipitating flux comparison
@info "elfin 1/6 @ $(elfin_prec_010621[1][1]) kev"
@info "elfin 1/12 @ $(elfin_prec_011221[1][1]) kev"
@info "sim @ $(E_bins[13]) kev"
# normalizer = mean([elfin_prec_010621[2][1], elfin_prec_010621[2][1]]) / results_array[1][2][13] # normalize max sim to measurements at 63 keV

# plot2 = plot(E_bins, normalizer*results_array[i][2], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(50,2000))
# plot!(elfin_prec_010621[1][1:9], elfin_prec_010621[2][1:9], label="elfin-b")    
# plot!(elfin_prec_010621[1][1:9], elfin_prec_010621[2][1:9], label="elfin-b")    
# plot!(ylim =(1e3,1e8))

ratio_normalizer = mean([elfin_ratio_010621[1], elfin_ratio_011221[1]]) ./ results_array[i][1][13]
plot2 = plot(E_bins, ratio_normalizer*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(50,2000), ylim=(1e-2, 2))
plot!(elfin_prec_010621[1][1:end-5], elfin_ratio_010621[1:end-5], label="elfin-b 1/6")
plot!(elfin_prec_011221[1][1:end-5], elfin_ratio_011221[1:end-5], label="elfin-a 1/12")
plot!(title = "[1/6,1/12] prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)")    
plot2 = plot!(dpi = 500,size=(1000,450), margin=5mm, bottom_margin=0mm)
savefig(plot2, "images/scenario1_ratiocomparison.png")


############
scenario = 2 # applies to elfin 1/5 and 1/8
############

elfin_prec_010521, elfin_error_010521 = extract_idl_csv("010521_time.csv", "010521_prec.csv",
                                                "010521_precerror.csv", "010521_precerror_time.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2021,1,5,14,57,44),DateTime(2021,1,5,14,57,49)); # time to sample from ELFIN measurements                                                                

elfin_trap_010521, elfin_error_010521 = extract_idl_csv("010521_time.csv", "010521_perp.csv",
                                                "010521_precerror.csv", "010521_precerror_time.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2021,1,5,14,57,44),DateTime(2021,1,5,14,57,49)); # time to sample from ELFIN measurements                                                                

elfin_prec_010821, elfin_error_010821 = extract_idl_csv("010821_time.csv", "010821_prec.csv",
                                                "010821_precerror.csv", "010821_precerror_time.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2021,1,8,6,46,49),DateTime(2021,1,08,6,46,56)); # time to sample from ELFIN measurements                                                                

elfin_trap_010821, elfin_error_010821 = extract_idl_csv("010821_time.csv", "010821_perp.csv",
                                                "010821_precerror.csv", "010821_precerror_time.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2021,1,8,6,46,49),DateTime(2021,1,08,6,46,56)); # time to sample from ELFIN measurements                                                                


#ratio comparison
elfin_ratio_010521 = elfin_prec_010521[2]./elfin_trap_010521[2]
elfin_ratio_010821 = elfin_prec_010821[2]./elfin_trap_010821[2]
i = scenario
plot1 = plot(E_bins, results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(50,2000), ylim=(1e-2, 15))
plot!(elfin_prec_010521[1][1:end-5], elfin_ratio_010521[1:end-5], label="elfin-a 1/5")    
plot!(elfin_prec_010821[1][1:end-5], elfin_ratio_010821[1:end-5], label="elfin-b 1/8")    

#precipitating flux comparison
# normalizer = mean([elfin_prec_010621[2][1], elfin_prec_010621[2][1]]) / results_array[1][2][13] # normalize max sim to measurements at 63 keV

# plot2 = plot(E_bins, normalizer*results_array[i][2], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(50,2000))
# plot!(elfin_prec_010621[1][1:9], elfin_prec_010621[2][1:9], label="elfin-b")    
# plot!(elfin_prec_010621[1][1:9], elfin_prec_010621[2][1:9], label="elfin-b")    
# plot!(ylim =(1e3,1e8))

ratio_normalizer = mean([elfin_ratio_010521[1], elfin_ratio_010821[1]]) ./ results_array[i][1][13]
plot2 = plot(E_bins, ratio_normalizer*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(50,2000), ylim=(1e-2, 2))
plot!(elfin_prec_010521[1][1:end-5], elfin_ratio_010521[1:end-5], label="elfin-a 1/5")
plot!(elfin_prec_010821[1][1:end-5], elfin_ratio_010821[1:end-5], label="elfin-b 1/8")
plot!(title = "[1/5,1/8] prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)", legend=:outerbottom)    
plot2 = plot!(dpi = 500,size=(1000,450), margin=5mm, bottom_margin=0mm)
savefig(plot2, "images/scenario2_ratiocomparison.png")




############
scenario = 3 # applies to elfin 1/11
############

elfin_prec_011121, elfin_error_011121 = extract_idl_csv("011121_time.csv", "011121_prec.csv",
                                                "011121_precerror.csv", "011121_precerror_time.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2021,1,11,17,50,50),DateTime(2021,1,11,17,50,58)); # time to sample from ELFIN measurements                                                                

elfin_trap_011121, elfin_error_011121 = extract_idl_csv("011121_time.csv", "011121_perp.csv",
                                                "011121_precerror.csv", "011121_precerror_time.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                DateTime(2021,1,11,17,50,50),DateTime(2021,1,11,17,50,58)); # time to sample from ELFIN measurements                                                                

#ratio comparison
elfin_ratio_011121 = elfin_prec_011121[2]./elfin_trap_011121[2]
i = scenario
plot1 = plot(E_bins, results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(50,2000), ylim=(1e-2, 15))
plot!(elfin_prec_011121[1][1:end-5], elfin_ratio_011121[1:end-5], label="elfin-a 1/11")    

#precipitating flux comparison
# normalizer = mean([elfin_prec_010621[2][1], elfin_prec_010621[2][1]]) / results_array[1][2][13] # normalize max sim to measurements at 63 keV

# plot2 = plot(E_bins, normalizer*results_array[i][2], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(50,2000))
# plot!(elfin_prec_010621[1][1:9], elfin_prec_010621[2][1:9], label="elfin-b")    
# plot!(elfin_prec_010621[1][1:9], elfin_prec_010621[2][1:9], label="elfin-b")    
# plot!(ylim =(1e3,1e8))

ratio_normalizer = elfin_ratio_011121[1] ./ results_array[i][1][13]
plot2 = plot(E_bins, ratio_normalizer*results_array[i][1], label=rm_array[i].label, xscale=:log10, yscale=:log10, xlim=(50,2000), ylim=(1e-2, 2))
plot!(elfin_prec_011121[1][1:end-5], elfin_ratio_011121[1:end-5], label="elfin-a 1/11")
plot!(title = "[1/11] prec/trap flux ratio comparison", ylabel="j_parallel/j_perp", xlabel = "Energy (keV)", legend=:outerbottom)    
plot2 = plot!(dpi = 500,size=(1000,450), margin=5mm, bottom_margin=0mm)
savefig(plot2, "images/scenario3_ratiocomparison.png")
