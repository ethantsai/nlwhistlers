function prec_to_trap_ratio(rm)
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
    # plot(PA_min:1:(PA_max-1), final_PA_dist.weights, label=false);
  
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

function extract_idyymmddHHMM(st::String)
    return st[3], st[8:9], st[10:11], st[12:13], st[15:16], st[17:18]
end



using CSV, DataFrames

function extract_idl_ratio(
    time_name::String,
    data_name::String,
    error_name::String,
    ebin_name::String,
    start::DateTime, stop::DateTime,
    indices_to_remove=[])

    time_csv_name = "idl_csvs/"*time_name
    data_csv_name = "idl_csvs/"*data_name
    ebins_csv_name = "idl_csvs/"*ebin_name
    error_csv_name = "idl_csvs/"*error_name

    times_df =  CSV.File(time_csv_name; header=false, delim=',', types=Float64) |> DataFrame
    time = unix2datetime.(times_df.Column1)
    indices = findall((time.>start).&(time.<stop)) # these are the indices corresponding to the time range to sum over
    for i in indices_to_remove
        deleteat!(indices, findall(indices.==i))
    end
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

function extract_elfin_p2t_ratio(datestring, tstart, tend, indices_to_remove=[])
    @info "Loading csv data from $datestring"
    elfin_p2t, elfin_p2t_error = extract_idl_ratio(datestring*"_time.csv", datestring*"_p2t.csv",
                                                    datestring*"_p2t_err.csv", "ebins.csv", # csvs containing ELFIN measurements
                                                    tstart, tend, indices_to_remove); # time to sample from ELFIN measurements                                                                
    #need to check if data-error is below zero
    for i in eachindex(elfin_p2t_error)
        if elfin_p2t[2][i] < elfin_p2t_error[i]
            elfin_p2t_error[i] = elfin_p2t[2][i]-1e-15
        end
    end 
    return elfin_p2t, elfin_p2t_error
end

function obtain_elfin_ratio(elfin_prec, elfin_trap)
    return elfin_prec[2]./elfin_trap[2]
end

function normalize_to_elfin(elfin_ratio, sim_ratio)
    # goal is to normalize to first two channels of ELFIN 
    # mean[63 and 97] keV is 80.61 keV
    # sim has 83 keV at index 6
    return mean(elfin_ratio[1:2]) / sim_ratio[6]
end

function make_elfin_error_bars(elfin_data, elfin_error)
    lo = elfin_data[2] .- elfin_error
    hi = elfin_data[2] .+ elfin_error
    return lo, hi
end

#################
## Plot Colors ##;
#################

using Plots.PlotMeasures;
function hexcolor(r::UInt8, g::UInt8, b::UInt8)
    return RGB(Int(r)/255,Int(g)/255,Int(b)/255)
end
bipride_pink = RGB(234/255, 2/255, 112/255);
bipride_orange = RGB(255/255, 79/255, 0/255);
bipride_lavender = RGB(155/255, 79/255, 150/255);
bipride_blue = RGB(0/255, 56/255, 168/255);
c0 = hexcolor(0xf5,0x97,0x00); # golden
c1 = hexcolor(0xff,0x4f,0x00); # orange
c2 = hexcolor(0xf8,0x00,0x4d); # reddishpink
c3 = hexcolor(0xd3,0x00,0x7d); # magenta
c4 = hexcolor(0x8f,0x19,0x9f); # purple
c5 = hexcolor(0x00,0x38,0xa8); # blue
c6 = hexcolor(0x00,0xA6,0xFF); # light blue
c7 = hexcolor(0x00,0xAB,0x12); # green
c8 = hexcolor(0x00,0xF7,0x1A); # light green


using DirectConvolution

# https://github.com/vincent-picaud/DirectConvolution.jl
# https://pixorblog.wordpress.com/2016/07/13/savitzky-golay-filters-julia/
function smooth(signal, filterwidth::Int, polydegree::Int)
    s = Float64[i for i in signal]
    sg = SG_Filter(halfWidth=filterwidth,degree=polydegree)
    ss = apply_SG_filter(s, sg)
    return ss # 1d savitzky-golay smoothed signal
end