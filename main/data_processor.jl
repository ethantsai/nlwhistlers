using CSV, DataFrames, BasicInterpolators

"""
    prec_to_trap_ratio(rm)

Takes in a Resultant_Matrix and obtains the precipitating to trapped ratio.
This is the ratio between precipitation flux and (flux in LC<α<3*LC)/2.
"""
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

function extract_idl_ratio(
    time_name::String,
    data_name::String,
    error_name::String,
    ebin_name::String,
    start::DateTime, stop::DateTime,
    indices_to_remove=[])

    time_csv_name = "external_data/idl_csvs/"*time_name
    data_csv_name = "external_data/idl_csvs/"*data_name
    ebins_csv_name = "external_data/idl_csvs/"*ebin_name
    error_csv_name = "external_data/idl_csvs/"*error_name

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

function obtain_diffusion_results(lo_or_hi, dawn_dusk_nite, wnX, freq_const_or_vary, FpeFceX)
    #=
    Chorus_Daa_Precip_[Frequency Model]_[fpe/fce ratio]_[L shell]_[MLT]_[Wave normal angle model]
    
    2 Frequency models
    FreConst: constant wave frequency along field line
    FreVary: decreasing wave frequency along field line
             if(lat<20) omega_m = 0.4-0.2*lat/20
             if(lat>=20) omega_m = 0.2
        
    a bunch of fpe/fce ratios
        FpeFceL: fpe/fce = L
        FpeFce2: fpe/fce = 2
        FpeFce3: fpe/fce = 3
        FpeFce3.5: fpe/fce = 3.5
        FpeFce4: fpe/fce = 4
        FpeFce4.5: fpe/fce = 4.5
        FpeFce5: fpe/fce = 5
        FpeFce5.5: fpe/fce = 5.5
        FpeFce6.5: fpe/fce = 6.5

    2 L shells:
        L4.50_LO: L = 4.5
        L6.50_HI: L = 6.5
        
    3 MLTs:
        NITE: nightside
        DAWN: dawnside
        DUSK: duskside
        
    4 wave normal angle models: °
        WN1: quasi-field aligned, theta_m = 0°, theta_w = 30°
        WN2: theta_m = (1-exp(-abs(lat/10°)))*theta_res, 
            theta_w = (theta_res-theta_m)/2
        WN3: theta_m = theta_gendrin+(theta_res-theta_gendrin)*(1-exp(-abs(lat/10°))), 
            theta_w = (theta_res-theta_m)/2
        WN4: theta_m = theta_res*lat/(lat+15.0)
            theta_w = (theta_res-theta_m)/2; if(theta_w > 10.0) theta_w=10.0
    =#
    if lo_or_hi ∈ ["lo", "Lo", "LO"]
        L_shell = "4.50"
    elseif lo_or_hi ∈ ["hi", "Hi", "HI"]
        L_shell = "6.50"
    else
        println("lo_or_hi should either be LO or HI")
    end

    if dawn_dusk_nite ∉ ["DAWN", "DUSK", "NITE"]
        println("dawn_dusk_nite should either be DAWN, DUSK, or NITE")
    end

    if wnX ∉ [1,2,3,4]
        println("wnX should either be 1, 2, 3 or 4")
    end

    if freq_const_or_vary ∈ ["const", "Const", "CONST"]
        freq_model = "Const"
    elseif freq_const_or_vary ∈ ["vary", "Vary", "VARY"]
        freq_model = "Vary"
    end

    if FpeFceX ∉ ["L", 2, 3, 3.5, 4, 4.5, 5, 5.5, 6.5]
        println("FpeFceX should either be \"L\", 2, 3, 3.5, 4, 4.5, 5, 5.5, or 6.5")
    end
    
    prefix = "external_data/diffusion_code_results/Chorus_Daa_Precip_Fre"
    name = "$(prefix)$(freq_model)_FpeFce$(FpeFceX)_L$(L_shell)_$(lo_or_hi)_$(dawn_dusk_nite)_WN$wnX.txt"
    @info "Loading diffusion code results from $name"
    data =  CSV.File(name; header=true, delim=',', types=Float64) |> DataFrame
    E = data.E
    Daa = data."<Daa>"
    prec_ratio = data.Prec_Ratio
    
    return E, Daa, prec_ratio
end


function obtain_diffusion_results(name::String)
    prefix = "external_data/diffusion_code_results/"
    @info "Loading diffusion code results from $name"
    data =  CSV.File(prefix*name; header=true, delim=',', types=Float64) |> DataFrame
    E = data.E
    Daa = data."<Daa>"
    prec_ratio = data.Prec_Ratio
    return E, Daa, prec_ratio
end


function obtain_elfin_ratio(elfin_prec, elfin_trap)
    return elfin_prec[2]./elfin_trap[2]
end

function normalize_to(normalization_energy, energy_1, energy_2, ratio_1, ratio_2)
    # given the energy to normalize at in keV
    # it will interpolate both to that energy then return the ratio of 1 to 2
    one = LinearInterpolator(energy_1, ratio_1)
    two = LinearInterpolator(energy_2, ratio_2)
    return one(normalization_energy) / two(normalization_energy)
end

function make_elfin_error_bars(elfin_data, elfin_error)
    lo = elfin_data[2] .- elfin_error
    hi = elfin_data[2] .+ elfin_error
    return lo, hi
end

calcb!(b::Vector{Float64}, lambda) = @. b = sqrt(1+3*sin(lambda)^2)/(cos(lambda)^6)
calcGamma!(gamma::Vector{Float64}, pz, mu, b::Vector{Float64}) = @. gamma = sqrt(1 + pz^2 + 2*mu*b)
calcAlpha!(alpha::Vector{Float64}, mu, gamma::Vector{Float64}) = @. alpha = rad2deg(asin(sqrt((2*mu)/(gamma^2 - 1))))

function extract(sol::EnsembleSolution)
    allZ = Vector{Vector{Float64}}();
    allPZ = Vector{Vector{Float64}}();
    allQ = Vector{Vector{Float64}}();
    allE = Vector{Vector{Float64}}();
    allZeta = Vector{Vector{Float64}}();
    allPhi = Vector{Vector{Float64}}();
    allPA = Vector{Vector{Float64}}();
    allT = Vector{Vector{Float64}}();
    allLambda = Vector{Vector{Float64}}();
    allBw = Vector{Vector{Float64}}();
    allK = Vector{Vector{Float64}}();
    for traj in sol      
        vars = Array(traj');
        timesteps = length(traj.t);
        b = zeros(timesteps);
        gamma = zeros(timesteps);
        Alpha = zeros(timesteps);

        @views calcb!(b,vars[:,5]);
        @views calcGamma!(gamma,vars[:,2],0.5.*vars[:,4].^2,b);
        @views calcAlpha!(Alpha,0.5.*vars[:,4].^2,gamma);
        @views push!(allT, traj.t);
        @views push!(allZ, vars[:,1]);
        @views push!(allPZ, vars[:,2]);
        @views push!(allQ, vars[:,4]);
        @views push!(allZeta, vars[:,3]);
        @views push!(allPhi, vars[:,6]);
        @views push!(allPA, Alpha);
        @views push!(allE, @. (511*(gamma - 1)));
        @views push!(allLambda, vars[:,5]);
        @views push!(allBw, vars[:,7]);
        @views push!(allK, vars[:,8]);
    end
    @info "$(length(sol)) particles loaded in."
    return allT, allZ, allPZ, allQ, allZeta, allPhi, allE, allPA, allLambda, allBw, allK;
end

function postProcessor(allT::Vector{Vector{Float64}}, allZ::Vector{Vector{Float64}}, allPZ::Vector{Vector{Float64}}, allE::Vector{Vector{Float64}}, allPA::Vector{Vector{Float64}})
    #=
    This function will take the output of the model and convert them into usable m x n matrices
    where m is max number of timesteps for the longest trajectory and N is number of particles
    Arrays that are not m long are filled with NaNs
    =#
    N = length(allT); # num particles from data
    @time tVec = allT[findall(i->i==maximum(length.(allT)),length.(allT))[1]]; # turns all time vectors into a single time vector spanning over the longest trajectory
    timeseriesLength = length(tVec); # all vectors must be this tall to ride
    Zmatrix = fill(NaN,timeseriesLength,N); 
    PZmatrix = fill(NaN,timeseriesLength,N); 
    Ematrix = fill(NaN,timeseriesLength,N); 
    PAmatrix = fill(NaN,timeseriesLength,N); 
    # iterate over each matrix column and fill with each vector
    @time for i = 1:N
        @views Zmatrix[1:length(allT[i]),i] = allZ[i]
        @views PZmatrix[1:length(allT[i]),i] = allPZ[i]
        @views Ematrix[1:length(allT[i]),i] = allE[i]
        @views PAmatrix[1:length(allT[i]),i] = allPA[i]
    end
    @info "Matrices generated."
    return tVec, Zmatrix, PZmatrix, Ematrix, PAmatrix
end

function countLostParticles(allT::Vector{Vector{Float64}}, endTime::Float64)
    #=
    Based on time vectors, counts which ones were lost
    and at what time. Returns a Nx2 array where the first
    column is the time at which the particle was lost, and
    the 2nd column denotes the number lost at that point.
    =#
    lossCounter = []; # initialize vector

    for ntime in allT # loop thru each particle
        if maximum(ntime) != endTime # particle lost if ended early
            push!(lossCounter, maximum(ntime)); # the final entry in the time vector is when the particle got lost
        end
    end
    if isempty(lossCounter) # if particle wasn't lost, then throw in a NaN
        push!(lossCounter, NaN) 
        @warn "All particles trapped!" # since all particles trapped
    end
    lostParticles = [sort(lossCounter) collect(1:length(lossCounter))] # sort it by time, so you can see when particles got lost

    if maximum(lostParticles[:,1]) != endTime # this adds in a final hline from last particle lost to end of simulation
        lostParticles = vcat(lostParticles, [endTime (maximum(lostParticles[:,2]))]);
    end 

    @info "Total of $(lostParticles[end,end]) particles lost during sim"
    return lostParticles
end
