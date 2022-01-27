using InteractiveUtils
using PlutoUI
using TickTock
using ConfParser
using Dates
using Random
using StaticArrays
using Distributed
using OrdinaryDiffEq
using LinearAlgebra
using JLD2
using Plots
using LoopVectorization
using BenchmarkTools
using StatsPlots
using DataFrames
using CSV
using LaTeXStrings
using Plots.PlotMeasures
plot(); # dummy plot just to get the dumb pkg to compile

#######################
## Constants n stuff ##
#######################

const Re   = 6370e3;        # Earth radius, f64
const c    = 3e8;           # speedo lite, f64
const Beq  = 3.e-5;         # B field at equator (T), f64

##################
## Loading data ##
##################

struct Resultant_Matrix
    label::String
    numParticles::Int64
    endTime::Float64
    allZ::Vector{Vector{Float64}}
    allPZ::Vector{Vector{Float64}}
    allT::Vector{Vector{Float64}}
    allPA::Vector{Vector{Float64}}
    allE::Vector{Vector{Float64}}
    lostParticles::Matrix{Float64}
    tVec::Vector{Float64}
    Zmatrix::Matrix{Float64}
    PZmatrix::Matrix{Float64}
    Ematrix::Matrix{Float64}
    PAmatrix::Matrix{Float64}
end

struct Precipitating_Particles
    label::String
    precipitating_fluxes_mean::Vector{Float64}
    precipitating_fluxes_plus::Vector{Float64}
    precipitating_fluxes_minus::Vector{Float64}
end

function parse_conf_file(directoryname::String, conffile::String)
    # Parsing the conf file
    conf = ConfParse(directoryname*"/"*conffile)
    parse_conf!(conf)
    numParticles = parse(Int64, retrieve(conf, "numberOfParticles"));
    startTime = parse(Float64, retrieve(conf, "startTime"));
    endTime = parse(Float64, retrieve(conf, "endTime"));
    tspan = (startTime, endTime); # integration time
    lossConeAngle = parse(Float64, retrieve(conf, "lossConeAngle"));
    saveDecimation = parse(Float64, retrieve(conf, "saveDecimation"));
    L = parse(Float64, retrieve(conf, "L"));
    omegam = parse(Float64, retrieve(conf, "omegam"));
    Omegape = parse(Float64, retrieve(conf, "Omegape"));
    z0 = parse(Float64, retrieve(conf, "z0"));
    lambda0 = parse(Float64, retrieve(conf, "lambda0"));
    waveAmplitudeModifier = parse(Float64, retrieve(conf, "waveAmplitudeModifier"));
    ELo = parse(Float64, retrieve(conf, "ELo"));
    EHi = parse(Float64, retrieve(conf, "EHi"));
    Esteps = parse(Float64, retrieve(conf, "Esteps"));
    PALo = parse(Float64, retrieve(conf, "PALo"));
    PAHi = parse(Float64, retrieve(conf, "PAHi"));
    PAsteps = parse(Float64, retrieve(conf, "PAsteps"));
    ICrange = [ELo, EHi, Esteps, PALo, PAHi, PAsteps];
    batches = parse(Int64, retrieve(conf, "batches"));
    numThreads = parse(Int64, retrieve(conf, "numberOfThreads"))
    @info "Parsed Config file: $conffile"
    return numParticles, endTime
end

function load_resultant_matrix(label::String, directoryname::String, basename::String, conffile::String, num_batches::Int64)
    # Parse conf file
    @time numParticles, endTime = parse_conf_file(directoryname, conffile)    
    # Load in data
    @time allZ, allPZ, allT, allPA, allE = loadData(directoryname, basename, num_batches);
    # Count lost particles
    @time lostParticles = countLostParticles(allT, endTime);
    # Convert into a workable matrix
    @time tVec, Zmatrix, PZmatrix, Ematrix, PAmatrix = postProcessor(allT, allZ, allPZ, allE, allPA);
    return Resultant_Matrix(label, numParticles, endTime, allZ, allPZ, allT, allPA, allE, lostParticles,tVec, Zmatrix, PZmatrix, Ematrix, PAmatrix)
end

function export_results(label::String, precipitating_flux_timeseries)
    ###
    # extract data and save it into 
    ###
    prec_flux_matrix = hcat(precipitating_flux_timeseries...)
    prec_flux_mean = Vector{Float64}()
    prec_flux_plus = Vector{Float64}()
    prec_flux_minus = Vector{Float64}()
    for energy_bin in eachrow(prec_flux_matrix)
        # remove zeros from each row
        row = filter(!iszero,energy_bin)
        if isempty(row)
            push!(prec_flux_mean, 1)
            push!(prec_flux_plus, 1)
            push!(prec_flux_minus, 1)
        elseif length(row) <= 2
            avg_val = mean(row)
            push!(prec_flux_mean, avg_val)
            push!(prec_flux_plus, avg_val*.25)
            push!(prec_flux_minus, avg_val*.25)
        else
            avg_val = mean(row)
            push!(prec_flux_mean, avg_val)
            push!(prec_flux_plus, maximum(row)-avg_val)
            push!(prec_flux_minus, avg_val-minimum(row))
        end
    end
    # test plot
    # plot(Egrid, prec_flux_mean, yerror=(prec_flux_minus, prec_flux_plus), ylim =(1e2,1e9), xlim=(50,800), yscale=:log10)
    
    return Precipitating_Particles(label, prec_flux_mean, prec_flux_plus, prec_flux_minus)
end
# themis_lolat = export_results("210429_themis_lolat", equatorial_fluxes_042921, prec_flux_timeseries_042921)
# @save "210429_data_storage.jld2" themis_lolat themis_hilat


###################
## PSD Functions ##
###################
# define dist function here

# generic
f0 = ((C::Float64, E::Float64, PA::Float64) -> ((C*(E/70)^-3)*(sin(deg2rad(PA)))))
# mms fit for 9/22/20
f0_092220 = function (E::Float64, PA)
    A = 4
    B0 = .25
    B1 = .3
    C = 1e5
    E0 = 300
    E1 = 70
    if E < E1
        B = B0
    elseif E <= E0
        B = B0 + (B1-B0)*(E-E0)/(E1-E0)
    elseif E > E0
        B = B1
    end
    return ((C*(E/70)^-A)*(sin(deg2rad(PA))-sin(deg2rad(8)))^B)
end
# mms fit for 10/27/20
f0_102720 = function (E::Float64, PA)
    A = 3
    B0 = .2
    B1 = .3
    C = 1.5e5
    E0 = 300
    E1 = 70
    if E < 70
        B = B0
    elseif E <= 300
        B = B0 + (B1-B0)*(E-E0)/(E1-E0)
    elseif E > 300
        B = B1
    end
    return ((C*(E/70)^-A)*(sin(deg2rad(PA))-sin(deg2rad(8)))^B)
end
# themis fit for 4/29/21
f0_042921 = function (E::Float64, PA)
    a1 = 5e-13
    a2 = 8e-16
    a3 = 2e-19
    b1 = .2
    b2 = 9
    b3 = 200
    if E < 5
        b = 0.5
    elseif E <= 20
        b = 0.5
    elseif E > 40
        b = 0
    end
    Epsd = (a1 * (1 + E/b1)^-5) + (a2 * (1 + E/b2)^-5) + (a3 * (1 + E/b3)^-3)
    PApsd = (sin(deg2rad(PA))-sin(deg2rad(8)))^b
    # km --> cm and 1/v^4 --> 1/E^2 conversion
    # since v = 438 * sqrt(1836 KeV)
    return (1e5*438^4*1836^2)*(Epsd*PApsd)
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


###############################
## Post Processing functions ##
###############################

calcb!(b::Vector{Float64}, lambda::SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}) = @. b = sqrt(1+3*sin(lambda)^2)/(cos(lambda)^6)
calcGamma!(gamma::Vector{Float64}, pz::SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, mu::SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, b::Vector{Float64}) = @. gamma = sqrt(1 + pz^2 + 2*mu*b)
calcAlpha!(alpha::Vector{Float64}, mu::SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, gamma::Vector{Float64}) = @. alpha = rad2deg(asin(sqrt((2*mu)/(gamma^2 - 1))))

function loadData(directory::String, basename::String, num_batches::Int64)
    #=
    Loads in the JLD2 file and creates timeseries of all the data.
    Takes in filename in string, returns 5 TxN vectors for time,
    z, pz, pitch angle, and energy timeseries for all N particles.
    =#

    allZ = Vector{Vector{Float64}}();
    allPZ = Vector{Vector{Float64}}();
    allE = Vector{Vector{Float64}}();
    allPA = Vector{Vector{Float64}}();
    allT = Vector{Vector{Float64}}();

    for i in 1:num_batches
        JLD2.@load directory*"/"*basename*"_$i.jld2" sol
        for traj in sol
            vars = Array(traj')
            timesteps = length(traj.t)
            b = zeros(timesteps)
            gamma = zeros(timesteps)
            Alpha = zeros(timesteps)
    
            @views calcb!(b,vars[:,5])
            @views calcGamma!(gamma,vars[:,2],vars[:,4],b)
            @views calcAlpha!(Alpha,vars[:,4],gamma)
            @views push!(allT, traj.t);
            @views push!(allZ, vars[:,1]);
            @views push!(allPZ, vars[:,2]);
            @views push!(allPA, Alpha);
            @views push!(allE, @. (511*(gamma - 1)));
        end
        @info "$(length(sol)) particles loaded in from $(basename*"_$i.jld2")"
        @info "Total $(length(allT)) particles loaded so far..."
    end
    @info "Loaded total $(length(allT)) particles!"
    return allZ, allPZ, allT, allPA, allE
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

function postProcessor(allT::Vector{Vector{Float64}}, allZ::Vector{Vector{Float64}}, allPZ::Vector{Vector{Float64}}, allE::Vector{Vector{Float64}}, allPA::Vector{Vector{Float64}})
    #=
    This function will take the output of the model and convert them into usable m x n matrices
    where m is max number of timesteps for the longest trajectory and N is number of particles
    Arrays that are not m long are filled with NaNs
    =#
    N = length(allT); # num particles from data
    tVec = allT[findall(i->i==maximum(length.(allT)),length.(allT))[1]]; # turns all time vectors into a single time vector spanning over the longest trajectory
    timeseriesLength = length(tVec); # all vectors must be this tall to ride
    Zmatrix = fill(NaN,timeseriesLength,N); 
    PZmatrix = fill(NaN,timeseriesLength,N); 
    Ematrix = fill(NaN,timeseriesLength,N); 
    PAmatrix = fill(NaN,timeseriesLength,N); 
    # iterate over each matrix column and fill with each vector
    for i = 1:N
        @views Zmatrix[1:length(allT[i]),i] = allZ[i]
        @views PZmatrix[1:length(allT[i]),i] = allPZ[i]
        @views Ematrix[1:length(allT[i]),i] = allE[i]
        @views PAmatrix[1:length(allT[i]),i] = allPA[i]
    end

    return tVec, Zmatrix, PZmatrix, Ematrix, PAmatrix
end

function animatePAD(gifFileName="PADanimation.gif", maxParticles=1000, binwidth=5)
    pyplot()
    animDec = 10; # make a png for animation every 10 points
    animScale = 10; # i.e. animscale = 10 means every 10 seconds in animation is 1 second of simulation time (increase for longer animation)
    anim = @animate for i=1:length(tVec)
        histogram(PAmatrix[i,:],bins=(0:binwidth:90), ylim = [0, 1.1*maxParticles], #xminorticks = 4, yminorticks = 2,
            legend = false, dpi = 100,
            xlabel="Pitch Angle (deg)", ylabel="Number of Particles", title="Time evolution of Pitch Angle Distribution");
        vline!([lossConeAngle]);
        annotate!(lossConeAngle+.5, 1.05*maxParticles, text("Lost Particle Count:", :left));
        annotate!(lossConeAngle+.5, .95*maxParticles, text((count(i->isnan(i), Ematrix[i,:])), :left));
        annotate!(80, 1.05*maxParticles, text("t = $(round(tVec[i]*R*L/(c),digits=4)) s"), :right)
    end every animDec
    savename = string("results/",gifFileName)
    gif(anim, savename, fps = (length(tVec)/animDec)/(animScale*endTime*Re*L/(c)))
end

function animateESD(gifFileName="ESDanimation.gif", maxParticles=1000, binwidth=50)
    pyplot()
    animDec = 10; # make a png for animation every 100 points
    animScale = 5; # i.e. animscale = 10 means every 10 seconds in animation is 1 second of simulation time (increase for longer animation)
    anim = @animate for i=1:length(tVec)
        histogram(Ematrix[i,:],bins=(0:binwidth:2000), ylim = [1, 1.1*maxParticles], #xminorticks = 4, yminorticks = 2,
            legend = false, dpi = 100, xscale=:log10, yscale =:log10, xticks = ([10,100,1000],["10", "100", "1000"]), xlim = [10,2000],
            xlabel="Energy (KeV)", ylabel="Number of Particles", title="Time evolution of Energy Spectra");
        annotate!(20, 1.05*maxParticles, text("t = $(round(tVec[i]*Re*L/(c),digits=4)) s"), :left)
    end every animDec
    savename = string("results/",gifFileName)
    gif(anim, savename, fps = (length(tVec)/animDec)/(animScale*endTime*Re*L/(c)))
end

function animatePSD(gifFileName="PSDanimation.gif", PAbinwidth=3, Ebinwidth=100, maxParticles=50)
    pyplot()
    animDec = 10; # make a png for animation every 10 points
    animScale = 1; # i.e. animscale = 10 means every 10 seconds in animation is 1 second of simulation time (increase for longer animation)
    anim = @animate for i=1:length(tVec)
        histogram2d(PAmatrix[i,:], Ematrix[i,:], bins=(0:PAbinwidth:90,0:Ebinwidth:2000),
            xlabel="Pitch Angle (deg)",ylabel="Energy (keV)",title="Time evolution of Phase Space Density",
            yscale=:log10, clims = (1,maxParticles), colorbar_scale=:log10, minorgrid=true,
            yticks = ([1, 10,100,1000,10000],["1", "10", "100", "1000", "10000"]));
        annotate!(lossConeAngle, 5, text("Lost Particle Count:", :left));
        annotate!(lossConeAngle, 3, text((count(i->isnan(i), Ematrix[i,:])), :left))
        annotate!(70, 5, text("t = $(round(tVec[i]*Re*L/(c),digits=4)) s"), :right)
        # count(i->isnan(i), Ematrix[1,:]) # counts number of NaNs (i.e particles lost) at certain time
    end every animDec
    savename = string("results/",gifFileName)
    gif(anim, savename, fps = (length(tVec)/animDec)/(animScale*endTime*Re*L/(c)))
end
            
function recalcDistFunc(Ematrix::Array{Float64,2},PAmatrix::Array{Float64,2},initial::Int64,final::Int64,distFunc,
    Egrid::Vector{Float64}, PAgrid::StepRange{Int64,Int64}, whistler_occurence_rate::Float64)
    #=
    Takes in matrix of Energy and PA, initial and final indices, and grid values as stepranges.
    Recalculates out a new PSD given a a distribution function at the initial and final indices.
    =#
    N = @views length(Ematrix[1,:]) # num particles
    EPA0i = @views round.([Ematrix[1,:] PAmatrix[1,:]])
    Nkl = @views N/(length(unique(EPA0i[:,1]))*length(unique(EPA0i[:,2]))) # obtain number of particles/bin

    EPAinitial = @views @. round([Ematrix[initial,:] PAmatrix[initial,:]])
    EPAfinal = @views @. round([Ematrix[final,:] PAmatrix[final,:]])
    EPAfinal_next = @views @. round([Ematrix[final+1,:] PAmatrix[final+1,:]])

    # initialize matrices to be used based on grid size
    f = zeros(length(Egrid), length(PAgrid));
    psd_init = zeros(length(Egrid), length(PAgrid));
    psd_final = zeros(length(Egrid), length(PAgrid));
    # initialize vectors to be filled in
    f0Vec = Vector{Float64}();
    indices = Vector{Tuple{Int64,Int64}}();
    psd_prec = zeros(length(Egrid));
    excludedParticles = 0;
    
    ### HANDLE NANS HERE ###
    # diff = @. abs(EPAinitial - EPAfinal)
    valid_rows = @. ~isnan(EPAfinal)[:,1] & ~isnan(EPAfinal[:,2]); # non-lost particles from final state
    valid_rows_next = @. ~isnan(EPAfinal_next)[:,1] & ~isnan(EPAfinal_next[:,2]); # non-lost particles from state at next timestep
    # if valid_rows!=valid_rows_next # this means that at least one of the particles will precipitate
        prec_rows = valid_rows .!= valid_rows_next # these are the indices of particles who are about to precipitate
        prec_indices = (ones(N).*vec(prec_rows))[vec(valid_rows)]
        # @info "$(length(ones(N)[prec_rows])) particles are about to precipitate"
    # end
    EPAinitial = EPAinitial[vec(valid_rows), :]; # new vec with no NaNs
    EPAfinal = EPAfinal[vec(valid_rows), :]; # new vec with no NaNs
    
    lostParticles = length(valid_rows) - length(EPAfinal[:,1]);
    # @info "$lostParticles particles lost at index $final."

    @inbounds for EPAi in eachrow(EPAinitial)
        k,l = 1, 1;
        while Egrid[k] < EPAi[1]; k+=1; end # energy
        while PAgrid[l] < EPAi[2]; l+=1; end # PA
        push!(indices, (k-1,l-1))
        push!(f0Vec, distFunc(EPAi[1], EPAi[2]))
        f[k-1,l-1] += 1;
    end

    psdVec = [(f0Vec[i]/f[indices[i][1],indices[i][2]]) for i in eachindex(f0Vec[:,1])]

    @inbounds for i in eachindex(EPAinitial[:,1]) # i iterates over each particle
        if ~(EPAfinal[i,1]>maximum(Egrid) || EPAfinal[i,2]>maximum(PAgrid)) # skip loop if data is outside of range
            k,l = 1, 1;
            while Egrid[k] < EPAinitial[i,1]; k+=1; end
            while PAgrid[l] < EPAinitial[i,2]; l+=1; end
            if k!=1 && l!=1; psd_init[k-1,l-1] += psdVec[i]; end
            k,l = 1, 1;
            while Egrid[k] < EPAfinal[i,1]; k+=1; end
            while PAgrid[l] < EPAfinal[i,2]; l+=1; end
            if k!=1 && l!=1; psd_final[k-1,l-1] += psdVec[i]; end
            if k!=1 && prec_indices[i] == 1.0;  # handle precipitaitng particles here
                psd_prec[k-1] += psdVec[i]*whistler_occurence_rate;
            end
        else
            excludedParticles += 1;
        end
    end
    # @info "Excluded $excludedParticles particles due to out of range."

    return f, psd_init, psd_final, psd_prec
end

function make_psd_timeseries(Ematrix,PAmatrix,tVec, dist_func, Egrid, PAgrid, whistler_occurence_rate)
    # recalc psd at every 
    psd_timeseries = Vector{Matrix{Float64}}()
    f_timeseries = copy(psd_timeseries)
    psd_prec_timeseries = Vector{Vector{Float64}}()
    @inbounds for time_index in eachindex(tVec[1:end-1])
        f, _, psd_final, psd_prec = recalcDistFunc(Ematrix, PAmatrix, 1, time_index, dist_func, Egrid, PAgrid, whistler_occurence_rate);
        push!(f_timeseries, f)
        push!(psd_timeseries, psd_final)
        push!(psd_prec_timeseries, psd_prec)
    end
    return f_timeseries, psd_timeseries, psd_prec_timeseries
end

function bin_psd_prec_timeseries(psd_prec_timeseries, indexArray)
    # given indexArray (an array with the start indices of each time bin)
    # iterate through and sum PSD to produce PSD of given time bins
    binned_psd_prec_timeseries = Vector{Vector{Float64}}() 
    index = 0
    # loop through each time bin
    @inbounds for time_bin_index in indexArray
        binned_psd = zeros(length(psd_prec_timeseries[1]))
        index += 1
        while index < time_bin_index # at every index, add in the 
            @. binned_psd += psd_prec_timeseries[index]
            index += 1
        end
        push!(binned_psd_prec_timeseries, binned_psd) 
    end
    return binned_psd_prec_timeseries
end

function calc_equatorial_fluxes(result_matrix::Resultant_Matrix, dist_func)
    # initial electron PSD is PSD_0i = flux(E0)/E0 of the recalced distribution
    _, psd_init, _, _ = recalcDistFunc(result_matrix.Ematrix,result_matrix.PAmatrix, 1, 1,dist_func, Egrid, PAgrid, 1.);
    @info "Equatorial fluxes calculated from $(result_matrix.label)."
    return [sum(Erow) for Erow in eachrow(psd_init)]  .* Egrid .* 1000 ./ (length(PAgrid))
end

function calc_prec_flux(result_matrix::Resultant_Matrix,timeBin::Int64,dist_func::Function,whistler_occurence_rate::Float64,count_threshold::Int64)
    # calc only precipitaitng fluxes
    allPrecip, indexArray, allPrecipInitial = precipitatingParticles(result_matrix, timeBin);
    f_timeseries, psd_timeseries, psd_prec_timeseries = make_psd_timeseries(result_matrix.Ematrix,result_matrix.PAmatrix,result_matrix.tVec, dist_func, Egrid, PAgrid, whistler_occurence_rate);
    binned_psd_prec_timeseries = bin_psd_prec_timeseries(psd_prec_timeseries, indexArray);
    @info "Precipitating fluxes calculated from $(result_matrix.label)."

    # remove all precipitating particles with counts per energy channel < count_threshold
    num_particles_per_channel = [length(findall(x -> Egrid[i]<=x<=Egrid[i+1],vcat(allPrecip...))) for i in 1:(length(Egrid)-1)];
    remove_low_counts = vcat([0],(@. (sign(num_particles_per_channel - count_threshold) + 1) / 2));
    binned_psd_prec_timeseries_clean = [list .* remove_low_counts for list in binned_psd_prec_timeseries]

    return calc_precipitating_flux_timeseries(binned_psd_prec_timeseries_clean);
end

function calc_precipitating_flux_timeseries(binned_psd_prec_timeseries)
    # convert prec flux into units of 1/cm^2/s/MeV where [s] is actually the width of timebin
    return [prec_psd.*Egrid.*1000 for prec_psd in binned_psd_prec_timeseries]
end

function animate_flux_comparison(gifFileName::String, equatorial_flux, elfin_measurements, thing::Vector{Vector{Float64}}, minY, maxY)
    pyplot()
    animDec = 1; # make a png for animation every 10 points
    animScale = 50; # i.e. animscale = 10 means every 10 seconds in animation is 1 second of simulation time (increase for longer animation)
    maxEnergy=1000
    anim = @animate for i in eachindex(thing[1:end-1])
        plot(Egrid, equatorial_flux)
        plot!(elfin_measurements)
        plot!(Egrid, thing[1:end-1], color = :gray, alpha = .5);
        plot!(Egrid,thing[i], color = :orange);
        plot!(ylim =(minY,maxY), xlim=(10,maxEnergy), yscale=:log10, legend=false);
        plot!(xlabel="Energy (keV)", ylabel="Flux (1/cm^2/s/sr/MeV)", title="Energy Fluxes of Precipitating Particles");
        annotate!(300, 0.1*maxY, text("t = $(round(tVec[indexArray[i]]*Re*L/(c),digits=3)) s"), :left)
    end every animDec
    savename = string("results/",gifFileName)
    gif(anim, savename, fps = (length(tVec)/animDec)/(animScale*endTime*Re*L/(c)))
    return anim
end

function animate_a_thing(gifFileName::String, thing::Vector{Matrix{Float64}})
    pyplot()
    animDec = 10; # make a png for animation every 10 points
    animScale = 1; # i.e. animscale = 10 means every 10 seconds in animation is 1 second of simulation time (increase for longer animation)
    initial, final = 1, 1
    anim = @animate for i in eachindex(thing)
        heatmap(PAgrid, Egrid, log10.(thing[i]), fc = :plasma, colorbar = false, 
            xlabel="Pitch Angle (deg)",ylabel="Energy (keV)", yscale=:log10, title = "Recalculated PSD at t = $(round(tVec[i]*Re*L/(c),digits=2)) s")
    end every animDec
    savename = string("results/",gifFileName)
    gif(anim, savename, fps = (length(tVec)/animDec)/(animScale*endTime*Re*L/(c)))
    return anim
end

function animate_a_thing(gifFileName::String, thing::Vector{Vector{Float64}}, minY, maxY)
    pyplot()
    animDec = 1; # make a png for animation every 10 points
    animScale = 50; # i.e. animscale = 10 means every 10 seconds in animation is 1 second of simulation time (increase for longer animation)
    maxEnergy=1000
    anim = @animate for i in eachindex(thing[1:end-1])
        plot(Egrid, thing[1:end-1], color = :gray, alpha = .5);
        plot!(Egrid,thing[i], color = :orange);
        plot!(ylim =(minY,maxY), xlim=(10,maxEnergy), yscale=:log10, xscale=:log10, legend=false);
        plot!(xlabel="Energy (keV)", ylabel="PSD", title="Energy Spectra of Precipitating Particles");
        annotate!(300, 0.1*maxY, text("t = $(round(tVec[indexArray[i]]*Re*L/(c),digits=3)) s"), :left)
    end every animDec
    savename = string("results/",gifFileName)
    gif(anim, savename, fps = (length(tVec)/animDec)/(animScale*endTime*Re*L/(c)))
    return anim
end

@userplot MarginalHist
@recipe function f(h::MarginalHist)
    if length(h.args) != 2 || !(typeof(h.args[1]) <: AbstractVector) ||
        !(typeof(h.args[2]) <: AbstractVector)
        error("Marginal Histograms should be given two vectors.  Got: $(typeof(h.args))")
    end
    x, y = h.args

    # set up the subplots
    legend := false
    link := :both
    framestyle := [:none :axes :none]
    grid := false
    layout := @layout [tophist           _
                       hist2d{0.9w,0.9h} righthist]

    # main histogram2d
    @series begin
        
        seriestype := :histogram2d
        subplot := 2
        x, y
    end

    # these are common to both marginal histograms
    fillcolor := :black
    fillalpha := 0.3
    linealpha := 0.3
    seriestype := :histogram

    # upper histogram
    @series begin
        subplot := 1
        x
    end

    # right histogram
    @series begin
        orientation := :h
        subplot := 3
        y
    end
end

function checkDistFunc(f::Array{Float64,2}, psd_init::Array{Float64,2}, psd_final::Array{Float64,2}, initial::Int64, final::Int64, Egrid::Vector{Float64}, PAgrid::StepRange{Int64,Int64})
    origPlot =  heatmap(PAgrid, Egrid, f, fc = :plasma,
    ylim = (35,1000), yscale=:log10, xlim = (8,90),
    xlabel="Pitch Angle (deg)",ylabel="Energy (keV)", title = "Original Flat PSD at t = 0");

    initPlot =  heatmap(PAgrid, Egrid, psd_init, fc = :plasma, colorbar_scale=:log10, clim = (1e3,1e4),
    ylim = (35,1000), yscale=:log10, xlim = (8,90),
    xlabel="Pitch Angle (deg)",ylabel="Energy (keV)", title = "Recalculated PSD at t = $(round(tVec[initial]*Re*L/(c),digits=3)) s");

    # finalPlot = heatmap(PAgrid, Egrid, log10.(psd_final), fc = :plasma, colorbar = false,
    #     xlabel="Pitch Angle (deg)",ylabel="Energy (keV)", title = "Recalculated PSD at t = $(round(tVec[final]*Re*L/(c),digits=2)) s");

    # origPlotAnime = marginalhist(round.(PAmatrix[final,:]), log10.(round.(Ematrix[final,:])),
    # xlim = (0,90), yscale = :identity ,yticks = ([1, 2,3,3.70],["10", "100", "1000", "5000"]),
    # fc =:plasma, bins= 30);

    initial_dist = plot(origPlot, initPlot, layout = (1,2), size = (1500,600)) # want this for just first two plots
    # savefig(initial_dist, string("results/plots/initial_recalced_dist3.png"))

    # plot(origPlot, initPlot, origPlotAnime, finalPlot, layout = (2,2), size = (1200,1200)) # want this for all plots
    # plot(finalPlot) # want this for just one plot
    
end


function animatePSD(gifFileName, psd_timeseries)
    pyplot()
    animDec = 5; # make a png for animation every 10 points
    animScale = 10; # i.e. animscale = 10 means every 10 seconds in animation is 1 second of simulation time (increase for longer animation)
    maxval=minimum(maximum.(psd_timeseries))
    anim = @inbounds @animate for i in eachindex(psd_timeseries)
        heatmap(PAgrid, Egrid,(psd_timeseries[i]), fc = :plasma, colorbar = true,
        colorbar_scale=:log10, clims=(1e3, maxval), ylim = (35,1000), yscale=:log10,
        xlim = (8,90),
        xlabel="Pitch Angle (deg)",ylabel="Energy (keV)", title = "Recalculated PSD at t = $(round(tVec[i]*Re*L/(c),digits=2)) s")
    end every animDec
    savename = string("results/",gifFileName)
    gif(anim, savename, fps = (length(tVec)/animDec)/(animScale*endTime*Re*L/(c)))
end

function get_Nloss_Ntotal_per_Energy(tVec, Ematrix)
    # get Nloss and Ntotal for each energy at each timestep 
    # should return a length(tVec) long vector of vectors, each w/ dimension length(Egrid)
    for index in eachindex(tVec)
        Ematrix()
        #wip
    end
end

function precipitatingParticles(rm::Resultant_Matrix, timeBin=10)
    unitSimTime = mean(diff(rm.tVec[begin:end-1]));
    simTimeBin = convert(Int64, ceil(timeBin/unitSimTime, digits=0))

    # produce a vector with indices demarcating time bin boundaries
    indexArray = [minimum([i*simTimeBin length(rm.tVec)]) for i in 1:convert(Int64,floor(rm.endTime/timeBin))]
    if length(rm.tVec)!=indexArray[end] push!(indexArray, length(rm.tVec)-1) end #guarantee last index corresponds with proper end
    # remove all particles that didn't get lost
    Ematrixprime = @views Matrix(rm.Ematrix[1:(end-1),findall(isnan, rm.Ematrix[indexArray[end]-1,:])]) # this is a new matrix that only has particles that have precipitated    # replace last value with second to last value to prevent off by one problems
    if indexArray[end] == length(rm.tVec)
        indexArray[end] = indexArray[end]-1
    end

    allPrecipInitial = Vector{Vector{Float64}}()
    allPrecip = Vector{Vector{Float64}}()
    for i in indexArray

        @debug "$(length(Ematrixprime[1,:])) remaining particles" size(Ematrixprime)
        if simTimeBin > (length(Ematrixprime[1,:]))
            @info "sumting wong"
            break
        end
            
        precipitatingIndices = @views findall(isnan, Ematrixprime[simTimeBin,:])
        notPrecipitatingIndices = @views findall(!isnan, Ematrixprime[simTimeBin,:])
        @debug "$(length(precipitatingIndices)) particles precipitated by index $i" precipitatingIndices
        # push!(allPrecip,[@views mean(Ematrixprime[:,j][findall(!isnan, Ematrixprime[:,j])]) for j in precipitatingIndices]) # provides mean of time bin
        push!(allPrecip,[Ematrixprime[:,j][last_index_before_nan(Ematrixprime[:,j])] for j in precipitatingIndices]) # provides last energy before precipitating
        push!(allPrecipInitial,[Ematrixprime[:,j][1] for j in precipitatingIndices]) # provides the initial energy of particle about to precipitate


        @debug "added stuff" size(Ematrixprime)
        if !isempty(notPrecipitatingIndices)
            Ematrixprime = @views Matrix(Ematrixprime[simTimeBin:end,notPrecipitatingIndices]) # purposely including one extra index in case particle is lost right after
        else
            break
        end
    end

    return allPrecip, indexArray, allPrecipInitial

end


function precipitating_initial_state_analyzer(rm::Resultant_Matrix,Echannel_start::Int64,Echannel_end::Int64)

    # remove all particles that didn't get lost
    Ematrixprime = @views Matrix(rm.Ematrix[:,findall(isnan, rm.Ematrix[length(rm.Ematrix[:,1])-1,:])]) # this is a new matrix that only has particles that have precipitated
    PAmatrixprime = @views Matrix(rm.PAmatrix[:,findall(isnan, rm.PAmatrix[length(rm.PAmatrix[:,1])-1,:])]) # this is a new matrix that only has particles that have precipitated

    E_allPrecipFinal = [c[last_index_before_nan(Vector(c))] for c in eachcol(Ematrixprime)] # a list of all final energies
    E_allPrecipInit = [c[1] for c in eachcol(Ematrixprime)] # a list of all initial energies of precipitating particles
    # PA_allPrecipFinal = [c[last_index_before_nan(Vector(c))] for c in eachcol(PAmatrixprime)]
    PA_allPrecipInit = [c[1] for c in eachcol(PAmatrixprime)]

    energies_to_plot = Vector{Vector{Float64}}()
    pitchangles_to_plot = Vector{Vector{Float64}}()
    plot_labels = Vector{String}()
    for i in Echannel_start:Echannel_end-1
        channel = findall(x -> Egrid[i]<=x<=Egrid[i+1],E_allPrecipFinal) # find energies that end in specified e channel
        push!(plot_labels, string(round(Egrid[i]))*" to "*string(round(Egrid[i+1]))*" keV")
        push!(energies_to_plot,E_allPrecipInit[channel])
        push!(pitchangles_to_plot,PA_allPrecipInit[channel])
    end

    return scatter(pitchangles_to_plot, energies_to_plot,
        xlim = [5,20], ylim = [100,500],
        palette = :seaborn_colorblind6,
        title=rm.label,
        xlabel = "Initial Pitch Angle",
        ylabel = "Initial Energy",
        label=plot_labels)
end

function animatePrecipitatingParticles(gifFileName, rm::Resultant_Matrix, timeBin=10, maxFraction=0.004, maxEnergy=1000)
    allPrecip, indexArray = precipitatingParticles(rm, timeBin)
    animDec = 1; # make a png for animation every 10 points
    animScale = 50; # i.e. animscale = 10 means every 10 seconds in animation is 1 second of simulation time (increase for longer animation)
    anim = @animate for i in eachindex(allPrecip)
        density(allPrecip, color = :gray, alpha = .5);
        density!(allPrecip[i], color = :orange);
        plot!(ylim =(0.,maxFraction), yscale=:identity, xlim=(20,maxEnergy), xscale=:log10, legend=false);
        plot!(xlabel="Energy (keV)", ylabel="Fraction of $(rm.numParticles) total particles", title="Energy Spectra of Precipitating Particles");
        annotate!(40, 0.1*maxFraction, text("t = $(round(rm.tVec[indexArray[i]]*Re*L/(c),digits=3)) s"), :left)
    end every animDec
    savename = string("results/",gifFileName)
    gif(anim, savename, fps = (length(rm.tVec)/animDec)/(animScale*rm.endTime*Re*L/(c)))
end

function trajectoryChecking(rm::Resultant_Matrix, index::Int64)
    converT = Re*L/(c)
    @info "starting PA: $(rm.PAmatrix[1,index]), startin E: $(rm.Ematrix[1,index])"
    plot(converT*rm.tVec, rm.PAmatrix[:,index]);
    plot!(converT*rm.tVec, rm.Ematrix[:,index])
    # plot(Zmatrix[:,index], PZmatrix[:,index])
end

function trajectoryChecking(rm::Resultant_Matrix, trajectories::Array{Int64}, plot_title::String)
    converT = Re*L/(c)
    PAplot = plot(xlim=(0,33.124), ylim = (0,90), title=plot_title, ylabel="Pitch Angle ("*L"\degree"*")");
    PAplot = plot!(converT*rm.tVec, rm.PAmatrix[:,trajectories], yminorticks=6, xminorticks=5, label = "", linewidth=2, palette = :seaborn_colorblind6);

    Eplot = plot(xlim=(0,33.124), ylim = (40,1000), yscale = :log10, xlabel = "Time (s)", ylabel="Energy (keV)");
    Eplot = plot!(converT*rm.tVec, rm.Ematrix[:,trajectories], yminorticks=10, xminorticks=5,  label = "", bottom_margin = 8mm, linewidth=2, palette = :seaborn_colorblind6);

    bigplot = plot(PAplot,Eplot, dpi = 96, layout = (2,1), left_margin = 8mm, size=(1000,500))
end


function trajectoryTracing(rm::Resultant_Matrix, trajectories, animDec)
    converT = Re*L/(c)

    # PAanim = @animate for i in eachindex(rm.tVec)
    #     plot(xlim=(0,12), ylim = (0,90));
    #     plot!(converT*rm.tVec, rm.PAmatrix[:,trajectories], color=:gray, alpha = .5, label = "");
    #     for j in trajectories
    #         plot!(converT*rm.tVec[1:i], rm.PAmatrix[:,j][1:i], label = "")
    #         scatter!((converT*rm.tVec[i], rm.PAmatrix[:,j][i]), label = "")
    #     end
    # end every animDec

    # savename1 = string("images/",rm.label,"_PAgif.gif")
    # gif(PAanim, savename1, fps = 20)
    
    
    # Eanim = @animate for i in eachindex(rm.tVec)
    #     Eplot = plot(xlim=(0,12), ylim = (10,5000), yscale = :log10)
    #     plot!(converT*rm.tVec, rm.Ematrix[:,trajectories], color=:gray, alpha = .5, label = "");
    #     for j in trajectories
    #         Eplot = plot!(converT*rm.tVec[1:i], rm.Ematrix[:,j][1:i], label = "")
    #         Eplot = scatter!((converT*rm.tVec[i], rm.Ematrix[:,j][i]), label = "")
    #     end
    # end every animDec

    # savename2 = string("images/",rm.label,"_Egif.gif")
    # gif(Eanim, savename2, fps = 20)
            
    # ZPZanim = @animate for i in eachindex(rm.tVec)
    #     ZPZplot = plot(xlim=(-1,1), ylim = (-7,7))
    #     for j in trajectories
    #         ZPZplot = plot!(rm.Zmatrix[:,j][1:i], rm.PZmatrix[:,j][1:i], alpha = max.((1:i) .+ 100 .- i, 0) / 100, label = "")
    #         ZPZplot = scatter!((rm.Zmatrix[:,j][i], rm.PZmatrix[:,j][i]), label = "")
    #     end
    # end every animDec

    # savename3 = string("images/",rm.label,"_ZPZgif.gif")
    # gif(ZPZanim, savename3, fps = 20)

    # combine all 3
    
    AllAnim = @animate for i in eachindex(rm.tVec)
        PAplot = plot(xlim=(0,33.124), ylim = (0,90), title="Pitch Angle Trajectory", xlabel = "Time(s)", ylabel="Pitch Angle");
        PAplot = plot!(converT*rm.tVec, rm.PAmatrix[:,trajectories], color=:gray, alpha = .5, label = "");
        for j in trajectories
            PAplot = plot!(converT*rm.tVec[1:i], rm.PAmatrix[:,j][1:i], label = "");
            PAplot = scatter!((converT*rm.tVec[i], rm.PAmatrix[:,j][i]), label = "");
        end
        Eplot = plot(xlim=(0,33.124), ylim = (10,1000), yscale = :log10, title="Energy Trajectory", xlabel = "Time(s)", ylabel="Energy (keV)");
        plot!(converT*rm.tVec, rm.Ematrix[:,trajectories], color=:gray, alpha = .5, label = "");
        for j in trajectories
            Eplot = plot!(converT*rm.tVec[1:i], rm.Ematrix[:,j][1:i], label = "");
            Eplot = scatter!((converT*rm.tVec[i], rm.Ematrix[:,j][i]), label = "");
        end
        ZPZplot = plot(xlim=(-1,1), ylim = (-7,7), title="Z-Pz phase space", xlabel="Z", ylabel="Pz");
        for j in trajectories
            ZPZplot = plot!(rm.Zmatrix[:,j][1:i], rm.PZmatrix[:,j][1:i], alpha = max.((1:i) .+ 100 .- i, 0) / 100, label = "");
            ZPZplot = scatter!((rm.Zmatrix[:,j][i], rm.PZmatrix[:,j][i]), label = "");
        end
        bigplot = plot(ZPZplot,PAplot,Eplot, dpi = 100, layout = (3,1), left_margin = 8mm, size=(500,1000))
    end every animDec

    savename4 = string("images/",rm.label,"_ALLgif.gif")
    gif(AllAnim, savename4, fps = 20)


end

function last_index_before_nan(x::Vector{Float64})
    k = 1
    @inbounds for i in eachindex(x)
        if isnan(x[i])
            k = i-1
            return k
        end
    end
    return
 end




function find_lost_particles(row1::Vector{Float64}, row2::Vector{Float64})
    # compares two rows and returns the index of all particles from row 1 that will be lost in row 2
    
    # valid_rows = @. ~isnan(EPAfinal)[:,1] & ~isnan(EPAfinal[:,2]); # non-lost particles from final state
    # valid_rows_next = @. ~isnan(EPAfinal_next)[:,1] & ~isnan(EPAfinal_next[:,2]); # non-lost particles from state at next timestep
    # # if valid_rows!=valid_rows_next # this means that at least one of the particles will precipitate
    # prec_rows = valid_rows .!= valid_rows_next # these are the indices of particles who are about to precipitate
    # prec_indices = (ones(N).*vec(prec_rows))[vec(valid_rows)]
    # @info "$(length(ones(N)[prec_rows])) particles are about to precipitate"

end

function extract_idl_csv(
    time_name::String,
    data_name::String,
    ebin_name::String,
    start::DateTime, stop::DateTime)

    time_csv_name = "idl_csvs/"*time_name
    data_csv_name = "idl_csvs/"*data_name
    ebins_csv_name = "idl_csvs/"*ebin_name

    times_df =  CSV.File(time_csv_name; header=false, delim=',', types=Float64) |> DataFrame
    time = unix2datetime.(times_df.Column1)
    indices = findall((time.>start).&(time.<stop)) # these are the indices corresponding to the time range to sum over
    time_of_interest = time[indices]
    @info "Summing over $(time_of_interest[end]-time_of_interest[1])"

    # import particle flux data 
    data_df  =  CSV.File(data_csv_name; header=false, delim=',', types=Float64) |> DataFrame
    data = [[data_df[row,col] for col in 1:length(data_df[1,:])] for row in 1:16]
    # get rid of NaNs and Infs
    data_of_interest = data
    for i = 1:16
        data_of_interest[i][findall(.!isfinite.(data[i]))] .= 0.0
    end
    flux = [sum(data_of_interest[energy][indices]) for energy in 1:16]

    ebins_df = CSV.File(ebins_csv_name; header=false, delim=',', types=Float64) |> DataFrame
    ebins = ebins_df.Column1

    return (ebins, flux)
end

# time_name = "092220_time.csv"
# data_name = "092220_prec.csv"
# error_name = "092220_precerror.csv"
# error_time_name = "092220_precerror_time.csv"
# ebin_name = "ebins.csv"
# start = DateTime(2020,9,22,9,16,38)
# stop = DateTime(2020,9,22,9,16,45)

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

function generate_flux_comparison(result_matrix::Resultant_Matrix, timeBin, dist_func, whistler_occurence_rate, Egrid, PAgrid,time_name, data_name, ebin_name, start, stop)
    allPrecip, indexArray, allPrecipInitial = precipitatingParticles(result_matrix.tVec, result_matrix.Ematrix, result_matrix.endTime, timeBin);
    f_timeseries, psd_timeseries, psd_prec_timeseries = make_psd_timeseries(result_matrix.Ematrix,result_matrix.PAmatrix,result_matrix.tVec, dist_func, Egrid, PAgrid, whistler_occurence_rate);
    binned_psd_prec_timeseries = bin_psd_prec_timeseries(psd_prec_timeseries, indexArray);
    equatorial_fluxes = calc_equatorial_fluxes(result_matrix.Ematrix,result_matrix.PAmatrix, dist_func, Egrid, PAgrid);
    prec_flux_timeseries = calc_precipitating_flux_timeseries(binned_psd_prec_timeseries);
    elfin_measurements = extract_idl_csv(time_name, data_name, ebin_name, start, stop);
    return equatorial_fluxes, elfin_measurements, prec_flux_timeseries
end

function generate_trapped_psd(dist_func, whistler_occurence_rate::Float64)
    f_timeseries, psd_timeseries, psd_prec_timeseries = make_psd_timeseries(Ematrix, PAmatrix, tVec, dist_func, Egrid, 8:2:90, whistler_occurence_rate)
    return calc_precipitating_flux_timeseries(bin_psd_prec_timeseries([psd_timeseries[particle][:,1] for particle in 1:length(tVec)-1], indexArray))
end

# Useful helpers
logrange(x1, x2, n::Int64) = [10^y for y in range(log10(x1), log10(x2), length=n)]

# df = DataFrame([allPrecip], :auto)
# CSV.write("for_james.csv",df)