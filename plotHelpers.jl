using InteractiveUtils
using PlutoUI
using TickTock
using ConfParser
using Dates
using Random
using StaticArrays
using Distributed
using OrdinaryDiffEq
using JLD2
using Plots
using LoopVectorization
using BenchmarkTools
using StatsPlots

#######################
## Constants n stuff ##
#######################

const Re   = 6370e3;        # Earth radius, f64
const c    = 3e8;           # speedo lite, f64
const Beq  = 3.e-5;         # B field at equator (T), f64

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

function countLostParticles(allT::Vector{Vector{Float64}})
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

function postProcessor2(allT::Vector{Vector{Float64}},
                        # allZ::Vector{Vector{Float64}},
                        # allPZ::Vector{Vector{Float64}},
                        allE::Vector{Vector{Float64}},
                        allPA::Vector{Vector{Float64}}
                        )
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
        # @views Zmatrix[1:length(allT[i]),i] = allZ[i]
        # @views PZmatrix[1:length(allT[i]),i] = allPZ[i]
        @views Ematrix[1:length(allT[i]),i] = allE[i]
        @views PAmatrix[1:length(allT[i]),i] = allPA[i]
    end

    return tVec, Ematrix, PAmatrix
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
    Egrid::Vector{Float64}, PAgrid::StepRange{Int64,Int64})
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
        @info "$(length(ones(N)[prec_rows])) particles are about to precipitate"
    # end
    EPAinitial = EPAinitial[vec(valid_rows), :]; # new vec with no NaNs
    EPAfinal = EPAfinal[vec(valid_rows), :]; # new vec with no NaNs
    
    lostParticles = length(valid_rows) - length(EPAfinal[:,1]);
    @info "$lostParticles particles lost at index $final."

    for EPAi in eachrow(EPAinitial)
        k,l = 1, 1;
        while Egrid[k] < EPAi[1]; k+=1; end # energy
        while PAgrid[l] < EPAi[2]; l+=1; end # PA
        push!(indices, (k-1,l-1))
        push!(f0Vec, distFunc(1., EPAi[1], EPAi[2]))
        f[k-1,l-1] += 1;
    end

    psdVec = [(f0Vec[i]/f[indices[i][1],indices[i][2]]) for i in eachindex(f0Vec[:,1])]

    for i in eachindex(EPAinitial[:,1]) # i iterates over each particle
        if ~(EPAfinal[i,1]>maximum(Egrid) || EPAfinal[i,2]>maximum(PAgrid)) # skip loop if data is outside of range
            k,l = 1, 1;
            while Egrid[k] < EPAinitial[i,1]; k+=1; end
            while PAgrid[l] < EPAinitial[i,2]; l+=1; end
            psd_init[k-1,l-1] += psdVec[i]
            k,l = 1, 1;
            while Egrid[k] < EPAfinal[i,1]; k+=1; end
            while PAgrid[l] < EPAfinal[i,2]; l+=1; end
            psd_final[k-1,l-1] += psdVec[i]

            # handle precipitaitng particles here
            if prec_indices[i] == 1.0
                psd_prec[k-1] += psdVec[i]
            end
        else
            excludedParticles += 1;
        end
    end
    @info "Excluded $excludedParticles particles due to out of range."

    return f, psd_init, psd_final, psd_prec
end

function make_psd_timeseries(Ematrix,PAmatrix,tVec, dist_func, Egrid, PAgrid)
    psd_timeseries = Vector{Matrix{Float64}}()
    f_timeseries = copy(psd_timeseries)
    psd_prec_timeseries = Vector{Vector{Float64}}()
    @inbounds for time_index in eachindex(tVec[1:end-1])
        f, _, psd_final, psd_prec = recalcDistFunc(Ematrix, PAmatrix, 1, time_index, dist_func, Egrid, PAgrid);
        push!(f_timeseries, f)
        push!(psd_timeseries, psd_final)
        push!(psd_prec_timeseries, psd_prec)
    end
    return f_timeseries, psd_timeseries, psd_prec_timeseries
end

function bin_psd_prec_timeseries(psd_prec_timeseries, indexArray)
    binned_psd_prec_timeseries = Vector{Vector{Float64}}() 
    index = 0
    # loop through each time bin
    for time_bin_index in indexArray
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

function animate_a_thing(gifFileName::String, thing::Vector{Matrix{Float64}})
    pyplot()
    animDec = 10; # make a png for animation every 10 points
    animScale = 1; # i.e. animscale = 10 means every 10 seconds in animation is 1 second of simulation time (increase for longer animation)
    initial, final = 1, 1
    anim = @animate for i in eachindex(thing)
        heatmap(PAgrid, Egrid, log10.(thing[i]), fc = :plasma, colorbar = false, 
            xlabel="Pitch Angle (deg)",ylabel="Energy (keV)", title = "Recalculated PSD at t = $(round(tVec[i]*Re*L/(c),digits=2)) s")
    end every animDec
    savename = string("results/",gifFileName)
    gif(anim, savename, fps = (length(tVec)/animDec)/(animScale*endTime*Re*L/(c)))
    return anim
end

function animate_a_thing(gifFileName::String, thing::Vector{Vector{Float64}})
    pyplot()
    animDec = 1; # make a png for animation every 10 points
    animScale = 50; # i.e. animscale = 10 means every 10 seconds in animation is 1 second of simulation time (increase for longer animation)
    maxEnergy=1000
    anim = @animate for i in eachindex(thing[1:end-1])
        plot(Egrid, thing[1:end-1], color = :gray, alpha = .5);
        plot!(Egrid,thing[i], color = :orange);
        plot!(ylim =(10e-10,10), xlim=(20,maxEnergy), yscale=:log10, xscale=:log10, legend=false);
        plot!(xlabel="Energy (keV)", ylabel="PSD", title="Energy Spectra of Precipitating Particles");
        annotate!(40, 0.1*maxPSD, text("t = $(round(tVec[indexArray[i]]*Re*L/(c),digits=3)) s"), :left)
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

function checkDistFunc(f::Array{Float64,2}, psd_init::Array{Float64,2}, psd_final::Array{Float64,2}, initial::Int64, final::Int64, Egrid::StepRange{Int64,Int64}, PAgrid::StepRange{Int64,Int64})
    origPlot =  heatmap(PAgrid, Egrid, f, fc = :plasma,
        xlabel="Pitch Angle (deg)",ylabel="Energy (keV)", title = "Original Flat PSD at t = 0");

    initPlot =  heatmap(PAgrid, Egrid, log10.(psd_init),clims=(-5.5,-2), fc = :plasma,
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


function animateNewPSD(gifFileName, Egrid::StepRange{Int64,Int64}, PAgrid::StepRange{Int64,Int64})
    pyplot()
    animDec = 1; # make a png for animation every 10 points
    animScale = 10; # i.e. animscale = 10 means every 10 seconds in animation is 1 second of simulation time (increase for longer animation)
    initial, final = 1, 1
    f, psd_init, psd_final = recalcDistFunc(Ematrix,PAmatrix,initial,final,f0,Egrid,PAgrid);
    anim = @animate for i in eachindex(tVec)
        _,_,psd_final = recalcDistFunc(Ematrix,PAmatrix,initial,i,f0,Egrid,PAgrid);
        # checkDistFunc(f, psd_init, psd_final, initial, i,Egrid, PAgrid)
        heatmap(PAgrid, Egrid, log10.(psd_final), fc = :plasma, colorbar = false, clims=(-5.5,-2),
            xlabel="Pitch Angle (deg)",ylabel="Energy (keV)", title = "Recalculated PSD at t = $(round(tVec[i]*Re*L/(c),digits=2)) s")
    end every animDec
    savename = string("results/",gifFileName)
    gif(anim, savename, fps = (length(tVec)/animDec)/(animScale*endTime*Re*L/(c)))
end

function precipitatingParticles(tVec, Ematrix, timeBin=10)
    unitSimTime = mean(diff(tVec[begin:end-1]));
    simTimeBin = convert(Int64, round(timeBin/unitSimTime, digits=0))

    # produce a vector with indices demarcating time bin boundaries
    indexArray = [minimum([i*simTimeBin length(tVec)]) for i in 1:convert(Int64,floor(endTime/timeBin))]
    if length(tVec)!=indexArray[end] push!(indexArray, length(tVec)-1) end #guarantee last index corresponds with proper end
    # remove all particles that didn't get lost
    Ematrixprime = @views Matrix(Ematrix[:,findall(isnan, Ematrix[indexArray[end-1],:])]) # this is a new matrix that only has particles that have precipitated
    # replace last value with second to last value to prevent off by one problems
    if indexArray[end] == length(tVec)
        indexArray[end] = indexArray[end]-1
    end

    allPrecipInitial = Vector{Vector{Float64}}()
    allPrecip = Vector{Vector{Float64}}()
    for i in indexArray[1:end-1]
        @debug "$(length(Ematrixprime[1,:])) remaining particles" size(Ematrixprime)
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



function animatePrecipitatingParticles(gifFileName, allPrecip, indexArray, maxFraction=0.004, maxEnergy=1000)
    animDec = 1; # make a png for animation every 10 points
    animScale = 50; # i.e. animscale = 10 means every 10 seconds in animation is 1 second of simulation time (increase for longer animation)
    anim = @animate for i in eachindex(allPrecip)
        density(allPrecip, color = :gray, alpha = .5);
        density!(allPrecip[i], color = :orange);
        plot!(ylim =(0.,maxFraction), yscale=:identity, xlim=(20,maxEnergy), xscale=:log10, legend=false);
        plot!(xlabel="Energy (keV)", ylabel="Fraction of $numParticles total particles", title="Energy Spectra of Precipitating Particles");
        annotate!(40, 0.1*maxFraction, text("t = $(round(tVec[indexArray[i]]*Re*L/(c),digits=3)) s"), :left)
    end every animDec
    savename = string("results/",gifFileName)
    gif(anim, savename, fps = (length(tVec)/animDec)/(animScale*endTime*Re*L/(c)))
end

function trajectoryChecking(index)
    converT = Re*L/(c)
    @info "starting PA: $(PAmatrix[1,index]), startin E: $(Ematrix[1,index])"
    plot(converT*tVec, PAmatrix[:,index]);
    plot!(converT*tVec, Ematrix[:,index])
    # plot(Zmatrix[:,index], PZmatrix[:,index])
end


function trajectoryTracing(trajectories, animDec)
    converT = Re*L/(c)

    PAanim = @animate for i in eachindex(tVec)
        plot(xlim=(0,12), ylim = (0,90));
        plot!(converT*tVec, PAmatrix[:,trajectories], color=:gray, alpha = .5, label = "");
        for j in trajectories
            plot!(converT*tVec[1:i], PAmatrix[:,j][1:i], label = "")
            scatter!((converT*tVec[i], PAmatrix[:,j][i]), label = "")
        end
    end every animDec

    savename1 = string("results/plots/","PAgif.gif")
    gif(PAanim, savename1, fps = 20)
    
    
    Eanim = @animate for i in eachindex(tVec)
        Eplot = plot(xlim=(0,12), ylim = (10,5000), yscale = :log10)
        plot!(converT*tVec, Ematrix[:,trajectories], color=:gray, alpha = .5, label = "");
        for j in trajectories
            Eplot = plot!(converT*tVec[1:i], Ematrix[:,j][1:i], label = "")
            Eplot = scatter!((converT*tVec[i], Ematrix[:,j][i]), label = "")
        end
    end every animDec

    savename2 = string("results/plots/","Egif.gif")
    gif(Eanim, savename2, fps = 20)
            
    ZPZanim = @animate for i in eachindex(tVec)
        ZPZplot = plot(xlim=(-1,1), ylim = (-7,7))
        for j in trajectories
            ZPZplot = plot!(Zmatrix[:,j][1:i], PZmatrix[:,j][1:i], alpha = max.((1:i) .+ 100 .- i, 0) / 100, label = "")
            ZPZplot = scatter!((Zmatrix[:,j][i], PZmatrix[:,j][i]), label = "")
        end
    end every animDec

    savename3 = string("results/plots/","ZPZgif.gif")
    gif(ZPZanim, savename3, fps = 20)

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


# df = DataFrame([allPrecip], :auto)
# CSV.write("for_james.csv",df)