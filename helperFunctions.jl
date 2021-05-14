@info "Compiling packages for running model..."
#meta
using TickTock
using ConfParser
using Profile
#sim
using Dates
using Random
using StaticArrays
using Distributed
@everywhere using OrdinaryDiffEq
using JLD2
using Plots
@info "Packges compiled, running model"

#######################
## Constants n stuff ##
#######################

const Re   = 6370e3;        # Earth radius, f64
const c    = 3e8;           # speedo lite, f64
const Beq  = 3.e-5;         # B field at equator (T), f64
const maxWave = deg2rad(1); # shoreeeeee

# Parsing the conf file
conf = ConfParse(conffile)
parse_conf!(conf)
basename = retrieve(conf, "basename");
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

####################
## Initialization ##
####################

function generateFlatParticleDistribution(numParticles::Int64, ICrange, z0=0::Float64, lambda0=0::Float64)
    ELo, EHi, Esteps, PALo, PAHi, PAsteps = ICrange
    @info "Generating a flat particle distribution with"
    @info "Energy from $ELo KeV to $EHi KeV in $Esteps KeV increments"
    @info "PA from $PALo deg to $PAHi deg in $PAsteps deg increments"

    nBins = ((EHi-ELo)÷Esteps + 1) * ((PAHi-PALo)÷PAsteps + 1)
    N = numParticles ÷ nBins # num of particles per bin
    @info "Flat distribution with $N particles/bin in $nBins bins"

    if numParticles%nBins != 0 @warn "Truncating $(numParticles%nBins) particles for an even distribution"
    end
    if iszero(N)
        N = 1;
        @warn "Use higher number of particles next time. Simulating 1 trajectory/bin."
        @warn "Minimum number of particles to simulate is 1 particles/bin."
    end

    @views f0 = [[(E+511.)/511. deg2rad(PA)] for PA in PALo:PAsteps:PAHi for E in ELo:Esteps:EHi for i in 1:N] # creates a 2xN array with initial PA and Energy
    #####       [[z0 pz0                          zeta0       mu0                          lambda0]             ]
    @views h0 = [[z0 sqrt(IC[1]^2 - 1)*cos(IC[2]) rand()*2*pi .5*(IC[1]^2-1)*sin(IC[2])^2  lambda0] for IC in f0] # creates a 5xN array with inital h0 terms
    f0 = vcat(f0...) # convert Array{Array{Float64,2},1} to Array{Float64,2}
    h0 = vcat(h0...) # since i used list comprehension it is now a nested list

    # Other ICs that are important
    # Define basic ICs and parameters
    B0          = Beq*sqrt(1. +3. *sin(lambda0)^2.)/L^3.;     # starting B field at eq
    Omegace0    = (1.6e-19*B0)/(9.11e-31);                    # electron gyrofreq @ the equator
    # todo, make lambda IC a funtion of z and L, get rid of dep on Beq
    eta         = Omegace0*L*Re/c;              # should be like 10^3
    epsilon     = waveAmplitudeModifier/eta;    # normalized wave large amplitude, .1 for small, 15 for large
    resolution  = .1/eta;                       # determines max step size of the integrator
    @info @views "Created Initial Conditions for $(length(h0[:,1])) particles"
    
    return h0, f0, eta, epsilon, resolution;
end

@everywhere function generateModifiableFunction(batches)
    ```
    Takes in the initial condition and splits them into batches.
    This way, we can feed each one in and get a percent completeness during sim.
    ```
    probGeneratorList = [];
    nPerBatch = numParticles÷batches;
    for j in 0:batches-1
        truncatedIC =  h0[nPerBatch*j+1:(nPerBatch*j+nPerBatch),:]
        push!(probGeneratorList, ((prob,i,repeat) -> remake(prob, u0 = truncatedIC[i,:], p = @SVector [eta, epsilon, Omegape, omegam])))
    end
    percentage = (round(100/batches))
    @info "Each batch will simulate $nPerBatch particles and correspond with $percentage%"
    return probGeneratorList, nPerBatch, percentage
end



##################
## Math n stuff ##
##################

function eom!(dH,H,p::SVector{4, Float64},t::Float64)
    ```
    These equations define the motion.
    ```
    @fastmath u = @views .5*(tanh(H[5]/deg2rad(1))+1); # single sided wave
    # u=1; # or not
    
    @fastmath sinlambda = @views sin(H[5]);
    @fastmath coslambda = @views cos(H[5]);
    @fastmath coszeta = @views cos(H[3]);

    @fastmath b = sqrt(1+3*sinlambda^2)/(coslambda^6);
    @fastmath db = @views (3*(27*sinlambda-5*sin(3*H[5])))/(coslambda^8*(4+12*sinlambda^2));
    @fastmath gamma = @views sqrt(1 + H[2]^2 + 2*H[4]*b);
    @fastmath K = @views (p[3] * (coslambda^-(5/2)))/sqrt(b/p[4] - 1);
    # eta = Omegace*L*Re/c;
    
    @fastmath etaepsilonusqrt2mubovergamma = @views p[1]*p[2]*u*sqrt(2*H[4]*b)/gamma

    @fastmath dH .= @views SizedVector{5}([ H[2]/gamma,
                    -(H[4]*db)/gamma - (etaepsilonusqrt2mubovergamma*coszeta),
                    p[1]*(K*dH[1] - p[4] + b/gamma) + (etaepsilonusqrt2mubovergamma*sin(H[3]))/(2*H[4]*K),
                    -(etaepsilonusqrt2mubovergamma*coszeta)/(K),
                    H[2]/(gamma*coslambda*sqrt(1+3*sinlambda^2)) ])
end

# Simple calcs
calcb(lambda::Float64)  = sqrt(1+3*sin(lambda)^2)/(cos(lambda)^6)
calcdb(lambda::Float64) = (3*(27*sin(lambda)-5*sin(3*lambda)))/(cos(lambda)^8*(4+12*sin(lambda)^2))
calcGamma(pz::Float64,mu::Float64,b::Float64) = sqrt(1 + pz^2 + 2*mu*b)
calcK(b::Float64,lambda::Float64) = (Omegape * (cos(lambda)^-(5/2)))/sqrt(b/omegam - 1)
calcAlpha(mu::Float64, gamma::Float64) = rad2deg(asin(sqrt((2*mu)/(gamma^2 - 1))))


###############################
## Post Processing functions ##
###############################

function loadSingleSolution(inFileName::String)
    ```
    Loads in the JLD2 file and creates timeseries of all the data.
    Takes in filename in string, returns 5 TxN vectors for time,
    z, pz, pitch angle, and energy timeseries for all N particles.
    ```
    JLD2.@load jldname = inFileName sol # load desired variables

    allZ = Vector{Vector{Float64}}();
    allPZ = Vector{Vector{Float64}}();
    allE = Vector{Vector{Float64}}();
    allPA = Vector{Vector{Float64}}();
    allT = Vector{Vector{Float64}}();
    
    for traj in sol # TODO make this multithreaded
        vars = hcat(traj.u...) # pulls out the canonical position/momentum
        # z        = vars[1,:]; 
        # pz       = vars[2,:];
        # zeta     = vars[3,:];
        # mu       = vars[4,:];
        # lambda   = vars[5,:];
        @views gamma = calcGamma.(vars[2,:],vars[4,:],calcb.(vars[5,:]));
        @views push!(allZ, vars[1,:]);
        @views push!(allPZ, vars[2,:]);
        push!(allT, traj.t);
        @views push!(allPA, calcAlpha.(vars[4,:], gamma));
        push!(allE, (511*(gamma .- 1)));
    end
    @info "All variables loaded in"
    return allZ, allPZ, allT, allPA, allE
end

function loadSolutions(inFileBaseName::String, batches::Int64)
    ```
    Loads in the JLD2 file and creates timeseries of all the data.
    Takes in filename in string, returns 5 TxN vectors for time,
    z, pz, pitch angle, and energy timeseries for all N particles.
    ```
    allZ = Vector{Vector{Float64}}();
    allPZ = Vector{Vector{Float64}}();
    allE = Vector{Vector{Float64}}();
    allPA = Vector{Vector{Float64}}();
    allT = Vector{Vector{Float64}}();
    for i in 1:batches
        JLD2.@load jldname = string("data/10000particles/",inFileBaseName,"$i.jld2") sol
        
        for traj in sol # TODO make this multithreaded
            vars = hcat(traj.u...) # pulls out the canonical position/momentum
            # z        = vars[1,:]; 
            # pz       = vars[2,:];
            # zeta     = vars[3,:];
            # mu       = vars[4,:];
            # lambda   = vars[5,:];
            @views gamma = calcGamma.(vars[2,:],vars[4,:],calcb.(vars[5,:]));
            @views push!(allZ, vars[1,:]);
            @views push!(allPZ, vars[2,:]);
            push!(allT, traj.t);
            @views push!(allPA, calcAlpha.(vars[4,:], gamma));
            push!(allE, (511*(gamma .- 1)));
        end
    end
    @info "All variables loaded in"
    return allZ, allPZ, allT, allPA, allE
end

function countLostParticles(allT)
    ```
    Based on time vectors, counts which ones were lost
    and at what time. Returns a Nx2 array where the first
    column is the time at which the particle was lost, and
    the 2nd column denotes the number lost at that point.
    ```
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
    #lostParticles[1:end .!= 6,: ] # clever 1 liner to rid a row
    @info "Total of $(lostParticles[end,end]) particles lost during sim"
    return lostParticles
end

function postProcessor(allT, allZ, allPZ, allE, allPA)
    ```
    This function will take the output of the model and convert them into usable m x n matrices
    where m is max number of timesteps for the longest trajectory and N is number of particles
    Arrays that are not m long are filled with NaNs
    ```
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
        annotate!(80, 1.05*maxParticles, text("t = $(round(tVec[i]*Re/(c),digits=4)) s"), :right)
    end every animDec
    savename = string("plots/",gifFileName)
    gif(anim, savename, fps = (length(tVec)/animDec)/(animScale*endTime*Re/(c)))
end

function animateESD(gifFileName="ESDanimation.gif", maxParticles=1500, binwidth=50)
    pyplot()
    animDec = 10; # make a png for animation every 10 points
    animScale = 10; # i.e. animscale = 10 means every 10 seconds in animation is 1 second of simulation time (increase for longer animation)
    anim = @animate for i=1:length(tVec)
        histogram(Ematrix[i,:],bins=(0:binwidth:2000), ylim = [0, 1.1*maxParticles], #xminorticks = 4, yminorticks = 2,
            legend = false, dpi = 100, xscale=:log10, xticks = ([10,100,1000],["10", "100", "1000"]), xlim = [10,2000],
            xlabel="Energy (KeV)", ylabel="Number of Particles", title="Time evolution of Energy Spectra");
        annotate!(20, 1.05*maxParticles, text("t = $(round(tVec[i]*Re/(c),digits=4)) s"), :left)
    end every animDec
    savename = string("plots/",gifFileName)
    gif(anim, savename, fps = (length(tVec)/animDec)/(animScale*endTime*Re/(c)))
end

function animatePSD(gifFileName="PSDanimation.gif", PAbinwidth=5, Ebinwidth=50, maxParticles=150)
    pyplot()
    animDec = 10; # make a png for animation every 10 points
    animScale = 10; # i.e. animscale = 10 means every 10 seconds in animation is 1 second of simulation time (increase for longer animation)
    anim = @animate for i=1:length(tVec)
        histogram2d(PAmatrix[i,:], Ematrix[i,:], bins=(0:PAbinwidth:90,0:Ebinwidth:2000),
            xlabel="Pitch Angle (deg)",ylabel="Energy (keV)",title="Time evolution of Phase Space Density",
            yscale=:log10, clims = (0,maxParticles), minorgrid=true,
            yticks = ([1, 10,100,1000],["1", "10", "100", "1000"]));
        annotate!(lossConeAngle, 5, text("Lost Particle Count:", :left));
        annotate!(lossConeAngle, 3.5, text((count(i->isnan(i), Ematrix[i,:])), :left))
        annotate!(70, 5, text("t = $(round(tVec[i]*Re/(c),digits=4)) s"), :right)
        # count(i->isnan(i), Ematrix[1,:]) # counts number of NaNs (i.e particles lost) at certain time
    end every animDec
    savename = string("plots/",gifFileName)
    gif(anim, savename, fps = (length(tVec)/animDec)/(animScale*endTime*Re/(c)))
end
            
function recalcDistFunc(Ematrix::Array{Float64,2},PAmatrix::Array{Float64,2},initial::Int64,final::Int64,distFunc,
    Egrid::StepRange{Int64,Int64}, PAgrid::StepRange{Int64,Int64})
    ```
    Takes in matrix of Energy and PA, initial and final indices, and grid values as stepranges.
    Recalculates out a new PSD given a a distribution function at the initial and final indices.
    ```
    N = @views length(Ematrix[1,:]) # num particles
    EPA0i = @views round.([Ematrix[1,:] PAmatrix[1,:]])
    Nkl = @views N/(length(unique(EPA0i[:,1]))*length(unique(EPA0i[:,2]))) # obtain number of particles/bin

    EPAinitial = @views @. round([Ematrix[initial,:] PAmatrix[initial,:]])
    EPAfinal = @views @. round([Ematrix[final,:] PAmatrix[final,:]])

    # initialize matrices to be used based on grid size
    f = zeros(length(Egrid), length(PAgrid));
    psd_init = zeros(length(Egrid), length(PAgrid));
    psd_final = zeros(length(Egrid), length(PAgrid));
    # initialize vectors to be filled in
    f0Vec = Vector{Float64}();
    indices = Vector{Tuple{Int64,Int64}}();
    excludedParticles = 0;
    
    ### HANDLE NANS HERE ###
    # diff = @. abs(EPAinitial - EPAfinal)
    valid_rows = @. ~isnan(EPAfinal)[:,1] & ~isnan(EPAfinal[:,2]); # non-lost particles from final state
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

    psdVec = [(f0Vec[i]/f[indices[i][1],indices[i][2]]) for i in 1:length(f0Vec[:,1])]

    for i in 1:length(EPAinitial[:,1])
        if ~(EPAfinal[i,1]>maximum(Egrid) || EPAfinal[i,2]>maximum(PAgrid)) # skip loop if data is outside of range
            k,l = 1, 1;
            while Egrid[k] < EPAinitial[i,1]; k+=1; end
            while PAgrid[l] < EPAinitial[i,2]; l+=1; end
            psd_init[k-1,l-1] += psdVec[i]
            k,l = 1, 1;
            while Egrid[k] < EPAfinal[i,1]; k+=1; end
            while PAgrid[l] < EPAfinal[i,2]; l+=1; end
            psd_final[k-1,l-1] += psdVec[i]
        else
            excludedParticles += 1;
        end
    end
    @info "Excluded $excludedParticles particles due to out of range."
    return f, psd_init, psd_final
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

function checkDistFunc(f::Array{Float64,2}, psd_init::Array{Float64,2}, psd_final::Array{Float64,2}, initial::Int64, final::Int64)
    # origPlot =  heatmap(0:5:90, 0:50:1000, f, fc = :plasma,
    #     xlabel="Pitch Angle (deg)",ylabel="Energy (keV)", title = "Original Flat PSD at t = 0");

    # initPlot =  heatmap(0:5:90, 0:50:1000, log10.(psd_init),clims=(-5.5,-2), fc = :plasma,
    #     xlabel="Pitch Angle (deg)",ylabel="Energy (keV)", title = "Recalculated PSD at t = $(round(tVec[initial]*Re/(c),digits=3)) s");

    finalPlot = heatmap(0:5:90, 0:50:1000, log10.(psd_final), fc = :plasma, colorbar = false,
        xlabel="Pitch Angle (deg)",ylabel="Energy (keV)", title = "Recalculated PSD at t = $(round(tVec[final]*Re/(c),digits=3)) s");

    # origPlotAnime = marginalhist(round.(PAmatrix[final,:]), log10.(round.(Ematrix[final,:])),
    # xlim = (0,90), yscale = :lin ,yticks = ([1, 2,3,3.70],["10", "100", "1000", "5000"]),
    # fc =:plasma, bins= 30);

    plot(finalPlot)

    # plot(origPlot, initPlot, origPlotAnime, finalPlot, layout = (2,2), size = (1200,1200))
end


function animateNewPSD(gifFileName="NewPSDanimation.gif")
    pyplot()
    animDec = 10; # make a png for animation every 10 points
    animScale = 10; # i.e. animscale = 10 means every 10 seconds in animation is 1 second of simulation time (increase for longer animation)
    initial, final = 1, 1
    f, psd_init, psd_final = recalcDistFunc(Ematrix,PAmatrix,initial,final,f0,0:50:1000,0:5:90);
    anim = @animate for i=1:length(tVec)
        _,_,psd_final = recalcDistFunc(Ematrix,PAmatrix,initial,i,f0,0:50:1000,0:5:90);
        checkDistFunc(f, psd_init, psd_final, initial, i)
        # heatmap(0:5:90, 0:50:1000, psd_final, clim = (0,.01), title = "PSD at t = $(round(tVec[i]*Re/(c),digits=3)) s")
    end every animDec
    savename = string("plots/",gifFileName)
    gif(anim, savename, fps = (length(tVec)/animDec)/(animScale*endTime*Re/(c)))
end




# ```
# Useful links
# Community stuff: https://discourse.julialang.org/c/domain/models/21
# Docs for ensemble stuff: https://diffeq.sciml.ai/dev/features/ensemble/
# Marginal Histograms: http://docs.juliaplots.org/latest/recipes/#Marginal-Histograms
# Dumb colorbar ticks issue: https://github.com/JuliaPlots/Plots.jl/issues/2308
# ```