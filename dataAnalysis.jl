### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 9da308cd-880e-45ca-9c9f-79315cb575c6
using InteractiveUtils, PlutoUI, TickTock, ConfParser, Dates, Random, StaticArrays, Distributed, OrdinaryDiffEq, JLD2, Plots, LoopVectorization, BenchmarkTools

# ╔═╡ fe4cd7b2-b5d2-11eb-29f5-81ec8d2a9a9e
md"
## Particle Tracing Simulation of Non-Linear Resonant Electron Scattering by Chorus Waves
#### General notes

 - default ingestion of data will be via `.jld2` files in a particular subdirectory
 
#### There's plenty of plots

Lost electrons are cool.

"

# ╔═╡ 9834a489-ea0f-45b7-a299-74454b4b3ace
md"
conf file = $(@bind conffile TextField())

directory = $(@bind dirname TextField())

basename  = $(@bind basename TextField())

batches   = $(@bind num_batches NumberField(1:100; default=1))
"

# ╔═╡ 8fa0edf6-5276-436f-bb76-e4e90429d0f2
begin
	const Re   = 6370e3;        # Earth radius, f64
	const c    = 3e8;           # speedo lite, f64
	const Beq  = 3.e-5;         # B field at equator (T), f64

	# Parsing the conf file
	conf = ConfParse(dirname*"/"*conffile)
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
end

# ╔═╡ 59c01c3b-effc-4260-be8f-981feda851c8
begin
	
	calcb(lambda::Float64)  = sqrt(1+3*sin(lambda)^2)/(cos(lambda)^6)
	
	calcdb(lambda::Float64) = (3*(27*sin(lambda)-5*sin(3*lambda))) / (cos(lambda)^8*(4+12*sin(lambda)^2))
	
	calcGamma(pz::Float64,mu::Float64,b::Float64) = sqrt(1 + pz^2 + 2*mu*b)
	
	calcK(b::Float64,lambda::Float64) = (Omegape * (cos(lambda)^-(5/2)))/sqrt(b/omegam - 1)
	
	calcAlpha(mu::Float64, gamma::Float64) = rad2deg(asin(sqrt((2*mu)/(gamma^2 - 1))))

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

		Threads.@threads for i in 1:num_batches
			JLD2.@load jldname = directory*"/"*basename*"_$i.jld2" sol

			for traj in sol # TODO make this multithreaded
				vars = hcat(traj.u...) # pulls out the canonical position/momentum
				# z        = vars[1,:]; 
				# pz       = vars[2,:];
				# zeta     = vars[3,:];
				# mu       = vars[4,:];
				# lambda   = vars[5,:];
				b = @views @avx @. sqrt(1+3*sin(vars[5,:])^2)/(cos(vars[5,:])^6)
				gamma = @views @avx @. sqrt(1 + vars[2,:]^2 + 2*vars[4,:]*b)
								
				@views push!(allT, traj.t);
				@views push!(allZ, vars[1,:]);
				@views push!(allPZ, vars[2,:]);
				@views push!(allPA, @. rad2deg( @avx asin(sqrt((2*vars[4,:])/(gamma^2 - 1)))));
				@views push!(allE, @. (511*(gamma - 1)));
			end
			@info "$(length(sol)) particles loaded in from $(basename*"_$i.jld2")"
		end
		@info "Loaded total $(length(allT)) particles!"
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

end

# ╔═╡ de04d879-e805-4d27-8383-6398fb44872d
allZ, allPZ, allT, allPA, allE = loadData(dirname, basename, num_batches);

# ╔═╡ 7af09807-4012-46e2-bb2e-d307fa092b28
# @time countLostParticles(allT)

# ╔═╡ dbd59420-c97c-48e0-81e4-57bd72a5e52b


# ╔═╡ Cell order:
# ╟─fe4cd7b2-b5d2-11eb-29f5-81ec8d2a9a9e
# ╠═9da308cd-880e-45ca-9c9f-79315cb575c6
# ╟─9834a489-ea0f-45b7-a299-74454b4b3ace
# ╠═8fa0edf6-5276-436f-bb76-e4e90429d0f2
# ╠═59c01c3b-effc-4260-be8f-981feda851c8
# ╠═de04d879-e805-4d27-8383-6398fb44872d
# ╠═7af09807-4012-46e2-bb2e-d307fa092b28
# ╠═dbd59420-c97c-48e0-81e4-57bd72a5e52b
# ╠═27be5594-8ad4-48e4-8a3d-c06c9e712311
