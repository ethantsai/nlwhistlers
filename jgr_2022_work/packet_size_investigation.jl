

idxno = [rand(1:themis_nopkt.numParticles) for i in 1:10]
p1 = plot((5.2*themis_nopkt.allT[idxno]*Re)./(c), themis_nopkt.allPA[idxno], title = "Long packet (dPhi = 1)",
legend=false,xlabel="Time (ms)", ylabel="Energy (keV)",
xlim = [0,12], ylim = [0,90])
idx30 = [rand(1:themis_hilat.numParticles) for i in 1:10]
p2 = plot((5.2*themis_hilat.allT[idxno]*Re)./(c), themis_hilat.allPA[idxno], title = "Short Packet (dPhi = 30)",
legend=false,xlabel="Time (ms)", ylabel="Energy (keV)",
xlim = [0,12], ylim = [0,90])
idx100 = [rand(1:themis_100pkt.numParticles) for i in 1:10]
p3 = plot((5.2*themis_100pkt.allT[idxno]*Re)./(c), themis_100pkt.allPA[idxno], title = "Short Packet (dPhi = 100)",
legend=false,xlabel="Time (ms)", ylabel="Energy (keV)",
xlim = [0,12], ylim = [0,90])
bigplot = plot(p1,p2,p3, dpi = 300, layout = (3,1))


ylimmax = maximum([themis_nopkt.lostParticles[end,2], themis_hilat.lostParticles[end,2], themis_100pkt.lostParticles[end,2]])
p1 = plot(L*themis_nopkt.lostParticles[:,1]*Re/(c), themis_nopkt.lostParticles[:,2], legend=false,
    xlim = [0,L*maximum(maximum(themis_nopkt.allT))*Re/(c)], ylim = [0,ylimmax],
    title = "Number of Lost Particles (dphi = 1)", xlabel="Time (s)", ylabel="Number of Particles")
p2 = plot(L*themis_hilat.lostParticles[:,1]*Re/(c), themis_hilat.lostParticles[:,2], legend=false,
    xlim = [0,L*maximum(maximum(themis_hilat.allT))*Re/(c)], ylim = [0,ylimmax],
    title = "Number of Lost Particles (dphi = 30)", xlabel="Time (s)", ylabel="Number of Particles")
p3 = plot(L*themis_100pkt.lostParticles[:,1]*Re/(c), themis_100pkt.lostParticles[:,2], legend=false,
    xlim = [0,L*maximum(maximum(themis_100pkt.allT))*Re/(c)], ylim = [0,ylimmax],
    title = "Number of Lost Particles (dphi = 100)", xlabel="Time (s)", ylabel="Number of Particles")
bigplot = plot(p1,p2,p3, dpi = 300, layout = (3,1))

p1 = histogram(diff(one.lostParticles[:,2]),L*one.lostParticles[1:end-1,1]*Re/(c), bins= 100,ylim = (0,250))
p2 = histogram(diff(themis_hilat.lostParticles[:,2]),L*themis_hilat.lostParticles[1:end-1,1]*Re/(c), bins= 100,ylim = (0,250))
p3 = histogram(diff(themis_100pkt.lostParticles[:,2]),L*themis_100pkt.lostParticles[1:end-1,1]*Re/(c), bins= 100,ylim = (0,250))
bigplot = plot(p1,p2,p3, dpi = 300, layout = (3,1))
plot!(xlim  = [2.5,10], ylim = [0,75])

idxno = [rand(1:one.numParticles) for i in 1:10]
p1 = plot((5.2*one.allT*Re)./(c), one.allE, title = "Long packet (dPhi = 1)",
legend=false,xlabel="Time (ms)", ylabel="Energy (keV)",
xlim = [0,1], ylim = [0,1000])
idx30 = [rand(1:hund.numParticles) for i in 1:10]
p2 = plot((5.2*hund.allT*Re)./(c), hund.allE, title = "Short Packet (dPhi = 100)",
legend=false,xlabel="Time (ms)", ylabel="Energy (keV)",
xlim = [0,12], ylim = [0,1000])
Eplot = plot(p1,p2, dpi = 300, layout = (2,1))
p3 = plot((5.2*one.allT*Re)./(c), one.allPA, title = "Long packet (dPhi = 1)",
legend=false,xlabel="Time (ms)", ylabel="PA (deg)",
xlim = [0,12], ylim = [0,90])
idx30 = [rand(1:hund.numParticles) for i in 1:10]
p4 = plot((5.2*hund.allT*Re)./(c), hund.allPA, title = "Short Packet (dPhi = 100)",
legend=false,xlabel="Time (ms)", ylabel="PA (deg)",
xlim = [0,12], ylim = [0,90])
PAPlot = plot(p3,p4, dpi = 300, layout = (2,1))




p1 = plot((5.2*one.allT*Re)./(c), one.allPA, title = "Short packet (dPhi = 1)",
legend=false,xlabel="Time (ms)", ylabel="PA (deg)",
xlim = [0,0.5], ylim = [0,15])
p2 = plot((5.2*thi.allT*Re)./(c), thi.allPA, title = "Normal Packet (dPhi = 30)",
legend=false,xlabel="Time (ms)", ylabel="PA (deg)",
xlim = [0,0.5], ylim = [0,15])
p3 = plot((5.2*hun.allT*Re)./(c), hun.allPA, title = "Long Packet (dPhi = 100)",
legend=false,xlabel="Time (ms)", ylabel="PA (deg)",
xlim = [0,0.5], ylim = [0,15])
bigplot = plot(p1,p2,p3, dpi = 300, layout = (3,1))

one = load_resultant_matrix("test1", "jgr_2022_work/results/jld2_211124_11", "PA_TEST_1_100", "setupasrun.conf", 1);
thi = load_resultant_matrix("test30", "jgr_2022_work/results/jld2_211124_09", "PA_TEST_30_100", "setupasrun.conf", 1);
hun = load_resultant_matrix("test100", "jgr_2022_work/results/jld2_211124_11", "PA_TEST_100_100", "setupasrun.conf", 1);

rm = hun
ipa = PAmatrix[1,:]
epa = Vector{Float64}()
for col in eachcol(PAmatrix)
    push!(epa,col[findall(!isnan, col)][end])
end
scatter!(ipa, epa)





directory="jgr_2022_work/results/jld2_211124_22"
basename = "PA_TEST_100_n_100"
allZ = Vector{Vector{Float64}}();
allPZ = Vector{Vector{Float64}}();
allE = Vector{Vector{Float64}}();
allPA = Vector{Vector{Float64}}();
allT = Vector{Vector{Float64}}();
allG = Vector{Vector{Float64}}();
# allLambda = Vector{Vector{Float64}}();
# allPsinou = Vector{Vector{Float64}}();
for i in 1:1
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
        # @views push!(allG, vars[:,7])
        # @views push!(allLambda, vars[:,5])
        # @views push!(allPsinou, vars[:,8])
    end
end
@time tVec, Zmatrix, PZmatrix, Ematrix, PAmatrix = postProcessor(allT, allZ, allPZ, allE, allPA);
ipa = PAmatrix[1,:]
epa = Vector{Float64}()
for col in eachcol(PAmatrix)
    push!(epa,col[findall(!isnan, col)][end])
end
scatter!(ipa, epa, legend=false)




p1 = plot((5.2*one.allT*Re)./(c), one.allPA, title = "Short packet (dPhi = 1)",
legend=false,xlabel="Time (ms)", ylabel="PA (deg)",
xlim = [0,0.5], ylim = [0,15])
p2 = plot((5.2*thi.allT*Re)./(c), thi.allPA, title = "Normal Packet (dPhi = 30)",
legend=false,xlabel="Time (ms)", ylabel="PA (deg)",
xlim = [0,0.5], ylim = [0,15])
p3 = plot(allLambda, allPsi, title = "Long Packet (dPhi = 100)",
legend=false,xlabel="Time (ms)", ylabel="Psi (rad)",
xlim = [-pi/4, pi/4], ylim = [-5,5])
bigplot = plot(p1,p2,p3, dpi = 300, layout = (3,1))