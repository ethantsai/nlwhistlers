include("agapitovHelpers.jl")
include("agapitovmodel.jl")

@info "Loading constants..."
save_dir = "results_ducting/"
folder = "run16/"

test_cases = [
    "a7_dphi3_PA3-15_constbw300pt_test3_5_41600.jld2"
]

rm_array = Vector{Resultant_Matrix}()
# loading loop™
for case in test_cases
    loadname = save_dir*folder*case
    @info "Loading solution from $loadname..."
    @time @load loadname sol
    label = case
    push!(rm_array, sol2rm(sol, label));
    @info "Loaded $label."
end


# for case #1, omega_m = 0.2, L = 5; contains both 63 and 138 keV
case = 1
rm_case = rm_array[case]
# first, obtain what B_w looks like:
full_length = length(rm_case.tVec)
time = rm_case.tVec
i=1
while full_length != length(rm_case.allBw[i])
    i+=1
end
b_w = rm_case.allBw[i]
lambda = rad2deg.(rm_case.allLambda[i])
plot_bw = plot(time, b_w, xlabel="time", ylabel="B_w", title="Wave amplitude as function of time", legend=false)
plot_bw = plot(lambda, b_w, xlabel="Latitude (deg)", ylabel="B_w", title="Wave amplitude as function of Latitude", legend=false)

# next, isolate just those with starting energy of 63 and 138 keV
initial_E = [rm_case.allE[i][1] for i = 1:length(rm_case.allT)]
final_E = [rm_case.allE[i][end] for i = 1:length(rm_case.allT)]
initial_PA = [rm_case.allPA[i][1] for i = 1:length(rm_case.allT)]
final_PA = [rm_case.allPA[i][end] for i = 1:length(rm_case.allT)]
PA_ROI = findall(x->4.7<x<=5.3, initial_PA) 
E_ROI = findall(x->9<x<30, initial_E) 
ROI = intersect(PA_ROI, E_ROI)

const ELo = 10;
const EHi = 1000;
const Esteps = 32; # double ELFIN E bins
const PALo = 3;
const PAHi = 15;
const PAsteps = 13;

E_bins = logrange(ELo,EHi, Int64(Esteps))
PA_bins = range(PALo, PAHi, 13)

PA_range = [findall(x -> floor(PA-1e-5) < x < ceil(PA+1e-5), initial_PA) for PA in PA_bins]
E_range = [findall(x -> floor(E-1e-5) < x < ceil(E+1e-5), initial_E) for E in E_bins]
permutedims(string.(floor.(E_bins)) .* " keV")

plot1 = scatter(xlabel="Initial Pitch Angle (deg)", ylabel="Final Pitch Angle (deg)", title="L = 5, omega_m = 0.35, dphi = 3", legend=:topleft)
plot1 = scatter!([initial_PA[range] for range in E_range],
                 [final_PA[range] for range in E_range],
                 label = permutedims(string.(floor.(E_bins)) .* " keV"),
                 xlim=[3,16], ylim=[0,30], ms=1, msw=0, legend = :outertopright)
plot1 = plot!(0:30, 0:30, label=false, ls=:dash, lc=:red)


# plot g(phi)

selection = 1:40:432
gofphi = [@. exp(-7 * (cos(PHI/(2*π*3))^2)) for PHI in rm_array[1].allPhi[ROI[selection]]]
BW_ROI = rm_array[1].allBw[ROI[selection]]
interaction = [gofphi[i] .* BW_ROI[i] for i in eachindex(BW_ROI)]
plot(rm_array[1].allT[ROI[selection]], interaction,legend=false)






# for case #1, omega_m = 0.35, L = 5; contains both 63 and 138 keV
case = 2
rm_case = rm_array[case]
# first, obtain what B_w looks like:
full_length = length(rm_case.tVec)
time = rm_case.tVec
i=1
while full_length != length(rm_case.allBw[i])
    i+=1
end
b_w = rm_case.allBw[i]
lambda = rad2deg.(rm_case.allLambda[i])
plot(time, b_w, xlabel="time", ylabel="B_w", title="Wave amplitude as function of time", legend=false)
plot(lambda, b_w, xlabel="Latitude (deg)", ylabel="B_w", title="Wave amplitude as function of Latitude", legend=false)

# next, isolate just those with starting energy of 63 and 138 keV
initial_E = [rm_case.allE[i][1] for i = 1:length(rm_case.allT)]
final_E = [rm_case.allE[i][end] for i = 1:length(rm_case.allT)]
initial_PA = [rm_case.allPA[i][1] for i = 1:length(rm_case.allT)]
final_PA = [rm_case.allPA[i][end] for i = 1:length(rm_case.allT)]
E63_range = findall(x->62<x<64, initial_E)
E138_range = findall(x->137<x<139, initial_E)

plot2 = scatter(xlabel="Initial Pitch Angle (deg)", ylabel="Final Pitch Angle (deg)", title="L = 5, omega_m = 0.35", legend=:topleft)
plot2 = scatter!([initial_PA[E63_range], initial_PA[E138_range]], [final_PA[E63_range], final_PA[E138_range]], label = ["63 keV" "138 keV"], xlim=[0,30], ylim=[0,90], ms=1, msw=0)
plot2 = plot!(0:30, 0:30, label=false, ls=:dash, lc=:red)

histogram2d(initial_PA[E63_range], final_PA[E63_range], xlim=[0,30], ylim=[0,90], nbins=(30,90), clim = (1,1000), cscale=log10)



# for case #3, omega_m = 0.2, L = 6; contains both 63 and 138 keV
case = 3
rm_case = rm_array[case]
# first, obtain what B_w looks like:
full_length = length(rm_case.tVec)
time = rm_case.tVec
i=1
while full_length != length(rm_case.allBw[i])
    i+=1
end
b_w = rm_case.allBw[i]
lambda = rad2deg.(rm_case.allLambda[i])
plot(time, b_w, xlabel="time", ylabel="B_w", title="Wave amplitude as function of time", legend=false)
plot(lambda, b_w, xlabel="Latitude (deg)", ylabel="B_w", title="Wave amplitude as function of Latitude", legend=false)

# next, isolate just those with starting energy of 63 and 138 keV
initial_E = [rm_case.allE[i][1] for i = 1:length(rm_case.allT)]
final_E = [rm_case.allE[i][end] for i = 1:length(rm_case.allT)]
initial_PA = [rm_case.allPA[i][1] for i = 1:length(rm_case.allT)]
final_PA = [rm_case.allPA[i][end] for i = 1:length(rm_case.allT)]
E63_range = findall(x->62<x<64, initial_E)
E138_range = findall(x->137<x<139, initial_E)

plot3 = scatter(xlabel="Initial Pitch Angle (deg)", ylabel="Final Pitch Angle (deg)", title="L = 6, omega_m = 0.2", legend=:topleft)
plot3 = scatter!([initial_PA[E63_range], initial_PA[E138_range]], [final_PA[E63_range], final_PA[E138_range]], label = ["63 keV" "138 keV"], xlim=[0,30], ylim=[0,90], ms=1, msw=0)
plot3 = plot!(0:30, 0:30, label=false, ls=:dash, lc=:red)

histogram2d(initial_PA[E63_range], final_PA[E63_range], xlim=[0,30], ylim=[0,90], nbins=(30,90), clim = (1,1000), cscale=log10)

# for case #4, omega_m = 0.35, L = 6; contains both 63 and 138 keV
case = 4
rm_case = rm_array[case]
# first, obtain what B_w looks like:
full_length = length(rm_case.tVec)
time = rm_case.tVec
i=1
while full_length != length(rm_case.allBw[i])
    i+=1
end
b_w = rm_case.allBw[i]
lambda = rad2deg.(rm_case.allLambda[i])
plot(time, b_w, xlabel="time", ylabel="B_w", title="Wave amplitude as function of time", legend=false)
plot(lambda, b_w, xlabel="Latitude (deg)", ylabel="B_w", title="Wave amplitude as function of Latitude", legend=false)

# next, isolate just those with starting energy of 63 and 138 keV
initial_E = [rm_case.allE[i][1] for i = 1:length(rm_case.allT)]
final_E = [rm_case.allE[i][end] for i = 1:length(rm_case.allT)]
initial_PA = [rm_case.allPA[i][1] for i = 1:length(rm_case.allT)]
final_PA = [rm_case.allPA[i][end] for i = 1:length(rm_case.allT)]
E63_range = findall(x->62<x<64, initial_E)
E138_range = findall(x->137<x<139, initial_E)

plot4 = scatter(xlabel="Initial Pitch Angle (deg)", ylabel="Final Pitch Angle (deg)", title="L = 6, omega_m = 0.35", legend=:topleft)
plot4 = scatter!([initial_PA[E63_range], initial_PA[E138_range]], [final_PA[E63_range], final_PA[E138_range]], label = ["63 keV" "138 keV"], xlim=[0,30], ylim=[0,90], ms=1, msw=0)
plot4 = plot!(0:30, 0:30, label=false, ls=:dash, lc=:red)

histogram2d(initial_PA[E63_range], final_PA[E63_range], xlim=[0,30], ylim=[0,90], nbins=(30,90), clim = (1,1000), cscale=log10)

savefig(plot_bw, "images/loss_cone_analysis/bw.png")
savefig(plot1, "images/loss_cone_analysis/l5_m2.png")
savefig(plot2, "images/loss_cone_analysis/l5_m5.png")
savefig(plot3, "images/loss_cone_analysis/l6_m2.png")
savefig(plot4, "images/loss_cone_analysis/l6_m5.png")