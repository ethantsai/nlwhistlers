include("agapitovHelpers.jl")
include("agapitovmodel.jl")

@info "Loading constants..."
save_dir = "results_ducting/"
folder = "run11/"
test_cases = [5 22  3  "L5";
              6 22  3  "L6"];
omega_m_cases = [0.2, 0.35] # these are the different frequencies to test
L_array = test_cases[:,1]

const numParticles = 2*2700*2;
const startTime = 0;
const endTime = 10;
tspan = (startTime, endTime); # integration time

const ELo = 63;
const EHi = 138;
const Esteps = 2;
const PALo = 4;
const PAHi = 30;
const PAsteps = 2700;
ICrange = [ELo, EHi, Esteps, PALo, PAHi, PAsteps];

const z0 = 0; # start at eq
const λ0 = 0; # start at eq

const lossConeAngle = 4;

const Bw = 1000;  # pT
const a = 3;     # exp(-a * (cos(Φ/dΦ)^2))
const dPhi = 300; # exp(-a * (cos(Φ/dΦ)^2)) number of waves in each packet

const Re   = 6370e3;        # Earth radius, f64
const c    = 3e8;           # speedo lite, f64
const Beq  = 3.e-5;         # B field at equator (T), f64

const saveDecimation = 10000; # really only need first and last point
@info "Done."

rm_array = Vector{Resultant_Matrix}()
# loading loop™
for case_index in eachindex(L_array)
    for omega_m in omega_m_cases
        loadname = save_dir*folder*test_cases[:,end][case_index]*"_"*string(omega_m)[end]*"_$numParticles.jld2"
        @info "Loading solution from $loadname..."
        @time @load loadname sol
        label = test_cases[:,end][case_index]*"_"*string(omega_m)[end]
        push!(rm_array, sol2rm(sol, label));
        @info "Loaded $label."
    end
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
prec_range = findall(x->4<=x<=4.1, final_PA) 
E63_range = findall(x->62<x<64, initial_E)
E138_range = findall(x->137<x<139, initial_E)

plot1 = scatter(xlabel="Initial Pitch Angle (deg)", ylabel="Final Pitch Angle (deg)", title="L = 5, omega_m = 0.2", legend=:topleft)
plot1 = scatter!([initial_PA[E63_range], initial_PA[E138_range]], [final_PA[E63_range], final_PA[E138_range]], label = ["63 keV" "138 keV"], xlim=[0,30], ylim=[0,90], ms=1, msw=0)
plot1 = plot!(0:30, 0:30, label=false, ls=:dash, lc=:red)

histogram2d(initial_PA[E63_range], final_PA[E63_range], xlim=[0,30], ylim=[0,90], nbins=(30,90), clim = (1,1000), cscale=log10)

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