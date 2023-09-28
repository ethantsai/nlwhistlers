# run this in julia to get import everything necessary
# and have an interactive repl for either running 
# particle tracing code or doing data analysis
include("util.jl")
include("constants.jl") # <-- adjust simulation parameters in here
include("sim_setup.jl")
include("agapitov_setup.jl")
include("particle_tracing.jl")
include("data_processor.jl")

save_dir = "main/data/"
folder = "test/"
mkpath(save_dir*folder)

            #N   L   MLT  Kp  ω_m  model Name
test_case = [100 6.5 23.0 3.0 0.35 "bw"  "testing"]
num_particles, L, MLT, Kp, ω_m, model, scenario_name = test_case

sol = @time run_model(num_particles, [50,500,10,3,5,2], L, MLT, Kp, ω_m, saveDecimation, 1, model);
sol_rm = @time sol2rm(sol, scenario_name);

savename = save_dir*folder*test_case[:,end][1]*"_$(model)_$(string(ω_m)[3:end])_$num_particles.jld2"
@info "Saving solution to $savename"
@save savename sol_rm

