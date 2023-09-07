include("agapitovHelpers.jl")
include("agapitovmodel.jl")
include("agapitovPlotHelpers.jl")

test_cases = [4.5 23   3  "LO_NITE_MODEL";
              4.5 16.5 3  "LO_DUSK_MODEL";
              4.5 8.0  3  "LO_DAWN_MODEL";
              6.5 23   3  "HI_NITE_MODEL";
              6.5 16.5 3  "HI_DUSK_MODEL";
              6.5 8.0  3  "HI_DAWN_MODEL";
              ]

wave_model_array, wave_model_coeff_array, wave_normalizer, wave_shifter_array = setup_wave_model(test_cases)

degree = 0:0.01:90
matrix = hcat(degree)

bw_plot = plot()
for case_index in eachindex(test_cases[:,1])
    wave_model_coeffs = wave_model_coeff_array[case_index]
    wave_shifter = wave_shifter_array[case_index]
    B_w(λ) = wave_normalizer * ((10 ^ abs( wave_model_coeffs[1] * (abs(λ) - wave_model_coeffs[4]) * exp(-abs(λ) * wave_model_coeffs[3] - wave_model_coeffs[2]))) - wave_shifter) * tanh(deg2rad(λ)/deg2rad(1))
    matrix = hcat(matrix, B_w.(degree))
    bw_plot = plot!(B_w, 0, 90, label = test_cases[case_index,4], xlabel="MLAT (deg)", ylabel="B_w")
end
bw_plot = plot!()
savefig(bw_plot, "bw_profiles.png")

df = DataFrame(matrix, :auto)
df = rename(df,String.(vcat("MLAT (deg)", test_cases[:,4])))

CSV.write("bw_profiles.csv", df)
