include("plotHelpers.jl")
Egrid, PAgrid = logrange(10,1000,21), 6:4:90
@load "210429_data_storage.jld2" themis_lolat themis_hilat
plot(Egrid, themis_lolat.precipitating_fluxes_mean,
        yerror=(themis_lolat.precipitating_fluxes_plus, themis_lolat.precipitating_fluxes_minus),
        ylim =(1e2,1e9), xlim=(50,800), yscale=:log10)
plot!(Egrid, themis_hilat.precipitating_fluxes_mean,
yerror=(themis_hilat.precipitating_fluxes_plus, themis_hilat.precipitating_fluxes_minus),
ylim =(1e2,1e9), xlim=(50,800), yscale=:log10)
