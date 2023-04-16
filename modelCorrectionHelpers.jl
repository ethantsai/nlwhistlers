using Dierckx

# import constants
const c      = 3e8; # speedo lite in m/s
const q_e    = 1.6e-19; # charge of electron in C
const m_e    = 9.11e-31; # mass of electron in kg
const m_er   = 510.999; # rest mass of electron in eV
B_eq = 3.e-5; # equatorial magnetic field
Omegace0 = (q_e*B_eq)/(m_e);  # electron gyrofreq @ the equator in kHz

# import energy resonance functions
ω_pe(L) = L
α_lc(L) = asin(sqrt(1/(L^3 * sqrt(4- (3/L)))))
B_lat(λ) = sqrt(1 + 3*sin(λ)^2) / cos(λ)^6
A(ω_m, L, λ) = ((ω_pe(L)/ω_m) * cos(λ)^-(5/2) * sqrt(1 - (sin(α_lc(L))^2)*B_lat(λ))) / sqrt((B_lat(λ)/ω_m) - 1)
B(ω_m, L, λ) = B_lat(λ)/ω_m
γ(ω_m, L, λ) = (-B(ω_m, L, λ) + (A(ω_m, L, λ) * sqrt( B(ω_m, L, λ)^2 + A(ω_m, L, λ)^2 - 1)) ) / (A(ω_m, L, λ)^2 -1)
E(γ) = m_er * (γ-1)
γ(E) = 1 + E/m_er

# plotting to check functions
# x = 0.1:0.001:0.5
# y = 0:0.1:50
# plot(y, @. E(γ(0.1, 6, deg2rad(y))))
# plot!(y, @. E(γ(0.2, 6, deg2rad(y))))
# plot!(y, @. E(γ(0.3, 6, deg2rad(y))))
# plot!(y, @. E(γ(0.4, 6, deg2rad(y))))
# plot!(xlabel = "Latitude", yscale = :log10, yminorticks=10, ylim=(1, 10000))

# z = @. E(γ(x', 5, deg2rad(y)))
# tv = [0,1,10,50,100,250,500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000]
# contour(x,y,z, xlabel = "f/f_ce", ylabel = "λ (deg)", fill=(true,cgrad(:seaborn_icefire_gradient, scale = :log10)), colorbar_scale = :log10, levels=tv, clabels=true, clim = (1,5000))

test_cases = [7.1 8.4  3  "ELB_SA_210106T1154"; # ELB SA 01/06 11:54
              6.5 19.8 3  "ELB_ND_210108T0646"; # ELB ND 01/08 06:46
              4.8 19.0 3  "ELA_SD_210111T1750"; # ELA SD 01/11 17:50
              6   8.4  3  "ELA_NA_210112T0226"; # ELA NA 01/12 02:26
              6.5 3.8  3  "ELA_ND_200904T0112"; # ELA ND 09/04 01:12
              4.8 2.6  4  "ELB_ND_200926T0101"; # ELB ND 09/26 01:01
              5.1 20.2 3  "ELA_SD_210203T0902"; # ELA SD 02/03 09:02
              6.6 19.3 3  "ELA_SD_210203T1342"; # ELA SD 02/03 13:42
              6.2 5.8  3  "ELA_NA_210203T2046"; # ELA NA 02/03 20:46
              7.5 5.7  3  "ELA_NA_210203T2047"; # ELA NA 02/03 20:47
              6.4 10.2 4  "ELA_SD_211001T0459"; # ELA SD 10/01 05:01
              6.6 8.8  3  "ELA_SD_211001T0810"; # ELA SD 10/01 08:10
              6.1 13.1 3  "ELA_SA_211101T0424"; # ELA SA 11/01 04:24
              4.5 20.8 3  "ELA_SA_211102T2218"] # ELA SA 11/02 22:18


function extract_elfin_comparison(scenario, start, stop)
    @time @load "result_matrix/"*scenario*".jld2" rm
    id, yy, mm, dd, HH, MM = extract_idyymmddHHMM(scenario)
    sim_ratio = prec_to_trap_ratio(rm)
    sim_ratio_sm = smooth(sim_ratio[1], 6, 5)
    elfin_p2t, elfin_p2t_error = extract_elfin_p2t_ratio(
        mm*dd*yy*"_"*HH, start, stop
        )
    @info "Model E_max = $(elfin_p2t[1][last_index_before_zero(elfin_p2t[2])]) keV for scenario $scenario"
    normalizer = normalize_to_elfin(elfin_p2t[2], sim_ratio_sm)
    index = findfirst(x->x==scenario, test_cases[:,4])
    L, MLT, Kp = test_cases[index,1:3]

    start_string = Dates.format(start, DateFormat("yyyy-mm-dd HH:MM:SS"))
    stop_string = Dates.format(stop, DateFormat("HH:MM:SS"))
    ELFIN_label = "ELFIN-$id $start_string-$stop_string"


    model_bins = E_bins
    model_ratio = normalizer*sim_ratio_sm
    ELFIN_bins, ELFIN_ratio = elfin_p2t
    return model_bins, model_ratio, ELFIN_bins, ELFIN_ratio, elfin_p2t_error, sim_ratio_sm, normalizer, L, MLT, Kp, ELFIN_label, index
end

function obtain_K_E(ELFIN_bins, ELFIN_ratio, model_bins, model_ratio)
    E_max_index = last_index_before_zero(ELFIN_ratio)
    while ELFIN_ratio[E_max_index-1] / ELFIN_ratio[E_max_index] > 2
        @info "reducing E_max_index; E_max = $(ELFIN_bins[E_max_index]) -> $(ELFIN_bins[E_max_index-1])"
        E_max_index -= 1
    end
    E_max = ELFIN_bins[E_max_index]
    @info "E_max = $E_max"

    # interpolate to 4x ELFIN bins over shorter energy range
    E_new = logrange(50,E_max,64)

    spl_model = Spline1D(model_bins, model_ratio)
    spl_elfin = Spline1D(ELFIN_bins, ELFIN_ratio)

    k = spl_elfin(E_new) ./ spl_model(E_new) 

    return E_new, spl_model, spl_elfin, k
end

function obtain_E2λ(k, E_new, model_bins, ELFIN_bins, resolution, omega_m, L)
    degrees = 0:resolution:60
    E2λ = Spline1D(degrees, @. E(γ(omega_m, L, deg2rad(degrees))))
    # plot(E2λ(degrees),degrees)
    
    # use E(λ) to find for each energy, what is the corresponding latitude
    E_of_λ = [degrees[findfirst(x->x>energy, E2λ(degrees))] for energy in E_new]
    E_of_λ_model_bins = [degrees[findfirst(x->x>energy, E2λ(degrees))] for energy in model_bins]
    E_of_λ_ELFIN_bins = [degrees[findfirst(x->x>energy, E2λ(degrees))] for energy in ELFIN_bins]

    # remove K < 1
    index_of_interest = findfirst(x->x>1, k)
    k_mult = Spline1D(E_of_λ[index_of_interest:end], k[index_of_interest:end])

    return k_mult, degrees, E2λ, E_of_λ, E_of_λ_model_bins, E_of_λ_ELFIN_bins
end

function obtain_B_w_mod(test_cases, scenario, k_mult, window, order, piecewise_location)
    wave_model_array, wave_model_coeff_array, wave_normalizer, wave_shifter_array = setup_wave_model(test_cases)
    case_index = findfirst(x->x==scenario, test_cases[:,4])
    wave_model_coeffs = wave_model_coeff_array[case_index]
    wave_shifter = wave_shifter_array[case_index]
    B_w_og(λ) = wave_normalizer * ((10 ^ abs( wave_model_coeffs[1] * (abs(λ) - wave_model_coeffs[4]) * exp(-abs(λ) * wave_model_coeffs[3] - wave_model_coeffs[2]))) - wave_shifter) * tanh(deg2rad(λ)/deg2rad(1))
    B_w_mod(λ) = k_mult(λ) * B_w_og(λ)


    deg = 0:0.1:90
    hyper_tan = @. tanh(deg2rad(deg)/deg2rad(1))
    B_w_mod_sm = hyper_tan .* smooth(B_w_mod.(deg), window, order)

    # here we splice the original B_w with the smoothed correction
    B_w_mod_pw = zeros(length(deg))
    for i in eachindex(deg)
        if deg[i] < piecewise_location
            B_w_mod_pw[i] = B_w_mod(deg[i])
        elseif deg[i] >= piecewise_location
            B_w_mod_pw[i] = B_w_mod_sm[i]
        end
    end

    B_w_mod_ip = Spline1D(deg, B_w_mod_pc)
    B_w_mod_ip_rad = Spline1D(deg2rad.(deg), B_w_mod_pc)

    return B_w_og, B_w_mod, B_w_mod_ip, B_w_mod_ip_rad
end


function eom_mod!(dH,H,p::SVector{7},t::Float64)
    # These equations define the motion.
    #                  lambda in radians                     
    # z, pz, zeta, mu, lambda, phi = H
    # p[1] p[2]     p[3]     p[4]    p[5] p[6]  p[7]
    # eta, epsilon, Omegape, omegam, a,   dPhi, B_w = p

    sinλ = sin(H[5]);
    cosλ = cos(H[5]);
    g = exp(-p[5] * (cos(H[6]/(2*π*p[6]))^2))
    dg = (p[5] / (2 * π * p[6])) * sin(H[6]/(π * p[6]))
    sinζ = sin(H[3])*g;
    cosζ = cos(H[3])*g;
    
    # helper variables
    b = sqrt(1+3*sinλ^2)/(cosλ^6);
    db = (3*(27*sinλ-5*sin(3*H[5])))/(cosλ^8*(4+12*sinλ^2));
    γ = sqrt(1 + H[2]^2 + H[4]^2*b);
    K = copysign(1, H[5]) * (p[3] * (cosλ^(-5/2)))/sqrt(b/p[4] - 1);

    # B_w
    H[7] = p[7](H[5])
    psi = p[1] * p[2]    * H[7] * sqrt(b)/γ;  

    # actual integration vars
    dH1 = H[2]/γ;
    dH2 = -(0.5 * H[4]^2 * db)/γ - (H[4] * psi * cosζ) - (H[4] * psi * cosζ * dg);
    dH3 = p[1]*(K*dH1 - p[4] + b/γ) + (psi*sinζ)/(H[4]*K);
    dH4 = -(psi*cosζ)/K;
    dH5 = H[2]/(γ*cosλ*sqrt(1+3*sinλ^2));
    dH6 = p[1]*(K*dH1 - p[4]);

    dH .= SizedVector{7}([ dH1, dH2, dH3, dH4, dH5, dH6, 0 ]);
    
end