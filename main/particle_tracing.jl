#= 
Functions for particle tracing are found here
eom_og! -> original equations of motion with only dphi and dlambda mods
eom_Bw! -> equations of motion with modified wave amplitude as a function of latitude
eom_Bwωm! -> eom with modified wave amplitude and wave frequency as a function of latitude
eom_wna1! -> eom w/ modded Bw with a slightly oblique wave normal angle model
eom_wna1! -> eom w/ modded Bw with a moderately oblique wave normal angle model
eom_wna3! -> eom w/ modded Bw with a very oblique wave normal angle model
=#

################################################
# generic particle tracing w/ dphi and dlambda #
################################################
#=
dphi determines wave packet sizes
dlambda determines an exponential decay
of wave power as function of latitude
This simulation was used as the basis
for the Tsai+ JGR 2022 publication,
but modified for compatibility
=#
function eom_og!(dH,H,p,t::Float64)
    #                       lambda in radians
    # H[1] H[2] H[3]  H[4], H[5],   H[6], H[7], H[8]                     
    # z,   pz,  zeta, mu,   lambda, phi,  u,    K = H
    # p[1] p[2]     p[3]     p[4]    p[5] p[6]  p[7]      p[8]      p[9]
    # eta, epsilon, Omegape, omegam, a,   dPhi, dLambda1, dLambda2, B_w_normalizer = p
    
    sinλ = sin(H[5]);
    cosλ = cos(H[5]);
    g = exp(-p[5] * (cos(H[6]/(2*π*p[6]))^2))
    dg = (p[5] / (2 * π * p[6])) * sin(H[6]/(π * p[6]))
    sinζ = sin(H[3])*g;
    cosζ = cos(H[3])*g;
    
    # helper variables
    b = sqrt(1+3*sinλ^2)/(cosλ^6);
    db = (3*(27*sinλ-5*sin(3*H[5])))/(cosλ^8*(4+12*sinλ^2));
    γ = sqrt(1 + H[2]^2 + 2*H[4]*b);
    H[8] = copysign(1, H[5]) * (p[3] * (cosλ^(-5/2)))/sqrt(b/p[4] - 1);

    # double sided wave, grows to max at dLambda1 deg, dissipates by dLambda2 deg
    H[7] = p[9]*tanh((H[5]/(deg2rad(p[7])))) * (exp(-(H[5]/(deg2rad(p[8])))^2)); 

    # new_psi = sqrt(2mu) * old_psi = q * old_psi 
    #     eta  * epsilon * B_w  * sqrt(b)/gamma
    psi = p[1] * p[2]    * H[7] * sqrt(b)/γ;
    
    # actual integration vars
    dH1 = H[2]/γ;
    dH2 = -(0.5 * H[4]^2 * db)/γ - (H[4] * psi * cosζ) - (H[4] * psi * cosζ * dg);
    dH3 = p[1]*(H[8]*dH1 - p[4] + b/γ) + (psi*sinζ)/(H[4]*H[8]);
    dH4 = -(psi*cosζ)/H[8];
    dH5 = H[2]/(γ*cosλ*sqrt(1+3*sinλ^2));
    dH6 = p[1]*(H[8]*dH1 - p[4]);
    dH .= SizedVector{8}([ dH1, dH2, dH3, dH4, dH5, dH6, 0, 0 ]);
end  
    
##################################################################
# particle tracing using modified B_w and ω_m from Agapitov+2018 #
##################################################################
#=
This uses the Bw profiles seen in Fig. 5 and 6 of Agapitov+ 2018.
It has been modified to be zero at λ=0°.
There is also a shifter (p[9]) which is usually zero (i.e. no changes).
However, if it is set in `setup_wave_model` then it will shift the
model such that it is zero at the desired latitude and higher.
This simulation was used as the basis
for the Tsai+ JGR 2023 publication
=#
function eom_Bw!(dH,H,p,t::Float64)
    #                       lambda in radians
    # H[1] H[2] H[3]  H[4], H[5],   H[6], H[7], H[8]                     
    # z,   pz,  zeta, mu,   lambda, phi,  u,    K = H
    # p[1] p[2]     p[3]     p[4]    p[5] p[6]  p[7] p[8]           p[9]
    # eta, epsilon, Omegape, omegam, a,   dPhi, B_w, B_w_normalizer B_w_shifter = p
    sinλ = sin(H[5]);
    cosλ = cos(H[5]);
    g = exp(-p[5] * (cos(H[6]/(2*π*p[6]))^2))
    dg = (p[5] / (2 * π * p[6])) * sin(H[6]/(π * p[6]))
    sinζ = sin(H[3])*g;
    cosζ = cos(H[3])*g;
    
    # helper variables
    b = sqrt(1+3*sinλ^2)/(cosλ^6);
    db = (3*(27*sinλ-5*sin(3*H[5])))/(cosλ^8*(4+12*sinλ^2));
    γ = sqrt(1 + H[2]^2 + 2*H[4]*b);
    H[8] = copysign(1, H[5]) * (p[3] * (cosλ^(-5/2)))/sqrt(b/p[4] - 1);

    # B_w
  #  normalizer * ( agapitov model                                                                                    ) - shifter * tanh
    H[7] = p[8] * ((10 ^ abs( p[7][1] * (abs(rad2deg(H[5])) - p[7][4]) * exp(-abs(rad2deg(H[5])) * p[7][3] - p[7][2]))) - p[9]) * tanh(rad2deg(H[5]/deg2rad(1)))
    
    #     eta  * epsilon * B_w  * sqrt(b)/gamma
    psi = p[1] * p[2]    * H[7] * sqrt(b)/γ;
    
    # actual integration vars
    dH1 = H[2]/γ;
    dH2 = -(0.5 * H[4]^2 * db)/γ - (H[4] * psi * cosζ) - (H[4] * psi * cosζ * dg);
    dH3 = p[1]*(H[8]*dH1 - p[4] + b/γ) + (psi*sinζ)/(H[4]*H[8]);
    dH4 = -(psi*cosζ)/H[8];
    dH5 = H[2]/(γ*cosλ*sqrt(1+3*sinλ^2));
    dH6 = p[1]*(H[8]*dH1 - p[4]);

    dH .= SizedVector{8}([ dH1, dH2, dH3, dH4, dH5, dH6, 0, 0 ]);
end

#=
this utilizes the wave frequency dependence as a
function of latitude seen in Agapitov+ 2018 Fig. 3
along with the Bw profiles seen in Fig. 5 and 6
=#
function eom_Bwωm!(dH,H,p,t::Float64)
    #                       lambda in radians
    # H[1] H[2] H[3]  H[4], H[5],   H[6], H[7], H[8]                     
    # z,   pz,  zeta, mu,   lambda, phi,  u,    K = H
    # p[1] p[2]     p[3]     p[4]    p[5] p[6]  p[7] p[8]           p[9]
    # eta, epsilon, Omegape, omegam, a,   dPhi, B_w, B_w_normalizer B_w_shifter = p
    sinλ = sin(H[5]);
    cosλ = cos(H[5]);
    g = exp(-p[5] * (cos(H[6]/(2*π*p[6]))^2))
    dg = (p[5] / (2 * π * p[6])) * sin(H[6]/(π * p[6]))
    sinζ = sin(H[3])*g;
    cosζ = cos(H[3])*g;
    
    # helper variables
    b = sqrt(1+3*sinλ^2)/(cosλ^6);
    db = (3*(27*sinλ-5*sin(3*H[5])))/(cosλ^8*(4+12*sinλ^2));
    γ = sqrt(1 + H[2]^2 + 2*H[4]*b);
    H[8] = copysign(1, H[5]) * (p[3] * (cosλ^(-5/2)))/sqrt(b/p[4] - 1);

    # B_w
  #  normalizer * ( agapitov model                                                                                    ) - shifter * tanh
    H[7] = p[8] * ((10 ^ abs( p[7][1] * (abs(rad2deg(H[5])) - p[7][4]) * exp(-abs(rad2deg(H[5])) * p[7][3] - p[7][2]))) - p[9]) * tanh(rad2deg(H[5]/deg2rad(1)))
    
    #     eta  * epsilon * B_w  * sqrt(b)/gamma
    psi = p[1] * p[2]    * H[7] * sqrt(b)/γ;
    

    # actual integration vars
    dH1 = H[2]/γ;
    
    # incorporate omega_m(lambda) from Agapitov+, 2018
    # deg2rad(16.8) = 0.2932153143350474
    if H[5] < 0.2932153143350474
        wave_num = H[8]*dH1 - (0.41 - 0.0125*rad2deg(H[5])) 
    else
        wave_num = H[8]*dH1 - (0.2)
    end

    dH2 = -(0.5 * H[4]^2 * db)/γ - (H[4] * psi * cosζ) - (H[4] * psi * cosζ * dg);
    dH3 = p[1]*(wave_num + b/γ) + (psi*sinζ)/(H[4]*H[8]);
    dH4 = -(psi*cosζ)/H[8];
    dH5 = H[2]/(γ*cosλ*sqrt(1+3*sinλ^2));
    dH6 = p[1]*(wave_num);

    dH .= SizedVector{8}([ dH1, dH2, dH3, dH4, dH5, dH6, 0, 0 ]);
end

#########################################################
# particle tracing including three different WNA models #
#########################################################

# WNA1: θ_g * (λ/15°)/((1+λ)/15°)
function eom_wna1!(dH,H,p,t::Float64)
    #                       lambda in radians
    # H[1] H[2] H[3]  H[4], H[5],   H[6], H[7], H[8]                     
    # z,   pz,  zeta, mu,   lambda, phi,  u,    K = H
    # p[1] p[2]     p[3]     p[4]    p[5] p[6]  p[7] p[8]           p[9]
    # eta, epsilon, Omegape, omegam, a,   dPhi, B_w, B_w_normalizer B_w_shifter = p
    sinλ = sin(H[5]);
    cosλ = cos(H[5]);
    g = exp(-p[5] * (cos(H[6]/(2*π*p[6]))^2))
    dg = (p[5] / (2 * π * p[6])) * sin(H[6]/(π * p[6]))
    sinζ = sin(H[3])*g;
    cosζ = cos(H[3])*g;
    
    # helper variables
    b = sqrt(1+3*sinλ^2)/(cosλ^6);
    db = (3*(27*sinλ-5*sin(3*H[5])))/(cosλ^8*(4+12*sinλ^2));
    γ = sqrt(1 + H[2]^2 + 2*H[4]*b);
    
    # model 1: slightly oblique
  # θ_g = acos(2*p[4]/b)                           deg2rad(15) = 0.2617993877991494;
    wna = acos(2*p[4]/b) * (H[5]/0.2617993877991494) / (1 + H[5]/0.2617993877991494);
    H[8] = copysign(1, H[5]) * (p[3] * (cosλ^(-5/2)))/sqrt((cos(wna)*b)/p[4] - 1);

    # B_w
   # normalizer * ( agapitov model                                                                                    ) - shifter * tanh
    H[7] = p[8] * ((10 ^ abs( p[7][1] * (abs(rad2deg(H[5])) - p[7][4]) * exp(-abs(rad2deg(H[5])) * p[7][3] - p[7][2]))) - p[9])   * tanh(rad2deg(H[5]/deg2rad(1)))

    #     eta  * epsilon * B_w  * sqrt(b)/gamma
    psi = p[1] * p[2]    * H[7] * sqrt(b)/γ;
    
    # actual integration vars
    dH1 = H[2]/γ;
    dH2 = -(0.5 * H[4]^2 * db)/γ - (H[4] * psi * cosζ) - (H[4] * psi * cosζ * dg);
    dH3 = p[1]*(H[8]*dH1 - p[4] + b/γ) + (psi*sinζ)/(H[4]*H[8]);
    dH4 = -(psi*cosζ)/H[8];
    dH5 = H[2]/(γ*cosλ*sqrt(1+3*sinλ^2));
    dH6 = p[1]*(H[8]*dH1 - p[4]);

    dH .= SizedVector{8}([ dH1, dH2, dH3, dH4, dH5, dH6, 0, 0 ]);
end

# WNA2: θ_r * (λ/15°)/((1+λ)/15°)
function eom_wna1!(dH,H,p,t::Float64)
    #                       lambda in radians
    # H[1] H[2] H[3]  H[4], H[5],   H[6], H[7], H[8]                     
    # z,   pz,  zeta, mu,   lambda, phi,  u,    K = H
    # p[1] p[2]     p[3]     p[4]    p[5] p[6]  p[7] p[8]           p[9]
    # eta, epsilon, Omegape, omegam, a,   dPhi, B_w, B_w_normalizer B_w_shifter = p
    sinλ = sin(H[5]);
    cosλ = cos(H[5]);
    g = exp(-p[5] * (cos(H[6]/(2*π*p[6]))^2))
    dg = (p[5] / (2 * π * p[6])) * sin(H[6]/(π * p[6]))
    sinζ = sin(H[3])*g;
    cosζ = cos(H[3])*g;
    
    # helper variables
    b = sqrt(1+3*sinλ^2)/(cosλ^6);
    db = (3*(27*sinλ-5*sin(3*H[5])))/(cosλ^8*(4+12*sinλ^2));
    γ = sqrt(1 + H[2]^2 + 2*H[4]*b);

    # model 2: moderately oblique
  # θ_r = acos(p[4]/b)                           deg2rad(15) = 0.2617993877991494;
    wna = acos(p[4]/b) * (H[5]/0.2617993877991494) / (1 + H[5]/0.2617993877991494);
    H[8] = copysign(1, H[5]) * (p[3] * (cosλ^(-5/2)))/sqrt((cos(wna)*b)/p[4] - 1);

    # B_w
   # normalizer * ( agapitov model                                                                                    ) - shifter * tanh
    H[7] = p[8] * ((10 ^ abs( p[7][1] * (abs(rad2deg(H[5])) - p[7][4]) * exp(-abs(rad2deg(H[5])) * p[7][3] - p[7][2]))) - p[9])   * tanh(rad2deg(H[5]/deg2rad(1)))

    #     eta  * epsilon * B_w  * sqrt(b)/gamma
    psi = p[1] * p[2]    * H[7] * sqrt(b)/γ;
    
    # actual integration vars
    dH1 = H[2]/γ;
    dH2 = -(0.5 * H[4]^2 * db)/γ - (H[4] * psi * cosζ) - (H[4] * psi * cosζ * dg);
    dH3 = p[1]*(H[8]*dH1 - p[4] + b/γ) + (psi*sinζ)/(H[4]*H[8]);
    dH4 = -(psi*cosζ)/H[8];
    dH5 = H[2]/(γ*cosλ*sqrt(1+3*sinλ^2));
    dH6 = p[1]*(H[8]*dH1 - p[4]);

    dH .= SizedVector{8}([ dH1, dH2, dH3, dH4, dH5, dH6, 0, 0 ]);
end

# WNA3: θ_r - 2°
function eom_wna3!(dH,H,p,t::Float64)
    #                       lambda in radians
    # H[1] H[2] H[3]  H[4], H[5],   H[6], H[7], H[8]                     
    # z,   pz,  zeta, mu,   lambda, phi,  u,    K = H
    # p[1] p[2]     p[3]     p[4]    p[5] p[6]  p[7] p[8]           p[9]
    # eta, epsilon, Omegape, omegam, a,   dPhi, B_w, B_w_normalizer B_w_shifter = p
    sinλ = sin(H[5]);
    cosλ = cos(H[5]);
    g = exp(-p[5] * (cos(H[6]/(2*π*p[6]))^2))
    dg = (p[5] / (2 * π * p[6])) * sin(H[6]/(π * p[6]))
    sinζ = sin(H[3])*g;
    cosζ = cos(H[3])*g;
    
    # helper variables
    b = sqrt(1+3*sinλ^2)/(cosλ^6);
    db = (3*(27*sinλ-5*sin(3*H[5])))/(cosλ^8*(4+12*sinλ^2));
    γ = sqrt(1 + H[2]^2 + 2*H[4]*b);
    
    # model 3: very oblique
  # θ_r = acos(p[4]/b);
          # deg2rad(2) = 0.03490658503988659;
    wna = acos(p[4]/b) - 0.03490658503988659;
    H[8] = copysign(1, H[5]) * (p[3] * (cosλ^(-5/2)))/sqrt((cos(wna)*b)/p[4] - 1);

    # B_w
   # normalizer * ( agapitov model                                                                                    ) - shifter * tanh
    H[7] = p[8] * ((10 ^ abs( p[7][1] * (abs(rad2deg(H[5])) - p[7][4]) * exp(-abs(rad2deg(H[5])) * p[7][3] - p[7][2]))) - p[9])   * tanh(rad2deg(H[5]/deg2rad(1)))

    #     eta  * epsilon * B_w  * sqrt(b)/gamma
    psi = p[1] * p[2]    * H[7] * sqrt(b)/γ;
    
    # actual integration vars
    dH1 = H[2]/γ;
    dH2 = -(0.5 * H[4]^2 * db)/γ - (H[4] * psi * cosζ) - (H[4] * psi * cosζ * dg);
    dH3 = p[1]*(H[8]*dH1 - p[4] + b/γ) + (psi*sinζ)/(H[4]*H[8]);
    dH4 = -(psi*cosζ)/H[8];
    dH5 = H[2]/(γ*cosλ*sqrt(1+3*sinλ^2));
    dH6 = p[1]*(H[8]*dH1 - p[4]);

    dH .= SizedVector{8}([ dH1, dH2, dH3, dH4, dH5, dH6, 0, 0 ]);
end

function palostcondition(H,t,integrator)
    # condition: if particle enters loss cone
    b = sqrt(1+3*sin(H[5])^2)/(cos(H[5])^6);
    γ = sqrt(1 + H[2]^2 + 2*H[4]*b);
    return (rad2deg(asin(sqrt( (2*H[4])/(γ^2 -1) )))) < (lossConeAngle/2)
end

function ixlostcondition(H,t,integrator)
    # condition: if I_x approaches 0
    #      2*mu  * b 
    return 2*H[4]*sqrt(1+3*sin(H[5])^2)/(cos(H[5])^6) < 3e-5
end

function eqlostcondition(H,t,integrator)
    # condition: if particle crosses eq in negative direction
    return H[1]<=0 && H[2]<0
end

affect!(integrator) = terminate!(integrator); # terminate if condition reached
cb1 = DiscreteCallback(palostcondition,affect!);
cb2 = DiscreteCallback(ixlostcondition,affect!);
cb3 = DiscreteCallback(eqlostcondition,affect!);