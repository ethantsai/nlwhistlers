const Re       = 6370e3;                      # Earth radius, f64
const c        = 3e8;                         # speedo lite, f64
eta(L)         = Omegace0*L*Re/c;

B0 = 188;                     # nT, measured from THEMIS A at L of 5.2 (conjunction w/ ELFIN)
L = 5.2;
Ew = [21.358639,36.114532,25.813985,48.348956,32.235131,17.927961,23.999402]; # mV/m, measured from THEMIS E bursts

Omegace0 = (1e-9*1.6e-19*B0)/(9.11e-31);  # electron gyrofreq @ the equator
Omegape = 8; # from THEMIS A density measurements
omegam = 0.22; # from THEMIS A at L of 5.2

N = (Omegape/omegam)*((1/omegam)-1)^(-1/2)

Bw = (Ew .* N ./ 300 ) # nT

wave_amp_modifier = ( Bw ./ B0 )*eta(L)
avg_wave_amp_mod = sum(wave_amp_modifier)/(2*length(wave_amp_modifier))









using Plots, LaTeXStrings

u(lambda) = .5*(tanh(deg2rad(lambda)/deg2rad(1))+1)
u20(lambda) = 0.929745945424938^-1*tanh((deg2rad(lambda)/(deg2rad(2)))) * (exp(-(deg2rad(lambda)/(deg2rad(20)))^2)); 
u30(lambda) = 0.9597545683166291^-1*tanh((deg2rad(lambda)/(deg2rad(2)))) * (exp(-(deg2rad(lambda)/(deg2rad(30)))^2));
u40(lambda) = 0.9733653674717542^-1*tanh((deg2rad(lambda)/(deg2rad(2)))) * (exp(-(deg2rad(lambda)/(deg2rad(40)))^2)); 

plot(u, 0,90,label=L"0.5\tanh(\lambda)/1 + 1")
plot!(u20, 0,90,label=L"\delta\lambda_2 = 20\degree")
plot!(u30, 0,90,label=L"\delta\lambda_2 = 30\degree")
plot!(u40, 0,90,label=L"\delta\lambda_2 = 40\degree")
plot!(xlabel="Latitude (degrees)", ylabel="B_w modifier", title="Latitude distribution of chorus waves")
plot!(yscale=:log10, ylim = (0.1,1.))


maximum(u20.(0:.1:90)) # 0.929745945424938^-1*
maximum(u30.(0:.1:90)) # 0.9597545683166291^-1*
maximum(u40.(0:.1:90)) # 0.9733653674717542^-1*