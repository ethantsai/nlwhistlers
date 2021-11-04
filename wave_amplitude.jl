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









using Plots

u(lambda) = .5*(tanh(deg2rad(lambda)/deg2rad(1))+1)
u10(lambda) = tanh((deg2rad(lambda)/(deg2rad(2)))) * (exp(-(deg2rad(lambda)/(deg2rad(10)))^2)); 
u15(lambda) = tanh((deg2rad(lambda)/(deg2rad(2)))) * (exp(-(deg2rad(lambda)/(deg2rad(15)))^2));
u40(lambda) = tanh((deg2rad(lambda)/(deg2rad(2)))) * (exp(-(deg2rad(lambda)/(deg2rad(40)))^2)); 

plot(u, 0,90)
plot!(u10, 0,90)
plot!(u15, 0,90)
plot!(u40, 0,90)
plot!(yscale=:log10, ylim = (1e-4,1.1))
