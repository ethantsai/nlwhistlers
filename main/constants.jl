# importing packages and constants
@info "Compiling packages for simulation..."
#meta
using TickTock
using ConfParser
using Logging
using Profile
#sim
using Dates
using Random
using OrdinaryDiffEq
using StaticArrays
#process
using JLD2
using Plots
using Plots.PlotMeasures
using StatsBase
using DirectConvolution
@info "Packages compiled."

#######################
## Constants n stuff ##
#######################
@info "Loading constants..."
save_dir = "results_ducting/"
folder = "run28/"
mkpath(save_dir*folder)

test_cases = [6.5 23   3  "HI_NITE_MODEL_small"]

omega_m_cases = [0.15, 0.45] # these are the different frequencies to test
L_array = test_cases[:,1]

const numParticles = 13;
const startTime = 0;
const endTime = 15;
tspan = (startTime, endTime); # integration time

const ELo = 52;
const EHi = 1000;
const Esteps = 32; # double ELFIN E bins
const PALo = 3;
const PAHi = 15;
const PAsteps = 1300; # only used for flat particle distribution
const factor = 40; #only used for skewed particle distribution
# num particles in highest energy bin = factor * num particles in lowest energy bin
ICrange = [ELo, EHi, Esteps, PALo, PAHi, PAsteps];
E_bins = logrange(ELo,EHi, Int64(Esteps))

const z0 = 0; # start at eq
const λ0 = 0; # start at eq

const lossConeAngle = 3;

const Bw = 300;  # pT
const a = 7;     # exp(-a * (cos(Φ/dΦ)^2))
const dPhi = 3; # exp(-a * (cos(Φ/dΦ)^2)) number of waves in each packet

const Re   = 6370e3;        # Earth radius, f64
const c    = 3e8;           # speedo lite, f64
const Beq  = 3.e-5;         # B field at equator (T), f64

saveDecimation = 40000; # really only need first and last point
@info "Done."


################
## Plot Setup ##;
################
function hexcolor(r::UInt8, g::UInt8, b::UInt8)
    return RGB(Int(r)/255,Int(g)/255,Int(b)/255)
end
bipride_pink = RGB(234/255, 2/255, 112/255);
bipride_orange = RGB(255/255, 79/255, 0/255);
bipride_lavender = RGB(155/255, 79/255, 150/255);
bipride_blue = RGB(0/255, 56/255, 168/255);
c0 = hexcolor(0xf5,0x97,0x00); # golden
c1 = hexcolor(0xff,0x4f,0x00); # orange
c2 = hexcolor(0xf8,0x00,0x4d); # reddishpink
c3 = hexcolor(0xd3,0x00,0x7d); # magenta
c4 = hexcolor(0x8f,0x19,0x9f); # purple
c5 = hexcolor(0x00,0x38,0xa8); # blue
c6 = hexcolor(0x00,0xA6,0xFF); # light blue
c7 = hexcolor(0x00,0xAB,0x12); # green
c8 = hexcolor(0x00,0xF7,0x1A); # light green

# colorblind friendlier colors
cb1 = hexcolor(0x00,0x00,0x00); # black
cb2 = hexcolor(0xE6,0x9F,0x00); # orange
cb3 = hexcolor(0x56,0xB4,0xE9); # light blue
cb4 = hexcolor(0x00,0x9E,0x73); # green
cb5 = hexcolor(0xF0,0xE4,0x42); # yellow
cb6 = hexcolor(0x00,0x72,0xB2); # blue
cb7 = hexcolor(0xD5,0x5E,0x00); # red
cb8 = hexcolor(0xCC,0x79,0xA7); # purple
