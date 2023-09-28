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
const numParticles = 1000000;
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
const ICrange = [ELo, EHi, Esteps, PALo, PAHi, PAsteps];
const E_bins = logrange(ELo,EHi, Int64(Esteps))

const z0 = 0; # start at eq
const λ0 = 0; # start at eq

const lossConeAngle = 3;

const Bw = 300;  # pT
const a = 7;     # exp(-a * (cos(Φ/dΦ)^2))
const dPhi = 3; # exp(-a * (cos(Φ/dΦ)^2)) number of waves in each packet

const Re   = 6370e3;        # Earth radius, f64
const c    = 3e8;           # speedo lite, f64
const Beq  = 3.e-5;         # B field at equator (T), f64

const saveDecimation = 40000; # really only need first and last point
@info "Done."


################
## Plot Setup ##;
################
function hexcolor(r::UInt8, g::UInt8, b::UInt8)
    return RGB(Int(r)/255,Int(g)/255,Int(b)/255)
end

# colorblind friendlier colors
c1 = hexcolor(0x00,0x00,0x00); # black
c2 = hexcolor(0xE6,0x9F,0x00); # orange
c3 = hexcolor(0x56,0xB4,0xE9); # light blue
c4 = hexcolor(0x00,0x9E,0x73); # green
c5 = hexcolor(0xF0,0xE4,0x42); # yellow
c6 = hexcolor(0x00,0x72,0xB2); # blue
c7 = hexcolor(0xD5,0x5E,0x00); # red
c8 = hexcolor(0xCC,0x79,0xA7); # purple
