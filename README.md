[![DOI](https://zenodo.org/badge/367134831.svg)](https://zenodo.org/badge/latestdoi/367134831)

# Setup/Update Instructions

## Install and set up Julia for Windows

  1. Download the 64 bit installer from the website: https://julialang.org/downloads/
  2. Once installed, edit path variables.
     -  `JULIA_HOME` variable should be set to `C:\Users\ethan\AppData\Local\Programs\Julia-X.X.X\bin` if default install location
     - `NUMBER_OF_PROCESSORS` should be set to whatever your max number of threads are e.g. `16`
  3. Install dependencies.
     - add TickTock, ConfParser, Profile, Dates, Random, StaticArrays, Distributed, OrdinaryDiffEq, JLD2, Plots

## For Mac

  1. remove any existing symlink: `rm -f /usr/local/bin/julia`
  2. make new symlink: `sudo ln -s /Applications/Julia-1.9.app/Contents/Resources/julia/bin/julia /usr/local/bin/julia`


## For Ubuntu

  1. In downloads, dl latest binaries, where X.X is gonna be version numbers:

  ```
  wget https://julialang-s3.julialang.org/bin/linux/x64/X.X/julia-X.X.X-linux-x86_64.tar.gz
  ```

  2. Extract with `tar -xvzf julia-X.X.X-linux-x86_64.tar.gz`.

  3. Move it to `/usr/lib/`

```
sudo cp -r julia-X.X.X /usr/lib/
```

  4. create a symlink in `/usr/local/bin`

```
sudo ln -s /usr/lib/julia-X.X.X/bin/julia /usr/local/bin/julia
```

  5. Ensure `bashrc` has the num_threads variable set correctly. I.e. `export JULIA_NUM_THREADS=16`
  
## Install packages

  1. Open terminal/powershell, and run `Julia`. Verify it's the correct version. You can do so by typing in `VERSION` into Julia.

  2. Install packages by hitting `]` and then entering:
```
add TickTock, ConfParser, Profile, Dates, Random, StaticArrays, Distributed, OrdinaryDiffEq, JLD2, Plots, PyPlot, InteractiveUtils, PlutoUI, LoopVectorization, BenchmarkTools, StatsPlots, DataFrames, CSV, LaTeXStrings, StatsBase, DirectConvolution
```

## Run code

  1. Modify the setup.conf to your liking. Make sure the basename is correct too.
  
  2. On any machine, you can just do `julia runEnsemble.jl` to run it.


## Profiling Code

Using Profile and PProf
  - Library: ```using Profile, PProf```
  - Syntax: `@profile [function]`
  - View: `pprof(;webport=58599)`

Make sure graphviz is installed. The arrows a -> b are a calls b and the numbers provide run time information, with second percentage being percent total run time and first percentage being how much itself takes (as compared to its subroutines).

```
ls|while read file; do mv "$file" "$(echo "$file"|sed -e 's/[() ]//g')";done
```

where stuff in square brackets is stuff you want to remove, and stuff between `//` is what you'll replace it with.
