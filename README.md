# Setup/Update Instructions

## Install and set up Julia for Windows

  1. Download the 64 bit installer from the website: https://julialang.org/downloads/
  2. Once installed, edit path variables.
     -  `JULIA_HOME` variable should be set to `C:\Users\ethan\AppData\Local\Programs\Julia-X.X.X\bin` if default install location
     - `NUMBER_OF_PROCESSORS` should be set to whatever your max number of threads are e.g. `16`
  3. Install packages:  
     - Open powershell, and run `Julia`. Verify it's the correct version.
     - Install packages by hitting `]` and then entering:
```
add TickTock, ConfParser, Distributed, DifferentialEquations, JLD2, Plots, Dates
```

## Run code

  1. Modify the setup.conf to your liking. Make sure the basename is correct too.
  2. On any machine, you can just do `julia runEnsemble.jl` to run it. However, there's a problem with windows machines and running in powershell that I haven't determined the cause yet. If you want to run in windows, you have to just open the repl and run stuff line by line in there.
  
