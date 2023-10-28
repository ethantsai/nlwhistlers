[![DOI](https://zenodo.org/badge/367134831.svg)](https://zenodo.org/badge/latestdoi/367134831)


# Code Organization

- Tsai, E., et al., Relativistic electron precipitation driven by non-linear resonance with whistler-mode waves. Journal of Geophysical Research: Space Physics, 127 (2022), doi: https://doi.org/10.1029/2022JA030338
    - Code for this paper is found in the `jgr_2022_work` directory
    - There was only a cursory check for filepath correctness when code got reorganized and it has not been validated to run without errors, especially since it relies on old versions of Julia/old packages
    - However, resulting images (some figures from paper and more) along with the resulting simulation data (stored in HDF5/JLD2 format) are saved in here too
    - Code here was run by setting up parameters in the `setup.conf` file then running `julia runEnsemble.jl`
- Tsai, E., et al., Investigating whistler-mode wave intensity along field lines using electron precipitation measurements, Journal of Geophysical Research Space Physics, 128 (2023), doi: https://doi.org/10.1029/2023JA031578
    - Code for this paper is found in the `jgr_2023_1_work` directory
    - There was only a cursory check for filepath correctness when code got reorganized and it has not been validated to run without errors
    - The saved data is too large, but some of the plots beyond what is published are saved here too
    - A lot of the fundamental code for this is used and incorporated in the `main` source code
    - Code here was run by setting up parameters in the `agapitovHelpers.jl` file, then running `julia run_ducting_analysis.jl`
- `external_data` hosts:
    - Processed ELFIN measurements using SPEDAS V4.1. This code is not publicly available, but available upon request.
    - Statistical output from ELFIN observations in `stats_csvs`.
    - Diffusion code results from the UCLA Full Diffusion Code (Ma, Q. 2018, https://doi.org/10.1002/2017JA025114)
- Current active work is being done in the `main` folder
    - Data here is output into a `data` directory in the root of the folder
    - Because simulation results are so large, they are not included in the repo. Please find and download them instead from here: https://ucla.box.com/s/stj0ilwbf3xk5nhk8qqp4ync3x71w1um. You can put this data inside of the results folder and run analysis plots on it.

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

## Packages

  1. Once Julia is set up, you should be able to enter the package manager via `]`.
  2. Type `activate .` to activate the manifest.
  3. Type `instantiate` to install all missing packages. When complete, you are ready to use the Julia code as long as you're in the proper virtual environment

To always activate the virtual environment of whatever project directory you're working in (i.e. if you'd like to run Julia directly from command line), you should add `export JULIA_PROJECT=@.` in the `bashrc` file.