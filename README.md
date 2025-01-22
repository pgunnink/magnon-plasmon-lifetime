# magnon-plasmon-lifetime

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project which produces all the plots neccesary for ...

To (locally) reproduce this project, do the following:

1. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
2. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate 
```
which auto-activate the project and enable local path handling from DrWatson.

# Reproducing figures
1. ``cuprate_def.jl`` and ``MnF2_DOS.jl`` contain the definitions and parameters for the cuprate (quasi-2D) and Rutile-like (3D) compounds.
2. Running ``cuprate_Gamma.jl`` and ``MnF2_Gamma.jl`` generates the plasmon self-energy plots 
3. Running ``cuprate_DOS.jl`` and ``MnF2_DOS.jl`` generates the magnon 2DOS plots
