#  First steps: 1. go into Reemergence project directory. 2. go into julia REPL package mode with `]`. 3. `activate .` 4. Backspace to return to command mode in REPL.
# Modeled after "filter_interannual.jl"
# Solve for linear trends at all points.
using Revise
using Reemergence
using Statistics, PyPlot, Distributions, LinearAlgebra, StatsBase
using MeshArrays, MITgcmTools

workdir = pwd()
push!(LOAD_PATH, workdir)
cd(workdir)

path_grid="../inputs/GRID_LLC90/"
γ = setupLLCgrid(path_grid)
lat,lon = latlon(γ)

# get standard levels of MITgcm
z = depthlevels(γ)
nz = length(z)

# get experiments on poseidon/batou
# list of experiments on poseidon
exppath = "/poseidon/ECCOv4r4/MITgcm/exps/"
runpath,diagpath = listexperiments(exppath);

# print output here
path_out = "/home/gebbie/julia/outputs/"

# abbreviations for each experiment for labels, etc.
shortnames = expnames()

## SELECT EXPERIMENTS TO ANALYZE #################################
# manually choose from available experiments listed above.
#exps = ("iter129_bulkformula","nointerannual")

# to do all experiments:
exps = keys(shortnames)
#################################################################

nexps = length(exps) # number of experiments

# assume monthly averages, ECCOv4r4
tstart = 1992 + 1/24
tend = 2018
tecco = range(tstart,step=1/12,stop=2018)
nt = length(tecco)

# get weight matrix for finding trends
E,F = trend_matrices(tecco)
path_out = "/home/gebbie/julia/outputs/"

# pre-allocate β, linear trends
β = MeshArray(γ,Float32,nz) # some nans here

# cycle through all chosen experiments
for exp in exps

    fill!(β,0.0f0) # initialize β
    trend_theta!(β,diagpath[exp],tecco,γ,F)

    # save β for each experiment
    # make an output directory for each experiment
    pathoutexp = path_out*exp
    !isdir(pathoutexp) ? mkdir(pathoutexp) : nothing;
    Toutname = path_out*"DthetaDt_"*exp*".data"
    # save to file before overwritten next time step.
    γ.write(Toutname,β)
end
