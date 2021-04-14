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

# get correct files
Troot = "state_3d_set1" # 1: θ, 2: S
path_out = "/home/gebbie/julia/outputs/"

# pre-allocate β, linear trends
β = MeshArray(γ,Float32,nz) # some nans here

# cycle through all chosen experiments
for exp in exps
    # name of file inside diagspath
    # Look at /poseidon ... exps/run/data.diagnostics for this info.
    filelist = searchdir(diagpath[exp],Troot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    # make an output directory for each experiment
    pathoutexp = path_out*exp
    !isdir(pathoutexp) ? mkdir(pathoutexp) : nothing;
    
    nt2 = length(datafilelist)
    if nt2 != nt
        error("incompatible t and data files")
    end

    # initialize β
    fill!(β,0.0f0)

    global tt = 0
    for Tname in datafilelist
        tt += 1
        println("year ",floor(tecco[tt])," month ",((tt-1)%12)+1)

        # read θ for timestep
        @time θ = γ.read(diagpath[exp]*Tname,MeshArray(γ,Float32,nz))

        # multiply by the correct part of F matrix
        # be sure to handle all points at once
        # add to existing solution for β
        # equal to matrix multiplication w/ columns (θ) times weighting F
        β += F[2,tt] * θ
    end

    # save β for each experiment
    Toutname = path_out*"DthetaDt_"*exp*".data"
    # save to file before overwritten next time step.
    γ.write(Toutname,β)
end

# make spatial plots. Use interpolation to regular grid.
# rectangular grid
longrid = -179.:2.0:179.; latgrid = -89.:2.0:89.;
f,i,j,w = prereginterp(latgrid,longrid,γ)
nx = length(longrid); ny = length(latgrid);

figure(101)
for zz = 1:nz
    βz = β[:,zz]
    βzreg = reginterp(βz,nx,ny,f,i,j,w)
    clf()
    contourf(longrid,latgrid,βzreg)

    # save it
    
end


