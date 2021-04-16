#  First steps: 1. go into Reemergence project directory. 2. go into julia REPL package mode with `]`. 3. `activate .` 4. Backspace to return to command mode in REPL.

# make spatial plots of trends computed from trends.jl.

using Revise
using Reemergence
using Statistics, PyPlot, Distributions
using LinearAlgebra, StatsBase
using MeshArrays, MITgcmTools

# get names of exps. 
workdir = pwd()
push!(LOAD_PATH, workdir)
cd(workdir)

path_grid="../inputs/GRID_LLC90/"
γ = setupLLCgrid(path_grid)
lat,lon = latlon(γ)
z = depthlevels(γ); nz = length(z)

# get experiments on poseidon/batou
exppath = "/poseidon/ECCOv4r4/MITgcm/exps/"
runpath,diagpath = listexperiments(exppath);
# abbreviations for each experiment for labels, etc.
shortnames = expnames()

# print output here
path_out = "/home/gebbie/julia/outputs/"

## SELECT EXPERIMENTS TO ANALYZE #################################
# manually choose from available experiments listed above.
#exps = ("iter129_bulkformula")

# to do all experiments:
exps = keys(shortnames)
##############################################################
nexps = length(exps) # number of experiments

# assume monthly averages, ECCOv4r4
tstart = 1992 + 1/24
tend = 2018
tecco = range(tstart,step=1/12,stop=2018)
nt = length(tecco)

#Use interpolation to regular grid.
# rectangular grid
longrid = -179.:2.0:179.; latgrid = -89.:2.0:89.;
f,i,j,w = prereginterp(latgrid,longrid,γ)
nx = length(longrid); ny = length(latgrid);
cmap_seismic =get_cmap("seismic")
for exp in exps
    dTname = path_out*"DthetaDt_"*exp*".data"
    @time β = γ.read(dTname,MeshArray(γ,Float32,nz))

    # read trends from file. 
    #β = γ.read(dTname)

    figure(101)
    for zz = 1:nz
        βz = 1.0f4*β[:,zz] # units K/yr -> cK/century
        βzreg = reginterp(βz,nx,ny,f,i,j,w)
        mx = maximum(βz,NaN32) # filter out NaN32
        mn = minimum(βz,NaN32) 
        extrm = max(mx,-mn)
        clf()
        lims = range(-extrm,step=2extrm/30,stop=extrm)
        contourf(longrid,latgrid,βzreg,lims,cmap=cmap_seismic)
        colorbar(label="cK/century",orientation="vertical",ticks=lims)
        contour(longrid,latgrid,βzreg,lims,colors="k")

        depthlbl = string(round(Int,-z[zz]))*" m"
        depthlbl2 = string(round(Int,-z[zz]))*"m"
        titlelbl = exp*", "*" "*depthlbl
        outfname = path_out*exp*"/DthetaDt_"*shortnames[exp]*"_"*depthlbl2*".eps"
        xlbl = "longitude "*L"[\degree E]"
        ylbl = "latitude "*L"[\degree N]"
        title(titlelbl)
        xlabel(xlbl)
        ylabel(ylbl)
        savefig(outfname)
    end
end
