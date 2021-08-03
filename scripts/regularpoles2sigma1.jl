# map regularpoles output to sigma 1 surfaces.
# ggebbie, 25-May-2021

using Revise 
using MITgcmTools, MeshArrays, Statistics
using Reemergence, JLD2, Dierckx, Interpolations

##########################################
# list of experiments on poseidon
# THIS IS USER INPUT
exppath = "/batou/ECCOv4r4/MITgcm/exps/"
runpath,diagpath,regpolespath = listexperiments(exppath);

## SELECT EXPERIMENTS TO COMPARE #################################
# manually choose from available experiments listed above.
expt = "iter129_bulkformula"
# print output here
path_out = regpolespath[expt]

## DEFINE THE LIST OF SIGMA1 VALUES.
sig1grid = sigma1grid()

# the state_3d monthly-average diagnostic output on regularpoles grid
varroot = ("EVEL","EVELMASS","EVELSTAR","NVEL","NVELMASS","NVELSTAR","PHIHYD","RHOAnoma","SALT","THETA")
nv = length(varroot)

# if splorder is larger than the number of points in a profile,
# then it will default to linear interpolation
#splorder = 3 # spline order
splorder = 100 # spline order

################################################################
# get MITgcm / ECCOv4r4 LLC grid and depth information. Store in γ.
path_grid="../inputs/GRID_LLC90/"
γ = setupLLCgrid(path_grid)
nf = length(γ.fSize)

# get standard levels of MITgcm
z = depthlevels(γ)
pstdz = pressurelevels(z)
p₀ = 1000 # dbar

# name of file inside diagspath
# Look at /batou ... exps/run/data.diagnostics for this info.

# first filter for state_3d_set1
regpolesroot = regpolespath[expt]

# get specific file names, one for each variable
fileroots = Dict{String,String}()

for tt = 1:312 # must be a better way
    for fldname in varroot
        push!(fileroots,fldname => regpolesroot*fldname*"/"*fldname)
    end

    #print timestamp
    year,month = timestamp_monthly_v4r4(tt)
    if month < 10
        tstamp = "_"*string(year)*"_0"*string(month)*".nc"
    else
        tstamp = "_"*string(year)*"_"*string(month)*".nc"
    end

    # add the tstamp
    ncfilenames = fileroots
    for (kk,vv) in ncfilenames
        ncfilenames[kk] = vv.*tstamp
    end

    # Read from filelist, map to sigma-1, write to file
    netcdf2sigma1(regpolesroot,regpolesroot,ncfilenames,γ,pstdz,sig1grid,splorder)

end
