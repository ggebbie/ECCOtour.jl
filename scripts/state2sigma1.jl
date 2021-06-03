# map θ, S, p to sigma 1 surfaces.
# This is a wrapper routine to read files on poseidon.
# ggebbie, 1-Apr-2021

using Revise 
using MITgcmTools, MeshArrays, Statistics
using Reemergence, JLD2, Dierckx, Interpolations

##########################################
# list of experiments on poseidon
# THIS IS USER INPUT
exppath = "/batou/ECCOv4r4/MITgcm/exps/"
runpath,diagpath = listexperiments(exppath);

## SELECT EXPERIMENTS TO COMPARE #################################
# manually choose from available experiments listed above.
expt = "iter129_bulkformula"
# print output here
path_out = "/batou/ECCOv4r4/MITgcm/exps/"*expt*"/run/sigma1/"

## DEFINE THE LIST OF SIGMA1 VALUES.
sig1grid = sigma1grid()

# the state_3d monthly-average diagnostic output
TSroot = "state_3d_set1" # 1: θ, 2: S
RProot = "state_3d_set2" # 1:rhoanoma, 2 phihyd

splorder = 3 # spline order

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
filelist = searchdir(diagpath[expt],TSroot)
 # second filter for "data"
datafilelist  = filter(x -> occursin("data",x),filelist)

filelist2 = searchdir(diagpath[expt],RProot) 
datafilelist2  = filter(x -> occursin("data",x),filelist)

# make an output directory for each expteriment
pathoutexpt = path_out*expt
!isdir(pathoutexpt) ? mkdir(pathoutexpt) : nothing;
nt = length(datafilelist)
    
global tt = 0
for datafile in datafilelist
    tt += 1

    #print timestamp
    year,month = timestamp_monthly_v4r4(tt)

    # eliminate suffix
    fileroot = rstrip(datafile,['.','d','a','t','a'])
    fileroot2 = RProot*fileroot[14:end] # a better way?
    fileroots = (fileroot,fileroot2)
    
    # Read from filelist, map to sigma-1, write to file
    files2sigma1(diagpath[expt],path_out,fileroots,γ,pstdz,sig1grid,splorder)
end
