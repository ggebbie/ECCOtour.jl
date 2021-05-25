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
path_out = "/home/gebbie/julia/outputs/"
## DEFINE THE LIST OF SIGMA1 VALUES.
sig1grid = sigma1grid()

# SHOULD ALLOW USER TO CHOOSE THIS SOMEHOW,
# TSROOT REQUIRED, OTHERS NOT REQUIRED
TSroot = "state_3d_set1" # 1: θ, 2: S
RProot = "state_3d_set2" # 1:rhoanoma, 2 phihyd

splorder = 3 # spline order

################################################################
# get MITgcm / ECCOv4r4 LLC grid and depth information. Store in γ.
path_grid="../inputs/GRID_LLC90/"
γ = setupLLCgrid(path_grid)
nf = length(γ.fSize)
# get standard levels of MITgcm
fileZ = "RC"
z = read_mdsio(path_grid,fileZ)
z = vec(z)
nz = length(z)

# ECCOv4r4 uses approximation for pressure without any horizontal deviations.
# Can precompute pressure for each depth level.
ρ₀ = 1029 # from "data" text file in run directory
g  = 9.81 # from "data" 
Pa2dbar = 1/10000 # standard pressures via hydrostatic balance
pstdz = -ρ₀ .*g .* Pa2dbar .* z # 10000 to change Pa to dbar
p₀ = 1000 # dbar

# name of file inside diagspath
# Look at /poseidon ... exps/run/data.diagnostics for this info.
filelist = searchdir(diagpath[expt],TSroot) # first filter for state_3d_set1
datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"

filelist2 = searchdir(diagpath[expt],RProot) # first filter for state_3d_set1
datafilelist2  = filter(x -> occursin("data",x),filelist) # second filter for "data"

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
        
    # Read All Variables And Puts Them Into "Vars" Dictionary
    vars = read_state_3d(diagpath[expt],fileroot,γ)
    vars2 = read_state_3d(diagpath[expt],fileroot2,γ)
    merge!(vars,vars2)
    
    # solve for sigma1 on depth levels.
    @time varsσ = vars2sigma1(vars,pstdz,sig1grid,γ,splorder)

    fileprefix = pathoutexpt*"/"
    filesuffix = "_on_sigma"*fileroot[14:end]*".data"
    write_vars(varsσ,fileprefix,filesuffix)
        
end

# extra scratch work, debugging, etc.
##     xx = 168; yy = 86; ff = 5
#     figure(1)
#     clf()
#     plot([vars["S"][5,iii][168,86] for iii = 1:50],-z,"x")
#     figure(2)
#     clf()
#     plot([varsσ["S"][ff,iii][xx,yy] for iii = 1:81],sig1grid,"x")
#     tmp = varsσ["S"]

#     faststats(tmp[:,25])
#     dryval = NaN32
#     [println(MeshArrays.maximum(tmp[:,40],dryval)) 
#         [println(MeshArrays.maximum(tmp[:,iii],dryval)) for iii = 1:81]
# MeshArrays.maximum(tmp[:,40],dryval) 
#         [println(MeshArrays.maximum(tmp[:,iii],dryval)) for iii = 1:81]
#     typeof(tmp)
#      methods(maximum)
