# map θ, S, p to sigma 1 surfaces.
# This is a wrapper routine to read files on poseidon.
# ggebbie, 1-Apr-2021

using Revise 
using MITgcmTools, MeshArrays, Statistics, Reemergence, JLD2

# save output to file?
output2file = true

# get MITgcm / ECCOv4r4 LLC grid and depth information. Store in γ.
path_grid="../inputs/GRID_LLC90/"
γ = setupLLCgrid(path_grid)

# get standard levels of MITgcm
fileZ = "RC"
z = read_mdsio(path_grid,fileZ)
z = vec(z)
nz = length(z)

# list of experiments on poseidon
exppath = "/poseidon/ECCOv4r4/MITgcm/exps/"
runpath,diagpath = listexperiments(exppath);

# abbreviations for each experiment for labels, etc.
shortnames = expnames()

## SELECT EXPERIMENTS TO COMPARE #################################
# manually choose from available experiments listed above.
exps = ("iter129_bulkformula","nointerannual")

# to do all experiments:
# exps = keys(shortnames)
#################################################################

nexps = length(exps) # number of experiments

## DEFINE THE LIST OF SIGMA1 VALUES.
sigma1 = 25:0.1:28.4 # sample values.

tecco= 1992+1/24:1/12:2018 # ecco years

TSroot = "state_3d_set1" # 1: θ, 2: S
RProot = "state_3d_set2" # 1:rhoanoma, 2 phihyd
UVWroot = "trsp_3d_set1" # 1: uvelmass, 2: vvelmass, 3:wvelmass

p₀ = 1000 # reference pressure = 1000 dbar

# ECCOv4r4 uses approximation for pressure without any horizontal deviations.
# Can precompute pressure for each depth level.
ρ₀ = 1029 # from "data" text file in run directory
g  = 9.81 # from "data" 

# standard pressures via hydrostatic balance
Pa2dbar = 1/10000
#pstdz = convert(Array{Float32,1},-ρ₀ .*g .* Pa2dbar .* z) # 10000 to change Pa to dbar
pstdz = -ρ₀ .*g .* Pa2dbar .* z # 10000 to change Pa to dbar
#p₀ = 1000f0 # dbar
p₀ = 1000 # dbar

# cycle through all chosen experiments
for exp in exps
    # name of file inside diagspath
    # Look at /poseidon ... exps/run/data.diagnostics for this info.
    filelist = searchdir(diagpath[exp],TSroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nt = length(datafilelist)
    
    # Improve code here to read meta file, make variable selection transparent.
    #nc = 2 # do first two (θ,S) of three properties in state_3d_set1

    # pre-allocate
    TS = Array{Float32, 2}(undef, nt, nz*2)
    UVW = Array{Float32, 2}(undef, nt, nz*3)
    #RP = Array{Float32, 2}(undef, nt, nz*2)
    
    tt = 0
    for TSname in datafilelist
        global tt += 1
        println("time index ",tecco[tt])

        # get θ, S
        @time TS = γ.read(diagpath[exp]*TSname,MeshArray(γ,Float32,nz*2))

        nstart = length(TSroot)+1;
        UVWname = UVWroot*TSname[nstart:end]
        # get velocity:
        @time UVW = γ.read(diagpath[exp]*UVWname,MeshArray(γ,Float32,nz*3))

        # get density, doesn't do any good because it is sigma0
        #Rname = RProot*TSname[nstart:end]
        #@time R = γ.read(diagpath[exp]*RPname,MeshArray(γ,Float32,nz*1))

        # sigma works in column space.
        # consider doing something similar for vertical interpolation.
        @time σ₁ = sigma(TS[:,1:nz],TS[:,nz+1:2*nz],pstdz,p₀)

        # put into variable names
        # solve for sigma1 on depth levels.
        #map_allvars(sigma1,pstdz,TS[:,1:50],TS[:,51:100],UVW[:,1:50],UVW[:,51:100],UVW[101:150])

        # save to file before analysis is overwritten the next time step.

    end
end

# save output as JLD2
if output2file
    outpath = "../outputs/"
    outfile = outpath*"faststats_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*".jld2"
    @save outfile absxbar xbar σx xmax xmin z 
end

