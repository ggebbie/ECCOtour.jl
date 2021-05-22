# map θ, S, p to sigma 1 surfaces.
# This is a wrapper routine to read files on poseidon.
# ggebbie, 1-Apr-2021

using Revise 
using MITgcmTools, MeshArrays, Statistics
using Reemergence, JLD2, Dierckx, Interpolations

##########################################
# list of experiments on poseidon
# THIS IS USER INPUT
exppath = "/poseidon/ECCOv4r4/MITgcm/exps/"
runpath,diagpath = listexperiments(exppath);

## SELECT EXPERIMENTS TO COMPARE #################################
# manually choose from available experiments listed above.
# EVENTUALLY THIS IS A FUNCTION, DOES ONE EXPERIMENT. LOOP OF EXPS IN
# DRIVER
expt = "iter129_bulkformula"
# print output here
path_out = "/home/gebbie/julia/outputs/"
## DEFINE THE LIST OF SIGMA1 VALUES.
sig1grid = sigma1grid()

# SHOULD ALLOW USER TO CHOOSE THIS SOMEHOW,
# TSROOT REQUIRED, OTHERS NOT REQUIRED
TSroot = "state_3d_set1" # 1: θ, 2: S
RProot = "state_3d_set2" # 1:rhoanoma, 2 phihyd
UVWroot = "trsp_3d_set1" # 1: uvelmass, 2: vvelmass, 3:wvelmass

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

# standard pressures via hydrostatic balance
Pa2dbar = 1/10000
pstdz = -ρ₀ .*g .* Pa2dbar .* z # 10000 to change Pa to dbar
p₀ = 1000 # dbar

# name of file inside diagspath
# Look at /poseidon ... exps/run/data.diagnostics for this info.
filelist = searchdir(diagpath[expt],TSroot) # first filter for state_3d_set1
datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    # make an output directory for each expteriment
pathoutexpt = path_out*expt
!isdir(pathoutexpt) ? mkdir(pathoutexpt) : nothing;
nt = length(datafilelist)
    
# Improve code here to read meta file, make variable selection transparent.
#nc = 2 # do first two (θ,S) of three properties in state_3d_set1
# READ STANDARD ROOT LIST, LOOKING FOR USER SUPPLIED VARIABLE STRINGS, READ WHAT IS NECESSARY

# WORKS FOR REGULARPOLES GRID TOO?

global tt = 0
tt =1
TSname = datafilelist[1]
for TSname in datafilelist
    tt += 1

    #print timestamp
    year,month = timestamp_monthly_v4r4(tt)
    # get θ, S
    @time TS = γ.read(diagpath[expt]*TSname,MeshArray(γ,Float32,nz*2))

    vars = Dict("θ" => TS[:,1:nz])
    push!(vars,"S" => TS[:,nz+1:2nz])

    #nstart = length(TSroot)+1;
    #    UVWname = UVWroot*TSname[nstart:end]
        # get velocity:
    #    @time UVW = γ.read(diagpath[expt]*UVWname,MeshArray(γ,Float32,nz*3))

        # put into variable names
        # solve for sigma1 on depth levels.
    @time varsσ = all2sigma1(vars,pstdz,sig1grid,γ,splorder)

    xx = 168; yy = 86; ff = 5
    figure(1)
    clf()
    plot([vars["S"][5,iii][168,86] for iii = 1:50],-z,"x")
    figure(2)
    clf()
    plot([varsσ["S"][ff,iii][xx,yy] for iii = 1:81],sig1grid,"x")
    tmp = varsσ["S"]

    faststats(tmp[:,50])
    dryval = NaN32
    [println(MeshArrays.maximum(tmp[:,40],dryval)) 
        [println(MeshArrays.maximum(tmp[:,iii],dryval)) for iii = 1:81]
MeshArrays.maximum(tmp[:,40],dryval) 
        [println(MeshArrays.maximum(tmp[:,iii],dryval)) for iii = 1:81]
    typeof(tmp)
     methods(maximum)
    # REPLACE WITH FUNCTION WRITE VARS_ON_SIGMA
        Toutname = pathoutexpt*"/theta_on_sigma1"*TSname[14:end]
        γ.write(Toutname,θσ)

        Soutname = pathoutexpt*"/S_on_sigma1"*TSname[14:end]
        γ.write(Soutname,Sσ)

        Poutname = pathoutexpt*"/p_on_sigma1"*TSname[14:end]
        γ.write(Poutname,pσ)

        Uoutname = pathoutexpt*"/u_on_sigma1"*TSname[14:end]
        γ.write(Uoutname,uσ)

        Voutname = pathoutexpt*"/v_on_sigma1"*TSname[14:end]
        γ.write(Voutname,vσ)

        Woutname = pathoutexpt*"/w_on_sigma1"*TSname[14:end]
        γ.write(Woutname,wσ)
        
    end
end
