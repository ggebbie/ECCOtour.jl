 ggebbie, 16-Apr-2021

using Revise 
using MITgcmTools, MeshArrays, Statistics, Reemergence

# get MITgcm / ECCOv4r4 LLC grid and depth information. Store in γ.
path_grid="../inputs/GRID_LLC90/"
γ = setupLLCgrid(path_grid)

# list of experiments on poseidon
exppath = "/poseidon/ECCOv4r4/MITgcm/exps/"
runpath,diagpath = listexperiments(exppath);

# print output here
path_out = "/home/gebbie/julia/outputs/"

# abbreviations for each experiment for labels, etc.
shortnames = expnames()

## SELECT EXPERIMENTS TO ANALYZE ##############
# manually choose from available experiments listed above.
exps = ("iter129_bulkformula","nointerannual")

# to do all experiments:
# exps = keys(shortnames)
###############################################

nexps = length(exps) # number of experiments

## DEFINE THE LIST OF SIGMA1 VALUES.
sig1grida = 24:0.05:31
sig1gridb = 31.02:0.02:33
sig1grid = vcat(sig1grida,sig1gridb)
sig1grid = sig1grid[1:3:end]

tecco= 1992+1/24:1/12:2018 # ecco years

TSroot = "state_3d_set1" # 1: θ, 2: S

# cycle through all chosen experiments
for exp in exps
    # name of file inside diagspath
    # Look at /poseidon ... exps/run/data.diagnostics for this info.
    filelist = searchdir(diagpath[exp],TSroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    # make an output directory for each experiment
    pathoutexp = path_out*exp
    !isdir(pathoutexp) ? mkdir(pathoutexp) : nothing;
    
    nt = length(datafilelist)
    
    # Improve code here to read meta file, make variable selection transparent.
    #nc = 2 # do first two (θ,S) of three properties in state_3d_set1

    # pre-allocate
    TS = Array{Float32, 2}(undef, nt, nz*2)
    UVW = Array{Float32, 2}(undef, nt, nz*3)
    #RP = Array{Float32, 2}(undef, nt, nz*2)
    
    global tt = 0
    for TSname in datafilelist
        tt += 1
        println("time index ",tecco[tt])

        # get θ, S
        @time TS = γ.read(diagpath[exp]*TSname,MeshArray(γ,Float32,nz*2))

        nstart = length(TSroot)+1;
        UVWname = UVWroot*TSname[nstart:end]
        # get velocity:
        @time UVW = γ.read(diagpath[exp]*UVWname,MeshArray(γ,Float32,nz*3))

        # Not needed in this program
        # sigma works column-by-column
        # consider doing something similar for vertical interpolation.
        #@time σ₁=sigma(TS[:,1:nz],TS[:,nz+1:2*nz],pstdz,p₀)

        # put into variable names
        # solve for sigma1 on depth levels.
        @time θσ,Sσ,pσ,uσ,vσ,wσ=all2sigma1(TS[:,1:nz],TS[:,nz+1:2*nz],pstdz,UVW[:,1:nz],UVW[:,nz+1:2*nz],UVW[:,2*nz+1:3*nz],sig1grid,γ,splorder)

        Toutname = pathoutexp*"/theta_on_sigma1"*TSname[14:end]
        γ.write(Toutname,θσ)

        Soutname = pathoutexp*"/S_on_sigma1"*TSname[14:end]
        γ.write(Soutname,Sσ)

        Poutname = pathoutexp*"/p_on_sigma1"*TSname[14:end]
        γ.write(Poutname,pσ)

        Uoutname = pathoutexp*"/u_on_sigma1"*TSname[14:end]
        γ.write(Uoutname,uσ)

        Voutname = pathoutexp*"/v_on_sigma1"*TSname[14:end]
        γ.write(Voutname,vσ)

        Woutname = pathoutexp*"/w_on_sigma1"*TSname[14:end]
        γ.write(Woutname,wσ)
        
    end
end
