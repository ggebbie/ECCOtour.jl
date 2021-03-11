# interactive analysis of divergence of runs, originally 10-Feb-2021
# recreate MATLAB difference plots
# use files directly on poseidon through cifs mount

using Revise # not sure if Revise is necessary
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
exppath = "/poseidon/ecco/ECCOv4r4/MITgcm/exps/"
runpath,diagpath = listexperiments(exppath);

# abbreviations for each experiment for labels, etc.
shortnames = expnames()

## SELECT EXPERIMENTS TO COMPARE #################################
# baseline experiment
expbase = "iter129_bulkformula"
# comparison experiment(s)
#expcompare = "iter129_fluxforced"
#expcompare = "iter0_bulkformula"
expcompare = "nointerannual"
#################################################################

# Improve code here to compare to multiple experiments.
#nexp = length(expcompare) # number of experiments

# name of file inside diagspath
fileroot = "state_3d_set1"
filelist = searchdir(diagpath[expbase],fileroot) # first filter for state_3d_set1
datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
nt = length(datafilelist)

# Improve code here to read meta file, make variable selection transparent.
nc = 2 # do first two of three properties in state_3d_set1

# pre-allocate
σx = Array{Float32, 3}(undef, nt, nz, nc)
# xmed = similar(σ); median was slow so I dropped it.
xbar = similar(σx); xmax = similar(σx); xmin = similar(σx)

tt = 0
for fname in datafilelist
    global tt += 1
    println("time index ",tt)

    # need to get all three tracers read.
    @time xbase = γ.read(diagpath[expbase]*fname,MeshArray(γ,Float32,nz*nc))
    @time x = γ.read(diagpath[expcompare]*fname,MeshArray(γ,Float32,nz*nc))

    x = x - xbase # turn into a difference 
    ntiles = size(x,1)
    
    # get a vector list of nonzero elements for each level.
    for zz = 1:nz
        for cc = 1:nc
            xtmp = x[:,(cc-1)*nz+zz]

            xcount = [sum(count(!iszero,xtmp[i])) for i = 1:ntiles]
            if sum(xcount) > 0

                xmax[tt,zz,cc] = maximum([maximum(filter(!iszero,xtmp[i])) for i = 1:ntiles if xcount[i] > 0])
                xmin[tt,zz,cc] = minimum([minimum(filter(!iszero,xtmp[i])) for i = 1:ntiles if xcount[i] > 0])

                # compute mean the old fashioned way
                xsum = sum([sum(xtmp[i]) for i = 1:ntiles if xcount[i] > 0]) # works b.c. 0 on land
                xbar[tt,zz,cc] = xsum/sum(xcount)

                # compute standard deviation
                x²sum = sum([sum((filter(!iszero,xtmp[i]).-xbar[tt,zz,cc]).^2) for i=1:ntiles if xcount[i]>0])
                σx[tt,zz,cc] = sqrt(x²sum/(sum(xcount)-1))
                
            else
                xbar[tt,zz,cc] = missing
                xmax[tt,zz,cc] = missing
                xmin[tt,zz,cc] = missing
                σx[tt,zz,cc] = missing
            end
        end
    end
end

# save output as JLD2?
if output2file
    outpath = "../outputs/"
    outfile = outpath*"divergence_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*".jld2"
    @save outfile xbar σx xmax xmin z 
end

# Work in progress
#include("plot_divergence.jl")
