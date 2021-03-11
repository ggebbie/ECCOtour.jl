# Example: read monthly average temperature/salinity fields on poseidon
# Here compare reading with read_mdsio versus MeshArrays built-in method.
# At bottom, also comparing statistics computed via two types of files.
using Revise
using Reemergence, MeshArrays, MITgcmTools, Statistics

# other useful packages in Reemergence project
# using Statistics, PyPlot, Distributions, FFTW, LinearAlgebra, StatsBase
path_grid="../inputs/GRID_LLC90/"
γ = setupLLCgrid(path_grid)

path = "/poseidon/ecco/ECCOv4r4/MITgcm/exps/iter129_bulkformula/run/diags/"
fileroot = "state_3d_set1" # contains θ, S, dρ/dz, see accompanying meta file

# get a list of all files in directory.
filelist = searchdir(path,fileroot) # first filter for state_3d_set1

# just save the list of ".data" files, eliminate ".meta"
datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"

# should read this from the meta file.
nz = 50  # number of depths
nc = 3   # number of variables

# read the first one into a 4D array: lon, lat, depth, variable type
@time state1 = read_mdsio(path,datafilelist[1]);
@time state2 = γ.read(path*datafilelist[1],MeshArray(γ,Float32,nz*nc));

#compare salinity at depth level 20.
level = 48
propertyno = 2;
ntiles = size(state2,1)
S1 = state1[:,:,level,propertyno]
S2 = state2[:,(propertyno-1)*nz+level]

# compare maxima
xcount = [sum(count(!iszero,S2[i])) for i = 1:ntiles]
@time maxtest = maximum([maximum(filter(!iszero,S2[i])) for i = 1:ntiles if xcount[i] > 0])
@assert maxtest == maximum(filter(!iszero,S1)) # filter out zeroes

# compare minimum
@time mintest = minimum([minimum(filter(!iszero,S2[i])) for i = 1:ntiles if xcount[i] > 0])
@assert mintest == minimum(filter(!iszero,S1)) # filter out zeroes

# compute mean of S1, S2
xsum = sum([sum(S2[i]) for i = 1:ntiles if xcount[i] > 0]) # works b.c. 0 on land
xbar = xsum/sum(xcount)
@assert xbar == mean(filter(!iszero,S1))

x²sum = sum([sum((filter(!iszero,S2[i]).-xbar).^2) for i=1:ntiles if xcount[i]>0])
σ = sqrt(x²sum/(sum(xcount)-1))
@assert σ == std(filter(!iszero,S1)) 
