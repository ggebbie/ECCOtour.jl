# Example: read monthly average temperature/salinity fields on poseidon

using Revise
using Reemergence, MeshArrays, MITgcmTools

# other useful packages in Reemergence project
# using Statistics, PyPlot, Distributions, FFTW, LinearAlgebra, StatsBase

path = "/poseidon/ecco/ECCOv4r4/MITgcm/exps/iter129_bulkformula/run/diags/"
fileroot = "state_3d_set1" # contains θ, S, dρ/dz, see accompanying meta file

# get a list of all files in directory.
filelist = searchdir(path,fileroot) # first filter for state_3d_set1

# just save the list of ".data" files, eliminate ".meta"
datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"

# read the first one into a 4D array: lon, lat, depth, variable type
state = read_mdsio(path,datafilelist[1])

θ = state[:,:,:,1]
S = state[:,:,:,2]

maximum(θ)
maximum(S)
minimum(θ)
minimum(S)

