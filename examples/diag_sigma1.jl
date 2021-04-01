# Example: Use monthly average temperature/salinity fields on poseidon to compute sigma1
using Revise
using Reemergence, MITgcmTools, MeshArrays

# other useful packages in Reemergence project
# using Statistics, PyPlot, Distributions, FFTW, LinearAlgebra, StatsBase

# get monthly average temperature and salinity
include("examples/read_TS_monavg_poseidon.jl")

p₀ = 1000 # reference pressure = 1000 dbar

# ECCOv4r4 uses approximation for pressure without any horizontal deviations.
# Can precompute pressure for each depth level.
ρ₀ = 1029 # from "data" text file in run directory
g  = 9.81 # from "data" 

# depths of MITgcm levels
pathLLC = "../inputs/GRID_LLC90/" # need to input grid earlier: see example
fileZ = "RC" # "R" = radius, "C"= center of grid cell
Z = read_mdsio(pathLLC,fileZ)
Z = vec(Z)

# standard pressures via hydrostatic balance
Pa2dbar = 1/10000
p = -ρ₀ .*g .* Pa2dbar .* Z # 10000 to change Pa to dbar

# mitgcmtools.seawaterdensity is set up for a column.
# random choice of a water column.
ix = 50
iy = 400

# filter out zeroes.
θz = filter(!iszero,θ[ix,iy,:])
Sz = filter(!iszero,S[ix,iy,:])
nz = length(θz)
pz = p[1:nz]
p₀ = 1000

# ρI is in-situ density, ρR is density referenced to p0.     
(ρP,ρI,ρR) = SeaWaterDensity(θz,Sz,pz,p₀)
σ₁ = ρR .- 1000

# could add loops for the entire domain.

# to do: plot it


