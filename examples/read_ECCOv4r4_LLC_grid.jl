# read the ECCOv4r4 lat-lon-cap (LLC) grid, return grid in γ

using Revise
using Reemergence, MeshArrays, MITgcmTools

# other useful packages in Reemergence project
# using Statistics, PyPlot, Distributions, FFTW, LinearAlgebra, StatsBase

  
workdir = pwd()
push!(LOAD_PATH, workdir)
cd(workdir)

pth="../inputs/GRID_LLC90/"

http="https://github.com/gaelforget/GRID_LLC90"
!isdir(pth) ? run(`git clone $http $pth`) : nothing;

γ=GridSpec("LatLonCap",pth)

D=γ.read(γ.path*"Depth.data",MeshArray(γ,Float64))
tmp1=write(D); tmp2=read(tmp1,D)
show(D)

LAT=γ.read(γ.path*"YC.data",MeshArray(γ,Float64))
tmp1=write(LAT); tmp2=read(tmp1,LAT)
show(LAT)

LON=γ.read(γ.path*"XC.data",MeshArray(γ,Float64))
tmp1=write(LON); tmp2=read(tmp1,LON)
show(LON)
