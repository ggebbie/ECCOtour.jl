# map θ, S, p to sigma 1 surfaces. Test for one ocean column.
# This is a wrapper routine to read files on poseidon.
# ggebbie, 7-Apr-2021

using Revise 
using Reemergence, Dierckx, PyPlot

pz = 0:500:4000 # pressure levels
nz = length(pz)
p₀ = 1000 # reference pressure = 1000 dbar
θz = range(20,stop=10,length=nz)
Sz = range(36,stop=35,length=nz)
sig1grida = 24:0.05:31
sig1gridb = 31.02:0.02:33
sig1grid = vcat(sig1grida,sig1gridb)
sig1grid = sig1grid[1:3:end]
splorder = 3

@time θσ,Sσ,pσ=TSP2sigma1column(collect(θz),collect(Sz),collect(pz),sig1grid,splorder)
@time θσ,Sσ,pσ=TSP2sigma1column(θz,Sz,pz,sig1grid,splorder)

figure(100)
plot(θz,-pz)
plot(θσ,-pσ,"o")

figure(101)
plot(Sz,-pz)
plot(Sσ,-pσ,"o")
