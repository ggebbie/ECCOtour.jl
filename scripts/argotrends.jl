
#= I’ve copied my Argo mean and trends into matlab files  to batou:/hpcshare1/shared/argoclim
there are individual variable files - both for sigma gridded and pressure gridded data.
Let me know if you can parse these. There are readme and note fields that I hope explain the variables. =#

using Revise
using Reemergence
using MAT
#using Statistics, PyPlot, Distributions, LinearAlgebra, StatsBase
using MeshArrays, MITgcmTools
using  PyPlot

argoclimpath = "/batou/argoclim"
cmap_seismic =get_cmap("seismic")

clim_sigma1_ptmp = matread("clim_sigma1_ptmp.mat")
clim_pres_ptmp = matread("clim_pres_ptmp.mat")

zz = 35
field = clim_pres_ptmp["tr"][:,:,zz]
lonargo = vec(clim_sigma1_ptmp["xi"])
latargo = vec(clim_sigma1_ptmp["yi"])
    mx = maximum(filter(!isnan,field)) # filter out NaN
mn = minimum(filter(!isnan,field)) # filter out NaN
extrm = max(mx,-mn)
lims = range(-extrm,step=2extrm/30,stop=extrm)

figure(10)
clf()
contourf(lonargo,latargo,field',lims,cmap=cmap_seismic)
colorbar()
