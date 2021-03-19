#  Instead use this to remove the interannual frequency energy in surface forcing fields.
#  Diagnostic plots have been removed from this version.
#  First steps: 1. go into Reemergence project directory. 2. go into julia REPL package mode with `]`. 3. `activate .` 4. Backspace to return to command mode in REPL.

using Revise
using Reemergence, MeshArrays, MITgcmTools
using PyPlot

workdir = pwd()
push!(LOAD_PATH, workdir)
cd(workdir)

# get MITgcm / ECCOv4r4 LLC grid and depth information. Store in γ.
path_grid="../inputs/GRID_LLC90/"
γ = setupLLCgrid(path_grid)

LAT=γ.read(γ.path*"YC.data",MeshArray(γ,Float64))
tmp1=write(LAT); tmp2=read(tmp1,LAT)
show(LAT)

LON=γ.read(γ.path*"XC.data",MeshArray(γ,Float64))
tmp1=write(LON); tmp2=read(tmp1,LON)
show(LON)

#
exppath = "/poseidon/ecco/ECCOv4r4/MITgcm/exps/"
runpath,diagpath = listexperiments(exppath);

# abbreviations for each experiment for labels, etc.
shortnames = expnames()

# find center of nino34 region. 145 W, equator
yv = 75
xv = 106
fv = 4
tmplat  = LAT[fv]; lat_point = tmplat[xv,yv]
tmplon  = LON[fv]; lon_point = tmplon[xv,yv]

fileroot = "state_3d_set1"

nino34 = zeros(6,312)
expno = 0
for (keys,values) in shortnames

    expname = keys
    expno += 1
    println(expname)
    println(expno)
    filelist = searchdir(diagpath[expname],fileroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    nt = length(datafilelist)
    tt = 0
    for fname in datafilelist
        tt += 1
        println("time index ",tt)

        x = γ.read(diagpath[expname]*fname,MeshArray(γ,Float32,1))

        nino34[expno,tt] = x[fv][xv,yv]
    end
    
end


# to do: remove seasonal climatology from model or from reality
# make figure with labels, legend
# get nino34 with mean of MeshArray, a patch function for this rectangle, area weighting, could involve updating "mean" function in Reemergence/MeshArrays 
