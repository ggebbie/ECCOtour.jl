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

(lat,lon) = readlatlon(γ)
area = readarea(γ)

exppath = "/poseidon/ECCOv4r4/MITgcm/exps/"
runpath,diagpath = listexperiments(exppath);

# abbreviations for each experiment for labels, etc.
shortnames = expnames()

# find center of nino34 region. 145 W, equator
#yv = 75
#xv = 106
#fv = 4
#tmplat  = LAT[fv]; lat_point = tmplat[xv,yv]
#tmplon  = LON[fv]; lon_point = tmplon[xv,yv]

fileroot = "state_3d_set1"

nino34monthly = zeros(6,312)
expno = 0
dryval = 0.0
iswet(x) = x != dryval # or was this defined in Reemergence.jl

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


        # area-weighted filtered mean.
        #nino34[expno,tt] = mean(x,area,isnino34)
        sst34monthly[expno,tt] = patchmean(x,area,ϕ,λ,isnino34,iswet)
        
        # old style for checking.
        #nino34old[expno,tt] = x[fv][xv,yv]
    end
    
end

# transfer to nino34 by removing seasonal climatology.

inputpath = "../inputs/"
@load inputpath*"nino34_1870-2020.jld2" nino34 nino34yr SSTclimatology

# to do: remove seasonal climatology from model or from reality
eccoyr = 1992+1/24:1/12:2018

# remove climatology
nino34monthly = similar(sst34monthly)
for nmon = 1:12
    nino34monthly[:,nmon:12:end] = sst34monthly[:,nmon:12:end] .- SSTclimatology[nmon]
end

nino34yr_1992 = nino34yr[findall(nino34yr.>1992)]
nino34_1992 = nino34[findall(nino34yr.>1992)]

# make figure with labels, legend
# get nino34 with mean of MeshArray, a patch function for this rectangle, area weighting, could involve updating "mean" function in Reemergence/MeshArrays 
ylbl  = "NINO3.4 "*L"[\degree C]"
figure()
plot(nino34yr_1992,nino34_1992,label="HADISST 1")
expno = 0
for (keys,values) in shortnames
    expno += 1
    plot(eccoyr,nino34monthly[expno,:],label = values)
end
grid(true)
legend()
ylabel(ylbl)
xlabel("calendar years")
savefig("nino34comparison.eps")

