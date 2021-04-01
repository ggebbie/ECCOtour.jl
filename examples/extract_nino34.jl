#  Instead use this to remove the interannual frequency energy in surface forcing fields.
#  Diagnostic plots have been removed from this version.
#  First steps: 1. go into Reemergence project directory. 2. go into julia REPL package mode with `]`. 3. `activate .` 4. Backspace to return to command mode in REPL.

using Revise
using Reemergence, MeshArrays, MITgcmTools
using PyPlot, JLD2

workdir = pwd()
push!(LOAD_PATH, workdir)
cd(workdir)

# get MITgcm / ECCOv4r4 LLC grid and depth information. Store in γ.
path_grid="../inputs/GRID_LLC90/"
γ = setupLLCgrid(path_grid)

# ϕ = lat, λ = lon
(ϕ,λ) = readlatlon(γ)
area = readarea(γ)

exppath = "/poseidon/ECCOv4r4/MITgcm/exps/"
runpath,diagpath = listexperiments(exppath);

# abbreviations for each experiment for labels, etc.
shortnames = expnames()
marks = expsymbols()
nexp = length(shortnames) # number of experiments
 
fileroot = "state_3d_set1"
dryval = 0.0
iswet(x) = x != dryval # or was this defined in Reemergence.jl

# transfer to nino34 by removing seasonal climatology. Read climatology.
inputpath = "../inputs/"
@load inputpath*"nino34_1870-2020.jld2" nino34 nino34yr SSTclimatology
nino34hadisst = nino34; thadisst = nino34yr; sst34hadisst = SSTclimatology;

# pre-define
nino34 = Dict{String,Array{Any,1}}() # don't forget trailing parentheses
nino34native = Dict{String,Array{Any,1}}() 
# to do: remove seasonal climatology from model or from reality
tecco = 1992+1/24:1/12:2018
monthsperyear = 12
fcycle = 1 # units: yr^{-1}
# for removing seasonal cycle from monthly averaged timeseries
Ecycle,Fcycle = seasonal_matrices(fcycle,tecco,overtones=4)

for (keys,values) in shortnames

    expname = keys
    println(expname)
    sst34 = extract_sst34(expname,diagpath,fileroot,γ,area,ϕ,λ,iswet)

    # remove climatology
    nino34[expname] = remove_climatology(sst34,sst34hadisst)
    nino34native[expname] = remove_seasonal(sst34,Ecycle,Fcycle)

end

# make figure with labels, legend
# get nino34 with mean of MeshArray, a patch function for this rectangle, area weighting, could involve updating "mean" function in Reemergence/MeshArrays 
ylbl  = "NINO3.4 "*L"[\degree C]"
figure(100); clf();
plot(thadisst,nino34hadisst,label="HADISST 1")
for (keys,values) in nino34
    plot(tecco,values,"-"*marks[keys],label = shortnames[keys])
end
grid(true)
legend()
ylabel(ylbl)
xlabel("calendar years")
axorig = axis()
axis((1992,2018,axorig[3],axorig[4]))
outputpath = "../outputs/"
savefig(outputpath*"nino34comparison_sameSSTscale.eps")

## make second figure with different SST baseline
ylbl  = "NINO3.4 relative to modeled SST "*L"[\degree C]"
figure(101); clf();
plot(thadisst,nino34hadisst,label="HADISST 1")
for (keys,values) in nino34native
    plot(tecco,values,"-"*marks[keys],label = shortnames[keys])
end
grid(true)
legend()
ylabel(ylbl)
xlabel("calendar years")
axorig = axis()
axis((1992,2018,axorig[3],axorig[4]))
outputpath = "../outputs/"
savefig(outputpath*"nino34comparison_nativeSSTscale.eps")

