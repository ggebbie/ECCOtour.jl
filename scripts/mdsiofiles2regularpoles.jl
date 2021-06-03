# transfer output from native to regularpoles grid
# 1. read monthly-average fields
# 2. interpolate to Cartesian grid
# 3. save to self-describing file.
# 4. repeat with all fields.

using Revise
using Reemergence
using MeshArrays, MITgcmTools
using NetCDF

# get names of exps. 
workdir = pwd()
push!(LOAD_PATH, workdir)
cd(workdir)

## SELECT EXPERIMENTS TO ANALYZE ##
#  manually choose from available experiments listed above.
expt = "iter129_bulkformula"
Froot = "state_3d_set1"
Froot = "state_3d_set2"
Froot = "state_2d_set1"
Froot = "state_2d_set2"

###################################

# get experiments on poseidon/batou
exppath = "/batou/ECCOv4r4/MITgcm/exps/"
runpath,diagpath,regpolespath = listexperiments(exppath);

# abbreviations for each experiment for labels, etc.
shortnames = expnames()

# to do all experiments:
exps = keys(shortnames)
nexps = length(exps) # number of experiments

# print output here
# pathout = "/batou/ECCOv4r4/MITgcm/exps/"
pathout = regpolespath[expt]
!isdir(path_out) ? mkdir(path_out) : nothing;

pathgrid="../inputs/GRID_LLC90/"
γ = setupLLCgrid(pathgrid)

lat,lon = latlon(γ)
# Set up Cartesian grid for interpolation.
#ϕG,ϕC = regularlatgrid(γ) # regular grid
ϕGarc,ϕCarc = latgridArctic(γ)
ϕGantarc,ϕCantarc = latgridAntarctic(γ) 
ϕGreg,ϕCreg = latgridRegular(γ) 
λG = -180.0:179.0
λC = -179.5:179.5

ϕG = vcat(ϕGantarc,ϕGreg,ϕGarc)
ϕC = vcat(ϕCantarc,ϕCreg,ϕCarc)

farc,iarc,jarc,warc = prereginterp(ϕCarc,λC,γ)
fantarc,iantarc,jantarc,wantarc = prereginterp(ϕCantarc,λC,γ)

# Fix this (?)
nx = length(λC)
ny = length(ϕC)
nyarc = length(ϕCarc)
nyantarc = length(ϕCantarc)

# get standard levels of MITgcm
z = depthlevels(γ)
nz = length(z)

tstart = 1992 + 1/24
tend = 2018
tecco = range(tstart,step=1/12,stop=2018)
nt = length(tecco)

lonatts = Dict("longname" => "Longitude", "units" => "degrees east")
latatts = Dict("longname" => "Latitude", "units" => "degrees north")
depthatts = Dict("longname" => "Depth", "units" => "m")

#for expt in exps
filelist = searchdir(diagpath[expt],Froot) 
datafilelist  = filter(x -> occursin("data",x),filelist)

global tt = 0

for Fname in datafilelist
    tt += 1
    println("filename ",Fname)

    year,month = timestamp_monthly_v4r4(tt)

    fileoutput = diagpath[expt]*Fname
    filelog = runpath[expt]*"available_diagnostics.log"

    if month < 10
        filesuffix = "_"*string(year)*"_0"*string(month)*".nc"
    else 
        filesuffix = "_"*string(year)*"_"*string(month)*".nc"
    end

    filein = Fname[1:end-5]
    pathin = diagpath[expt]

    @time varsregpoles =  mdsio2regularpoles(pathin,filein,γ,nx,ny,nyarc,farc,iarc,jarc,warc,nyantarc,fantarc,iantarc,jantarc,wantarc)

    @time writeregularpoles(varsregpoles,γ,pathout,filesuffix,filelog,λC,lonatts,ϕC,latatts,z,depthatts)

end

