# Special case: Take EVEL* NVEL* from nctiles_monthly and map onto regular grid. 

# 1. read monthly-average fields
# 2. interpolate to Cartesian grid
# 3. save to self-describing file.
# 4. repeat

using Revise
using Reemergence
using MeshArrays, MITgcmTools
using NetCDF


exppath = "/batou/ECCOv4r4/MITgcm/exps/"
runpath,diagpath,regpolespath = listexperiments(exppath);
varnames = ("EVEL","EVELMASS","EVELSTAR","NVEL","NVELMASS","NVELSTAR","WVELMASS","WVELSTAR","oceTAUN","oceTAUE")

# get names of exps. 
workdir = pwd()
push!(LOAD_PATH, workdir)
cd(workdir)

# Read nctiles_monthly from ECCO Drive
exppath = "/batou/eccodrive/nctiles_monthly/"
# print output here
pathout = "/batou/ECCOv4r4/MITgcm/exps/iter129_bulkformula/run/regularpoles/"

expt = "iter129_bulkformula"
filelog = runpath[expt]*"available_diagnostics.log"

path_grid="../inputs/GRID_LLC90/"
γ = setupLLCgrid(path_grid)

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

# Fix this
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

for varname in varnames

    println("varname ",varname)
    exppathvar = exppath*varname*"/"

    pathoutvar = pathout*varname*"/"
    !isdir(pathoutvar) ? mkdir(pathoutvar) : nothing;

    # make directory if doesn't exist
    for tt = 1:length(tecco)

        year,month = timestamp_monthly_v4r4(tt)

        if month < 10
            filesuffix = "_"*string(year)*"_0"*string(month)*".nc"
            nctilename = exppathvar*string(year)*"/"*varname*filesuffix
        else
            filesuffix = "_"*string(year)*"_"*string(month)*".nc"
            nctilename = exppathvar*string(year)*"/"*varname*filesuffix
        end
        
        @time varsregpoles =  netcdf2regularpoles(nctilename,varname,γ,nx,ny,nyarc,farc,iarc,jarc,warc,nyantarc,fantarc,iantarc,jantarc,wantarc)

        @time writeregularpoles(varsregpoles,γ,pathout,filesuffix,filelog,λC,lonatts,ϕC,latatts,z,depthatts)
        
    end
end

