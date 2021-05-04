# Special case: Take EVEL* NVEL* from nctiles_monthly and map onto regular grid. 

# 1. read monthly-average fields
# 2. interpolate to Cartesian grid
# 3. save to self-describing file.
# 4. repeat with all surface fields.

using Revise
using Reemergence
using MeshArrays, MITgcmTools
#using Statistics, PyPlot, Distributions
#using LinearAlgebra, StatsBase
using NetCDF

# get names of exps. 
workdir = pwd()
push!(LOAD_PATH, workdir)
cd(workdir)

exppath = "/batou/ECCOv4r4/MITgcm/exps/"
runpath,diagpath = listexperiments(exppath);

# Read nctiles_monthly from ECCO Drive
exppath = "/batou/eccodrive/nctiles_monthly/"
# print output here
pathout = "/batou/ECCOv4r4/MITgcm/exps/iter129_bulkformula/run/regularpoles/"

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
fileZ = "RC"
z = read_mdsio(path_grid,fileZ)
z = vec(z)
nz = length(z)

tstart = 1992 + 1/24
tend = 2018
tecco = range(tstart,step=1/12,stop=2018)
nt = length(tecco)

lonatts = Dict("longname" => "Longitude", "units" => "degrees east")
latatts = Dict("longname" => "Latitude", "units" => "degrees north")
depthatts = Dict("longname" => "Depth", "units" => "m")
varnames = ("EVEL","EVELMASS","EVELSTAR","NVEL","NVELMASS","NVELSTAR","WVELMASS","WVELSTAR")

for varname in varnames
    exppathvar = exppath*varname*"/"

    pathoutvar = pathout*varname*"/"
    !isdir(pathoutvar) ? mkdir(pathoutvar) : nothing;

    # make directory if doesn't exist
    for tt = 1:length(tecco)

        year = Int(floor(tecco[tt]))
        month = ((tt-1)%12)+1
        println("year ",year," month ",month)
        if month < 10
            nctilename = exppathvar*string(year)*"/"*varname*"_"*string(year)*"_0"*string(month)*".nc"
        else
            nctilename = exppathvar*string(year)*"/"*varname*"_"*string(year)*"_"*string(month)*".nc"
        end
        #output = ncread(nctilename,varname)
        #output1 = reshape(output[:,:,:,:,1],90,90*13,50)
        #θ = convert2gcmfaces(output1,γ)

        θ = read_netcdf(nctilename,varname,γ);
        
        if varname[1] == 'E'
            field = "UE_VEL_C"
        elseif varname[1] == 'N'
            field = "VN_VEL_C"
        elseif varname[1] == 'W'
            field = "WVEL"
        end
        fieldDict = read_available_diagnostics(field,filename=runpath["iter129_bulkformula"]*"available_diagnostics.log")
        levs = size(θ,2)
        global θinterp = Array{Float64, 3}(undef, nx, ny, levs)
rec = 1
        for rec = 1:levs
            θcrop =  LLCcropC(θ[:,rec],γ) # get regular grid by cropping
            θarc = reginterp(θ[:,rec],nx,nyarc,farc,iarc,jarc,warc) # interpolate to "LLCregular"
            θantarc = reginterp(θ[:,rec],nx,nyantarc,fantarc,iantarc,jantarc,wantarc)
            θinterp[:,:,rec]=hcat(θantarc',θcrop,θarc')

        end
        replace!(θinterp,0.0 => NaN)
        
        # get filename for this month.
        if month < 10
            fileout = pathoutvar*varname*"_"*string(year)*"_0"*string(month)*".nc"
        else
            fileout = pathoutvar*varname*"_"*string(year)*"_"*string(month)*".nc"
        end
        
        # save in a NetCDF file with info from fieldDict
        varatts = Dict("longname" => fieldDict["title"], "units" => fieldDict["units"])
         
        isfile(fileout) && rm(fileout)
        nccreate(
            fileout,
            varname,
            "lon",
            λC,
            lonatts,
            "lat",
            ϕC,
            latatts,
            "depth",
            -z,
            depthatts,
            atts = varatts,
        )
            ncwrite(θinterp, fileout, varname)
        end
    end
end

