#  First steps: 1. go into Reemergence project directory. 2. go into julia REPL package mode with `]`. 3. `activate .` 4. Backspace to return to command mode in REPL.

# 1. read monthly-average surface (2d) fields
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


pathout = "/batou/ECCOv4r4/MITgcm/exps/"

path_grid="../inputs/GRID_LLC90/"
γ = setupLLCgrid(path_grid)

lat,lon = latlon(γ)
# Set up Cartesian grid for interpolation.
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

lonatts = Dict("longname" => "Longitude", "units" => "degrees east")
latatts = Dict("longname" => "Latitude", "units" => "degrees north")
depthatts = Dict("longname" => "Depth", "units" => "m")

# make an output directory for each experiment
pathoutexp = pathout*ex*"/run/regularpoles/"
!isdir(pathoutexp) ? mkdir(pathoutexp) : nothing;


D=γ.read(γ.path*"Depth.data",MeshArray(γ,Float64))

# choice: replace zeros with NaN32?
replace!(D,0.0 => NaN)


θcrop =  LLCcropC(D,γ) # get regular grid by cropping
θarc = reginterp(D,nx,nyarc,farc,iarc,jarc,warc) # interpolate to "LLCregular"
θantarc = reginterp(D,nx,nyantarc,fantarc,iantarc,jantarc,wantarc)
θinterp=hcat(θantarc',θcrop,θarc')
            
        # make a directory for this output
pathoutexpvar = "/home/gebbie/julia/outputs/"
!isdir(pathoutexpvar) ? mkdir(pathoutexpvar) : nothing;
fileout = "Depth_regularpoles.nc"
        
# save in a NetCDF file with info from fieldDict
varatts = Dict("longname" => "Depth", "units" => "m")
         
# Now we can write values to the file
        isfile(fileout) && rm(fileout)
        nccreate(
            fileout,
            "Depth",
            "lon",
            λC,
            lonatts,
            "lat",
            ϕC,
            latatts,
            atts = varatts,
        )
        ncwrite(θinterp, fileout, "Depth")

