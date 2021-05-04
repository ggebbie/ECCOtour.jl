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

# get experiments on poseidon/batou
exppath = "/batou/ECCOv4r4/MITgcm/exps/"
runpath,diagpath = listexperiments(exppath);
# abbreviations for each experiment for labels, etc.
shortnames = expnames()

# print output here
pathout = "/batou/ECCOv4r4/MITgcm/exps/"

## SELECT EXPERIMENTS TO ANALYZE #################################
# manually choose from available experiments listed above.
#exps = ("iter129_bulkformula")

# to do all experiments:
exps = keys(shortnames)
##############################################################
nexps = length(exps) # number of experiments

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

tstart = 1992 + 1/24
tend = 2018
tecco = range(tstart,step=1/12,stop=2018)
nt = length(tecco)

lonatts = Dict("longname" => "Longitude", "units" => "degrees east")
latatts = Dict("longname" => "Latitude", "units" => "degrees north")
depthatts = Dict("longname" => "Depth", "units" => "m")

# oceqnet
Froot = "state_2d_set1"
Froot = "state_2d_set2"
#Froot = "budg2d_"

# set this up as an initial argument to the script.
ex = "iter129_bulkformula"
#for ex in exps
filelist = searchdir(diagpath[ex],Froot) 
datafilelist  = filter(x -> occursin("data",x),filelist)

# make an output directory for each experiment
pathoutexp = pathout*ex*"/run/regularpoles/"
!isdir(pathoutexp) ? mkdir(pathoutexp) : nothing;

global tt = 0

for Fname in datafilelist
    tt += 1
    println("filename ",Fname)

    year = Int(floor(tecco[tt]))
    month = ((tt-1)%12)+1
    println("year ",year," month ",month)

    #θ = read_bin(diagpath[ex]*Fname,γ)
    # had InexactError with above read. 
    θmeta = read_meta(diagpath[ex]*Fname)
    # hard coded Float32, should use θmeta.dataprec but type is wrong
    θ = γ.read(diagpath[ex]*Fname,MeshArray(γ,Float32,θmeta.nrecords))
    
    #loop over all the records in the file
    global rec = 0
    for field in θmeta.fldList
        # assume levs = 1    
        fieldDict = read_available_diagnostics(field,filename=runpath[ex]*"available_diagnostics.log")
        rec += 1
        #global θinterp = Array{Float64, 3}(undef, nx, ny)
        
        θcrop =  LLCcropC(θ[:,rec],γ) # get regular grid by cropping
        θarc = reginterp(θ[:,rec],nx,nyarc,farc,iarc,jarc,warc) # interpolate to "LLCregular"
        θantarc = reginterp(θ[:,rec],nx,nyantarc,fantarc,iantarc,jantarc,wantarc)
        θinterp=hcat(θantarc',θcrop,θarc')
        replace!(θinterp,0.0 => NaN)
            
        # make a directory for this output
        pathoutexpvar = pathoutexp*fieldDict["fldname"]*"/" 
        !isdir(pathoutexpvar) ? mkdir(pathoutexpvar) : nothing;

        # get filename for this month.
        if month < 10
            fileout = pathoutexpvar*fieldDict["fldname"]*"_"*string(year)*"_0"*string(month)*".nc"
        else 
            fileout = pathoutexpvar*fieldDict["fldname"]*"_"*string(year)*"_"*string(month)*".nc"
        end
        
        # save in a NetCDF file with info from fieldDict
        varatts = Dict("longname" => fieldDict["title"], "units" => fieldDict["units"])
         
        # Now we can write values to the file
        isfile(fileout) && rm(fileout)
        nccreate(
            fileout,
            fieldDict["fldname"],
            "lon",
            λC,
            lonatts,
            "lat",
            ϕC,
            latatts,
            atts = varatts,
        )
        ncwrite(θinterp, fileout, fieldDict["fldname"])
    end
end

end
