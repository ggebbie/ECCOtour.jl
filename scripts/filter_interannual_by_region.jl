#  Instead use this to remove the interannual frequency energy in surface forcing fields.
#  Diagnostic plots have been removed from this version.
#  This code extends filter_interannual.jl to allow for regional masks. 
#  First steps: 1. go into ECCOonPoseidon project directory. 2. go into julia REPL package mode with `]`. 3. `activate .` 4. Backspace to return to command mode in REPL.

using Revise
using Reemergence # this name may change
using Statistics, PyPlot, Distributions, FFTW, LinearAlgebra, StatsBase
using MeshArrays, MITgcmTools

workdir = pwd()
push!(LOAD_PATH, workdir)
cd(workdir)

# get MIT GCM native grid
path_grid="../inputs/GRID_LLC90/"
γ = setupLLCgrid(path_grid)

# depth [meters]
D=γ.read(γ.path*"Depth.data",MeshArray(γ,Float64))

# read lat, lon at center of grid cell
(ϕC,λC) = latlonC(γ)
# on the vector (Staggered) grid
(ϕG,λG) = latlonG(γ)

# Read one flux at one point to get information about number of time points.
# read directly from poseidon
inputdir = "/poseidon/ecco.jpl.nasa.gov/drive/files/Version4/Release4/other/flux-forced/forcing/"

# change the directory name, keep file names the same as previous versions
regionname = "interannual_southpac"
outputdir="../output_" * regionname * "/"
# ultimately want it here /poseidon/ecco.jpl.nasa.gov/drive/files/Version4/Release4/other/flux-forced/forcing_interannual_southpac
 
## Define region of interest.
latrect = (-90, -15) # immutable
lonrect = [150,-67] # mutable for wraparound 
dlat = 10
dlon = 10

if lonrect[1] > lonrect[2] # then handle wraparound
    lonrect[2] += 360  # shift one way, could've shifted the opposite way
end
# center the longitudinal coordinate system
lonmid =  (lonrect[1]+lonrect[2])/2

# now λC, λG are close to centered for all calcs
centerlon!(λC,lonmid)
centerlon!(λG,lonmid)

if !isdir(outputdir)
    mkdir(outputdir)
end
midname = "_6hourlyavg_"
varnames = ("atmPload","oceFWflx","oceQsw","oceSflux","oceSPflx","oceTAUE","oceTAUN","oceTAUX",
            "oceTAUY","sIceLoad","sIceLoadPatmPload","sIceLoadPatmPload_nopabar","TFLUX")

# sample calcs at one point. Get the length of timeseries in nseries.
yv = 40
xv = 180
fv = 5
tmplat  = ϕC[fv]; lat_point = tmplat[xv,yv]
tmplon  = λC[fv]; lon_point = tmplon[xv,yv]
frootsample = inputdir*varnames[end]*midname
years = 1992:2017
fluxsample_point,nseries = extract_timeseries(frootsample,years,γ,xv,yv,fv)

# take biweekly mean. Use triangular filter with break points at:
Δi14day = 4*14 # grid index range
nt6hr = length(fluxsample_point)
i14day = 2.5:Δi14day:nt6hr+56 # goes past end of time by 14 days to be sure

t6hr_start = 1/8 # 3Z Jan 1 1992
Δt6hr = 6/24 # 1/4 of a day
t6hr = range(t6hr_start,step=Δt6hr,length=nt6hr)

# Given a full-resolution timeseries, what values at the tiepoints best represent the timeseries?
# Determine tiepoints (in time) where fluxes are adjusted by ECCO in optimization procedure.
Δt14day = 14 # units: days
t14day_start = 1/2 # 12Z Jan 1 1992 by inspection of figure
t14day = range(t14day_start,step=Δt14day,stop=t6hr[end]+14)# goes past the end by one
nt14day = length(t14day)

# Values are added to the t14day tiepoints, then linearly interpolated to fill gaps.
# Careful not to store high-resolution data all at same time.
# Make a function that takes values at the tiepoints and then makes a full-resolution timeseries.
@time E14to6,F6to14 = get_filtermatrix(t6hr,t14day)

daysperyear = 365.25
fcycle = 1/(daysperyear) # units: day^{-1}
# for removing seasonal cycle from 14-day averaged timeseries
Ecycle,Fcycle = seasonal_matrices(fcycle,t14day)

# interannual filter is a Hann(ing) filter
Thann = 100.0 # days

# Is it possible to solve for regional mask before the variable loop? Yes, if the variables are on the same grid. (Double check that they all apply to the center of a grid cell.)
#maskN, maskS, maskW, maskE = regional_mask(ϕC,λC,latrect,lonrect,dlat,dlon)
# clash with `mask` name
msk = regional_mask(ϕC,λC,latrect,lonrect,dlat,dlon)

# check the spatial pattern of the mask
# Set up Cartesian grid for interpolation.
λCregpoles,λGregpoles,ϕCregpoles,ϕGregpoles,nx,ny,nyarc,nyantarc,farc,iarc,jarc,warc,fantarc,iantarc,jantarc,wantarc =
    factors4regularpoles(γ)
msk_regpoles =  var2regularpoles(msk,γ,nx,ny,nyarc,farc,iarc,jarc,warc,nyantarc,fantarc,iantarc,jantarc,wantarc)

figure()
clf()
cmap_seismic =get_cmap("seismic")
lims = range(0.0,step=0.05,stop=1.0)
contourf(λCregpoles,ϕCregpoles,msk_regpoles',lims,cmap=cmap_seismic)
colorbar(label="weight",orientation="vertical",ticks=lims)
outfname = outputdir*"southpac_mask.eps"
xlbl = "longitude "*L"[\degree E]"
ylbl = "latitude "*L"[\degree N]"
titlelbl = "Southpac Mask"
title(titlelbl)
xlabel(xlbl)
ylabel(ylbl)
savefig(outfname)

#vname = varnames[1] # for interactive use
for vname ∈ varnames
    filein = inputdir*vname*midname
    fileout = outputdir*vname*midname
    println(filein)

    # some tracers have a staggered grid.
    # will have to re-do the mask unfortunately to be sure.
    if vname == "oceTAUX"
        msk = regional_mask(ϕC,λG,latrect,lonrect,dlat,dlon)
    elseif vname == "oceTAUY"
        msk = regional_mask(ϕG,λC,latrect,lonrect,dlat,dlon)
    else
        msk = regional_mask(ϕC,λC,latrect,lonrect,dlat,dlon)
    end
    
    # use this F to decompose controllable/uncontrollable parts
    # i.e., 6 hourly to 14 day
    flux_14day = matrixfilter(F6to14,filein,years,γ)

    # remove seasonal cycle from 14-day averaged timeseries
    # solve for seasonal cycle parameters
    # use 14-day because it's computationally efficient
    βcycle = matmul(Fcycle,flux_14day,γ)

    # reconstruct the full seasonal cycle.
    flux_14day_seasonal = matmul(Ecycle,βcycle,γ)

    # remove seasonal signal from total signal
    flux_14day_noseasonal = flux_14day - flux_14day_seasonal

    # take the noseasonal timeseries and run a filter (Thann [days] cutoff)
    # hanning filter for all locations.
    # output is interannual signal.
    flux_14day_lopass = hannfilter(flux_14day_noseasonal,t14day,t14day,Thann,γ)

    #use function that takes forcing field and multiplies each entry by a number
    #from 0 to 1 (1 being removed entirely) based on a provided lat/lon box and sponge layer widths
    # pre-compute mask before looping over all variables.
    apply_regional_mask!(flux_14day_lopass,1.0 .- msk)

    # put tflux_14day_lopass on to 6hr
    # check for NaN's in output
    nancount_lopass = sum(nancount_gcmarray(flux_14day_lopass))
    nancount_14day = sum(nancount_gcmarray(flux_14day))
    nancount_14day_seasonal = sum(nancount_gcmarray(flux_14day_seasonal))
    nancount_14day_noseasonal = sum(nancount_gcmarray(flux_14day_noseasonal))

    if nancount_lopass + nancount_14day + nancount_14day_seasonal + nancount_14day_noseasonal > 0
        error("NaNs in the filtered output")
    end

    # write the hi-pass filtered tflux
    # need to get it on the 6hr timesteps
    # need to write it for every year
    matrixspray(E14to6,-flux_14day_lopass,filein,fileout,years,γ) #changed rmfield to _regional

end
