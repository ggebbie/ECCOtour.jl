#  Instead use this to remove the interannual frequency energy in surface forcing fields.
#  Diagnostic plots have been removed from this version.
#  First steps: 1. go into Reemergence project directory. 2. go into julia REPL package mode with `]`. 3. `activate .` 4. Backspace to return to command mode in REPL.

using Revise
using Reemergence
using Statistics, PyPlot, Distributions, FFTW, LinearAlgebra, StatsBase
using MeshArrays, MITgcmTools

workdir = pwd()
push!(LOAD_PATH, workdir)
cd(workdir)

pth="../inputs/GRID_LLC90/"

http="https://github.com/gaelforget/GRID_LLC90"
!isdir(pth) ? run(`git clone $http $pth`) : nothing;

γ=GridSpec("LatLonCap",pth)

D=γ.read(γ.path*"Depth.data",MeshArray(γ,Float64))
tmp1=write(D); tmp2=read(tmp1,D)
show(D)

LAT=γ.read(γ.path*"YC.data",MeshArray(γ,Float64))
tmp1=write(LAT); tmp2=read(tmp1,LAT)
show(LAT)

LON=γ.read(γ.path*"XC.data",MeshArray(γ,Float64))
tmp1=write(LON); tmp2=read(tmp1,LON)
show(LON)

# p=dirname(pathof(MeshArrays))
# include(joinpath(p,"../examples/Plots.jl"))
# heatmap(D,title="Ocean Depth",clims=(0.,6000.))
# size(D)
# eltype(D)
# view(D,:)

###################################################
# Determine tiepoints where fluxes are adjusted.
# Use standard forcing to determine these points.

# read directly from poseidon
inputdir = "/poseidon/ecco/ecco.jpl.nasa.gov/drive/files/Version4/Release4/other/flux-forced/forcing/"
outputdir="../outputs/"
#outputdir = pwd()*"/flux-forced-nointerannual-test/"

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
tmplat  = LAT[fv]; lat_point = tmplat[xv,yv]
tmplon  = LON[fv]; lon_point = tmplon[xv,yv]
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
Δt14day = 14 # units: days
t14day_start = 1/2 # 12Z Jan 1 1992 by inspection of figure
t14day = range(t14day_start,step=Δt14day,stop=t6hr[end]+14)# goes past the end by one
nt14day = length(t14day)

# Values are added to the t14day tiepoints, then linearly interpolated to fill gaps.
# careful not to store high-resolution all at same time.
# Make a function that takes values at the tiepoints and then makes a full-resolution timeseries.
@time E14to6,F6to14 = get_filtermatrix(t6hr,t14day)

daysperyear = 365.25
fcycle = 1/(daysperyear) # units: day^{-1}
# for removing seasonal cycle from 14-day averaged timeseries
Ecycle,Fcycle = seasonal_matrices(fcycle,t14day)

Lhann = 100.0 # days

for vname ∈ varnames
    filein = inputdir*vname*midname
    fileout = outputdir*vname*midname*"nointerannual_"
    println(filein)

    # use this F to decompose controllable/uncontrollable parts
    flux_14day = matrixfilter(F6to14,filein,years,γ)

    # remove seasonal cycle from 14-day averaged timeseries
    # solve for seasonal cycle parameters
    βcycle = matmul(Fcycle,flux_14day,γ)

    # reconstruct the full seasonal cycle.
    flux_14day_seasonal = matmul(Ecycle,βcycle,γ)

    # remove it from total signal
    flux_14day_noseasonal = flux_14day - flux_14day_seasonal

    # take the noseasonal timeseries and run a high-pass filter (9-month cutoff)
    # hanning filter for all locations.
    flux_14day_lopass = hannfilter(flux_14day_noseasonal,t14day,t14day,Lhann,γ)

    # put tflux_14day_lopass on to 6hr
    ## MOVE THIS TO THE STEP BEFORE WRITING TO FILE.
    ##  check for NaN's in output   ###################
    nancount_lopass = sum(nancount_gcmarray(flux_14day_lopass))
    nancount_14day = sum(nancount_gcmarray(flux_14day))
    nancount_14day_seasonal = sum(nancount_gcmarray(flux_14day_seasonal))
    nancount_14day_noseasonal = sum(nancount_gcmarray(flux_14day_noseasonal))

    if nancount_lopass + nancount_14day + nancount_14day_seasonal + nancount_14day_noseasonal > 0.5
        error("NaNs in the filtered output")
    end
    
    # write the hi-pass filtered tflux
    # need to get it on the 6hr timesteps
    # need to write it for every year
    matrixspray(E14to6,-flux_14day_lopass,filein,fileout,years,γ)

end



