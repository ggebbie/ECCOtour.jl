module ECCOtour
# write new functions and put them in this module.
# add them with text below, or create a new file in "src" and include it.

using MeshArrays, MITgcmTools, SigmaShift
using Statistics, PyPlot, Distributions, FFTW,
    LinearAlgebra, StatsBase, LaTeXStrings,
    Dierckx, DelimitedFiles, Interpolations, NetCDF

import Statistics.mean, Base.maximum, Base.minimum, Base.replace!, Statistics.std

export hanncoeffs, hannsum, hannsum!, hannfilter
export get_filtermatrix, matrixfilter, matrixspray, columnscale!
export seasonal_matrices, trend_matrices
export position_label, searchdir, setupLLCgrid
export listexperiments, expnames, expsymbols, time_label
export inrectangle, isnino34, issouthpac, isnino3, isnino4, isnino12
export latlon, latlonC, latlonG
export depthlevels, pressurelevels, readarea, patchmean
export nino34mean, nino3mean, nino4mean, nino12mean, extract_sst34
export remove_climatology, remove_seasonal
export prereginterp, reginterp, trend_theta!, trend_theta
export regularlatgrid, LLCcropC, LLCcropG, croplimitsLLC
export latgridAntarctic, latgridArctic, latgridRegular
export timestamp_monthly_v4r4, netcdf2dict, write_vars
export mdsio2dict
export mdsio2sigma1, ncwritefromtemplate
export netcdf2sigma1, mdsio2regularpoles
export writeregularpoles, vars2regularpoles, var2regularpoles
export netcdf2regularpoles, factors4regularpoles
export regional_mask, apply_regional_mask!, zero2one!, wrapdist
export centerlon!, read_netcdf
export extract_timeseries,matmul,position_label,nancount
export faststats, allstats, std, mean
export maximum, minimum, mean, std, replace!
export velocity2center, rotate_uv, rotate_velocity!

include("HannFilter.jl")
include("MatrixFilter.jl")
include("SeasonalCycle.jl")

"""
    function position_label(lon,lat)
# Arguments
- `lon`: longitude
- `lat`: latitude
# Output
- `lbl`: LaTex string used for figure labels
"""
function position_label(lon,lat)
    # produce label for a title by rounding to nearest whole integer.
    if lat >= 0
        latlbl = string(round(Integer,lat))* L"{\degree}N"
    else
        latlbl = string(round(Integer,-lat))* L"{\degree}S"
    end

    if lon >= 0
        lonlbl = string(round(Integer,lon))* L"{\degree}E"
    else
        lonlbl = string(round(Integer,-lon))* L"{\degree}W"
    end

    lbl = latlbl * " " * lonlbl
    return lbl
end

"""
    function searchdir(path,key)
    a useful one-line function
# Arguments
- `path`: directory to search for file
- `key`: expression to search for in path
"""
searchdir(path,key) = filter(x->occursin(key,x), readdir(path))

"""
    function inrectangle(lonin,latin,lons,lats)
    find all gridcells within Cartesian rectangle
# Arguments
- `lons`: tuple of longitude range
- `lats`: tuple of latitude range
# Output
- `rectangle`: boolean of type gcmarray
"""
function inrectangle(latpt,lonpt,latrect,lonrect)

    # load latitude and longitude on this grid.
    if lonrect[1] < lonrect[2]
        inrectangle = minimum(lonrect) <= lonpt <= maximum(lonrect) &&
            minimum(latrect) <= latpt <= maximum(latrect)
    else # wraps around the date line
        inrectangle = (lonpt <= minimum(lonrect) || #this OR was an AND before -- think it should be an OR?
            lonpt >= maximum(lonrect)) &&
            minimum(latrect) <= latpt <= maximum(latrect)
    end
    # watch out if lonrect[1] == lonrect[2]
end

issouthpac(lat,lon) = inrectangle(lat,lon,(-90,-15),(150,-67))
isnino34(lat,lon) = inrectangle(lat,lon,(-5,5),(-170,-90))
isnino3(lat,lon) = inrectangle(lat,lon,(-5,5),(-150,-90))
isnino4(lat,lon) = inrectangle(lat,lon,(-5,5),(-170,-120))
isnino12(lat,lon) = inrectangle(lat,lon,(-10,0),(-90,-80))

function latlon(γ)
    ϕ,λ = latlonC(γ)
return ϕ,λ
end

"""
    function latlonC(γ)
    Latitude-longitude of ECCOv4r4 "C" (tracer) grid
"""
function latlonC(γ)
    ϕ=γ.read(γ.path*"YC.data",MeshArray(γ,Float64))
    λ=γ.read(γ.path*"XC.data",MeshArray(γ,Float64))
return ϕ,λ
end

"""
    function latlongG(γ)
    Latitude-longitude of ECCOv4r4 "G" (velocity) grid
"""
function latlonG(γ)
    ϕ=γ.read(γ.path*"YG.data",MeshArray(γ,Float64))
    λ=γ.read(γ.path*"XG.data",MeshArray(γ,Float64))
return ϕ,λ
end

"""
    function depthlevels(γ)
    Depths of ECCO v4r4 grid (positive vals)
"""
function depthlevels(γ)
    z = read_mdsio(γ.path,"RC")
    z = -vec(z)
    return z
end

"""
    function pressurelevels(z)
    Standard pressures of ECCO v4r4 grid
"""
function pressurelevels(z)
    # ECCOv4r4 approximates pressure without any horizontal deviations in EOS.
    # Can precompute pressure for each depth level.
    ρ₀ = 1029 # from "data" text file in run directory
    g  = 9.81 # from "data"
    Pa2dbar = 1/10000 # standard pressures via hydrostatic balance
    pstdz = ρ₀ .*g .* Pa2dbar .* z # 10000 to change Pa to dbar
    return pstdz
end

"""
    function readarea(γ)
    area of ECCO v4r4 grid
"""
function readarea(γ)
    area=γ.read(γ.path*"RAC.data",MeshArray(γ,Float64))
return area
end

"""
     function patchmean(x,area,ϕ,λ,ispatch,iswet)
     get weight for rectangle region.
# Arguments
- `x`: variable of interest
- `area`: weighting
- `ϕ`: lat
- `λ`: lon
- `ispatch`: true in patch of interest
- `iswet` = function that is true if in ocean
# Output
- `xbar`: weighted filtered average of x
"""
function patchmean(x,area,ϕ,λ,ispatch,iswet)
    x = x*area
    xsum = 0.0
    asum = 0.0
    for i in eachindex(x)
        mask = ispatch.(ϕ[i],λ[i]) .* iswet.(x[i])
        for j in eachindex(ϕ[i])
            if mask[j]
                xsum += x[i][j]
                asum += area[i][j]
            end
        end
    end
    if asum > 0.0
        xbar = xsum/asum
    else
        xbar = NaN32
    end
    return xbar
end

nino34mean(x,area,ϕ,λ,iswet) = patchmean(x,area,ϕ,λ,isnino34,iswet)
nino3mean(x,area,ϕ,λ,iswet) = patchmean(x,area,ϕ,λ,isnino3,iswet)
nino4mean(x,area,ϕ,λ,iswet) = patchmean(x,area,ϕ,λ,isnino4,iswet)
nino12mean(x,area,ϕ,λ,iswet) = patchmean(x,area,ϕ,λ,isnino12,iswet)

isnotpositive(x) = (abs(x) == -x)


"""
    function time_label(index)
    Time label useful in plots. For ECCOv4r4.
# Input
- `index::Integer`: index number for monthly average field, number of months since Jan. 1992
# Output
- `tlbl`: label with month name (abbrev. with 3 letters) and calendar year CE
"""
function time_label(index::Integer)
    startyear = 1992 # could be an input argument
    monlbl = "Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"
    year = convert(Integer,floor(startyear + ( index./12)))
    mon  = index%12 + 1
    monthlbl = monlbl[mon]
    yearlbl = string(year)
    tlbl = monthlbl*" "*yearlbl
    return tlbl
end

"""
    function allstats(x)
    Compute fast statistics of gcmgrid type using function calls, same as faststats
# Input
- `x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}`: input of gcmarray type
# Output
- `xbar::Float32`: mean value after discarding dryval
- `xmax::Float32`: maximum value
- `xmin::Float32`: minimum value
- `σx::Float32`: standard deviation
- `absxbar::Float32`: mean of absolute value
"""
function allstats(x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}})

    dryval = 0.f0 # land points are zeroes, use NaN32 for NaN's

    xmax = maximum(x,dryval)
    xmin = minimum(x,dryval)
    xbar = mean(x,dryval)

    # why does it require an explicit call to MeshArrays?
    absxbar = MeshArrays.mean(abs.(x),dryval)
    σx   = std(x,xbar,dryval)

    return xbar, xmax, xmin, σx, absxbar
end

"""
    function faststats(x)
    Compute fast statistics of gcmgrid type using function calls, eliminate redundancy
# Input
- `x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}`: input of gcmarray type
# Output
- `xbar::Float32`: mean value after discarding dryval
- `xmax::Float32`: maximum value
- `xmin::Float32`: minimum value
- `σx::Float32`: standard deviation
- `absxbar::Float32`: mean of absolute value
"""
function faststats(x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}})

    #dryval = 0.f0 # land points are zeroes, use NaN32 for NaN's
    #dryval = NaN32
    #isdry(z) = (z == dryval)

    xcount = [sum(count(notnanorzero,x[i])) for i in eachindex(x)]

    if sum(xcount)>0
        xmax = maximum([maximum(filter(notnanorzero,x[i])) for i ∈ eachindex(x) if xcount[i] > 0])
        xmin = minimum([minimum(filter(notnanorzero,x[i])) for i ∈ eachindex(x) if xcount[i] > 0])

        # compute mean the old fashioned way
        xsum = sum([sum(filter(notnanorzero,x[i])) for i ∈ eachindex(x) if xcount[i] > 0])
        xbar = xsum/sum(xcount)

        # compute standard deviation
        x′ = x.-xbar
        x²sum = sum([sum(filter(notnanorzero,x′[i]).^2) for i ∈ eachindex(x′) if xcount[i]>0])
        σx = sqrt(x²sum/(sum(xcount)-1))

        absxsum = sum([sum(filter(notnanorzero,abs.(x[i]))) for i ∈ eachindex(x) if xcount[i] > 0]) # works b.c. 0 on land
        absxbar = absxsum/sum(xcount)
    else
        xbar = NaN
        xmax = NaN
        xmin = NaN
        σx = NaN
        absxbar = NaN
    end
    return xbar, xmax, xmin, σx, absxbar
end

"""
    function extract_sst34
    extract by reading multiple files
"""
function extract_sst34(expname,diagpath,fileroot,γ,area,ϕ,λ,iswet)
    filelist = searchdir(diagpath[expname],fileroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nt = length(datafilelist)
    tt = 0
    sst34 = Float32[]
    for fname in datafilelist
        tt += 1
        println("time index ",tt)

        x = γ.read(diagpath[expname]*fname,MeshArray(γ,Float32,1))

        # area-weighted filtered mean.
        sst34 = push!(sst34,nino34mean(x,area,ϕ,λ,iswet))
    end
    return sst34
end

"""
    function remove_climatology(x,xbar)
    `x` = long monthly timeseries, starting in Jan
    `xbar` = 12-month climatology, starting in Jan
    remove climatology
"""
function remove_climatology(x,xbar)
    xprime = similar(x)
    for nmon = 1:12
        xprime[nmon:12:end] = x[nmon:12:end] .- xbar[nmon]
    end
    return xprime
end

"""
    function remove_seasonal(x,Ecycle,Fcycle,γ) natively
    `x` = long monthly timeseries, starting in Jan

    remove seasonal cycle of this timeseries
"""
function remove_seasonal(x,Ecycle,Fcycle)
    # remove seasonal cycle from 14-day averaged timeseries
    # solve for seasonal cycle parameters
    βcycle = Fcycle*x
    xseasonal = Ecycle*βcycle

    # remove it from total signal
    x′ = x - xseasonal
    return x′
end

"""
    function remove_seasonal(x,Ecycle,Fcycle,γ) 
    remove seasonal cycle of this timeseries on the native gcmgrid

# Arguments
-  `x::gcmarray{T,N,Array{T,2}}` = gcmarray of a long timeseries where time is dimension N
- `Ecycle::Matrix`: matrix that operates on seasonal cycle parameters and return a seasonal cycle timeseries
- `Fcycle::Matrix`: precomputed pseudo-inverse of `Ecycle`
- `γ::gcmgrid`: MITgcm grid
# Output
- `x′::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}`: potential density referenced to p₀ minus 1000 kg/m³
"""
function remove_seasonal(x::MeshArrays.gcmarray{T,N,Array{T,2}},Ecycle,Fcycle,γ)::MeshArrays.gcmarray{T,N,Array{T,2}} where T <: AbstractFloat where N

    # remove seasonal cycle from 14-day averaged timeseries
    # solve for seasonal cycle parameters

    # for gcmarray input:
    βcycle = matmul(Fcycle,x,γ)

    # reconstruct the full seasonal cycle.
    xseasonal = matmul(Ecycle,βcycle,γ)
    x′ = x - xseasonal
    return x′
end

"""
    function factors4regularpoles(γ)
    Get interpolation factors for regularpoles grid in one place
"""
function factors4regularpoles(γ)
    ϕGarc,ϕCarc = latgridArctic(γ)
    ϕGantarc,ϕCantarc = latgridAntarctic(γ) 
    ϕGreg,ϕCreg = latgridRegular(γ) 
    λG = -180.0:179.0
    λC = -179.5:179.5
    ϕG = vcat(ϕGantarc,ϕGreg,ϕGarc)
    ϕC = vcat(ϕCantarc,ϕCreg,ϕCarc)
    λarc = prereginterp(ϕCarc,λC,γ)
    λantarc = prereginterp(ϕCantarc,λC,γ)
    nx = length(λC)
    ny = length(ϕC)
    nyarc = length(ϕCarc)
    nyantarc = length(ϕCantarc)

    return λC,λG,ϕC,ϕG,nx,ny,nyarc,nyantarc,λarc,λantarc
end

"""
    function prereginterp(latgrid,longrid,γ)
    prepare for regular interpolation
    regular = onto rectangular 2d map
# Arguments
- `latgrid`: 1d array of latitude
- `longrid`: 1d array of longitude
- `γ`: GCM grid
# Output
- `f,i,j,w`: interpolation factors
"""
function prereginterp(latgrid,longrid,γ)
    Γ = GridLoad(γ;option="full")
    lon2d=[i for i=longrid, j=latgrid]
    lat2d=[j for i=longrid, j=latgrid]
    @time (f,i,j,w)=InterpolationFactors(Γ,vec(lon2d),vec(lat2d))
    λ=(lon=lon2d,lat=lat2d,f=f,i=i,j=j,w=w);
    return λ
end

"""
    function reginterp(fldin,nx,ny,f,i,j,w)
    regular interpolation with precomputed factors
    regular = onto rectangular 2d map
# Arguments
- `fldin`: gcmarray field of input
- `nx,ny`: x,y dimension of output regular field
- `λ=(f,i,j,w)`: precomputed interpolation factors
# Output
- `fldout`: interpolated regular 2d field
"""
function reginterp(fldin,nx,ny,λ)
    fldout = Interpolate(fldin,λ.f,λ.i,λ.j,λ.w)

    # think about moving this before interpolation
    # replace!(fldout,0.0 => NaN)
    fldout = transpose(reshape(fldout,nx,ny))
    
    return fldout
end

"""
    function trend_theta!(β,diagpathexp,path_out,tecco,γ,F)
    get linear trend of potential temperature
    from monthly-average ECCOv4r4 gcmarray fields
# Arguments
- `β`: gcmarray field of trends
- `diagpathexp`: where to find files
- `tecco`: time stamps of monthly output
- `γ`: GCM grid information
- `F`: linear estimator of trend
# Output
- `β`: updated trends
"""
function trend_theta!(β,diagpathexp,tecco,γ,F)
    # name of file inside diagspath
    # Look at /poseidon ... exps/run/data.diagnostics for this info.
    Troot = "state_3d_set1" # hardcoded for ECCOv4r4

    fill!(β,0.0f0) # initialize β
    nz = size(β,2)
    filelist = searchdir(diagpathexp,Troot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    nt = length(tecco)
    nt2 = length(datafilelist)
    if nt2 != nt
        error("incompatible t and data files")
    end

    global tt = 0
    for Tname in datafilelist
        tt += 1
        #println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)

        # read θ for timestep
        θ = γ.read(diagpathexp*Tname,MeshArray(γ,Float32,nz))

        # multiply by the correct part of F matrix
        # be sure to handle all points at once
        # add to existing solution for β
        # equal to matrix multiplication w/ columns (θ) times weighting F
        β += F[2,tt] * θ
    end
end

"""
    function trend_theta(diagpathexp,path_out,tecco,γ,F)
    get linear trend of potential temperature
    from monthly-average ECCOv4r4 gcmarray fields
# Arguments
- `diagpathexp`: where to find files
- `tecco`: time stamps of monthly output
- `γ`: GCM grid information
- `F`: linear estimator of trend
# Output
- `β`: updated trends
"""
function trend_theta(diagpathexp,tecco,γ,F)
    # name of file inside diagspath
    # Look at /poseidon ... exps/run/data.diagnostics for this info.
    Troot = "state_3d_set1" # hardcoded for ECCOv4r4

    # 50 shouldn't be hard coded
    β = MeshArray(γ,Float32,50) # some nans here
    fill!(β,0.0f0) # initialize β
    nz = size(β,2)
    filelist = searchdir(diagpathexp,Troot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    nt = length(tecco)
    nt2 = length(datafilelist)
    if nt2 != nt
        error("incompatible t and data files")
    end

    global tt = 0
    for Tname in datafilelist
        tt += 1
        #println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)

        # read θ for timestep
        θ = γ.read(diagpathexp*Tname,MeshArray(γ,Float32,nz))

        # multiply by the correct part of F matrix
        # be sure to handle all points at once
        # add to existing solution for β
        # equal to matrix multiplication w/ columns (θ) times weighting F
        β += F[2,tt] * θ
    end
    return β
end

function croplimitsLLC(ϕ,λ)
    # There is a regular grid inside the LLC grid.
    # Find the limiting indices of that regular grid.

    nx,ny = size(ϕ[1])
    # find edge of antarctic distortion
    jantarc = 1
    while sum(abs.(diff(vec(ϕ[1][:,jantarc])))) > 0.01
        jantarc += 1
    end
    ϕantarc= ϕ[1][1,jantarc]

    # println("Southern edge of regular grid found")
    # println("at latitude ",ϕantarc,"°N")

    # find edge of arctic distortion
    jarc = ny
    while sum(abs.(diff(vec(ϕ[1][:,jarc])))) > 0.01
        jarc -= 1
    end
    ϕarc = ϕ[1][1,jarc]
    #println("Northern edge of regular grid found")
    #println("at latitude ",ϕarc,"°N")
    return jarc, jantarc, ϕarc, ϕantarc
end

function regularlatgrid(γ)
    ϕ,λ = latlonG(γ)
    jarc,jantarc,ϕarc,ϕantarc = croplimitsLLC(ϕ,λ)

    # get grid spacing for polar regions.
    Δϕarc= ϕ[1][1,jarc] - ϕ[1][1,jarc-1]
    Δϕantarc= ϕ[1][1,jantarc+1] - ϕ[1][1,jantarc]

    # get number of regularly-spaced polar gridcells
    ymax = 90; ymin = -82; # trim reg grid at 80 S
    narc = Int(ceil((ymax-ϕarc)/Δϕarc))
    nantarc = Int(floor((ϕantarc-ymin)/Δϕantarc))

    ϕG = collect(reverse(range(ϕantarc,length=nantarc,step=-Δϕantarc)))
    append!(ϕG,ϕ[1][1,jantarc+1:jarc-1])
    append!(ϕG,range(ϕarc,length=narc,step=Δϕarc))

    # Don't put ϕC halfway between ϕG points.
    # It doesn't match the MITgcm LLC grid.

    ϕ,λ = latlonC(γ)
    jarc -= 1 # update for C grid, subtract one here
    ϕarc = ϕ[1][1,jarc]
    ϕantarc= ϕ[1][1,jantarc]
    ϕC = collect(reverse(range(ϕantarc,length=nantarc,step=-Δϕantarc)))
    append!(ϕC,ϕ[1][1,jantarc+1:jarc-1])
    append!(ϕC,range(ϕarc,length=narc,step=Δϕarc))

    return ϕG, ϕC
end

function latgridRegular(γ)
    ϕ,λ = latlonG(γ)
    jarc,jantarc,ϕarc,ϕantarc = croplimitsLLC(ϕ,λ)

    ϕG = ϕ[1][1,jantarc:jarc]

    ϕ,λ = latlonC(γ)
    jarc -= 1 # update for C grid, subtract one here
    ϕC = ϕ[1][1,jantarc:jarc]

    return ϕG, ϕC
end

function latgridArctic(γ)
    ϕ,λ = latlonG(γ)
    jarc,jantarc,ϕarc,ϕantarc = croplimitsLLC(ϕ,λ)

    # get grid spacing for polar regions.
    Δϕarc= ϕ[1][1,jarc] - ϕ[1][1,jarc-1]

    # get number of regularly-spaced polar gridcells
    ymax = 90 # trim reg grid at 80 S
    narc = Int(ceil((ymax-ϕarc)/Δϕarc))

    ϕG=collect(range(ϕarc,length=narc,step=Δϕarc))
    popfirst!(ϕG)

    # same thing for centered (tracer) grid
    ϕ,λ = latlonC(γ)
    jarc -= 1 # update for C grid, subtract one here
    ϕarc = ϕ[1][1,jarc]
    ϕC=collect(range(ϕarc,length=narc,step=Δϕarc))
    popfirst!(ϕC)

    return ϕG, ϕC
end

function latgridAntarctic(γ)
    ϕ,λ = latlonG(γ)
    jarc,jantarc,ϕarc,ϕantarc = croplimitsLLC(ϕ,λ)

    # get grid spacing for polar regions.
    Δϕantarc= ϕ[1][1,jantarc+1] - ϕ[1][1,jantarc]

    # get number of regularly-spaced polar gridcells
    ymin = -82; # trim reg grid at 80 S
    nantarc = Int(floor((ϕantarc-ymin)/Δϕantarc))

    ϕG = collect(reverse(range(ϕantarc,length=nantarc,step=-Δϕantarc)))
    # eliminate northernmost grid point.
    # it's inside regular subgrid.
    pop!(ϕG)

    ϕ,λ = latlonC(γ)
    ϕantarc= ϕ[1][1,jantarc]
    ϕC = collect(reverse(range(ϕantarc,length=nantarc,step=-Δϕantarc)))
    pop!(ϕC)

    return ϕG, ϕC
end

function LLCcropC(gcmfield,γ)
    # There is a regular grid inside the LLC grid.
    # Subsample/crop just those points. Put them together correctly.
    ϕ,λ = latlonG(γ)
    jarc,jantarc,ϕarc,ϕantarc = croplimitsLLC(ϕ,λ)
    jarc -= 1 # update for C grid, subtract one here

    regfield = gcmfield[1][:,jantarc:jarc]
    regfield = vcat(regfield,gcmfield[2][:,jantarc:jarc])
    for i = 4:5
        tmp = reverse(transpose(gcmfield[i]),dims=2)
        regfield = vcat(regfield,tmp[:,jantarc:jarc])
    end

    # handle longitudinal wraparound by hand
    wrapval = 218
    xwrap = vcat(wrapval+1:size(regfield,1),1:wrapval)
    regfield = regfield[xwrap,:]
    return regfield
end

function LLCcropG(gcmfield,γ)
    # There is a regular grid inside the LLC grid.
    # Subsample/crop just those points. Put them together correctly.
    ϕ,λ = latlonG(γ)
    jarc,jantarc,ϕarc,ϕantarc = croplimitsLLC(ϕ,λ)

    regfield = gcmfield[1][:,jantarc:jarc]
    regfield = vcat(regfield,gcmfield[2][:,jantarc:jarc])
    for i = 4:5
        tmp = reverse(transpose(gcmfield[i]),dims=2)
        regfield = vcat(regfield,tmp[:,jantarc:jarc])
    end

    # handle longitudinal wraparound by hand
    wrapval = 218
    xwrap = vcat(wrapval+1:size(regfield,1),1:wrapval)
    regfield = regfield[xwrap,:]
    return regfield
end

"""
    function timestamp_monthly_v4r4(t)
    print year and month given time index
    assuming using ECCOv4r4
"""
function timestamp_monthly_v4r4(t)
    tstart = 1992 + 1/24
    tend = 2018
    tecco = range(tstart,step=1/12,stop=2018)
    year = Int(floor(tecco[t]))
    month = ((t-1)%12)+1
    println("year ",year," month ",month)
    return year,month
end

notnanorzero(z) = !iszero(z) && !isnan(z)

function mdsio2dict(pathin,filein,γ)

    metafile = filein*".meta"
    datafile = filein*".data"
    meta = read_meta(pathin,metafile);

    state =  read_bin(pathin*datafile,missing,missing,Float32,γ)
    ndimz = meta[1].nDims

    if ndimz == 3         # 3d
        nz = meta[1].dimList[ndimz,1]
        vars = Dict{String,MeshArrays.gcmarray{Float32,2,Array{Float32,2}}}()
    elseif ndimz == 2
        nz = 1
        vars = Dict{String,MeshArrays.gcmarray{Float32,1,Array{Float32,2}}}()
    end

    zlo = 1
    for fldname in meta[1].fldList

        # use a dictionary
        zhi = zlo + nz -1
        push!(vars,fldname => state[:,zlo:zhi])
        zlo = zhi + 1
    end
    return vars
end

"""
   function write_vars
   This version writes mdsio output on the native grid
   Missing: write the accompanying meta file
"""
function write_vars(vars::Dict{String,MeshArrays.gcmarray{Float32,2,Matrix{Float32}}},fileprefix::String,filesuffix::String)
    for (fldname,fldvals) in vars
        filename = fileprefix*fldname*filesuffix
        write(filename,fldvals)
    end
end

"""
   function ncwritefromtemplate
   This version writes NetCDF output on the regularpoles grid
"""
function ncwritefromtemplate(vars::Dict{String,Array{Float64,3}},fileprefix::String,filesuffixold::String,filesuffixnew::String,sig1grid)

    # 2 WAYS TO GO:
    #1. DO LIKE REGULARPOLES AND TRY TO OUTPUT.
    #2. READ NC, SAVE ATTRIBUTES, ADD NEW SIGMA1 FIELD.
    # choose method 2 here

    for (fldname,fldval) in vars
        println(fldname)
        println(fldval)
        
        # get a filename with _on_sigma1
        if fldname == "p"
            fileold = fileprefix*"THETA/THETA"*filesuffixold
        else
            fileold = fileprefix*fldname*"/"*fldname*filesuffixold
        end

        filenewpath = fileprefix*fldname*"/"
        filenew = filenewpath*fldname*filesuffixnew
        # make sure path exists, was a problem for pressure
        !isdir(filenewpath) ? mkpath(filenewpath) : nothing 
        
        # recover information from template/existing file
        lon = ncread(fileold,"lon")
        lonunits = ncgetatt(fileold, "lon", "units")
        lonlongname = ncgetatt(fileold, "lon", "longname")
        lonatts = Dict("longname" => lonlongname, "units" => lonunits)

        lat = ncread(fileold,"lat")
        latunits = ncgetatt(fileold, "lat", "units")
        latlongname = ncgetatt(fileold, "lat", "longname")
        latatts = Dict("longname" => latlongname, "units" => latunits)

        sigmaatts = Dict("longname" => "Sigma-1", "units" => "kg/m^3 - 1000")

        if fldname == "p"
            varatts = Dict("longname" => "pressure", "units" => "dbar")
        else
            vunits = ncgetatt(fileold, fldname, "units")
            vlongname = ncgetatt(fileold, fldname, "longname")
            varatts = Dict("longname" => vlongname, "units" => vunits)
        end

        isfile(filenew) && rm(filenew)
        nccreate(
            filenew,
            fldname,
            "lon",
            lon,
            lonatts,
            "lat",
            lat,
            latatts,
            "sigma1",
            sig1grid,
            sigmaatts,
            atts = varatts,
        )
            ncwrite(fldval, filenew, fldname)
    end
end

"""
    function mdsio2sigma1
    Take variables in a filelist, read, map to sigma1, write to file.
"""
function mdsio2sigma1(pathin,pathout,fileroots,γ,pstdz,sig1grid,splorder) 
    # Read All Variables And Puts Them Into "Vars" Dictionary

    # ideally would be more generic and Float64 would be covered.
    vars = Dict{String,MeshArrays.gcmarray{Float32,2,Matrix{Float32}}}()
    for fileroot in fileroots
        merge!(vars,mdsio2dict(pathin,fileroot,γ))
    end

    # check for fields on the velocity (staggered) grid
    # and rotate them if necessary
    Γ = GridLoad(γ;option="full")
    rotate_velocity!(vars,Γ)

    # solve for sigma1 on depth levels.
    @time varsσ = vars2sigma1(vars,pstdz,sig1grid,γ,splorder)

    fileprefix = pathout
    # use first filename to get timestamp
    filesuffix = "_on_sigma1"*fileroots[1][14:end]*".data"
    write_vars(varsσ,fileprefix,filesuffix)

    return varsσ
end

"""
    function netcdf2sigma1
    Take variables in a NetCDF filelist, read, map to sigma1, write to file.
"""
function netcdf2sigma1(pathin,pathout,ncfilenames,γ,pstdz,sig1grid,splorder)
    # Read All Variables And Puts Them Into "Vars" Dictionary
    println(ncfilenames)
    vars = Dict{String,Array{Float64,3}}()
    for (ncvarname,ncfilename) in ncfilenames
        merge!(vars,netcdf2dict(ncfilename,ncvarname))
    end

    # solve for sigma1 on depth levels.
    @time varsσ = vars2sigma1(vars,pstdz,sig1grid,splorder)

    fileprefix = pathout
    # use first filename to get timestamp
    filesuffixold = ncfilenames["THETA"][end-10:end]
    #    filesuffixnew = "_on_sigma1"*ncfilenames["THETA"][end-10:end]
    filesuffixnew = "_on_sigma1"*filesuffixold
    ncwritefromtemplate(varsσ,fileprefix,filesuffixold,filesuffixnew,sig1grid)
end

"""
netcdf2dict array output (?)
"""
function netcdf2dict(ncfilename,ncvarname)
    # 3d
    println(ncfilename,ncvarname)
    state = ncread(ncfilename,ncvarname);
    if size(state,3) == 1
        vars  = Dict(ncvarname => state[:,:,1])
    else
        vars  = Dict(ncvarname => state)
    end
    return vars
end

"""
netcdf2dict gcmarray output
"""
function netcdf2dict(ncfilename::String,ncvarname::String,γ::gcmgrid)
    # 3d
    println(ncfilename,ncvarname)
    state = read_netcdf(ncfilename,ncvarname,γ);
    if size(state,2) == 1
        vars  = Dict(ncvarname => state[:,1])
    else
        vars  = Dict(ncvarname => state)
    end

    return vars
end

"""
function replace!(f,a::MeshArrays.gcmarray{T,N,Array{T,2}}}
# Arguments
- `f::function` = replace a with f(a)
- `a::MeshArrays.gcmarray{T,N,Array{T,2}}`: gcmarray with variable type and time-dimension
"""
function replace!(f,a::MeshArrays.gcmarray{T,N,Array{T,2}}) where T<:AbstractFloat where N

    #nf = size(a,1)
    for ff = eachindex(a)
        Base.replace!(f,a[ff])
    end
    return a
end

"""
function replace!(a::MeshArrays.gcmarray{T,N,Array{T,2}}},b::Pair) where T
# Arguments
- `a::MeshArrays.gcmarray{T,N,Array{T,2}}`: gcmarray with variable type and time-dimension
- `b::Pair`: replace `a` elements of pair 1 with pair 2
"""
function replace!(a::MeshArrays.gcmarray{T,N,Array{T,2}},b::Pair) where T<:AbstractFloat where N

    # nf = size(a,1)
    for ff = eachindex(a)
        Base.replace!(a[ff],b)
    end
    return a
end

function netcdf2regularpoles(ncfilename,ncvarname,γ,nx,ny,nyarc,λarc,nyantarc,λantarc)

    vars = netcdf2dict(ncfilename,ncvarname,γ)
    varsregpoles = vars2regularpoles(vars,γ,nx,ny,nyarc,λarc,nyantarc,λantarc)

end

function mdsio2regularpoles(pathin,filein,γ,nx,ny,nyarc,λarc,nyantarc,λantarc)

    vars = mdsio2dict(pathin,filein,γ)
    Γ = GridLoad(γ;option="full")
    rotate_velocity!(vars,Γ)
    varsregpoles = vars2regularpoles(vars,γ,nx,ny,nyarc,λarc,nyantarc,λantarc)

end

function vars2regularpoles(vars::Dict{String,MeshArrays.gcmarray{T,2,Matrix{T}}},γ,nx,ny,nyarc,λarc,nyantarc,λantarc) where T<:AbstractFloat

    varsregpoles = Dict{String,Array{T,3}}()
    NaNT = zero(T)/zero(T)
    for (varname, varvals) in vars

        # remove contamination from land
        
        replace!(varvals, 0.0 => NaNT)
        nz = size(varvals,2)

        #pre-allocate dict
        varsregpoles[varname] = fill(NaNT,(nx,ny,nz))

        for zz = 1:nz
            # get regular grid by cropping
            varsregpoles[varname][:,:,zz]=var2regularpoles(varvals[:,zz],γ,nx,ny,nyarc,λarc,nyantarc,λantarc)
        end
    end
    return varsregpoles
end

"""
function vars2regularpoles(vars::Dict{String,MeshArrays.gcmarray{T,1,Matrix{T}}},γ,nx,ny,nyarc,λarc,nyantarc,λantarc)
     variables interpolated onto regularpoles grid
"""
function vars2regularpoles(vars::Dict{String,MeshArrays.gcmarray{T,1,Matrix{T}}},γ,nx,ny,nyarc,λarc,nyantarc,λantarc) where T<:AbstractFloat

    varsregpoles = Dict{String,Matrix{T}}()
    for (varname, varvals) in vars

        varsregpoles[varname] = var2regularpoles(varvals,γ,nx,ny,nyarc,λarc,nyantarc,λantarc)

    end
    return varsregpoles
end

"""
   var2regularpoles
   Take one gcmarray in memory, put on regularpoles grid
"""
function var2regularpoles(var,γ,nx,ny,nyarc,λarc,nyantarc,λantarc)

    T = eltype(var)
    NaNT = zero(T)/zero(T)
    # remove contamination from land
    replace!(var, 0.0 => NaNT)

    #pre-allocate output
    varregpoles = fill(NaNT,(nx,ny))

    # get regular grid by cropping
    θcrop =  LLCcropC(var,γ)

    # interpolate to "regularpoles"
    θarc = reginterp(var,nx,nyarc,λarc)
    θantarc = reginterp(var,nx,nyantarc,λantarc)
    varregpoles=hcat(θantarc',θcrop,θarc')

    return varregpoles
    
end

"""
function writeregularpoles(vars,γ,pathout,filesuffix,filelog,λC,lonatts,ϕC,latatts,z,depthatts)
"""
function writeregularpoles(vars::Dict{String,Array{Float32,3}},γ,pathout,filesuffix,filelog,λC,lonatts,ϕC,latatts,z,depthatts)

    for (varname,varvals) in vars
        println(varname)
        # Some names are a little inconsistent
        # Useful to pull meta-data from these locations
        if varname == "EVELMASS"
            field = "UE_VEL_C"
        elseif varname == "NVELMASS"
            field = "VN_VEL_C"
        elseif varname == "WVELMASS"
            field = "WVEL"
        elseif varname == "oceTAUN"
            field = "oceTAUY"
        elseif varname == "oceTAUE"
            field = "oceTAUX"
        elseif varname == "GM_PsiN"
            field = "GM_PsiY"
        elseif varname == "GM_PsiE"
            field = "GM_PsiX"
        else
            field = varname
        end

        if varname == "p"
            fieldDict = Dict("fldname" => "p","title" => "standard pressure", "units" => "dbar", "levs" => length(z))
        else
            fieldDict = read_available_diagnostics(field,filename=filelog)
            
        end

        if depthatts["longname"] == "Sigma-1"
            depthname = "sigma-1"
        else
            depthname = "depth"
        end
        
        # make a directory for this output
        pathoutdir = pathout*varname*"/"
        !isdir(pathoutdir) ? mkdir(pathoutdir) : nothing;

        # get filename for this month.
        fileout = pathoutdir*varname*filesuffix
        println(fileout)

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
            depthname,
            z,
            depthatts,
            atts = varatts,
        )
        ncwrite(varvals, fileout, varname)
    end
end

"""
function writeregularpoles(vars,pathout,filesuffix,filelog,λC,lonatts,ϕC,latatts)
"""
function writeregularpoles(vars::Dict{String,Array{Float32,2}},γ,pathout,filesuffix,filelog,λC,lonatts,ϕC,latatts,z,depthatts)

    for (varname,varvals) in vars
        if varname[1] == 'E'
            field = "UE_VEL_C"
        elseif varname[1] == 'N'
            field = "VN_VEL_C"
        elseif varname[1] == 'W'
            field = "WVEL"
        elseif varname[end] == 'N'
            field = "oceTAUY"
        elseif varname[end] == 'E'
            field = "oceTAUX"
        else
            field = varname
        end
        fieldDict = read_available_diagnostics(field,filename=filelog)

        # make a directory for this output
        pathoutdir = pathout*varname*"/"
        !isdir(pathoutdir) ? mkdir(pathoutdir) : nothing;

        # get filename for this month.
        fileout = pathoutdir*varname*filesuffix
        println(fileout)

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
            atts = varatts,
        )
        ncwrite(varvals, fileout, varname)
    end
end

"""
function writeregularpoles(vars,γ,pathout,filesuffix,filelog,λC,lonatts,ϕC,latatts,z,depthatts)
"""
function writeregularpoles(vars::Dict{String,Array{Float64,3}},γ,pathout,filesuffix,filelog,λC,lonatts,ϕC,latatts,z,depthatts)

    for (fldname,varvals) in vars
        println(fldname)
        fieldDict = read_available_diagnostics(fldname,filename=filelog)

        # make a directory for this output
        pathoutdir = pathout*fieldDict["fldname"]*"/"
        !isdir(pathoutdir) && mkdir(pathoutdir)

        # get filename for this month.
        fileout = pathoutdir*fieldDict["fldname"]*filesuffix
        println(fileout)

        # save in a NetCDF file with info from fieldDict
        varatts = Dict("longname" => fieldDict["title"], "units" => fieldDict["units"])

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
            "depth",
            z,
            depthatts,
            atts = varatts,
        )
        ncwrite(varvals, fileout, fieldDict["fldname"])
    end
end

"""
function writeregularpoles(vars,pathout,filesuffix,filelog,λC,lonatts,ϕC,latatts)
"""
function writeregularpoles(vars::Dict{String,Array{Float64,2}},γ,pathout,filesuffix,filelog,λC,lonatts,ϕC,latatts,z,depthatts)

    for (fldname,varvals) in vars
        fieldDict = read_available_diagnostics(fldname,filename=filelog)

        # make a directory for this output
        pathoutdir = pathout*fieldDict["fldname"]*"/"
        !isdir(pathoutdir) ? mkdir(pathoutdir) : nothing;

        # get filename for this month.
        fileout = pathoutdir*fieldDict["fldname"]*filesuffix
        println(fileout)

        # save in a NetCDF file with info from fieldDict
        varatts = Dict("longname" => fieldDict["title"], "units" => fieldDict["units"])

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
        ncwrite(varvals, fileout, fieldDict["fldname"])
    end
end

"""
    function regional_mask(latpt,lonpt,latrect,lonrect,dlat,dlon)
# Arguments
- `latpt`: latitude grid
- `lonpt`: longitude grid
- `latrect`: tuple of latitude range of central rectangle
- `lonrect`: tuple of longitude range of central rectangle
- `dlon`: width in degrees longitude of East/West fade zones
- 'dlat': width in degrees latitude of North/South fade zones
# Output
- 'mask': space and time field of surface forcing, value of zero inside
designated lat/lon rectangle and fading to 1 outside sponge zone on each edge. This is
because this field ends up being SUBTRACTED from the total forcing
"""
function regional_mask(latpt,lonpt,latrect,lonrect,dlat,dlon)

    # North
    # make meridional plane with
    # gradient determined by sponge layers
    mask = 1 .+ (latrect[2] .- latpt) / dlat
    
    # only keep values between 0 and 1
    zero2one!(mask)

    # south
    mask2 = 1 .+ (latpt .- latrect[1]) / dlat
    zero2one!(mask2)
    mask *= mask2

    # west
    mask2 = 1 .+ (lonrect[2] .- lonpt) / dlon
    zero2one!(mask2)
    mask *= mask2

    # east
    mask2 = 1 .+ (lonpt .- lonrect[1]) / dlon
    zero2one!(mask2)
    mask *= mask2

    return mask
    
end

"""
    function wrapdist
    Longitudinal distance in degrees
    Passing the date line may be the shortest distance.
"""
function wrapdist(lonpt,lonmid)

    dist = lonpt .- lonmid

    # greater than 180, then subtract 360
    # less than -180, then add 360
    wraphi = (x -> x > 180 ? x -= 360 : x)
    ECCOtour.replace!(wraphi,dist)
    wraplo = (x -> x < -180 ? x += 360 : x)
    ECCOtour.replace!(wraplo,dist)
    return dist
end

"""
    function centerlon!
    Make a longitudinal coordinate system
    centered at lonmid.
    Makes handling wraparound easy later on.
"""
function centerlon!(lonpt,lonmid)

    # greater than 180, then subtract 360
    # less than -180, then add 360
    wraphi = (x -> x > lonmid+180 ? x -= 360 : x)
    ECCOtour.replace!(wraphi,lonpt)
    wraplo = (x -> x < lonmid-180 ? x += 360 : x)
    ECCOtour.replace!(wraplo,lonpt)
    
end

function zero2one!(v)
    ECCOtour.replace!(x -> x > 1.0 ? 1.0 : x, v)
    ECCOtour.replace!(x -> x < 0.0 ? 0.0 : x, v)
end

"""
    function regional_mask(field,mask)
# Arguments
- `field`: input field that is mutated/modified
- `mask`: multiplicative factor, same spatial grid as field
"""
function apply_regional_mask!(field,mask)

    #once you have the mask, apply it to the input forcing field
    #regionalforcingfield = maskmultiplier.*forcingfield
    # Do this in a separate function that can be called every timestep.

    for i = 1:size(field,2)
        println(i)
        field[:,i] = field[:,i].*mask
    end
    
end

"""
    read_netcdf(fileName,fldName,mygrid)

Read model output from NetCDF files that are global and convert to MeshArray instance.
```
"""
function read_netcdf(fileName::String,fldName::String,mygrid::gcmgrid)
    if (mygrid.class!="LatLonCap")||(mygrid.ioSize!=[90 1170])
        error("non-llc90 cases not implemented yet")
    end

    fileIn= fileName #@sprintf("%s.%04d.nc",fileName,1)
    #    println(fileIn)
    xglobal = ncread(fileIn,fldName)
    s = [size(xglobal,i) for i in 1:ndims(xglobal)]
    n=length(size(xglobal))
    
    # manually, drop dimension 5.
    if n == 5
        xglobal = xglobal[:,:,:,:,1]
        n = ndims(xglobal)
    end
    f0=Array{Float64}(undef,90,0,s[4])
    f00=Array{Float64}(undef,0,90,s[4])
    f=[f0,f0,f0,f00,f00]

    # take the 13 tiles and put them into 5 faces
    #xarray = Array{Float64}(undef,90,1170,50)
    f[1] = cat(xglobal[:,:,1,:],xglobal[:,:,2,:],xglobal[:,:,3,:];dims=2)
    f[2] = cat(xglobal[:,:,4,:],xglobal[:,:,5,:],xglobal[:,:,6,:];dims=2)
    f[3] = xglobal[:,:,7,:];
    f[4] = cat(xglobal[:,:,8,:],xglobal[:,:,9,:],xglobal[:,:,10,:];dims=1)
    f[5] = cat(xglobal[:,:,11,:],xglobal[:,:,12,:],xglobal[:,:,13,:];dims=1)

    fld=MeshArray(mygrid,f)
    return fld
end


"""
    extract_timeseries(froot,years,γ,xval,yval,fval)
# Arguments
- `froot::String`: filename root
- `years::StepRange`: iterator for multiple files
- `γ::gcmarray`: GCM grid from MeshArrays.jl
- `xval::Integer`: x grid index
- `yval::Integer`: y grid index
- `fval::Integer`: face index
# Output    
- `tseries`: timeseries
- `nseries`: nseries = number of timeseries elements in each year
```
"""
function extract_timeseries(froot,years,γ,xval,yval,fval)
    tseries = []
    nseries = []
    nyr = length(years)
    for tt = 1:nyr
        fname = froot*string(years[tt])
        println(fname)
        field = read_bin(fname,Float32,γ)
        series = extract_timeseries(field,xval,yval,fval)
        push!(nseries,length(series))
        append!(tseries,series[:])
    end
    return tseries,nseries
end

################################################################################################
"""
    extract_timeseries(flux,xval,yval,fval)
# Arguments
- `flux::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}`: input array of arrays
- `years::StepRange`: iterator for multiple files
- `xval::Integer`: x grid index
- `yval::Integer`: y grid index
- `fval::Integer`: face index
# Output    
- `series`: timeseries
"""
function extract_timeseries(flux::MeshArrays.gcmarray{Float32,2,Array{Float32,2}},xval,yval,fval)
    nt = size(flux,2)
    series = Float32[]
    
    for tval ∈ 1:nt
        tmp = flux[fval,tval]
        push!(series,tmp[xval,yval])
    end
    return series
end

###############################################################################################
"""
    extract_timeseries(flux,xval,yval,fval)
# Arguments
- `flux::MeshArrays.gcmarray{Float64,2,Array{Float64,2}}`: input array of arrays
- `years::StepRange`: iterator for multiple files
- `xval::Integer`: x grid index
- `yval::Integer`: y grid index
- `fval::Integer`: face index
# Output    
- `series`: timeseries
"""
function extract_timeseries(flux::MeshArrays.gcmarray{Float64,2,Array{Float64,2}},xval,yval,fval)
    nt = size(flux,2)
    series = Float64[]
    
    for tval ∈ 1:nt
        tmp = flux[fval,tval]
        push!(series,tmp[xval,yval])
    end
    return series
end

#################################################################################################
# tell julia how to do matrix multiplication
function matmul(M::Array{Float64,2},flux::MeshArrays.gcmarray{Float64,2,Array{Float64,2}},γ) ::MeshArrays.gcmarray{Float64,2,Array{Float64,2}}

    nM = size(M,1)  # matrix size M
    nN = size(flux,2) # matrix size N
    nQ  = size(flux,1) # repeat matmul Q times
    
    # initialize product
    product = 0.0 .* MeshArray(γ,Float64,nM)
    
    for mm = 1:nM
        for qq = 1:nQ # inner product over faces
            for nn = 1:nN # inner product over time
                product[qq,mm] += flux[qq,nn] * M[mm,nn]
                #product[ff,pp] +=  transpose(F[pp,tt]'*flux[ff,tt]')
            end
        end
    end
    return product
end

function matmul(M::Array{Float32,2},flux::MeshArrays.gcmarray{Float32,2,Array{Float32,2}},γ::gcmgrid) ::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}

    nM = size(M,1)  # matrix size M
    nN = size(flux,2) # matrix size N
    nQ  = size(flux,1) # repeat matmul Q times
    
    # initialize product
    #product = 0.0f0 .* MeshArray(γ,Float32,nM)
    println("sub-optimal initialization")
    product = MeshArray(γ,Float32,nM) # some nans here
    tmp1=zeros(Float32,Tuple(γ.ioSize))
    for tt= 1:nM
        # initialize a sub-optimal way
        product[:,tt]=γ.read(tmp1,MeshArray(γ,Float32))
    end

    for mm = 1:nM
        for qq = 1:nQ # inner product over faces
            for nn = 1:nN # inner product over time
                product[qq,mm] += flux[qq,nn] * M[mm,nn]
                #product[ff,pp] +=  transpose(F[pp,tt]'*flux[ff,tt]')
            end
        end
    end
    return product
end

function matmul(M::Array{Float32,1},flux::MeshArrays.gcmarray{Float32,1,Array{Float32,2}},γ::gcmgrid) ::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}

    nM = size(M,1)  # matrix size M
    nQ  = size(flux,1) # repeat matmul Q times
    
    # initialize product
    #product = 0.0f0 .* MeshArray(γ,Float32,nM)
    println("sub-optimal initialization")
    product = MeshArray(γ,Float32,nM) # some nans here
    tmp1=zeros(Float32,Tuple(γ.ioSize))
    for tt= 1:nM
        # initialize a sub-optimal way
        product[:,tt]=γ.read(tmp1,MeshArray(γ,Float32))
    end

    for mm = 1:nM
        for qq = 1:nQ # inner product over faces
            product[qq,mm] += flux[qq] * M[mm]
                #product[ff,pp] +=  transpose(F[pp,tt]'*flux[ff,tt]')
        end
    end
    return product
end

function matmul(M::Array{Float64,1},flux::MeshArrays.gcmarray{Float64,1,Array{Float64,2}},γ) ::MeshArrays.gcmarray{Float64,2,Array{Float64,2}}

    nM = size(M,1)  # matrix size M
    nQ  = size(flux,1) # repeat matmul Q times
    
    # initialize product
    #product = 0.0 .* MeshArray(γ,Float64,nM)
    println("sub-optimal initialization")
    product = MeshArray(γ,Float64,nM) # some nans here
    tmp1=zeros(Float64,Tuple(γ.ioSize))
    for tt= 1:nM
        # initialize a sub-optimal way
        product[:,tt]=γ.read(tmp1,MeshArray(γ,Float64))
    end
    
    for mm = 1:nM
        for qq = 1:nQ # inner product over faces
            product[qq,mm] += flux[qq] * M[mm]
                #product[ff,pp] +=  transpose(F[pp,tt]'*flux[ff,tt]')
        end
    end
    return product
end

"""
    function nancount

    Count number of NaN's in gcmgrid type 
# Input
- `field::MeshArrays.gcmarray{T,N,Matrix{T}}`: input of gcmarray type, can have multiple records
# Output
- `nannum::Int`: number of NaNs
"""
function nancount(field::MeshArrays.gcmarray{T, N, Matrix{T}}) where T<:AbstractFloat where N

    nannum = zeros(Int64,size(field))
    for i1 = eachindex(field)
        nannum[i1] = sum(isnan,field[i1])
    end
    return nannum
end

"""
    function maximum

    Compute maximum value of gcmgrid type 
# Input
- `x::MeshArrays.gcmarray{T,N,Array{T,2}}`: input of gcmarray type and filter out `dryval`s
- `dryval::T`: land value to be eliminated in calculation
# Output
- `xmax::T`: maximum value of 2D field
"""
function maximum(x::MeshArrays.gcmarray{T,N,Array{T,2}},dryval::T) where T<:AbstractFloat where N

    # problem NaN == NaN is false
    # if isnan(z)
    #     isdry = isnan
    # else
    #     isdry(z) = (z == dryval)
    # end
    isnan(dryval) ? isdry = isnan : isdry(z) = (z == dryval)
    
    #  vector list of non-dry (wet) elements
    xcount = [sum(count(!isdry,x[i])) for i in eachindex(x)]

    if sum(xcount) > 0
        xmax = maximum([maximum(filter(!isdry,x[i])) for i in eachindex(x) if xcount[i] > 0])
    else
        xmax = convert(T,NaN)
    end
    return xmax
end

"""
    function minimum

    Compute minimum value of gcmgrid type 
    and filter out `dryval`s
# Input
- `x::MeshArrays.gcmarray{T,N,Array{T,2}}`: input of gcmarray type and filter out `dryval`s
- `dryval::T`: land value to be eliminated in calculatio# Output
- `xmin::T`: minimum value of 2D field
"""
function minimum(x::MeshArrays.gcmarray{T,N,Array{T,2}},dryval::T) where T<:AbstractFloat where N
    xmin = -maximum(-x,dryval)
    return xmin
end

"""
    function mean

    Compute mean of gcmgrid type using function calls, eliminate redundancy
    Eliminate all values = dryval
# Input
- `x::MeshArrays.gcmarray{T,N,Array{T,2}}`: input of gcmarray type
- `dryval::T`: land value (doesn't work for NaN32)
# Output
- `xbar::T`: mean value (unweighted)
"""
function mean(x::MeshArrays.gcmarray{T,N,Array{T,2}},dryval::T)::T where T<:AbstractFloat where N

    #isdry(z) = (z == dryval)
    isnan(dryval) ? isdry = isnan : isdry(z) = (z == dryval)

    #  vector list of nonzero elements
    xcount = [sum(count(!isdry,x[i])) for i in eachindex(x)]
    if sum(xcount) > 0
        # don't assume 0 on land
        xsum = sum([sum(filter(!isdry,x[i])) for i in eachindex(x) if xcount[i] > 0])
        xbar = xsum/sum(xcount)
    else
        #xbar = NaN
        xbar = convert(T,NaN)
    end
    return xbar
end

"""
    function mean
    Compute area-weighted mean of gcmgrid type using function calls, eliminate redundancy
    Area weighting = area
    Eliminate all values = dryval

# Input
- `x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}`: input of gcmarray type
- `weight::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}`: weighting variable of gcmarray type
- `dryval::Float32`: land value (doesn't work for NaN32)
# Output
- `xbar::Float32`: mean value (unweighted)
"""
function mean(x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}},weight::MeshArrays.gcmarray{Float64,1,Array{Float64,2}},dryval::Float64)

    isnan(dryval) ? isdry = isnan : isdry(z) = (z == dryval)
    #isdry(z) = (z == dryval)
    
    #  vector list of nonzero elements
    xcount = [sum(count(!isdry,x[i])) for i in eachindex(x)]
    if sum(xcount) > 0
        # don't assume 0 on land
        xsum = sum([sum(filter(!isdry,x[i].*weight[i])) for i in eachindex(x) if xcount[i] > 0])
        xdenom = sum([sum(filter(!isdry,weight[i])) for i in eachindex(weight) if xcount[i] > 0])
        xbar = xsum/xdenom
    else
        xbar = NaN32
    end
    return xbar
end

"""
    function mean
    Compute area-weighted mean of gcmgrid type using filtered with function isgood
    Area weighting = area
    Eliminate all values = isgood(x) -> true

# Input
- `x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}`: input of gcmarray type
- `weight::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}`: weighting variable of gcmarray type
- `isgood::Function`: returns true is a value to be used in the mean
# Output
- `xbar::Float32`: mean value (weighted and filtered)
"""
function mean(x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}},weight::MeshArrays.gcmarray{Float64,1,Array{Float64,2}},isgood)

    #  vector list of nonzero elements
    xcount = [sum(count(isgood,x[i])) for i in eachindex(x)]
    if sum(xcount) > 0
        # don't assume 0 on land
        xsum = sum([sum(filter(isgood,x[i].*weight[i])) for i in eachindex(x) if xcount[i] > 0])
        xdenom = sum([sum(filter(isgood,weight[i])) for i in eachindex(weight) if xcount[i] > 0])
        xbar = xsum/xdenom
    else
        xbar = NaN32
    end
    return xbar
end

"""
    function std
    Compute standard deviation of gcmgrid type using function calls, eliminate redundancy
    Eliminate all values = dryval

# Input
- `x::MeshArrays.gcmarray{T,N,Array{T,2}}`: input of gcmarray type
- `xbar::T`: mean value
- `dryval::T`: land value 
# Output
- `σx::T`: standard deviation 
"""
function std(x::MeshArrays.gcmarray{T,N,Array{T,2}},xbar::T,dryval::T) where T<:AbstractFloat where N

    isnan(dryval) ? isdry = isnan : isdry(z) = (z == dryval)

    # x prime = fluctuation
    x′ = x .- xbar
    
    #  vector list of nonzero elements
    xcount = [sum(count(!isdry,x′[i])) for i in eachindex(x′)]
    if sum(xcount) > 0
        # don't assume 0 on land
        x²sum = sum([sum(filter(!isdry,x′[i]).^2) for i in eachindex(x) if xcount[i] > 0])
        σx = sqrt(x²sum/(sum(xcount)-1))
    else
        xbar = convert(T,NaN)
    end
    return σx
end

function rotate_velocity!(vars,Γ)

    velvars = (("UVELMASS","VVELMASS"),("GM_PsiX","GM_PsiY"),("oceTAUX","oceTAUY"))
    velchange = Dict("UVELMASS" => "EVELMASS",
                     "VVELMASS" => "NVELMASS",
                     "GM_PsiX" => "GM_PsiE",
                     "GM_PsiY" => "GM_PsiN",
                     "oceTAUX" => "oceTAUE",
                     "oceTAUY" => "oceTAUN")

    for (kk,vv) in velvars
        if haskey(vars,kk) && haskey(vars,vv)
            evel = similar(vars[kk])
            nvel = similar(vars[vv])
            uC = similar(vars[kk][:,1]) # one level
            vC = similar(vars[vv][:,1]) # one level
            for zz = 1: size(vars[kk],2)

                #interpolate velocity to center of C-grid
                velocity2center!(uC,vC,vars[kk][:,zz],vars[vv][:,zz],Γ)
                
                evel[:,zz],nvel[:,zz] = rotate_uv(uC,vC,Γ);
            end
            push!(vars,velchange[kk] => evel)
            push!(vars,velchange[vv] => nvel)
            delete!(vars,kk)
            delete!(vars,vv)
        end
    end
end

"""
    velocity2center!(uC,vC,u,v,G)
  
    From Gael Forget, JuliaClimateNotebooks/Transport
    #1. Convert to `Sv` units and mask out land
    2. Interpolate `x/y` transport to grid cell center
"""
function velocity2center!(uC,vC,u,v,G)

    #G = GridLoad(γ;option="full")
    #u[findall(G.hFacW[:,1].==0)].=NaN
    #v[findall(G.hFacS[:,1].==0)].=NaN;

    # replace slow findall with this kludge
    # kludge makes routine 2x faster
    # tmp variables below may also be slow (not views)
    # warning: may introduce bugs
    T = eltype(u)
    replace!(u,zero(T) => zero(T)/zero(T)) # NaN
    replace!(v,zero(T) => zero(T)/zero(T)) # NaN
    
    nanmean(x) = mean(filter(!isnan,x))
    nanmean(x,y) = mapslices(nanmean,x,dims=y)
    (u,v)=exch_UV(u,v); 
    for iF=1:u.grid.nFaces
        tmp1=u[iF][1:end-1,:]; tmp2=u[iF][2:end,:]
        uC[iF]=reshape(nanmean([tmp1[:] tmp2[:]],2),size(tmp1))
        tmp1=v[iF][:,1:end-1]; tmp2=v[iF][:,2:end]
        vC[iF]=reshape(nanmean([tmp1[:] tmp2[:]],2),size(tmp1))
    end
end

"""
    velocity2center(u,v,G)
  
    From Gael Forget, JuliaClimateNotebooks/Transport
    #1. Convert to `Sv` units and mask out land
    2. Interpolate `x/y` transport to grid cell center
"""
function velocity2center(u,v,G)

    #G = GridLoad(γ;option="full")
    u[findall(G.hFacW[:,1].==0)].=NaN
    v[findall(G.hFacS[:,1].==0)].=NaN;

    nanmean(x) = mean(filter(!isnan,x))
    nanmean(x,y) = mapslices(nanmean,x,dims=y)
    (u,v)=exch_UV(u,v); uC=similar(u); vC=similar(v)
    for iF=1:u.grid.nFaces
        tmp1=u[iF][1:end-1,:]; tmp2=u[iF][2:end,:]
        uC[iF]=reshape(nanmean([tmp1[:] tmp2[:]],2),size(tmp1))
        tmp1=v[iF][:,1:end-1]; tmp2=v[iF][:,2:end]
        vC[iF]=reshape(nanmean([tmp1[:] tmp2[:]],2),size(tmp1))
    end

    return uC, vC
end

"""
function rotate_uv(uvel,vvel,G)
From Gael Forget, JuliaClimateNotebooks/Transport
    3. Convert to `Eastward/Northward` transport
    4. Display Subdomain Arrays (optional)
"""
function rotate_uv(uvel,vvel,G)
    cs=G.AngleCS
    sn=G.AngleSN
    evel=uvel.*cs-vvel.*sn
    nvel=uvel.*sn+vvel.*cs

    return evel,nvel
end

end #module
