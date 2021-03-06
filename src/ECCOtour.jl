module ECCOtour
# write new functions and put them in this module.
# add them with text below, or create a new file in "src" and include it.

using MeshArrays, MITgcmTools, IsopycnalSurfaces
using Statistics, PyPlot, Distributions, FFTW,
    LinearAlgebra, StatsBase, LaTeXStrings,
    Dierckx, DelimitedFiles, Interpolations, NetCDF

import Statistics.mean, Statistics.std,
       Base.maximum, Base.minimum, Base.replace!,
       IsopycnalSurfaces.vars2sigma1, IsopycnalSurfaces.sigma1grid

export hanncoeffs, hannsum, hannsum!, hannfilter
export get_filtermatrix, matrixfilter, matrixspray, columnscale!
export seasonal_matrices, trend_matrices
export position_label, searchdir
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
export mdsio2dict, mdsio2sigma1, ncwritefromtemplate
export netcdf2sigma1, mdsio2regularpoles
export writeregularpoles, vars2regularpoles, var2regularpoles
export netcdf2regularpoles, factors4regularpoles
export regional_mask, apply_regional_mask!, zero2one!, wrapdist
export centerlon!, read_netcdf
export extract_timeseries,matmul,position_label,nancount
export faststats, allstats, std, mean
export mean, std, replace!
export velocity2center, rotate_uv, rotate_velocity!
export vars2sigma1, sigma1grid

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

function latlon(??)
    ??,?? = latlonC(??)
return ??,??
end

"""
    function latlonC(??)
    Latitude-longitude of ECCOv4r4 "C" (tracer) grid
"""
function latlonC(??)
    ??=??.read(??.path*"YC.data",MeshArray(??,Float64))
    ??=??.read(??.path*"XC.data",MeshArray(??,Float64))
return ??,??
end

"""
    function latlongG(??)
    Latitude-longitude of ECCOv4r4 "G" (velocity) grid
"""
function latlonG(??)
    ??=??.read(??.path*"YG.data",MeshArray(??,Float64))
    ??=??.read(??.path*"XG.data",MeshArray(??,Float64))
return ??,??
end

"""
    function depthlevels(??)
    Depths of ECCO v4r4 grid (positive vals)
"""
function depthlevels(??)
    z = read_mdsio(??.path,"RC")
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
    ????? = 1029 # from "data" text file in run directory
    g  = 9.81 # from "data"
    Pa2dbar = 1/10000 # standard pressures via hydrostatic balance
    pstdz = ????? .*g .* Pa2dbar .* z # 10000 to change Pa to dbar
    return pstdz
end

"""
    function readarea(??)
    area of ECCO v4r4 grid
"""
function readarea(??)
    area=??.read(??.path*"RAC.data",MeshArray(??,Float64))
return area
end

"""
     function patchmean(x,area,??,??,ispatch,iswet)
     get weight for rectangle region.
# Arguments
- `x`: variable of interest
- `area`: weighting
- `??`: lat
- `??`: lon
- `ispatch`: true in patch of interest
- `iswet` = function that is true if in ocean
# Output
- `xbar`: weighted filtered average of x
"""
function patchmean(x,area,??,??,ispatch,iswet)
    x = x*area
    xsum = 0.0
    asum = 0.0
    for i in eachindex(x)
        mask = ispatch.(??[i],??[i]) .* iswet.(x[i])
        for j in eachindex(??[i])
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

nino34mean(x,area,??,??,iswet) = patchmean(x,area,??,??,isnino34,iswet)
nino3mean(x,area,??,??,iswet) = patchmean(x,area,??,??,isnino3,iswet)
nino4mean(x,area,??,??,iswet) = patchmean(x,area,??,??,isnino4,iswet)
nino12mean(x,area,??,??,iswet) = patchmean(x,area,??,??,isnino12,iswet)

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
    Compute all statistics of gcmgrid type using function calls, same as faststats
    Warning: assumes land is zero
# Input
- `x::MeshArrays.gcmarray{T,N,Matrix{T}}`: input of gcmarray type
# Output
- `xbar::T`: mean value after discarding dryval
- `xmax::T`: maximum value
- `xmin::T`: minimum value
- `??x::T`: standard deviation
- `absxbar::T`: mean of absolute value
"""
function allstats(x::MeshArrays.gcmarray{T,N,Matrix{T}}) where {N,T<:AbstractFloat}

    dryval = 0.f0 # land points are zeroes, use NaN32 for NaN's

    xmax = -Inf
    xmin = Inf
    for n = 1:size(x,2)
        xmax = maximum(mask(x[:,n],xmax))# some trickery
        xmin = minimum(mask(x[:,n],xmin))
    end
    
    #xmax = maximum(x,dryval)
    #xmin = minimum(x,dryval)
    xbar = mean(x,dryval)

    # why does it require an explicit call to MeshArrays?
    absxbar = MeshArrays.mean(abs.(x),dryval)
    ??x   = std(x,xbar,dryval)

    return xbar, xmax, xmin, ??x, absxbar
end

"""
    function faststats(x)
    Compute fast statistics of gcmgrid type using function calls, eliminate redundancy
# Input
- `x::MeshArrays.gcmarray{T,1,Matrix{T}}`: input of gcmarray type
# Output
- `xbar::T`: mean value after discarding dryval
- `xmax::T`: maximum value
- `xmin::T`: minimum value
- `??x::T`: standard deviation
- `absxbar::T`: mean of absolute value
"""
function faststats(x::MeshArrays.gcmarray{T,1,Matrix{T}}) where T<:AbstractFloat

    xcount = [sum(count(notnanorzero,x[i])) for i in eachindex(x)]

    xmax = maximum(MeshArrays.mask(x,-Inf))
    xmin = minimum(MeshArrays.mask(x,Inf))
    
    if sum(xcount)>0
        #xmax = maximum([maximum(filter(notnanorzero,x[i])) for i ??? eachindex(x) if xcount[i] > 0])
        #xmin = minimum([minimum(filter(notnanorzero,x[i])) for i ??? eachindex(x) if xcount[i] > 0])

        # compute mean the old fashioned way
        xsum = sum([sum(filter(notnanorzero,x[i])) for i ??? eachindex(x) if xcount[i] > 0])
        xbar = xsum/sum(xcount)

        # compute standard deviation
        x??? = x.-xbar
        x??sum = sum([sum(filter(notnanorzero,x???[i]).^2) for i ??? eachindex(x???) if xcount[i]>0])
        ??x = sqrt(x??sum/(sum(xcount)-1))

        absxsum = sum([sum(filter(notnanorzero,abs.(x[i]))) for i ??? eachindex(x) if xcount[i] > 0]) # works b.c. 0 on land
        absxbar = absxsum/sum(xcount)
    else
        xbar = zero(T)/zero(T)
        #xmax = zero(T)/zero(T)
        #xmin = zero(T)/zero(T)
        ??x = zero(T)/zero(T)
        absxbar = zero(T)/zero(T)
    end
    return xbar, xmax, xmin, ??x, absxbar
end

"""
    function extract_sst34
    extract by reading multiple files
"""
function extract_sst34(expname,diagpath,fileroot,??,area,??,??,iswet)
    filelist = searchdir(diagpath[expname],fileroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nt = length(datafilelist)
    tt = 0
    sst34 = Float32[]
    for fname in datafilelist
        tt += 1
        println("time index ",tt)

        x = ??.read(diagpath[expname]*fname,MeshArray(??,Float32,1))

        # area-weighted filtered mean.
        sst34 = push!(sst34,nino34mean(x,area,??,??,iswet))
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
    function remove_seasonal(x,Ecycle,Fcycle,??) natively
    `x` = long monthly timeseries, starting in Jan

    remove seasonal cycle of this timeseries
"""
function remove_seasonal(x,Ecycle,Fcycle)
    # remove seasonal cycle from 14-day averaged timeseries
    # solve for seasonal cycle parameters
    ??cycle = Fcycle*x
    xseasonal = Ecycle*??cycle

    # remove it from total signal
    x??? = x - xseasonal
    return x???
end

"""
    function remove_seasonal(x,Ecycle,Fcycle,??) 
    remove seasonal cycle of this timeseries on the native gcmgrid

# Arguments
-  `x::gcmarray{T,N,Array{T,2}}` = gcmarray of a long timeseries where time is dimension N
- `Ecycle::Matrix`: matrix that operates on seasonal cycle parameters and return a seasonal cycle timeseries
- `Fcycle::Matrix`: precomputed pseudo-inverse of `Ecycle`
- `??::gcmgrid`: MITgcm grid
# Output
- `x???::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}`: potential density referenced to p??? minus 1000 kg/m??
"""
function remove_seasonal(x::MeshArrays.gcmarray{T,N,Array{T,2}},Ecycle,Fcycle,??)::MeshArrays.gcmarray{T,N,Array{T,2}} where T <: AbstractFloat where N

    # remove seasonal cycle from 14-day averaged timeseries
    # solve for seasonal cycle parameters

    # for gcmarray input:
    ??cycle = matmul(Fcycle,x,??)

    # reconstruct the full seasonal cycle.
    xseasonal = matmul(Ecycle,??cycle,??)
    x??? = x - xseasonal
    return x???
end

"""
    function factors4regularpoles(??)
    Get interpolation factors for regularpoles grid in one place
"""
function factors4regularpoles(??)
    ??Garc,??Carc = latgridArctic(??)
    ??Gantarc,??Cantarc = latgridAntarctic(??) 
    ??Greg,??Creg = latgridRegular(??) 
    ??G = -180.0:179.0
    ??C = -179.5:179.5
    ??G = vcat(??Gantarc,??Greg,??Garc)
    ??C = vcat(??Cantarc,??Creg,??Carc)
    ??arc = prereginterp(??Carc,??C,??)
    ??antarc = prereginterp(??Cantarc,??C,??)
    nx = length(??C)
    ny = length(??C)
    nyarc = length(??Carc)
    nyantarc = length(??Cantarc)

    return ??C,??G,??C,??G,nx,ny,nyarc,nyantarc,??arc,??antarc
end

"""
    function prereginterp(latgrid,longrid,??)
    prepare for regular interpolation
    regular = onto rectangular 2d map
# Arguments
- `latgrid`: 1d array of latitude
- `longrid`: 1d array of longitude
- `??`: GCM grid
# Output
- `f,i,j,w`: interpolation factors
"""
function prereginterp(latgrid,longrid,??)
    ?? = GridLoad(??;option="full")
    lon2d=[i for i=longrid, j=latgrid]
    lat2d=[j for i=longrid, j=latgrid]
    @time (f,i,j,w)=InterpolationFactors(??,vec(lon2d),vec(lat2d))
    ??=(lon=lon2d,lat=lat2d,f=f,i=i,j=j,w=w);
    return ??
end

"""
    function reginterp(fldin,nx,ny,f,i,j,w)
    regular interpolation with precomputed factors
    regular = onto rectangular 2d map
# Arguments
- `fldin`: gcmarray field of input
- `nx,ny`: x,y dimension of output regular field
- `??=(f,i,j,w)`: precomputed interpolation factors
# Output
- `fldout`: interpolated regular 2d field
"""
function reginterp(fldin,nx,ny,??)
    fldout = Interpolate(fldin,??.f,??.i,??.j,??.w)

    # think about moving this before interpolation
    # replace!(fldout,0.0 => NaN)
    fldout = transpose(reshape(fldout,nx,ny))
    
    return fldout
end

"""
    function trend_theta!(??,diagpathexp,path_out,tecco,??,F)
    get linear trend of potential temperature
    from monthly-average ECCOv4r4 gcmarray fields
# Arguments
- `??`: gcmarray field of trends
- `diagpathexp`: where to find files
- `tecco`: time stamps of monthly output
- `??`: GCM grid information
- `F`: linear estimator of trend
# Output
- `??`: updated trends
"""
function trend_theta!(??,diagpathexp,tecco,??,F)
    # name of file inside diagspath
    # Look at /poseidon ... exps/run/data.diagnostics for this info.
    Troot = "state_3d_set1" # hardcoded for ECCOv4r4

    fill!(??,0.0f0) # initialize ??
    nz = size(??,2)
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

        # read ?? for timestep
        ?? = ??.read(diagpathexp*Tname,MeshArray(??,Float32,nz))

        # multiply by the correct part of F matrix
        # be sure to handle all points at once
        # add to existing solution for ??
        # equal to matrix multiplication w/ columns (??) times weighting F
        ?? += F[2,tt] * ??
    end
end

"""
    function trend_theta(diagpathexp,path_out,tecco,??,F)
    get linear trend of potential temperature
    from monthly-average ECCOv4r4 gcmarray fields
# Arguments
- `diagpathexp`: where to find files
- `tecco`: time stamps of monthly output
- `??`: GCM grid information
- `F`: linear estimator of trend
# Output
- `??`: updated trends
"""
function trend_theta(diagpathexp,tecco,??,F)
    # name of file inside diagspath
    # Look at /poseidon ... exps/run/data.diagnostics for this info.
    Troot = "state_3d_set1" # hardcoded for ECCOv4r4

    # 50 shouldn't be hard coded
    ?? = MeshArray(??,Float32,50) # some nans here
    fill!(??,0.0f0) # initialize ??
    nz = size(??,2)
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

        # read ?? for timestep
        ?? = ??.read(diagpathexp*Tname,MeshArray(??,Float32,nz))

        # multiply by the correct part of F matrix
        # be sure to handle all points at once
        # add to existing solution for ??
        # equal to matrix multiplication w/ columns (??) times weighting F
        ?? += F[2,tt] * ??
    end
    return ??
end

function croplimitsLLC(??,??)
    # There is a regular grid inside the LLC grid.
    # Find the limiting indices of that regular grid.

    nx,ny = size(??[1])
    # find edge of antarctic distortion
    jantarc = 1
    while sum(abs.(diff(vec(??[1][:,jantarc])))) > 0.01
        jantarc += 1
    end
    ??antarc= ??[1][1,jantarc]

    # println("Southern edge of regular grid found")
    # println("at latitude ",??antarc,"??N")

    # find edge of arctic distortion
    jarc = ny
    while sum(abs.(diff(vec(??[1][:,jarc])))) > 0.01
        jarc -= 1
    end
    ??arc = ??[1][1,jarc]
    #println("Northern edge of regular grid found")
    #println("at latitude ",??arc,"??N")
    return jarc, jantarc, ??arc, ??antarc
end

function regularlatgrid(??)
    ??,?? = latlonG(??)
    jarc,jantarc,??arc,??antarc = croplimitsLLC(??,??)

    # get grid spacing for polar regions.
    ????arc= ??[1][1,jarc] - ??[1][1,jarc-1]
    ????antarc= ??[1][1,jantarc+1] - ??[1][1,jantarc]

    # get number of regularly-spaced polar gridcells
    ymax = 90; ymin = -82; # trim reg grid at 80 S
    narc = Int(ceil((ymax-??arc)/????arc))
    nantarc = Int(floor((??antarc-ymin)/????antarc))

    ??G = collect(reverse(range(??antarc,length=nantarc,step=-????antarc)))
    append!(??G,??[1][1,jantarc+1:jarc-1])
    append!(??G,range(??arc,length=narc,step=????arc))

    # Don't put ??C halfway between ??G points.
    # It doesn't match the MITgcm LLC grid.

    ??,?? = latlonC(??)
    jarc -= 1 # update for C grid, subtract one here
    ??arc = ??[1][1,jarc]
    ??antarc= ??[1][1,jantarc]
    ??C = collect(reverse(range(??antarc,length=nantarc,step=-????antarc)))
    append!(??C,??[1][1,jantarc+1:jarc-1])
    append!(??C,range(??arc,length=narc,step=????arc))

    return ??G, ??C
end

function latgridRegular(??)
    ??,?? = latlonG(??)
    jarc,jantarc,??arc,??antarc = croplimitsLLC(??,??)

    ??G = ??[1][1,jantarc:jarc]

    ??,?? = latlonC(??)
    jarc -= 1 # update for C grid, subtract one here
    ??C = ??[1][1,jantarc:jarc]

    return ??G, ??C
end

function latgridArctic(??)
    ??,?? = latlonG(??)
    jarc,jantarc,??arc,??antarc = croplimitsLLC(??,??)

    # get grid spacing for polar regions.
    ????arc= ??[1][1,jarc] - ??[1][1,jarc-1]

    # get number of regularly-spaced polar gridcells
    ymax = 90 # trim reg grid at 80 S
    narc = Int(ceil((ymax-??arc)/????arc))

    ??G=collect(range(??arc,length=narc,step=????arc))
    popfirst!(??G)

    # same thing for centered (tracer) grid
    ??,?? = latlonC(??)
    jarc -= 1 # update for C grid, subtract one here
    ??arc = ??[1][1,jarc]
    ??C=collect(range(??arc,length=narc,step=????arc))
    popfirst!(??C)

    return ??G, ??C
end

function latgridAntarctic(??)
    ??,?? = latlonG(??)
    jarc,jantarc,??arc,??antarc = croplimitsLLC(??,??)

    # get grid spacing for polar regions.
    ????antarc= ??[1][1,jantarc+1] - ??[1][1,jantarc]

    # get number of regularly-spaced polar gridcells
    ymin = -82; # trim reg grid at 80 S
    nantarc = Int(floor((??antarc-ymin)/????antarc))

    ??G = collect(reverse(range(??antarc,length=nantarc,step=-????antarc)))
    # eliminate northernmost grid point.
    # it's inside regular subgrid.
    pop!(??G)

    ??,?? = latlonC(??)
    ??antarc= ??[1][1,jantarc]
    ??C = collect(reverse(range(??antarc,length=nantarc,step=-????antarc)))
    pop!(??C)

    return ??G, ??C
end

function LLCcropC(gcmfield,??)
    # There is a regular grid inside the LLC grid.
    # Subsample/crop just those points. Put them together correctly.
    ??,?? = latlonG(??)
    jarc,jantarc,??arc,??antarc = croplimitsLLC(??,??)
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

function LLCcropG(gcmfield,??)
    # There is a regular grid inside the LLC grid.
    # Subsample/crop just those points. Put them together correctly.
    ??,?? = latlonG(??)
    jarc,jantarc,??arc,??antarc = croplimitsLLC(??,??)

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

"""
function mdsio2dict(pathin,filein,??)

Read native (i.e., mdsio) MITgcm files and return a Dictionary.
"""
function mdsio2dict(pathin,filein,??)

    metafile = filein*".meta"
    datafile = filein*".data"
    println(joinpath(pathin,metafile))
    meta = read_meta(joinpath(pathin,metafile));
    println(joinpath(pathin,datafile))
    state =  read_bin(joinpath(pathin,datafile),missing,missing,Float32,??)

    # changed this but don't understand why
    # other changes below, check ECCOtour for errors
    #    ndimz = meta[1].nDims
    ndimz = meta.nDims

    # should be able to use multiple dispatch with "N" parametric type below
    if ndimz == 3         # 3d
        #nz = meta[1].dimList[ndimz,1]
        nz = meta.dimList[ndimz,1]
        vars = Dict{String,MeshArrays.gcmarray{Float32,2,Matrix{Float32}}}()
    elseif ndimz == 2
        nz = 1
        vars = Dict{String,MeshArrays.gcmarray{Float32,1,Matrix{Float32}}}()
    end

    zlo = 1
    #for fldname in meta[1].fldList
    for fldname in meta.fldList

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
        filename = joinpath(fileprefix,fldname*filesuffix)
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
    function sigma(??,S,pz,p???)
    sigma values from SeaWaterDensity for gcmarrays
# Arguments
- `??::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}`: potential temperature
- `S::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}`: practical salinity
- `pz::Array{Float64,1}` : vertical profile of standard pressures
- `p???::Int64` : reference pressure for sigma value
# Output
- `??::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}`: potential density referenced to p??? minus 1000 kg/m??
"""
function sigma(??::MeshArrays.gcmarray{T,N,Array{T,2}},S::MeshArrays.gcmarray{T,N,Array{T,2}},pz::Array{Float64,1},p???::Int64) where T<:AbstractFloat where N

    # loop over faces
    nf,nz = size(??)
    ????? = similar(??)
    for ff = 1:nf
        nx,ny = size(??[ff,1])
        for xx = 1:nx
            for yy = 1:ny
                ??z = [??[ff,zz][xx,yy] for zz = 1:nz]
                Sz = [S[ff,zz][xx,yy] for zz = 1:nz]
                ??1,??2,??3 = MITgcmTools.SeaWaterDensity(??z,Sz,pz,p???)
                [?????[ff,zz][xx,yy] = convert(T,??3[zz]) .- 1000.0 for zz = 1:nz]
            end
        end
    end
    return ?????
end

"""
    function vars2sigma1(vars,p,sig1grid,??,spline_order)
    map variables onto sigma1 surfaces for gcmarrays
# Arguments
- `vars::Dict{String,MeshArrays.gcmarray{T,N,Array{T,2}}}`: dict of gcmarrays
- `p::Array{Float64,1}` : vertical profile of standard pressures
- `sig1grid`: ????? surface values
- `??`: grid description needed for preallocation
- `splorder`: 1-5, order of spline
- `linearinterp`: optional logical
# Output
- `vars??::Dict{String,MeshArrays.gcmarray{T,N,Array{T,2}}}`: dict of gcmarrays of variables on sigma1 surfaces
"""
function vars2sigma1(vars::Dict{String,MeshArrays.gcmarray{T,2,Matrix{T}}},pressure::Vector{Float64},sig1grid::Vector{Float64},??::gcmgrid;splorder=3,linearinterp=false,eos="JMD95") where T<:AbstractFloat 

    # ?? and S must exist
    !haskey(vars,"THETA") && error("?? missing")
    !haskey(vars,"SALT") && error("S missing")

    # loop over faces
    nf,nz = size(vars["THETA"])
    n?? = length(sig1grid)

    # vcol = Dict with profile/column data
    # pre-allocate each key in vars
    vcol = Dict{String,Vector{T}}() # vars in a column
    vars?? = Dict{String,MeshArrays.gcmarray{T,2,Matrix{T}}}()

    for (key, value) in vars
        vcol[key] = fill(convert(T,NaN),nz)
        vars??[key] = MeshArray(??,T,n??); fill!(vars??[key],convert(T,NaN))
    end
    # allocate standard pressure by hand.
    vars??["p"] = MeshArray(??,T,n??); fill!(vars??["p"],convert(T,NaN))

    for ff = 1:nf
        nx,ny = size(vars["THETA"][ff,1])
        for xx = 1:nx
            for yy = 1:ny
                for (vcolname, vcolval) in vars
                    # vcol = Dict with profile/column data
                    vcol[vcolname] = [vcolval[ff,zz][xx,yy] for zz = 1:nz]
                end

                nw = count(notnanorzero,vcol["THETA"]) # number of wet points in column
                nwS = count(notnanorzero,vcol["SALT"])
                if nw != nwS
                    error("T,S zeroes inconsistent")
                end

                if nw > 3 #need >=n+1 points to do order-n interpolation
                    ?????=sigma1column(vcol["THETA"][1:nw],vcol["SALT"][1:nw],pressure[1:nw],eos)

                    for (vckey,vcval) in vcol
                        var?? = var2sigmacolumn(?????,vcval[1:nw],sig1grid,splorder=splorder,linearinterp=linearinterp)
                        [vars??[vckey][ff,ss][xx,yy] = convert(T,var??[ss]) for ss = 1:n??]
                    end

                    # do standard pressure by hand.
                    p?? = var2sigmacolumn(?????,pressure[1:nw],sig1grid,splorder=splorder,linearinterp=linearinterp)
                    [vars??["p"][ff,ss][xx,yy] = convert(T,p??[ss]) for ss = 1:n??]
                    
                end
            end
        end
    end
    return vars??
end

"""
    function mdsio2sigma1(pathin::String,pathout::String,fileroots::Vector{String},??,pstdz,sig1grid;splorder=3,linearinterp=false,eos="JMD95") 

    Take variables in a filelist, read, map to sigma1, write to file.

# Arguments
- `pathin::String`: path of input
- `pathout::String`: path of output
- `fileroots::String`: beginning of file names in pathin
- `??`: grid description needed for preallocation
- `splorder::Integer`: optional keyword argument, default = 3, order of spline
- `linearinterp::Logical`: optional keyword argument, do linear interpolation?, default = false
- `eos::String`: optional key argument for equation of state, default = "JMD95"
"""
function mdsio2sigma1(pathin::String,pathout::String,fileroots::Vector{String},??,pstdz,sig1grid;splorder=3,linearinterp=false,eos="JMD95") 
    # Read All Variables And Puts Them Into "Vars" Dictionary

    # ideally would be more generic and Float64 would be covered.
    vars = Dict{String,MeshArrays.gcmarray{Float32,2,Matrix{Float32}}}()
    for fileroot in fileroots
        merge!(vars,mdsio2dict(pathin,fileroot,??))
    end

    # check for fields on the velocity (staggered) grid
    # and rotate them if necessary
    ?? = GridLoad(??;option="full")
    rotate_velocity!(vars,??)

    # solve for sigma1 on depth levels.
    @time vars?? = vars2sigma1(vars,pstdz,sig1grid,??,splorder=splorder,linearinterp=linearinterp,eos=eos)

    fileprefix = pathout
    # use first filename to get timestamp
    filesuffix = "_on_sigma1"*fileroots[1][14:end]*".data"
    write_vars(vars??,fileprefix,filesuffix)

    return vars??
end

"""
    function netcdf2sigma1
    Take variables in a NetCDF filelist, read, map to sigma1, write to file.
"""
function netcdf2sigma1(pathin,pathout,ncfilenames,??,pstdz,sig1grid,splorder)
    # Read All Variables And Puts Them Into "Vars" Dictionary
    println(ncfilenames)
    vars = Dict{String,Array{Float64,3}}()
    for (ncvarname,ncfilename) in ncfilenames
        merge!(vars,netcdf2dict(ncfilename,ncvarname))
    end

    # solve for sigma1 on depth levels.
    @time vars?? = vars2sigma1(vars,pstdz,sig1grid,splorder)

    fileprefix = pathout
    # use first filename to get timestamp
    filesuffixold = ncfilenames["THETA"][end-10:end]
    #    filesuffixnew = "_on_sigma1"*ncfilenames["THETA"][end-10:end]
    filesuffixnew = "_on_sigma1"*filesuffixold
    ncwritefromtemplate(vars??,fileprefix,filesuffixold,filesuffixnew,sig1grid)
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
function netcdf2dict(ncfilename::String,ncvarname::String,??::gcmgrid)
    # 3d
    println(ncfilename,ncvarname)
    state = read_netcdf(ncfilename,ncvarname,??);
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
function replace!(f::Function,a::MeshArrays.gcmarray{T,N,Matrix{T}}) where {N,T<:AbstractFloat}

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

function netcdf2regularpoles(ncfilename,ncvarname,??,nx,ny,nyarc,??arc,nyantarc,??antarc)

    vars = netcdf2dict(ncfilename,ncvarname,??)
    varsregpoles = vars2regularpoles(vars,??,nx,ny,nyarc,??arc,nyantarc,??antarc)

end

"""
function mdsio2regularpoles(pathin,filein,??,nx,ny,nyarc,??arc,nyantarc,??antarc)
"""
function mdsio2regularpoles(pathin,filein,??,nx,ny,nyarc,??arc,nyantarc,??antarc)

    vars = mdsio2dict(pathin,filein,??)
    ?? = GridLoad(??;option="full")
    rotate_velocity!(vars,??)
    varsregpoles = vars2regularpoles(vars,??,nx,ny,nyarc,??arc,nyantarc,??antarc)

end

function vars2regularpoles(vars::Dict{String,MeshArrays.gcmarray{T,2,Matrix{T}}},??,nx,ny,nyarc,??arc,nyantarc,??antarc) where T<:AbstractFloat

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
            varsregpoles[varname][:,:,zz]=var2regularpoles(varvals[:,zz],??,nx,ny,nyarc,??arc,nyantarc,??antarc)
        end
    end
    return varsregpoles
end

"""
function vars2regularpoles(vars::Dict{String,MeshArrays.gcmarray{T,1,Matrix{T}}},??,nx,ny,nyarc,??arc,nyantarc,??antarc)
     variables interpolated onto regularpoles grid
"""
function vars2regularpoles(vars::Dict{String,MeshArrays.gcmarray{T,1,Matrix{T}}},??,nx,ny,nyarc,??arc,nyantarc,??antarc) where T<:AbstractFloat

    varsregpoles = Dict{String,Matrix{T}}()
    for (varname, varvals) in vars

        varsregpoles[varname] = var2regularpoles(varvals,??,nx,ny,nyarc,??arc,nyantarc,??antarc)

    end
    return varsregpoles
end

"""
   var2regularpoles
   Take one gcmarray in memory, put on regularpoles grid
"""
function var2regularpoles(var,??,nx,ny,nyarc,??arc,nyantarc,??antarc)

    T = eltype(var)
    NaNT = zero(T)/zero(T)
    # remove contamination from land
    replace!(var, 0.0 => NaNT)

    #pre-allocate output
    varregpoles = fill(NaNT,(nx,ny))

    # get regular grid by cropping
    ??crop =  LLCcropC(var,??)

    # interpolate to "regularpoles"
    ??arc = reginterp(var,nx,nyarc,??arc)
    ??antarc = reginterp(var,nx,nyantarc,??antarc)
    varregpoles=hcat(??antarc',??crop,??arc')

    return varregpoles
    
end

"""
function writeregularpoles(vars,??,pathout,filesuffix,filelog,??C,lonatts,??C,latatts,z,depthatts)
"""
function writeregularpoles(vars::Dict{String,Array{Float32,3}},??,pathout,filesuffix,filelog,??C,lonatts,??C,latatts,z,depthatts)

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
            ??C,
            lonatts,
            "lat",
            ??C,
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
function writeregularpoles(vars::Dict{String,Array{Float32,2}},??,pathout,filesuffix,filelog,??C,lonatts,??C,latatts,z,depthatts)
"""
function writeregularpoles(vars::Dict{String,Array{Float32,2}},??,pathout,filesuffix,filelog,??C,lonatts,??C,latatts,z,depthatts)

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
            ??C,
            lonatts,
            "lat",
            ??C,
            latatts,
            atts = varatts,
        )
        ncwrite(varvals, fileout, varname)
    end
end

"""
function writeregularpoles(vars,??,pathout,filesuffix,filelog,??C,lonatts,??C,latatts,z,depthatts)
"""
function writeregularpoles(vars::Dict{String,Array{Float64,3}},??,pathout,filesuffix,filelog,??C,lonatts,??C,latatts,z,depthatts)

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
            ??C,
            lonatts,
            "lat",
            ??C,
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
function writeregularpoles(vars,pathout,filesuffix,filelog,??C,lonatts,??C,latatts)
"""
function writeregularpoles(vars::Dict{String,Array{Float64,2}},??,pathout,filesuffix,filelog,??C,lonatts,??C,latatts,z,depthatts)

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
            ??C,
            lonatts,
            "lat",
            ??C,
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
    extract_timeseries(froot,years,??,xval,yval,fval)
# Arguments
- `froot::String`: filename root
- `years::StepRange`: iterator for multiple files
- `??::gcmarray`: GCM grid from MeshArrays.jl
- `xval::Integer`: x grid index
- `yval::Integer`: y grid index
- `fval::Integer`: face index
# Output    
- `tseries`: timeseries
- `nseries`: nseries = number of timeseries elements in each year
```
"""
function extract_timeseries(froot,years,??,xval,yval,fval)
    tseries = []
    nseries = []
    nyr = length(years)
    for tt = 1:nyr
        fname = froot*string(years[tt])
        println(fname)
        field = read_bin(fname,Float32,??)
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
    
    for tval ??? 1:nt
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
    
    for tval ??? 1:nt
        tmp = flux[fval,tval]
        push!(series,tmp[xval,yval])
    end
    return series
end

#################################################################################################
# tell julia how to do matrix multiplication
function matmul(M::Array{Float64,2},flux::MeshArrays.gcmarray{Float64,2,Array{Float64,2}},??) ::MeshArrays.gcmarray{Float64,2,Array{Float64,2}}

    nM = size(M,1)  # matrix size M
    nN = size(flux,2) # matrix size N
    nQ  = size(flux,1) # repeat matmul Q times
    
    # initialize product
    product = 0.0 .* MeshArray(??,Float64,nM)
    
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

function matmul(M::Array{Float32,2},flux::MeshArrays.gcmarray{Float32,2,Array{Float32,2}},??::gcmgrid) ::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}

    nM = size(M,1)  # matrix size M
    nN = size(flux,2) # matrix size N
    nQ  = size(flux,1) # repeat matmul Q times
    
    # initialize product
    #product = 0.0f0 .* MeshArray(??,Float32,nM)
    println("sub-optimal initialization")
    product = MeshArray(??,Float32,nM) # some nans here
    tmp1=zeros(Float32,Tuple(??.ioSize))
    for tt= 1:nM
        # initialize a sub-optimal way
        product[:,tt]=??.read(tmp1,MeshArray(??,Float32))
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

function matmul(M::Array{Float32,1},flux::MeshArrays.gcmarray{Float32,1,Array{Float32,2}},??::gcmgrid) ::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}

    nM = size(M,1)  # matrix size M
    nQ  = size(flux,1) # repeat matmul Q times
    
    # initialize product
    #product = 0.0f0 .* MeshArray(??,Float32,nM)
    println("sub-optimal initialization")
    product = MeshArray(??,Float32,nM) # some nans here
    tmp1=zeros(Float32,Tuple(??.ioSize))
    for tt= 1:nM
        # initialize a sub-optimal way
        product[:,tt]=??.read(tmp1,MeshArray(??,Float32))
    end

    for mm = 1:nM
        for qq = 1:nQ # inner product over faces
            product[qq,mm] += flux[qq] * M[mm]
                #product[ff,pp] +=  transpose(F[pp,tt]'*flux[ff,tt]')
        end
    end
    return product
end

function matmul(M::Array{Float64,1},flux::MeshArrays.gcmarray{Float64,1,Array{Float64,2}},??) ::MeshArrays.gcmarray{Float64,2,Array{Float64,2}}

    nM = size(M,1)  # matrix size M
    nQ  = size(flux,1) # repeat matmul Q times
    
    # initialize product
    #product = 0.0 .* MeshArray(??,Float64,nM)
    println("sub-optimal initialization")
    product = MeshArray(??,Float64,nM) # some nans here
    tmp1=zeros(Float64,Tuple(??.ioSize))
    for tt= 1:nM
        # initialize a sub-optimal way
        product[:,tt]=??.read(tmp1,MeshArray(??,Float64))
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
- `??x::T`: standard deviation 
"""
function std(x::MeshArrays.gcmarray{T,N,Array{T,2}},xbar::T,dryval::T) where T<:AbstractFloat where N

    isnan(dryval) ? isdry = isnan : isdry(z) = (z == dryval)

    # x prime = fluctuation
    x??? = x .- xbar
    
    #  vector list of nonzero elements
    xcount = [sum(count(!isdry,x???[i])) for i in eachindex(x???)]
    if sum(xcount) > 0
        # don't assume 0 on land
        x??sum = sum([sum(filter(!isdry,x???[i]).^2) for i in eachindex(x) if xcount[i] > 0])
        ??x = sqrt(x??sum/(sum(xcount)-1))
    else
        xbar = convert(T,NaN)
    end
    return ??x
end

function rotate_velocity!(vars,??)

    velvars = (("UVELMASS","VVELMASS"),("GM_PsiX","GM_PsiY"),("oceTAUX","oceTAUY"))
    velchange = Dict("UVELMASS" => "EVELMASS",
                     "VVELMASS" => "NVELMASS",
                     "GM_PsiX" => "GM_PsiE",
                     "GM_PsiY" => "GM_PsiN",
                     "oceTAUX" => "oceTAUE",
                     "oceTAUY" => "oceTAUN")

    for (kk,vv) in velvars
        if haskey(vars,kk) && haskey(vars,vv)

            println("rotate variable pair:",kk,"-",vv)
            evel = similar(vars[kk])
            nvel = similar(vars[vv])

            if ndims(vars[kk]) > 1

                uC = similar(vars[kk][:,1]) # one level
                vC = similar(vars[vv][:,1]) # one level
                for zz = 1: size(vars[kk],2)

                    #interpolate velocity to center of C-grid
                    velocity2center!(uC,vC,vars[kk][:,zz],vars[vv][:,zz],??)
                
                    evel[:,zz],nvel[:,zz] = rotate_uv(uC,vC,??);
                end

            elseif ndims(vars[kk]) == 1

                uC = similar(vars[kk]) # one level
                vC = similar(vars[vv]) # one level

                #interpolate velocity to center of C-grid
                velocity2center!(uC,vC,vars[kk],vars[vv],??)               
                evel,nvel = rotate_uv(uC,vC,??);

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

    #G = GridLoad(??;option="full")
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
function rotate_uv(uvel::MeshArrays.gcmarray{T,N,Matrix{T}},vvel::MeshArrays.gcmarray{T,N,Matrix{T}},G) where T<:AbstractFloat where N

    # this kludge must be very slow
    cs=convert.(T,G.AngleCS)
    sn=convert.(T,G.AngleSN)
    
    evel=uvel.*cs-vvel.*sn
    nvel=uvel.*sn+vvel.*cs

    return evel,nvel
end

end #module
