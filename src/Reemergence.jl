module Reemergence
# write new functions and put them in this module.
# add them with text below, or create a new file in "src" and include it.

using Statistics, PyPlot, Distributions, FFTW, LinearAlgebra, StatsBase
using MeshArrays, MITgcmTools, LaTeXStrings, Dierckx
using DelimitedFiles, Interpolations, NetCDF

export hanncoeffs, hannsum, hannsum!, hannfilter
export get_filtermatrix, matrixfilter, matrixspray, columnscale!
export seasonal_matrices, trend_matrices
export position_label, searchdir, setupLLCgrid
export listexperiments, expnames, expsymbols, time_label
export faststats, allstats, std, mean
export inrectangle, isnino34, isnino3, isnino4, isnino12
export latlon, latlonC, latlonG
export depthlevels, pressurelevels, readarea, patchmean
export nino34mean, nino3mean, nino4mean, nino12mean, extract_sst34
export remove_climatology, remove_seasonal, sigma, vars2sigma1
export historicalNino34, prereginterp, reginterp, trend_theta!
export regularlatgrid, LLCcropC, LLCcropG, croplimitsLLC
export latgridAntarctic, latgridArctic, latgridRegular
export timestamp_monthly_v4r4, sigma1column, var2sigmacolumn
export sigma1grid, mdsio2dict, netcdf2dict, write_vars
export mdsio2sigma1, ncwritefromtemplate
export netcdf2sigma1, replace!, mdsio2regularpoles
export writeregularpoles, vars2regularpoles
export netcdf2regularpoles, mixinversions!, dedup!
export infadezoneE, infadezoneW, infadezoneN, infadezoneS

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
    function setupLLCgrid(path_grid)
    Load the ECCOv4r4 grid at location path_grid into a MeshArray variable
# Arguments
- `path_grid`: directory to search for grid files
# Output
- `γ::gcmgrid`: MITgcm grid
"""
function setupLLCgrid(path_grid)
    http="https://github.com/gaelforget/GRID_LLC90"
    !isdir(path_grid) ? run(`git clone $http $path_grid`) : nothing;
    γ=GridSpec("LatLonCap",path_grid)
    return γ
end

"""
    function Nino34file()
    Get location of the historical Nino34 Google Drive file
    Download if necessary
# Output
- `fileloc`: local location of gdrive file
"""
function Nino34file()
    # exist on local machine?
    if isdir(expanduser("~/gdrive/"))
         filename = "Shared_drives/NSF_reemergence/ECCOv4r4/data/nino34.hadisst1.1870-2020.txt"
        fileloc = expanduser("~/gdrive/")*filename
    else
        inputdir = "../inputs"
        !isdir(inputdir) ? mkdir(inputdir) : nothing

        # download from google drive and save location
        run(`wget "https://drive.google.com/file/d/1kOtOnD6B3Y9SAI5W6ezP-o_tgNkonn5n/view?usp=sharing" -O ../inputs/nino34.hadisst1.1870-2020.txt`)
        fileloc = inputdir*"nino34.hadisst1.1870-2020.txt"
    end
    return fileloc
end

"""
    function readNino34()
    Get historical Nino3.4 data from HadISST 1
# Output
- `SST_nino34`: local location of gdrive file
"""
function readNino34()
    #println("start")
    filename = Nino34file()
    #println(typeof(filename))
    SST_nino34 = DelimitedFiles.readdlm(filename)
    return SST_nino34
end

"""
    function historicalNino34(baselineyears)
    Get historical Nino3.4 data and SST climatology from HadISST 1
# Argument
- `baselineyears`: for computation of SST climatology, i.e., (1960,2020)
# Output
- `nino34`: historical Nino3.4 index
- `tnino34`: time in years CE corresponding to index
- `SSTclimatology`: monthly mean values in Nino3.4 patch
"""
function historicalNino34(baselineyears)
    SSTnino34 = readNino34()
    nyr,nmon = size(SSTnino34)

    SSTclimatology = zeros(12) # 12 months
    count = 0
    for i = 1:nyr
        if baselineyears[1] <= SSTnino34[i,1] <= baselineyears[2]
            SSTclimatology += SSTnino34[i,2:end]  # skip the year in column 1
            count += 1
        end
    end
    SSTclimatology = SSTclimatology/count

    # get nino3.4 timeseries
    nino34 = []
    tnino34 = []
    monlist= 1/24:1/12:1
    for i  = 1:nyr
        append!(nino34,SSTnino34[i,2:end]-SSTclimatology)
        append!(tnino34,SSTnino34[i,1] .+ collect(monlist))
    end
    return nino34,tnino34,SSTclimatology
end

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

"""
    function infadezoneE(latpt,lonpt,latrect,lonrect,dlon)
    find all gridcells within East 'wing of box' around Cartesian rectangle
# Arguments
- `latpt`: latitude grid
- `lonpt`: longitude grid
- `latrect`: tuple of latitude range of central rectangle
- `lonrect`: tuple of longitude range of central rectangle
- `dlon`: width in degrees longitude of East fade zone
# Output
- `fadezone`: boolean of type gcmarray
"""
function infadezoneE(latpt,lonpt,latrect,lonrect,dlon)
    # load latitude and longitude on this grid.
    if lonrect[1] < lonrect[2] # if center rectangle doesn't wrap around dateline
        if lonrect[1] - dlon >= -180 && lonrect[2] + dlon <= 180 # if East fade zone does NOT wrap around date line
            infadezoneE = (maximum(lonrect) <= lonpt <= (maximum(lonrect) + dlon))
            && minimum(latrect) <= latpt <= maximum(latrect)
        elseif lonrect[2] + dlon > 180 # if East fade zone crosses over dateline
            infadezoneE = ((maximum(lonrect) <= lonpt || lonpt <= -360 + maximum(lonrect) + dlon)
            && minimum(latrect) <= latpt <= maximum(latrect)
        end
    else # central rectangle wraps around the date line (therefore east fade zone will not)
        infadezoneE =  minimum(lonrect) <= lonpt <= (minimum(lonrect) + dlon)
            && minimum(latrect) <= latpt <= maximum(latrect)
    end
    # watch out if lonrect[1] == lonrect[2]
end

"""
    function infadezoneW(latpt,lonpt,latrect,lonrect,dlon)
    find all gridcells within West 'wing of box' around Cartesian rectangle
# Arguments
- `latpt`: latitude grid
- `lonpt`: longitude grid
- `latrect`: tuple of latitude range
- `lonrect`: tuple of longitude range
- `dlon`: width in degrees longitude of West fade zone
# Output
- `fadezone`: boolean of type gcmarray
"""
function infadezoneW(latpt,lonpt,latrect,lonrect,dlon)
    # load latitude and longitude on this grid.
    if lonrect[1] < lonrect[2] # if center rectangle doesn't wrap around dateline
        if lonrect[1] - dlon >= -180 && lonrect[2] + dlon <= 180 # if West fade zone does NOT wrap around date line
            infadezoneW = ((minimum(lonrect)-dlon) <= lonpt <= minimum(lonrect))
            && minimum(latrect) <= latpt <= maximum(latrect)
        elseif lonrect[2] + dlon > 180 # if West fade zone crosses over dateline
            infadezoneW = ((minimum(lonrect) >= lonpt || lonpt >= 360 - minimum(lonrect) - dlon)
            && minimum(latrect) <= latpt <= maximum(latrect)
        end
    else # central rectangle wraps around the date line (therefore west fade zone will not)
        infadezoneW = (maximum(lonrect)-dlon) <= lonpt <= (maximum(lonrect))
            && minimum(latrect) <= latpt <= maximum(latrect)
    end
    # watch out if lonrect[1] == lonrect[2]
end

"""
    function infadezoneN(latpt,lonpt,latrect,lonrect,dlat)
    find all gridcells within North 'wing of box' around Cartesian rectangle
# Arguments
- `latpt`: latitude grid
- `lonpt`: longitude grid
- `latrect`: tuple of latitude range
- `lonrect`: tuple of longitude range
- `dlat`: width in degrees latitude of North fade zone
# Output
- `fadezone`: boolean of type gcmarray
"""
function infadezoneN(latpt,lonpt,latrect,lonrect,dlat)
    # load latitude and longitude on this grid.
    if lonrect[1] < lonrect[2] #if center rectangle crosses dateline
        infadezoneN = minimum(lonrect) <= lonpt <= maximum(lonrect) &&
            maximum(latrect) <= latpt <= (maximum(latrect) + dlat)
    else # center rectangle wraps around the date line
        infadezoneN = (lonpt <= minimum(lonrect) ||
            lonpt >= maximum(lonrect)) &&
            maximum(latrect) <= latpt <= (maximum(latrect) + dlat)
    end
    # watch out if lonrect[1] == lonrect[2]
end

"""
    function infadezoneS(latpt,lonpt,latrect,lonrect,dlat)
    find all gridcells within South 'wing of box' around Cartesian rectangle
# Arguments
- `latpt`: latitude grid
- `lonpt`: longitude grid
- `latrect`: tuple of latitude range
- `lonrect`: tuple of longitude range
- `dlat`: width in degrees latitude of South fade zone
# Output
- `rectangle`: boolean of type gcmarray
"""
function infadezoneS(latpt,lonpt,latrect,lonrect,dlat)
    # load latitude and longitude on this grid.
    if lonrect[1] < lonrect[2] #if center rectangle crosses dateline
        infadezoneS = minimum(lonrect) <= lonpt <= maximum(lonrect) &&
            minimum(latrect) >= latpt >= (minimum(latrect) - dlat)
    else # center rectangle wraps around the date line
        infadezoneS = (lonpt <= minimum(lonrect) ||
            lonpt >= maximum(lonrect)) &&
            minimum(latrect) >= latpt >= (minimum(latrect) - dlat)
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
    function listexperiments(exppath)
    get a dictionary of experiments and their locations
# Arguments
- `exppath`: directory of experiments
# Output
- `runpath`: dictionary with keys=experiments, values=run paths
- `diagpath`: dictionary with keys=experiments, values=diagnostic paths
"""
function listexperiments(exppath)
    dirlist = searchdir(exppath,"") # all files in directory
    explist  = filter(x -> !occursin("README",x),dirlist) # remove README to get explist

    runpath = Dict(explist .=> exppath.*explist.*"/run/")
    diagpath = Dict(explist .=> exppath.*explist.*"/run/diags/")
    regpolespath = Dict(explist .=> exppath.*explist.*"/run/regularpoles/")

    # print to screen all the available
    println("Available experiments are")
    for (key,value) in runpath
        println(key)
    end
    return runpath,diagpath,regpolespath
end

"""
    function expnames()
    Abbreviations/shortnames for some experiments, useful for labels
    Hand-coded and needs manual change with new experiments.
# Output
- `shortnames`: dictionary with keys=experiments, values=abbreviations
"""
function expnames()
    shortnames = Dict("iter129_bulkformula" => "129bf")
    push!(shortnames,"iter0_bulkformula" => "0bf")
    push!(shortnames,"iter129_fluxforced" => "129ff")
    push!(shortnames,"noinitadjust" => "noINIT")
    push!(shortnames,"nosfcadjust" => "noSFC")
    push!(shortnames,"nointerannual" => "noIA")
    return shortnames
end

"""
    function expsymbols()
    List of symbols for some experiments, useful for distinguishing plots
    Hand-coded and needs manual change with new experiments.
# Output
- `marks`: dictionary with keys=experiments, values=abbreviations
"""
function expsymbols()
    marks = Dict("iter129_bulkformula" => "o")
    push!(marks,"iter0_bulkformula" => "x")
    push!(marks,"iter129_fluxforced" => "*")
    push!(marks,"noinitadjust" => ".")
    push!(marks,"nosfcadjust" => "+")
    push!(marks,"nointerannual" => "s")
    return marks
end

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
    println(xcount)
    if sum(xcount)>0
        xmax = maximum([maximum(filter(notnanorzero,x[i])) for i ∈ eachindex(x) if xcount[i] > 0])
        xmin = minimum([minimum(filter(notnanorzero,x[i])) for i ∈ eachindex(x) if xcount[i] > 0])

        # compute mean the old fashioned way
        xsum = sum([sum(filter(notnanorzero,x[i])) for i ∈ eachindex(x) if xcount[i] > 0])
        println(xsum)

        xbar = xsum/sum(xcount)
        println(xbar)

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
    function std(x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}},xbar::Float32,dryval)
    Compute standard deviation of gcmgrid type using function calls, eliminate redundancy
    Eliminate all values = dryval
    Avoid a name clash with StatsBase and MeshArrays
"""
function std(x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}},xbar::Float32,dryval)
     MeshArrays.std(x,xbar,dryval)
end

# avoid another name clash.
"""
    function mean(x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}},dryval::Float32)
    Compute mean of gcmgrid type using function calls, eliminate redundancy
    Eliminate all values = dryval
    Avoid a name clash with StatsBase and MeshArrays
"""
function mean(x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}},dryval::Float32)
     MeshArrays.mean(x,dryval)
end

"""
    function mean(x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}},dryval::Float32)
    Compute area-weighted mean of gcmgrid type using function calls, eliminate redundancy
    Area weighting = area
    Eliminate all values = dryval
    Avoid a name clash with StatsBase and MeshArrays
"""
function mean(x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}},area::MeshArrays.gcmarray{Float64,1,Array{Float64,2}},dryval::Float64)
    MeshArrays.mean(x,area,dryval)
end

"""
    function mean(x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}},dryval::Float32)
    Compute area-weighted mean of gcmgrid type using filtered with function isgood
    Area weighting = area
    Eliminate all values = isgood(x) -> true
    Avoid a name clash with StatsBase and MeshArrays
"""
function mean(x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}},area::MeshArrays.gcmarray{Float64,1,Array{Float64,2}},isgood)
    MeshArrays.mean(x,area,isgood)
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
    xprime = x - xseasonal
    return xprime
end

"""
    function remove_seasonal(x,Ecycle,Fcycle,γ) natively
    `x` = gcmarray of long monthly timeseries, starting in Jan
    remove seasonal cycle of this timeseries
"""
function remove_seasonal(x,Ecycle,Fcycle,γ)

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
    function sigma(θ,S,pz,p₀)
    sigma values from SeaWaterDensity for gcmarrays
# Arguments
- `θ::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}`: potential temperature
- `S::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}`: practical salinity
- `pz::Array{Float64,1}` : vertical profile of standard pressures
- `p₀::Int64` : reference pressure for sigma value
# Output
- `σ::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}`: potential density referenced to p₀ minus 1000 kg/m³
"""
function sigma(θ::MeshArrays.gcmarray{Float32,2,Array{Float32,2}},S::MeshArrays.gcmarray{Float32,2,Array{Float32,2}},pz::Array{Float64,1},p₀::Int64)

    # loop over faces
    nf,nz = size(θ)
    σ₁ = similar(θ)
    for ff = 1:nf
        nx,ny = size(θ[ff,1])
        for xx = 1:nx
            for yy = 1:ny
                θz = [θ[ff,zz][xx,yy] for zz = 1:nz]
                Sz = [S[ff,zz][xx,yy] for zz = 1:nz]
                σ1,σ2,σ3 = MITgcmTools.SeaWaterDensity(θz,Sz,pz,p₀)
                [σ₁[ff,zz][xx,yy] = convert(Float32,σ3[zz]) .- 1000 for zz = 1:nz]
            end
        end
    end
    return σ₁
end

"""
    function vars2sigma1(vars,p,sig1grid,γ,spline_order)
    map variables onto sigma1 surfaces for gcmarrays
# Arguments
- `vars::Dict{String,MeshArrays.gcmarray{Float32,2,Array{Float32,2}}}`: dict of gcmarrays
- `p::Array{Float64,1}` : vertical profile of standard pressures
- `sig1grid`: σ₁ surface values
- `γ`: grid description needed for preallocation
- `splorder`: 1-5, order of spline
# Output
- `varsσ::Dict{String,MeshArrays.gcmarray{Float32,2,Array{Float32,2}}}`: dict of gcmarrays of variables on sigma1 surfaces
"""
function vars2sigma1(vars::Dict{String,MeshArrays.gcmarray{Float32,2,Array{Float32,2}}},pressure::Array{Float64,1},sig1grid::Array{Float64,1},γ::gcmgrid,splorder::Integer)

    # check that θ and S exist. They must.
    (haskey(vars,"THETA") && haskey(vars,"SALT")) ? nothing : error("Need θ and S in vars")

    # loop over faces
    nf,nz = size(vars["THETA"])
    nσ = length(sig1grid)

    # vcol = Dict with profile/column data
    # pre-allocate each key in vars
    vcol = Dict{String,Array{Float64,1}}() # vars in a column
    varsσ = Dict{String,MeshArrays.gcmarray{Float32,2,Array{Float32,2}}}()

    for (key, value) in vars
        vcol[key] = fill(NaN32,nz)
        varsσ[key] = MeshArray(γ,Float32,nσ); fill!(varsσ[key],NaN32)
    end
    # allocate standard pressure by hand.
    # CONSIDER ANOTHER FUNCTION TO DO PHIHYD.
    varsσ["p"] = MeshArray(γ,Float32,nσ); fill!(varsσ["p"],NaN32)

    for ff = 1:nf
        nx,ny = size(vars["THETA"][ff,1])
        for xx = 1:nx
            for yy = 1:ny
                for (vcolname, vcolval) in vars
                    # vcol = Dict with profile/column data
                    vcol[vcolname] = [vcolval[ff,zz][xx,yy] for zz = 1:nz]
                end

                # also need to filter dry values and to change zz
                # Consider using `isdry` function and dryval in future.
                nw = count(notnanorzero,vcol["THETA"]) # number of wet points in column
                nwS = count(notnanorzero,vcol["SALT"])
                if nw != nwS
                    error("T,S zeroes inconsistent")
                end

                # incurs error if splorder > number of points in column
                if nw > splorder #need >=n+1 points to do order-n interpolation
                    #println(length(vcol["θ"][1:nw]))
                    σ₁=sigma1column(vcol["THETA"][1:nw],vcol["SALT"][1:nw],pressure[1:nw])

                    # 1) no inversions or homogeneity
                    # 2) range of sig1, 3) no extrapolation
                    if sum(diff(σ₁).<0)==0 && count(minimum(σ₁).<=sig1grid.<=maximum(σ₁)) > 0

                        # eliminate any extrapolation
                        sgood = findall(minimum(σ₁).<=sig1grid.<=maximum(σ₁))
                        ngood = length(sgood)

                        for (vckey,vcval) in vcol
                            varσ = var2sigmacolumn(σ₁,vcval[1:nw],sig1grid[sgood],splorder)

                            [varsσ[vckey][ff,sgood[ss]][xx,yy] = convert(Float32,varσ[ss]) for ss = 1:ngood]
                        end

                        # do standard pressure by hand.
                        pσ = var2sigmacolumn(σ₁,pressure[1:nw],sig1grid[sgood],splorder)
                        [varsσ["p"][ff,sgood[ss]][xx,yy] = convert(Float32,pσ[ss]) for ss = 1:ngood]

                    end
                end
            end
        end
    end
    return varsσ
end

"""
    function vars2sigma1(vars,p,sig1grid,γ,spline_order)
    map variables from regularpoles grid onto sigma1 surfaces
# Arguments
- `vars::Dict{String,Array{Float64,3}}}`: dict of 3d arrays
- `p::Array{Float64,1}` : vertical profile of standard pressures
- `sig1grid`: σ₁ surface values
- `γ`: grid description needed for preallocation
- `splorder`: 1-5, order of spline
# Output
- `varsσ::Dict{String,MeshArrays.gcmarray{Float32,2,Array{Float32,2}}}`: dict of gcmarrays of variables on sigma1 surfaces
"""
function vars2sigma1(vars::Dict{String,Array{Float64,3}},pressure::Array{Float64,1},sig1grid::Array{Float64,1},splorder::Int64)

    # check that θ and S exist. They must.
    (haskey(vars,"THETA") && haskey(vars,"SALT")) ? nothing : error("Need θ and S in vars")

    # loop over faces
    nx,ny,nz = size(vars["THETA"])
    nσ = length(sig1grid)

    # vcol = Dict with profile/column data
    # pre-allocate each key in vars
    vcol = Dict{String,Array{Float64,1}}() # vars in a column
    varsσ = Dict{String,Array{Float64,3}}()

    for (key, value) in vars
        vcol[key] = fill(NaN,nz)
        varsσ[key] = fill(NaN,(nx,ny,nσ))
    end

    # allocate standard pressure by hand.
    varsσ["p"] = fill(NaN,(nx,ny,nσ))

    for xx = 1:nx
        for yy = 1:ny
            for (vcolname, vcolval) in vars
                # vcol = Dict with profile/column data
                vcol[vcolname] = vcolval[xx,yy,:]
            end

            # also need to filter dry values and to change zz
            # Consider using `isdry` function and dryval in future.
            nw = count(notnanorzero,vcol["THETA"]) # number of wet points in column
            nwS = count(notnanorzero,vcol["SALT"])
            if nw != nwS
                error("T,S zeroes inconsistent")
            end

            if nw > 3
                # incurs error if splorder > number of points in column
                # if nw > splorder #need >=n+1 points to do order-n interpolation
                σ₁=sigma1column(vcol["THETA"][1:nw],vcol["SALT"][1:nw],pressure[1:nw])

                for (vckey,vcval) in vcol
                    varσ = var2sigmacolumn(σ₁,vcval[1:nw],sig1grid,splorder)
                    # eliminate sgood here.
                    [varsσ[vckey][xx,yy,ss] = convert(Float32,varσ[ss]) for ss = 1:nσ]
                end

                    # do standard pressure by hand.
                pσ = var2sigmacolumn(σ₁,pressure[1:nw],sig1grid,splorder)

                # again eliminate sgood
                [varsσ["p"][xx,yy,ss] = convert(Float32,pσ[ss]) for ss = 1:nσ]

            end
        end
    end
    return varsσ
end

"""
    function sigmacolumn(θ,S,p,p0)
    σ for a water column
# Arguments
- `θz::Array{Float,1}}`: potential temperature
- `Sz::Array{Float,1}}`: practical salinity
- `pz::Array{Float,1}`: vertical profile of standard pressures
- `p0::Float`: reference pressure
# Output
- `σ`:  sigma for wet points in column
"""
function sigmacolumn(θz,Sz,pz,p0)
    nz = length(θz)
    σ = similar(θz)
    # hard coded for sigma1
    σa,σb,σc = MITgcmTools.SeaWaterDensity(θz,Sz,pz,p0)
    [σ[zz] = convert(Float32,σc[zz]) .- 1000 for zz = 1:nz]
    return σ
end

"""
    function sigma1column(θ,S,p)
    σ₁ for a water column
# Arguments
- `θz::Array{Float,1}}`: potential temperature
- `Sz::Array{Float,1}}`: practical salinity
- `pz::Array{Float,1}`: vertical profile of standard pressures
# Output
- `σ₁`:  sigma-1 for wet points in column
"""
sigma1column(θz,Sz,pz) = sigmacolumn(θz,Sz,pz,1000f0)

"""
    function var2sigmacolumn(σ,v,sig1grid,splorder)
    map θ,S, p onto σ₁ surfaces for a water column
# Arguments
- `σ`::Array{Float,1}}`: sigma values of input variable
- `v::Array{Float,1}}`: variable of interest
- `sig1`: σ surface values
- `splorder`: 1-5, order of spline
# Output
- `θonσ`: variable on sig1 sigma surfaces
"""
function var2sigmacolumn(σorig,v,σgrid,splorder)
    # choose a univariate spline with s = magic number
    #θspl = Spline1D(σ₁,θz;k=splorder,s=length(σ₁))

    linearinterp = false

    σ = copy(σorig) # make sure sigma-1 doesn't mutate and pass back

    nσout = length(σgrid)
    θonσ = fill(NaN,nσout)

    # eliminate homogeneities
    dedup!(σ,v)

    # mix inversions
    mixinversions!(σ,v)

    nσin = length(σ)

    # 1) no inversions or homogeneity (this constraint relaxed now because too many profiles thrown out)
    # 2) range of sig1, 3) no extrapolation
    # if sum(diff(σ₁).<0)==0 && count(minimum(σ₁).<=sig1grid.<=maximum(σ₁)) > 0
    if count(minimum(σ).<=σgrid.<=maximum(σ)) > 0

        # eliminate any extrapolation
        sgood = findall(minimum(σ).<=σgrid.<=maximum(σ))
        ngood = length(sgood)

        if nσin > splorder
            θspl = Spline1D(σ,v;k=splorder)
            prinln("doing spline")
            #println(size(θonσ))
#            println(ngood)
 #           println(size(σgrid))
            for ss in sgood
                θonσ[ss] = θspl(σgrid[ss])
            end

            # check for spline instability
            if maximum(θonσ[sgood]) - minimum(θonσ[sgood]) > 1.05 * (maximum(v) - minimum(v))
                linearinterp = true
                println("unstable spline")
            end
        else
            linearinterp = true
        end # spline interp

        if linearinterp
            interp_linear = LinearInterpolation(σ, v)
            for ss in sgood
                θonσ[ss] = interp_linear(σgrid[ss])
            end
        end # linearinterp
    end # any good points?

    return θonσ
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
    Γ = GridLoad(γ)
    lon2d=[i for i=longrid, j=latgrid]
    lat2d=[j for i=longrid, j=latgrid]
    @time (f,i,j,w)=InterpolationFactors(Γ,vec(lon2d),vec(lat2d))
    return f,i,j,w
end

"""
    function reginterp(fldin,nx,ny,f,i,j,w)
    regular interpolation with precomputed factors
    regular = onto rectangular 2d map
# Arguments
- `fldin`: gcmarray field of input
- `nx,ny`: x,y dimension of output regular field
- `f,i,j,w`: precomputed interpolation factors
# Output
- `fldout`: interpolated regular 2d field
"""
function reginterp(fldin,nx,ny,f,i,j,w)
    fldout = Interpolate(fldin,f,i,j,w)

    # think about moving this before interpolation
    # replace!(fldout,0.0 => NaN)
    fldout = transpose(reshape(fldout,nx,ny))
    return fldout
end

function mixinversions!(a,b)
    while sum(diff(a).<=0) > 0
        length(a) == 1 ? da = 1. : da = diff(a)
        ii = findfirst(isnotpositive,da)
        a[ii] = (a[ii] + a[ii+1])/2
        b[ii] = (b[ii] + b[ii+1])/2
        deleteat!(a,ii+1)
        deleteat!(b,ii+1)
    end
end

function dedup!(a,b)
    length(a) == 1 ? da = 1. : da = diff(a)
    while count(iszero,da) > 0
        dedupfirst!(a,b)
        length(a) ==1 ? da = 1. : da = diff(a)
    end
end

function dedupfirst!(a,b)
    counter = 1
    da = diff(a) # requires length of 2 or more
    ii = findfirst(iszero,da)
    while iszero(da[ii])
        b[ii+1] += b[ii]
        deleteat!(a,ii)
        deleteat!(b,ii)
        deleteat!(da,ii)
        counter += 1
        (length(da) < ii) ? break : nothing
    end
    b[ii] /= counter
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
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)

        # read θ for timestep
        @time θ = γ.read(diagpathexp*Tname,MeshArray(γ,Float32,nz))

        # multiply by the correct part of F matrix
        # be sure to handle all points at once
        # add to existing solution for β
        # equal to matrix multiplication w/ columns (θ) times weighting F
        β += F[2,tt] * θ
    end

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

""" function sigma1grid()
    Standard choice of sigma1 surfaces
"""
function sigma1grid()
    sig1grida = 24:0.05:31
    sig1gridb = 31.02:0.02:33
    sig1grid = vcat(sig1grida,sig1gridb)
    sig1grid = sig1grid[1:3:end]
    return sig1grid
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
"""
function write_vars(vars::Dict{String,MeshArrays.gcmarray{Float32,2,Array{Float32,2}}},fileprefix::String,filesuffix::String)
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

        # get a filename with _on_sigma1
        if fldname == "p"
            fileold = fileprefix*"THETA"*"/"*"THETA"*filesuffixold
        else
            fileold = fileprefix*fldname*"/"*fldname*filesuffixold
        end
        filenew = fileprefix*fldname*"/"*fldname*filesuffixnew

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
    function mdsio2sigma1()
    Take variables in a filelist, read, map to sigma1, write to file.
"""
function mdsio2sigma1(pathin,pathout,fileroots,γ,pstdz,sig1grid,splorder)
    # Read All Variables And Puts Them Into "Vars" Dictionary

    vars = Dict{String,MeshArrays.gcmarray{Float32,2,Array{Float32,2}}}()
    for fileroot in fileroots
        merge!(vars,mdsio2dict(pathin,fileroot,γ))
    end

    # solve for sigma1 on depth levels.
    @time varsσ = vars2sigma1(vars,pstdz,sig1grid,γ,splorder)

    fileprefix = pathout
    # use first filename to get timestamp
    filesuffix = "_on_sigma1"*fileroots[1][14:end]*".data"
    write_vars(varsσ,fileprefix,filesuffix)
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
    filesuffixnew = "_on_sigma1"*ncfilenames["THETA"][end-10:end]
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
function replace!(a::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}},b::Pair)
"""
function replace!(a::MeshArrays.gcmarray{Float32,1,Array{Float32,2}},b::Pair)

    nf = size(a,1)
    for ff = 1:nf
        Base.replace!(a[ff],b)
    end
    return a
end

"""
function replace!(a::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}},b::Pair)
"""
function replace!(a::MeshArrays.gcmarray{Float64,1,Array{Float64,2}},b::Pair)

    nf = size(a,1)
    for ff = 1:nf
        Base.replace!(a[ff],b)
    end
    return a
end

"""
function replace!(a::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}},b::Pair)
"""
function replace!(a::MeshArrays.gcmarray{Float32,2,Array{Float32,2}},b::Pair)

    # tried to recursively call replace! but failed.
    # instead do double nested case.
    nf = size(a,1)
    nrec = size(a,2)
    for ff = 1:nf
        for rr = 1:nrec
            Base.replace!(a[ff,rr],b)
        end
    end
    return a
end

"""
function replace!(a::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}},b::Pair)
"""
function replace!(a::MeshArrays.gcmarray{Float64,2,Array{Float64,2}},b::Pair)

    # tried to recursively call replace! but failed.
    # instead do double nested case.
    nf = size(a,1)
    nrec = size(a,2)
    for ff = 1:nf
        for rr = 1:nrec
            Base.replace!(a[ff,rr],b)
        end
    end
    return a
end

function netcdf2regularpoles(ncfilename,ncvarname,γ,nx,ny,nyarc,farc,iarc,jarc,warc,nyantarc,fantarc,iantarc,jantarc,wantarc)

    vars = netcdf2dict(ncfilename,ncvarname,γ)

    varsregpoles = vars2regularpoles(vars,γ,nx,ny,nyarc,farc,iarc,jarc,warc,nyantarc,fantarc,iantarc,jantarc,wantarc)

end

function mdsio2regularpoles(pathin,filein,γ,nx,ny,nyarc,farc,iarc,jarc,warc,nyantarc,fantarc,iantarc,jantarc,wantarc)

    vars = mdsio2dict(pathin,filein,γ)
    varsregpoles = vars2regularpoles(vars,γ,nx,ny,nyarc,farc,iarc,jarc,warc,nyantarc,fantarc,iantarc,jantarc,wantarc)

end

function vars2regularpoles(vars::Dict{String,MeshArrays.gcmarray{Float32,2,Array{Float32,2}}},γ,nx,ny,nyarc,farc,iarc,jarc,warc,nyantarc,fantarc,iantarc,jantarc,wantarc)

    varsregpoles = Dict{String,Array{Float32,3}}()

    for (varname, varvals) in vars

        # remove contamination from land
        replace!(varvals, 0.0 => NaN)

        nz = size(varvals,2)

        #pre-allocate dict
        varsregpoles[varname] = fill(NaN,(nx,ny,nz))

        for zz = 1:nz
            # get regular grid by cropping
            θcrop =  LLCcropC(varvals[:,zz],γ)

            # interpolate to "LLCregular"
            θarc = reginterp(varvals[:,zz],nx,nyarc,farc,iarc,jarc,warc)
            θantarc = reginterp(varvals[:,zz],nx,nyantarc,fantarc,iantarc,jantarc,wantarc)
            varsregpoles[varname][:,:,zz]=hcat(θantarc',θcrop,θarc')

        end
    end
    return varsregpoles
end

function vars2regularpoles(vars::Dict{String,MeshArrays.gcmarray{Float64,2,Array{Float64,2}}},γ,nx,ny,nyarc,farc,iarc,jarc,warc,nyantarc,fantarc,iantarc,jantarc,wantarc)

    varsregpoles = Dict{String,Array{Float32,3}}()

    for (varname, varvals) in vars

        # remove contamination from land
        replace!(varvals, 0.0 => NaN)

        nz = size(varvals,2)

        #pre-allocate dict
        varsregpoles[varname] = fill(NaN,(nx,ny,nz))

        for zz = 1:nz
            # get regular grid by cropping
            θcrop =  LLCcropC(varvals[:,zz],γ)

            # interpolate to "LLCregular"
            θarc = reginterp(varvals[:,zz],nx,nyarc,farc,iarc,jarc,warc)
            θantarc = reginterp(varvals[:,zz],nx,nyantarc,fantarc,iantarc,jantarc,wantarc)
            varsregpoles[varname][:,:,zz]=hcat(θantarc',θcrop,θarc')

        end
    end
    return varsregpoles
end

function vars2regularpoles(vars::Dict{String,MeshArrays.gcmarray{Float64,1,Array{Float64,2}}},γ,nx,ny,nyarc,farc,iarc,jarc,warc,nyantarc,fantarc,iantarc,jantarc,wantarc)

    varsregpoles = Dict{String,Array{Float32,2}}()
    for (varname, varvals) in vars

        # remove contamination from land
        replace!(varvals, 0.0 => NaN)

        #pre-allocate dict
        varsregpoles[varname] = fill(NaN,(nx,ny))

        # get regular grid by cropping
        θcrop =  LLCcropC(varvals,γ)

        # interpolate to "LLCregular"
        θarc = reginterp(varvals,nx,nyarc,farc,iarc,jarc,warc)
        θantarc = reginterp(varvals,nx,nyantarc,fantarc,iantarc,jantarc,wantarc)
        varsregpoles[varname]=hcat(θantarc',θcrop,θarc')
    end
    return varsregpoles
end

"""
vars 2 regularpoles for netcdf input
"""
function vars2regularpoles(vars::Dict{String,MeshArrays.gcmarray{Float32,1,Array{Float32,2}}},γ,nx,ny,nyarc,farc,iarc,jarc,warc,nyantarc,fantarc,iantarc,jantarc,wantarc)

    varsregpoles = Dict{String,Array{Float32,2}}()
    for (varname, varvals) in vars

        # remove contamination from land
        replace!(varvals, 0.0f0 => NaN32)

        #pre-allocate dict
        varsregpoles[varname] = fill(NaN32,(nx,ny))

        # get regular grid by cropping
        θcrop =  LLCcropC(varvals,γ)

        # interpolate to "LLCregular"
        θarc = reginterp(varvals,nx,nyarc,farc,iarc,jarc,warc)
        θantarc = reginterp(varvals,nx,nyantarc,fantarc,iantarc,jantarc,wantarc)
        varsregpoles[varname]=hcat(θantarc',θcrop,θarc')
    end
    return varsregpoles
end

"""
function writeregularpoles(vars,γ,pathout,filesuffix,filelog,λC,lonatts,ϕC,latatts,z,depthatts)
"""
function writeregularpoles(vars::Dict{String,Array{Float32,3}},γ,pathout,filesuffix,filelog,λC,lonatts,ϕC,latatts,z,depthatts)

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
            "depth",
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

end

"""
    function apply_regional_forcing(forcingfield,latpt,lonpt,latrect,lonrect,dlat,dlon)
# Arguments
- `forcingfield`: space and time field of surface forcing
- `latpt`: latitude grid
- `lonpt`: longitude grid
- `latrect`: tuple of latitude range of central rectangle
- `lonrect`: tuple of longitude range of central rectangle
- `dlon`: width in degrees longitude of East/West fade zones
- 'dlat': width in degrees latitude of North/South fade zones
# Output
- 'regionalforcingfield': space and time field of surface forcing, value of zero inside
designated lat/lon rectangle and fading to 1 outside sponge zone on each edge. This is
because this field ends up being SUBTRACTED from the total forcing
"""
function apply_regional_forcing(forcingfield,latpt,lonpt,latrect,lonrect,dlat,dlon)
    maskmultiplier = ones(size(forcingfield)) #set all to one for starters
    maskmultiplier[inrectangle(latpt,lonpt,latrect,lonrect)] = 0 #set mask to zero inside central rectangle
    zonalslope = 1/dlon
    meridionalslope = 1/dlat
    for xx = #run through longitude
        #one approach: for loop through longitude grid, calculate distance in degrees longitude
        #from both lonrect[1] and lonrect[2]. Will have to deal with date line
        #if distance is <= dlon from East edge of center rectangle (lonrect[2]), set value of
        #maskmultiplier at that longitude for all points infadezoneE to zonalslope*distance (distance in degrees lon)
        #same for West: if distance is <=dlon from West edge of rectangle (lonrect[1]), set value of
        #maskmultiplier at that longitude for all points infadezoneW to zonalslope*distance (distance in degrees lon)

        #this way maybe defeats the purpose of making all those infadezone functions?

    end

    for yy = #run through latitude grid
        #one approach: for loop through latitude grid, calculate distance in degrees latitude
        #from both latrect[1] and latrect[2]
        #if distance is <= dlat from North edge of center rectangle (latrect[2]), set value of
        #maskmultiplier at that latitude for all points infadezoneN to meridionalslope*distance (distance in degrees lat)
        #same for South: if distance is <=dlon from South edge of rectangle (latrect[1]), set value of
        #maskmultiplier at that longitude for all points infadezoneS to meridionalslope*distance (distance in degrees lat)
    end

    #once you have the mask, apply it to the input forcing field
    regionalforcingfield = maskmultiplier.*forcingield

end
