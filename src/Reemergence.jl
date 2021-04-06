module Reemergence
# write new functions and put them in this module.
# add them with text below, or create a new file in "src" and include it.

using Statistics, PyPlot, Distributions, FFTW, LinearAlgebra, StatsBase
using MeshArrays, MITgcmTools, LaTeXStrings

export hanncoeffs, hannsum, hannsum!, hannfilter
export get_filtermatrix, matrixfilter, matrixspray, columnscale!
export seasonal_matrices, position_label, searchdir, setupLLCgrid
export listexperiments, expnames, expsymbols, time_label
export faststats, allstats, std, mean
export inrectangle, isnino34, isnino3, isnino4, isnino12, readlatlon, readarea, patchmean
export nino34mean, nino3mean, nino4mean, nino12mean, extract_sst34
export remove_climatology, remove_seasonal, sigma, sigmacolumn

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
    # watch out for wraparound; didn't code a solution
    inrectangle = minimum(lonrect) <= lonpt <= maximum(lonrect) &&
        minimum(latrect) <= latpt <= maximum(latrect)

end

isnino34(lat,lon) = inrectangle(lat,lon,(-5,5),(-170,-90))
isnino3(lat,lon) = inrectangle(lat,lon,(-5,5),(-150,-90))
isnino4(lat,lon) = inrectangle(lat,lon,(-5,5),(-170,-120))
isnino12(lat,lon) = inrectangle(lat,lon,(-10,0),(-90,-80))

function readlatlon(γ)
    ϕ=γ.read(γ.path*"YC.data",MeshArray(γ,Float64))
    λ=γ.read(γ.path*"XC.data",MeshArray(γ,Float64))
return ϕ,λ
end

function readarea(γ)
    area=γ.read(γ.path*"RAC.data",MeshArray(γ,Float64))
return area
end


# get weight for rectangle region.
function patchmean(x,area,ϕ,λ,ispatch,iswet)

    # x = variable
    # area = weighting
    # ϕ = lat, λ = lon
    # ispatch = true in patch of interest
    # iswet = function that is true if in ocean

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

    # print to screen all the available
    println("Available experiments are")
    for (key,value) in runpath
        println(key)
    end

    return runpath,diagpath
    
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

    dryval = 0.f0 # land points are zeroes, use NaN32 for NaN's

    isdry(z) = (z == dryval)
    xcount = [sum(count(!isdry,x[i])) for i in eachindex(x)]
    if sum(xcount)>0
        xmax = maximum([maximum(filter(!iszero,x[i])) for i ∈ eachindex(x) if xcount[i] > 0])
        xmin = minimum([minimum(filter(!iszero,x[i])) for i ∈ eachindex(x) if xcount[i] > 0])

        # compute mean the old fashioned way
        xsum = sum([sum(filter(!isdry,x[i])) for i ∈ eachindex(x) if xcount[i] > 0]) 
        xbar = xsum/sum(xcount)

        # compute standard deviation
        x′ = x.-xbar
        x²sum = sum([sum(filter(!isdry,x′[i]).^2) for i ∈ eachindex(x′) if xcount[i]>0])
        σx = sqrt(x²sum/(sum(xcount)-1))

        absxsum = sum([sum(abs.(x[i])) for i ∈ eachindex(x) if xcount[i] > 0]) # works b.c. 0 on land
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

# """
#     sigma values from SeaWaterDensity for gcmarrays
# """
# function sigmacolumn(θ::Array{Float32,1},S::Array{Float32,1},pz::Array{Float64,1},p₀::Int64)

#     (ρP,ρI,ρR) = SeaWaterDensity(θ,S,pz,p₀)
    
#     σ = ρR .- 1000
#     return σ
# end

    
end
