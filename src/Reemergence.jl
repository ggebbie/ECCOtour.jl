module Reemergence
# write new functions and put them in this module.
# add them with text below, or create a new file in "src" and include it.

using Statistics, PyPlot, Distributions, FFTW, LinearAlgebra, StatsBase
using MeshArrays, MITgcmTools, LaTeXStrings

export hanncoeffs, hannsum, hannsum!, hannfilter
export get_filtermatrix, matrixfilter, matrixspray, columnscale!
export seasonal_matrices, position_label, searchdir, setupLLCgrid
export listexperiments, expnames, time_label, faststats, allstats, std, mean

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
    function function mean(x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}},dryval::Float32)
    Compute mean of gcmgrid type using function calls, eliminate redundancy
    Eliminate all values = dryval
    Avoid a name clash with StatsBase and MeshArrays
"""
function mean(x::MeshArrays.gcmarray{Float32,1,Array{Float32,2}},dryval::Float32)
     MeshArrays.mean(x,dryval)
end

end
