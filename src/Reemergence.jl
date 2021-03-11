module Reemergence
# write new functions and put them in this module.
# add them with text below, or create a new file in "src" and include it.

using Statistics, PyPlot, Distributions, FFTW, LinearAlgebra, StatsBase
using MeshArrays, MITgcmTools, LaTeXStrings

export hanncoeffs, hannsum, hannsum!, hannfilter
export get_filtermatrix, matrixfilter, matrixspray, columnscale!
export seasonal_matrices, position_label, searchdir, setupLLCgrid
export listexperiments

include("HannFilter.jl")
include("MatrixFilter.jl")
include("SeasonalCycle.jl")

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

# a useful one-line function
searchdir(path,key) = filter(x->occursin(key,x), readdir(path))

function setupLLCgrid(path_grid)
    http="https://github.com/gaelforget/GRID_LLC90"
    !isdir(path_grid) ? run(`git clone $http $path_grid`) : nothing;
    γ=GridSpec("LatLonCap",path_grid)
    return γ
end

function listexperiments(exppath)
    
    dirlist = searchdir(exppath,"") # all files in directory
    explist  = filter(x -> !occursin("README",x),dirlist) # remove README to get explist
    #runpath = 
    #diagpath = runpath.*"diags/"

    runpath = Dict(explist .=> exppath.*explist.*"/run/")
    diagpath = Dict(explist .=> exppath.*explist.*"/run/diags/")

    # print to screen all the available
    println("Available experiments are")
    for (key,value) in runpath
        println(key)
    end

    return runpath,diagpath
    
end       

end
