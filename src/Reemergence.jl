module Reemergence

# The Reemergence Project has already added these packages.
#import Pkg; Pkg.add("Distributions")
#import Pkg; Pkg.add("StatsBase")
#import Pkg; Pkg.add("LaTeXStrings")
#import Pkg; Pkg.add("FFTW")

using Statistics, PyPlot, Distributions, FFTW, LinearAlgebra, StatsBase
using MeshArrays, MITgcmTools, LaTeXStrings

export hanncoeffs, hannsum, hannsum!, hannfilter
export get_filtermatrix, matrixfilter, matrixspray, columnscale!
export seasonal_matrices, position_label

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

end

       
