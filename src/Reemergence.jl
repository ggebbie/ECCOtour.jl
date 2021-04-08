module Reemergence
# write new functions and put them in this module.
# add them with text below, or create a new file in "src" and include it.

using Statistics, PyPlot, Distributions, FFTW, LinearAlgebra, StatsBase
using MeshArrays, MITgcmTools, LaTeXStrings, Dierckx, DelimitedFiles

export hanncoeffs, hannsum, hannsum!, hannfilter
export get_filtermatrix, matrixfilter, matrixspray, columnscale!
export seasonal_matrices, position_label, searchdir, setupLLCgrid
export listexperiments, expnames, expsymbols, time_label
export faststats, allstats, std, mean
export inrectangle, isnino34, isnino3, isnino4, isnino12, readlatlon, readarea, patchmean
export nino34mean, nino3mean, nino4mean, nino12mean, extract_sst34
export remove_climatology, remove_seasonal, sigma, TSP2sigma1, all2sigma1
export historicalNino34

include("HannFilter.jl")
include("MatrixFilter.jl")
include("SeasonalCycle.jl")
include("interp_TSP_to_sigma1.jl")

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

"""
    function TSP2sigma1(θ,S,p,sig1grid,spline_order)
    sigma values from SeaWaterDensity for gcmarrays
# Arguments
- `θ::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}`: potential temperature
- `S::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}`: practical salinity
- `p::Array{Float64,1}` : vertical profile of standard pressures
- `sig1grid`: σ₁ surface values
- `γ`: grid description needed for preallocation
- `splorder`: 1-5, order of spline
# Output
- `θσ`: gcmarray of potential temperature on sigma-1 surfaces
- `Sσ`: gcmarray of practical salinity on sigma-1 surfaces
- `pσ`: gcmarray of pressure on sigma-1 surfaces
"""
function TSP2sigma1(θ::MeshArrays.gcmarray{Float32,2,Array{Float32,2}},S::MeshArrays.gcmarray{Float32,2,Array{Float32,2}},p::Array{Float64,1},sig1grid,γ::gcmgrid,splorder)
    # loop over faces
    nf,nz = size(θ)
    nσ = length(sig1grid)

    # pre-allocate θσ,Sσ,pσ
    #σ₁ = similar(θ)
    θσ = MeshArray(γ,Float32,nσ); fill!(θσ,NaN32) 
    Sσ = MeshArray(γ,Float32,nσ); fill!(Sσ,NaN32) 
    pσ = MeshArray(γ,Float32,nσ); fill!(pσ,NaN32) 
    
    for ff = 1:nf
        nx,ny = size(θ[ff,1])
        for xx = 1:nx
            for yy = 1:ny
                θz = [θ[ff,zz][xx,yy] for zz = 1:nz]
                Sz = [S[ff,zz][xx,yy] for zz = 1:nz]

                # also need to filter dry values and to change zz
                # Consider using `isdry` function and dryval in future.
                filter!(!iszero,θz)
                filter!(!iszero,Sz)
                if length(θz) != length(Sz)
                    error("T,S zeroes inconsistent")
                end
                nw = length(θz) # number of wet points in column
                # incurs error if splorder > number of points in column
                if nw > splorder #need at least n+1 points to do order-n interpolation
                    pz = p[1:nw]
                    if length(θz) < 3
                        println("Warming")
                        println(θz)
                        println(splorder)
                        println(nw)
                    end
                    
                    σ₁=sigma1column(θz,Sz,pz)

                    # proceed if no inversions
                    if sum(diff(σ₁).<0)==0
                    # eliminate any extrapolation
                        sgood = findall(minimum(σ₁).<=sig1grid.<=maximum(σ₁))
                        ngood = length(sgood)

                        if ngood > 0 # don't make the call if it doesn't span
                            # enough of sigma-1 space
                            θonσ=var2sigma1column(σ₁,θz,sig1grid[sgood],splorder)
                            [θσ[ff,sgood[ss]][xx,yy] = convert(Float32,θonσ[ss]) for ss = 1:ngood]
                            Sonσ=var2sigma1column(σ₁,Sz,sig1grid[sgood],splorder)
                            [Sσ[ff,sgood[ss]][xx,yy] = convert(Float32,Sonσ[ss]) for ss = 1:ngood]
                            ponσ=var2sigma1column(σ₁,pz,sig1grid[sgood],splorder)
                            [pσ[ff,sgood[ss]][xx,yy] = convert(Float32,ponσ[ss]) for ss = 1:ngood]
                        end
                    end
                end
            end
        end
    end
    return θσ,Sσ,pσ
end

"""
    function all2sigma1(θ,S,p,u,v,w,sig1grid,γ,spline_order)
    sigma values from SeaWaterDensity for gcmarrays
# Arguments
- `θ::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}`: potential temperature
- `S::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}`: practical salinity
- `p::Array{Float64,1}` : vertical profile of standard pressures
- `u::Array{Float64,1}` : vertical profile of zonal velocity
- `v::Array{Float64,1}` : vertical profile of meridional velocity
- `w::Array{Float64,1}` : vertical profile of vertical velocity
- `sig1grid`: σ₁ surface values
- `γ`: grid description needed for preallocation
- `splorder`: 1-5, order of spline
# Output
- `θσ`: gcmarray of potential temperature on sigma-1 surfaces
- `Sσ`: gcmarray of practical salinity on sigma-1 surfaces
- `pσ`: gcmarray of pressure on sigma-1 surfaces
- `uσ`: gcmarray of zonal velocity on sigma-1 surfaces
- `vσ`: gcmarray of meridional velocity on sigma-1 surfaces
- `wσ`: gcmarray of verticl velocity on sigma-1 surfaces
"""
function all2sigma1(θ::MeshArrays.gcmarray{Float32,2,Array{Float32,2}},S::MeshArrays.gcmarray{Float32,2,Array{Float32,2}},p::Array{Float64,1},u::MeshArrays.gcmarray{Float32,2,Array{Float32,2}},v::MeshArrays.gcmarray{Float32,2,Array{Float32,2}},w::MeshArrays.gcmarray{Float32,2,Array{Float32,2}},sig1grid,γ::gcmgrid,splorder)
    # loop over faces
    nf,nz = size(θ)
    nσ = length(sig1grid)

    # pre-allocate θσ,Sσ,pσ
    #σ₁ = similar(θ)
    θσ = MeshArray(γ,Float32,nσ); fill!(θσ,NaN32) 
    Sσ = MeshArray(γ,Float32,nσ); fill!(Sσ,NaN32) 
    pσ = MeshArray(γ,Float32,nσ); fill!(pσ,NaN32) 
    uσ = MeshArray(γ,Float32,nσ); fill!(uσ,NaN32) 
    vσ = MeshArray(γ,Float32,nσ); fill!(vσ,NaN32) 
    wσ = MeshArray(γ,Float32,nσ); fill!(wσ,NaN32) 
    
    for ff = 1:nf
        nx,ny = size(θ[ff,1])
        for xx = 1:nx
            for yy = 1:ny
                θz = [θ[ff,zz][xx,yy] for zz = 1:nz]
                Sz = [S[ff,zz][xx,yy] for zz = 1:nz]
                uz = [u[ff,zz][xx,yy] for zz = 1:nz]
                vz = [v[ff,zz][xx,yy] for zz = 1:nz]
                wz = [w[ff,zz][xx,yy] for zz = 1:nz]

                # also need to filter dry values and to change zz
                # Consider using `isdry` function and dryval in future.
                filter!(!iszero,θz)
                filter!(!iszero,Sz)
                if length(θz) != length(Sz)
                    error("T,S zeroes inconsistent")
                end
                nw = length(θz) # number of wet points in column
                # incurs error if splorder > number of points in column
                if nw > splorder #need at least n+1 points to do order-n interpolation
                    pz = p[1:nw]
                    σ₁=sigma1column(θz,Sz,pz)

                    # proceed if no inversions
                    if sum(diff(σ₁).<0)==0
                        # get velocities
                        uz = uz[1:nw]
                        vz = vz[1:nw]
                        wz = wz[1:nw]

                        # eliminate any extrapolation
                        sgood = findall(minimum(σ₁).<=sig1grid.<=maximum(σ₁))
                        ngood = length(sgood)

                        if ngood > 0 # don't make the call if it doesn't span
                            # enough of sigma-1 space
                            θonσ=var2sigma1column(σ₁,θz,sig1grid[sgood],splorder)
                            [θσ[ff,sgood[ss]][xx,yy] = convert(Float32,θonσ[ss]) for ss = 1:ngood]
                            Sonσ=var2sigma1column(σ₁,Sz,sig1grid[sgood],splorder)
                            [Sσ[ff,sgood[ss]][xx,yy] = convert(Float32,Sonσ[ss]) for ss = 1:ngood]
                            ponσ=var2sigma1column(σ₁,pz,sig1grid[sgood],splorder)
                            [pσ[ff,sgood[ss]][xx,yy] = convert(Float32,ponσ[ss]) for ss = 1:ngood]
                            uonσ=var2sigma1column(σ₁,uz,sig1grid[sgood],splorder)
                            [uσ[ff,sgood[ss]][xx,yy] = convert(Float32,uonσ[ss]) for ss = 1:ngood]
                            vonσ=var2sigma1column(σ₁,vz,sig1grid[sgood],splorder)
                            [vσ[ff,sgood[ss]][xx,yy] = convert(Float32,vonσ[ss]) for ss = 1:ngood]
                            wonσ=var2sigma1column(σ₁,wz,sig1grid[sgood],splorder)
                            [wσ[ff,sgood[ss]][xx,yy] = convert(Float32,wonσ[ss]) for ss = 1:ngood]
                        end
                    end
                end
            end
        end
    end
    return θσ,Sσ,pσ,uσ,vσ,wσ
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
function sigma1column(θz,Sz,pz)
    # p0 hard coded to 1000 dbar for sigma1
    nz = length(θz)
    σ₁ = similar(θz)
    # hard coded for sigma1
    σa,σb,σc = MITgcmTools.SeaWaterDensity(θz,Sz,pz,1000f0)
    [σ₁[zz] = convert(Float32,σc[zz]) .- 1000 for zz = 1:nz]
    return σ₁
end

"""
    function var2sigma1column(σ₁,θ,sig1grid,splorder)
    map θ,S, p onto σ₁ surfaces for a water column
# Arguments
- `σ₁`::Array{Float,1}}`: sigma-1
- `θ::Array{Float,1}}`: variable of interest
- `sig1`: σ₁ surface values
- `splorder`: 1-5, order of spline
# Output
- `θonσ`: variable on sigma-1 surfaces
"""
function var2sigma1column(σ₁,θz,sig1,splorder)
    θspl = Spline1D(σ₁,θz;k=splorder)
    #θspl = Spline1D(σ₁,θz)
    θonσ = θspl(sig1)
    return θonσ
end

# """
#     function interp_TSP_to_sigma1(sig1grid,ix,iy,filename,spline_order)
# # Arguments
# - `sig1grid`: sigma1 grid to interpolate to
# - `ix`: longitude index
# - 'iy': latitude index
# - 'it': month index
# - 'spline_order': 1-5, order of spline interpolation
# # Output
# - 'θonσ': potential temperature on sigma1 grid (vector length of sigma1 grid with NaNs outside of data)
# - 'Sonσ': salinity on sigma1 grid
# - 'ponσ': pressure on sigma1 grid
# """
# function TSP2sigma1(sig1grid,ix,iy,filename,spline_order)

#     p₀ = 1000 # reference pressure = 1000 dbar

#     # ECCOv4r4 uses approximation for pressure without any horizontal deviations.
#     # Can precompute pressure for each depth level.
#     ρ₀ = 1029 # from "data" text file in run directory
#     g  = 9.81 # from "data"

#     # depths of MITgcm levels
#     pathLLC = "../inputs/GRID_LLC90/" # need to input grid earlier: see example
#     fileZ = "RC" # "R" = radius, "C"= center of grid cell
#     Z = read_mdsio(pathLLC,fileZ)
#     Z = vec(Z)

#     # standard pressures via hydrostatic balance
#     Pa2dbar = 1/10000
#     p = -ρ₀ .*g .* Pa2dbar .* Z # 10000 to change Pa to dbar

#     # filter out zeroes.
#     θz = filter(!iszero,θ[ix,iy,:])
#     Sz = filter(!iszero,S[ix,iy,:])
#     nz = length(θz)
#     pz = p[1:nz]
#     p₀ = 1000

# end

end
