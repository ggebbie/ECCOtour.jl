var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = ECCOtour","category":"page"},{"location":"#ECCOtour","page":"Home","title":"ECCOtour","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [ECCOtour]","category":"page"},{"location":"#Base.replace!-Union{Tuple{T}, Tuple{N}, Tuple{Function, MeshArrays.gcmarray{T, N, Matrix{T}}}} where {N, T<:AbstractFloat}","page":"Home","title":"Base.replace!","text":"function replace!(f,a::MeshArrays.gcmarray{T,N,Array{T,2}}}\n\nArguments\n\nf::function = replace a with f(a)\na::MeshArrays.gcmarray{T,N,Array{T,2}}: gcmarray with variable type and time-dimension\n\n\n\n\n\n","category":"method"},{"location":"#Base.replace!-Union{Tuple{T}, Tuple{N}, Tuple{MeshArrays.gcmarray{T, N, Matrix{T}}, Pair}} where {N, T<:AbstractFloat}","page":"Home","title":"Base.replace!","text":"function replace!(a::MeshArrays.gcmarray{T,N,Array{T,2}}},b::Pair) where T\n\nArguments\n\na::MeshArrays.gcmarray{T,N,Array{T,2}}: gcmarray with variable type and time-dimension\nb::Pair: replace a elements of pair 1 with pair 2\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.allstats-Union{Tuple{MeshArrays.gcmarray{T, N, Matrix{T}}}, Tuple{T}, Tuple{N}} where {N, T<:AbstractFloat}","page":"Home","title":"ECCOtour.allstats","text":"function allstats(x)\nCompute all statistics of gcmgrid type using function calls, same as faststats\nWarning: assumes land is zero\n\nInput\n\nx::MeshArrays.gcmarray{T,N,Matrix{T}}: input of gcmarray type\n\nOutput\n\nxbar::T: mean value after discarding dryval\nxmax::T: maximum value\nxmin::T: minimum value\nσx::T: standard deviation\nabsxbar::T: mean of absolute value\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.annualcontrib!-Tuple{Any, Any, Any}","page":"Home","title":"ECCOtour.annualcontrib!","text":"function matrixfilter_annualcontrib!(θout,F,θin)\n\nArguments\n\nθout: filtered values\nF: filter function as a matrix\nθin: input values at all times t (e.g., θ(t))\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.apply_regional_mask!-Tuple{Any, Any}","page":"Home","title":"ECCOtour.apply_regional_mask!","text":"function regional_mask(field,mask)\n\nArguments\n\nfield: input field that is mutated/modified\nmask: multiplicative factor, same spatial grid as field\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.centerlon!-Tuple{Any, Any}","page":"Home","title":"ECCOtour.centerlon!","text":"function centerlon!\nMake a longitudinal coordinate system\ncentered at lonmid.\nMakes handling wraparound easy later on.\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.columnscale!-Tuple{MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}, Matrix{Float32}, MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}}","page":"Home","title":"ECCOtour.columnscale!","text":"function columnscale!(product,M,flux)\n\nproduct::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}\nM::Array{Float32,1}\nflux::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.columnscale!-Tuple{MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}, Vector{Float32}, MeshArrays.gcmarray{Float32, 1, Matrix{Float32}}}","page":"Home","title":"ECCOtour.columnscale!","text":"function columnscale!(product,M,flux)\n\nproduct::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}\nM::Array{Float32,1}\nflux::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.depthlevels-Tuple{Any}","page":"Home","title":"ECCOtour.depthlevels","text":"function depthlevels(γ)\nDepths of ECCO v4r4 grid (positive vals)\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.extract_sst34-NTuple{8, Any}","page":"Home","title":"ECCOtour.extract_sst34","text":"function extract_sst34\nextract by reading multiple files\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.extract_timeseries-NTuple{6, Any}","page":"Home","title":"ECCOtour.extract_timeseries","text":"extract_timeseries(froot,years,γ,xval,yval,fval)\n\nArguments\n\nfroot::String: filename root\nyears::StepRange: iterator for multiple files\nγ::gcmarray: GCM grid from MeshArrays.jl\nxval::Integer: x grid index\nyval::Integer: y grid index\nfval::Integer: face index\n\nOutput\n\ntseries: timeseries\nnseries: nseries = number of timeseries elements in each year\n\n```\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.extract_timeseries-Tuple{MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}, Any, Any, Any}","page":"Home","title":"ECCOtour.extract_timeseries","text":"extract_timeseries(flux,xval,yval,fval)\n\nArguments\n\nflux::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}: input array of arrays\nyears::StepRange: iterator for multiple files\nxval::Integer: x grid index\nyval::Integer: y grid index\nfval::Integer: face index\n\nOutput\n\nseries: timeseries\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.extract_timeseries-Tuple{MeshArrays.gcmarray{Float64, 2, Matrix{Float64}}, Any, Any, Any}","page":"Home","title":"ECCOtour.extract_timeseries","text":"extract_timeseries(flux,xval,yval,fval)\n\nArguments\n\nflux::MeshArrays.gcmarray{Float64,2,Array{Float64,2}}: input array of arrays\nyears::StepRange: iterator for multiple files\nxval::Integer: x grid index\nyval::Integer: y grid index\nfval::Integer: face index\n\nOutput\n\nseries: timeseries\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.factors4regularpoles-Tuple{Any}","page":"Home","title":"ECCOtour.factors4regularpoles","text":"function factors4regularpoles(γ)\nGet interpolation factors for regularpoles grid in one place\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.faststats-Union{Tuple{MeshArrays.gcmarray{T, 1, Matrix{T}}}, Tuple{T}} where T<:AbstractFloat","page":"Home","title":"ECCOtour.faststats","text":"function faststats(x)\nCompute fast statistics of gcmgrid type using function calls, eliminate redundancy\n\nInput\n\nx::MeshArrays.gcmarray{T,1,Matrix{T}}: input of gcmarray type\n\nOutput\n\nxbar::T: mean value after discarding dryval\nxmax::T: maximum value\nxmin::T: minimum value\nσx::T: standard deviation\nabsxbar::T: mean of absolute value\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.get_filtermatrix-Tuple{Any, Any}","page":"Home","title":"ECCOtour.get_filtermatrix","text":"get_filtermatrix(tin,tout)\n\nArguments\n\ntin: input times\ntout: output times\n\nOutput\n\nEout2in: matrix that maps from tout to tin\nFin2out: matrix that maps from tin to tout\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.hanncoeffs-Tuple{Any, Any, Any}","page":"Home","title":"ECCOtour.hanncoeffs","text":"hanncoeffs(tin,tout,L)\n\nArguments\n\ntin: input dataset time (vector valued)\ntout: time where filtered output will be created\nL: filter length \n\nOutput\n\nh: normalized coefficients of Hann(ing) filter\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.hannfilter-Tuple{MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}, Float64, MeshArrays.gcmgrid}","page":"Home","title":"ECCOtour.hannfilter","text":"function hannfilter(θin,tin,tout,L,γ)\n\nArguments\n\nθin::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}\ntin::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}\ntout::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}\nL::Float64\nγ::gcmgrid\n\nOutput\n\nθout\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.hannfilter-Tuple{Vector{Float32}, Any, Any, Any}","page":"Home","title":"ECCOtour.hannfilter","text":"function hannfilter(θin,tin,tout,L)\n\nArguments\n\nθin::Array{Float32,1}: input\ntin\ntout\nL\n\nOutput\n\nθout: output\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.hannsum!-NTuple{5, Any}","page":"Home","title":"ECCOtour.hannsum!","text":"hannsum!(θout,θin,tin,tout,L)\nuse this when a timeseries is available\n\nArguments\n\nθout: output value at time th (e.g., θh(th))\nθin: input value at all times t (e.g., θ(t))\ntin: input dataset time (vector valued)\ntout: time where filtered output will be created\nL: filter length \n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.hannsum!-Tuple{MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}, MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}, StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}, Integer, Float64}","page":"Home","title":"ECCOtour.hannsum!","text":"function hannsum!(θout,θin,tin,tout,toutindex,L::Float64)\n\n Call this version for mesharrays\n\nArguments\n\nθout::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}: output value at time tout (e.g., θout(tout))\nθin::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}: input value at all times t (e.g., θ(t))\ntin::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}: input dataset time (vector valued)\ntout::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}: time where filtered output will be created\ntoutindex::Integer: output times?\nL::Float64: filter length\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.hannsum-NTuple{4, Any}","page":"Home","title":"ECCOtour.hannsum","text":"hannsum(θin,tin,tout,L)\nuse this when a timeseries is available\n\nArguments\n\nθ: input value at all times t (e.g., θ(t))\ntout: time where filtered output will be created\ntin: input dataset time (vector valued)\nL: filter length \n\nOutput\n\nθout: output value at time th (e.g., θh(th))\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.inrectangle-NTuple{4, Any}","page":"Home","title":"ECCOtour.inrectangle","text":"function inrectangle(lonin,latin,lons,lats)\nfind all gridcells within Cartesian rectangle\n\nArguments\n\nlons: tuple of longitude range\nlats: tuple of latitude range\n\nOutput\n\nrectangle: boolean of type gcmarray\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.latlonC-Tuple{Any}","page":"Home","title":"ECCOtour.latlonC","text":"function latlonC(γ)\nLatitude-longitude of ECCOv4r4 \"C\" (tracer) grid\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.latlonG-Tuple{Any}","page":"Home","title":"ECCOtour.latlonG","text":"function latlongG(γ)\nLatitude-longitude of ECCOv4r4 \"G\" (velocity) grid\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.matrixfilter-NTuple{4, Any}","page":"Home","title":"ECCOtour.matrixfilter","text":"matrixfilter(F,froot,years,γ)\nwriting it in a funny way to save computation\nissue with timeseries being read in different files\n\nArguments\n\nF: filter in matrix form\nfroot: filename root\nyears: iterator for multiple files\nγ: GCM grid (meshArray type)\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.matrixspray-NTuple{6, Any}","page":"Home","title":"ECCOtour.matrixspray","text":"matrixspray(F,rmfield,frootin,frootout,years,γ)\n\nwriting it in a funny way to save computation\nissue with timeseries being read in different files\n\nArguments\n\nF: spray operator in matrix form\nrmfield: field to be removed from files\nfrootin: filename root of input field\nfrootout: filename root of output field\nyears: iterator for multiple files\nγ: GCM grid (meshArray type)\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.mdsio2dict-Tuple{Any, Any, Any}","page":"Home","title":"ECCOtour.mdsio2dict","text":"function mdsio2dict(pathin,filein,γ)\n\nRead native (i.e., mdsio) MITgcm files and return a Dictionary.\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.mdsio2regularpoles-NTuple{9, Any}","page":"Home","title":"ECCOtour.mdsio2regularpoles","text":"function mdsio2regularpoles(pathin,filein,γ,nx,ny,nyarc,λarc,nyantarc,λantarc)\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.mdsio2sigma1","page":"Home","title":"ECCOtour.mdsio2sigma1","text":"function mdsio2sigma1\nTake variables in a filelist, read, map to sigma1, write to file.\n\n\n\n\n\n","category":"function"},{"location":"#ECCOtour.nancount-Union{Tuple{MeshArrays.gcmarray{T, N, Matrix{T}}}, Tuple{T}, Tuple{N}} where {N, T<:AbstractFloat}","page":"Home","title":"ECCOtour.nancount","text":"function nancount\n\nCount number of NaN's in gcmgrid type\n\nInput\n\nfield::MeshArrays.gcmarray{T,N,Matrix{T}}: input of gcmarray type, can have multiple records\n\nOutput\n\nnannum::Int: number of NaNs\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.ncwritefromtemplate-Tuple{Dict{String, Array{Float64, 3}}, String, String, String, Any}","page":"Home","title":"ECCOtour.ncwritefromtemplate","text":"function ncwritefromtemplate    This version writes NetCDF output on the regularpoles grid\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.netcdf2dict-Tuple{Any, Any}","page":"Home","title":"ECCOtour.netcdf2dict","text":"netcdf2dict array output (?)\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.netcdf2dict-Tuple{String, String, MeshArrays.gcmgrid}","page":"Home","title":"ECCOtour.netcdf2dict","text":"netcdf2dict gcmarray output\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.netcdf2sigma1-NTuple{7, Any}","page":"Home","title":"ECCOtour.netcdf2sigma1","text":"function netcdf2sigma1\nTake variables in a NetCDF filelist, read, map to sigma1, write to file.\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.patchmean-NTuple{6, Any}","page":"Home","title":"ECCOtour.patchmean","text":" function patchmean(x,area,ϕ,λ,ispatch,iswet)\n get weight for rectangle region.\n\nArguments\n\nx: variable of interest\narea: weighting\nϕ: lat\nλ: lon\nispatch: true in patch of interest\niswet = function that is true if in ocean\n\nOutput\n\nxbar: weighted filtered average of x\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.position_label-Tuple{Any, Any}","page":"Home","title":"ECCOtour.position_label","text":"function position_label(lon,lat)\n\nArguments\n\nlon: longitude\nlat: latitude\n\nOutput\n\nlbl: LaTex string used for figure labels\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.prereginterp-Tuple{Any, Any, Any}","page":"Home","title":"ECCOtour.prereginterp","text":"function prereginterp(latgrid,longrid,γ)\nprepare for regular interpolation\nregular = onto rectangular 2d map\n\nArguments\n\nlatgrid: 1d array of latitude\nlongrid: 1d array of longitude\nγ: GCM grid\n\nOutput\n\nf,i,j,w: interpolation factors\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.pressurelevels-Tuple{Any}","page":"Home","title":"ECCOtour.pressurelevels","text":"function pressurelevels(z)\nStandard pressures of ECCO v4r4 grid\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.read_netcdf-Tuple{String, String, MeshArrays.gcmgrid}","page":"Home","title":"ECCOtour.read_netcdf","text":"read_netcdf(fileName,fldName,mygrid)\n\nRead model output from NetCDF files that are global and convert to MeshArray instance. ```\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.readarea-Tuple{Any}","page":"Home","title":"ECCOtour.readarea","text":"function readarea(γ)\narea of ECCO v4r4 grid\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.reginterp-NTuple{4, Any}","page":"Home","title":"ECCOtour.reginterp","text":"function reginterp(fldin,nx,ny,f,i,j,w)\nregular interpolation with precomputed factors\nregular = onto rectangular 2d map\n\nArguments\n\nfldin: gcmarray field of input\nnx,ny: x,y dimension of output regular field\nλ=(f,i,j,w): precomputed interpolation factors\n\nOutput\n\nfldout: interpolated regular 2d field\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.regional_mask-NTuple{6, Any}","page":"Home","title":"ECCOtour.regional_mask","text":"function regional_mask(latpt,lonpt,latrect,lonrect,dlat,dlon)\n\nArguments\n\nlatpt: latitude grid\nlonpt: longitude grid\nlatrect: tuple of latitude range of central rectangle\nlonrect: tuple of longitude range of central rectangle\ndlon: width in degrees longitude of East/West fade zones\n'dlat': width in degrees latitude of North/South fade zones\n\nOutput\n\n'mask': space and time field of surface forcing, value of zero inside\n\ndesignated lat/lon rectangle and fading to 1 outside sponge zone on each edge. This is because this field ends up being SUBTRACTED from the total forcing\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.remove_climatology-Tuple{Any, Any}","page":"Home","title":"ECCOtour.remove_climatology","text":"function remove_climatology(x,xbar)\n`x` = long monthly timeseries, starting in Jan\n`xbar` = 12-month climatology, starting in Jan\nremove climatology\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.remove_seasonal-Tuple{Any, Any, Any}","page":"Home","title":"ECCOtour.remove_seasonal","text":"function remove_seasonal(x,Ecycle,Fcycle,γ) natively\n`x` = long monthly timeseries, starting in Jan\n\nremove seasonal cycle of this timeseries\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.remove_seasonal-Union{Tuple{T}, Tuple{N}, Tuple{MeshArrays.gcmarray{T, N, Matrix{T}}, Any, Any, Any}} where {N, T<:AbstractFloat}","page":"Home","title":"ECCOtour.remove_seasonal","text":"function remove_seasonal(x,Ecycle,Fcycle,γ) \nremove seasonal cycle of this timeseries on the native gcmgrid\n\nArguments\n\nx::gcmarray{T,N,Array{T,2}} = gcmarray of a long timeseries where time is dimension N\nEcycle::Matrix: matrix that operates on seasonal cycle parameters and return a seasonal cycle timeseries\nFcycle::Matrix: precomputed pseudo-inverse of Ecycle\nγ::gcmgrid: MITgcm grid\n\nOutput\n\nx′::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}: potential density referenced to p₀ minus 1000 kg/m³\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.rotate_uv-Union{Tuple{T}, Tuple{N}, Tuple{MeshArrays.gcmarray{T, N, Matrix{T}}, MeshArrays.gcmarray{T, N, Matrix{T}}, Any}} where {N, T<:AbstractFloat}","page":"Home","title":"ECCOtour.rotate_uv","text":"function rotate_uv(uvel,vvel,G) From Gael Forget, JuliaClimateNotebooks/Transport     3. Convert to Eastward/Northward transport     4. Display Subdomain Arrays (optional)\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.searchdir-Tuple{Any, Any}","page":"Home","title":"ECCOtour.searchdir","text":"function searchdir(path,key)\na useful one-line function\n\nArguments\n\npath: directory to search for file\nkey: expression to search for in path\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.seasonal_matrices","page":"Home","title":"ECCOtour.seasonal_matrices","text":"function seasonal_matrices(fcycle,t,overtones=1)\n\nArguments\n\nfcycle: frequency of seasonal cycle\nt: time\novertones=1: optional argument for number of overtones\n\nOutput\n\nE: matrix that solves E*parameters= seasonal cycle\nF=E†: generalized inverse of E\n\n\n\n\n\n","category":"function"},{"location":"#ECCOtour.sigma-Union{Tuple{T}, Tuple{N}, Tuple{MeshArrays.gcmarray{T, N, Matrix{T}}, MeshArrays.gcmarray{T, N, Matrix{T}}, Vector{Float64}, Int64}} where {N, T<:AbstractFloat}","page":"Home","title":"ECCOtour.sigma","text":"function sigma(θ,S,pz,p₀)\nsigma values from SeaWaterDensity for gcmarrays\n\nArguments\n\nθ::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}: potential temperature\nS::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}: practical salinity\npz::Array{Float64,1} : vertical profile of standard pressures\np₀::Int64 : reference pressure for sigma value\n\nOutput\n\nσ::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}: potential density referenced to p₀ minus 1000 kg/m³\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.time_label-Tuple{Integer}","page":"Home","title":"ECCOtour.time_label","text":"function time_label(index)\nTime label useful in plots. For ECCOv4r4.\n\nInput\n\nindex::Integer: index number for monthly average field, number of months since Jan. 1992\n\nOutput\n\ntlbl: label with month name (abbrev. with 3 letters) and calendar year CE\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.timestamp_monthly_v4r4-Tuple{Any}","page":"Home","title":"ECCOtour.timestamp_monthly_v4r4","text":"function timestamp_monthly_v4r4(t)\nprint year and month given time index\nassuming using ECCOv4r4\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.trend_matrices-Tuple{Any}","page":"Home","title":"ECCOtour.trend_matrices","text":"function trend_matrices(t)\n\nArguments\n\nt: time\n\nOutput\n\nE: matrix that solves E*parameters= timeseries\nF=E†: generalized inverse of E\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.trend_theta!-NTuple{5, Any}","page":"Home","title":"ECCOtour.trend_theta!","text":"function trend_theta!(β,diagpathexp,path_out,tecco,γ,F)\nget linear trend of potential temperature\nfrom monthly-average ECCOv4r4 gcmarray fields\n\nArguments\n\nβ: gcmarray field of trends\ndiagpathexp: where to find files\ntecco: time stamps of monthly output\nγ: GCM grid information\nF: linear estimator of trend\n\nOutput\n\nβ: updated trends\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.trend_theta-NTuple{4, Any}","page":"Home","title":"ECCOtour.trend_theta","text":"function trend_theta(diagpathexp,path_out,tecco,γ,F)\nget linear trend of potential temperature\nfrom monthly-average ECCOv4r4 gcmarray fields\n\nArguments\n\ndiagpathexp: where to find files\ntecco: time stamps of monthly output\nγ: GCM grid information\nF: linear estimator of trend\n\nOutput\n\nβ: updated trends\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.var2regularpoles-NTuple{8, Any}","page":"Home","title":"ECCOtour.var2regularpoles","text":"var2regularpoles    Take one gcmarray in memory, put on regularpoles grid\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.vars2regularpoles-Union{Tuple{T}, Tuple{Dict{String, MeshArrays.gcmarray{T, 1, Matrix{T}}}, Any, Any, Any, Any, Any, Any, Any}} where T<:AbstractFloat","page":"Home","title":"ECCOtour.vars2regularpoles","text":"function vars2regularpoles(vars::Dict{String,MeshArrays.gcmarray{T,1,Matrix{T}}},γ,nx,ny,nyarc,λarc,nyantarc,λantarc)      variables interpolated onto regularpoles grid\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.velocity2center!-NTuple{5, Any}","page":"Home","title":"ECCOtour.velocity2center!","text":"velocity2center!(uC,vC,u,v,G)\n\nFrom Gael Forget, JuliaClimateNotebooks/Transport\n#1. Convert to `Sv` units and mask out land\n2. Interpolate `x/y` transport to grid cell center\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.velocity2center-Tuple{Any, Any, Any}","page":"Home","title":"ECCOtour.velocity2center","text":"velocity2center(u,v,G)\n\nFrom Gael Forget, JuliaClimateNotebooks/Transport\n#1. Convert to `Sv` units and mask out land\n2. Interpolate `x/y` transport to grid cell center\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.wrapdist-Tuple{Any, Any}","page":"Home","title":"ECCOtour.wrapdist","text":"function wrapdist\nLongitudinal distance in degrees\nPassing the date line may be the shortest distance.\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.write_vars-Tuple{Dict{String, MeshArrays.gcmarray{Float32, 2, Matrix{Float32}}}, String, String}","page":"Home","title":"ECCOtour.write_vars","text":"function write_vars    This version writes mdsio output on the native grid    Missing: write the accompanying meta file\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.writeregularpoles-Tuple{Dict{String, Array{Float32, 3}}, Any, Any, Any, Any, Any, Any, Any, Any, Any, Any}","page":"Home","title":"ECCOtour.writeregularpoles","text":"function writeregularpoles(vars,γ,pathout,filesuffix,filelog,λC,lonatts,ϕC,latatts,z,depthatts)\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.writeregularpoles-Tuple{Dict{String, Array{Float64, 3}}, Any, Any, Any, Any, Any, Any, Any, Any, Any, Any}","page":"Home","title":"ECCOtour.writeregularpoles","text":"function writeregularpoles(vars,γ,pathout,filesuffix,filelog,λC,lonatts,ϕC,latatts,z,depthatts)\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.writeregularpoles-Tuple{Dict{String, Matrix{Float32}}, Any, Any, Any, Any, Any, Any, Any, Any, Any, Any}","page":"Home","title":"ECCOtour.writeregularpoles","text":"function writeregularpoles(vars::Dict{String,Array{Float32,2}},γ,pathout,filesuffix,filelog,λC,lonatts,ϕC,latatts,z,depthatts)\n\n\n\n\n\n","category":"method"},{"location":"#ECCOtour.writeregularpoles-Tuple{Dict{String, Matrix{Float64}}, Any, Any, Any, Any, Any, Any, Any, Any, Any, Any}","page":"Home","title":"ECCOtour.writeregularpoles","text":"function writeregularpoles(vars,pathout,filesuffix,filelog,λC,lonatts,ϕC,latatts)\n\n\n\n\n\n","category":"method"},{"location":"#IsopycnalSurfaces.vars2sigma1-Union{Tuple{T}, Tuple{Dict{String, MeshArrays.gcmarray{T, 2, Matrix{T}}}, Vector{Float64}, Vector{Float64}, MeshArrays.gcmgrid, Integer}, Tuple{Dict{String, MeshArrays.gcmarray{T, 2, Matrix{T}}}, Vector{Float64}, Vector{Float64}, MeshArrays.gcmgrid, Integer, Any}} where T<:AbstractFloat","page":"Home","title":"IsopycnalSurfaces.vars2sigma1","text":"function vars2sigma1(vars,p,sig1grid,γ,spline_order)\nmap variables onto sigma1 surfaces for gcmarrays\n\nArguments\n\nvars::Dict{String,MeshArrays.gcmarray{T,N,Array{T,2}}}: dict of gcmarrays\np::Array{Float64,1} : vertical profile of standard pressures\nsig1grid: σ₁ surface values\nγ: grid description needed for preallocation\nsplorder: 1-5, order of spline\nlinearinterp: optional logical\n\nOutput\n\nvarsσ::Dict{String,MeshArrays.gcmarray{T,N,Array{T,2}}}: dict of gcmarrays of variables on sigma1 surfaces\n\n\n\n\n\n","category":"method"},{"location":"#Statistics.mean-Tuple{MeshArrays.gcmarray{Float32, 1, Matrix{Float32}}, MeshArrays.gcmarray{Float64, 1, Matrix{Float64}}, Any}","page":"Home","title":"Statistics.mean","text":"function mean\nCompute area-weighted mean of gcmgrid type using filtered with function isgood\nArea weighting = area\nEliminate all values = isgood(x) -> true\n\nInput\n\nx::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}: input of gcmarray type\nweight::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}: weighting variable of gcmarray type\nisgood::Function: returns true is a value to be used in the mean\n\nOutput\n\nxbar::Float32: mean value (weighted and filtered)\n\n\n\n\n\n","category":"method"},{"location":"#Statistics.mean-Tuple{MeshArrays.gcmarray{Float32, 1, Matrix{Float32}}, MeshArrays.gcmarray{Float64, 1, Matrix{Float64}}, Float64}","page":"Home","title":"Statistics.mean","text":"function mean\nCompute area-weighted mean of gcmgrid type using function calls, eliminate redundancy\nArea weighting = area\nEliminate all values = dryval\n\nInput\n\nx::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}: input of gcmarray type\nweight::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}: weighting variable of gcmarray type\ndryval::Float32: land value (doesn't work for NaN32)\n\nOutput\n\nxbar::Float32: mean value (unweighted)\n\n\n\n\n\n","category":"method"},{"location":"#Statistics.mean-Union{Tuple{T}, Tuple{N}, Tuple{MeshArrays.gcmarray{T, N, Matrix{T}}, T}} where {N, T<:AbstractFloat}","page":"Home","title":"Statistics.mean","text":"function mean\n\nCompute mean of gcmgrid type using function calls, eliminate redundancy\nEliminate all values = dryval\n\nInput\n\nx::MeshArrays.gcmarray{T,N,Array{T,2}}: input of gcmarray type\ndryval::T: land value (doesn't work for NaN32)\n\nOutput\n\nxbar::T: mean value (unweighted)\n\n\n\n\n\n","category":"method"},{"location":"#Statistics.std-Union{Tuple{T}, Tuple{N}, Tuple{MeshArrays.gcmarray{T, N, Matrix{T}}, T, T}} where {N, T<:AbstractFloat}","page":"Home","title":"Statistics.std","text":"function std\nCompute standard deviation of gcmgrid type using function calls, eliminate redundancy\nEliminate all values = dryval\n\nInput\n\nx::MeshArrays.gcmarray{T,N,Array{T,2}}: input of gcmarray type\nxbar::T: mean value\ndryval::T: land value \n\nOutput\n\nσx::T: standard deviation \n\n\n\n\n\n","category":"method"}]
}
