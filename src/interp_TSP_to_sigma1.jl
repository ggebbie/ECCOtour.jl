"""
    function interp_TSP_to_sigma1(sig1grid,ix,iy,filename,spline_order)
# Arguments
- `sig1grid`: sigma1 grid to interpolate to
- `ix`: longitude index
- 'iy': latitude index
- 'it': month index
- 'spline_order': 1-5, order of spline interpolation
# Output
- 'θonσ': potential temperature on sigma1 grid (vector length of sigma1 grid with NaNs outside of data)
- 'Sonσ': salinity on sigma1 grid
- 'ponσ': pressure on sigma1 grid
"""
function interp_TSP_to_sigma1(sig1grid,ix,iy,filename,spline_order)

    # read ECCO data file into a 4D array: lon, lat, depth, variable type
    @time state = read_mdsio(filename)
    #monthnum = string(it) #some function of it--how to convert month index to the number after the filename?
    #@time state = read_mdsio("state_3d_set1",monthnum,'.data')
    # @time state = read_mdsio("/Users/cora/Desktop/Reemergence.jl/state_3d_set1.0000014604.data")

    θ = state[:,:,:,1]
    S = state[:,:,:,2]

    #then, calculate sigma1 from the T,S,P data
    #Use monthly average temperature/salinity fields on poseidon to compute sigma1

    p₀ = 1000 # reference pressure = 1000 dbar

    # ECCOv4r4 uses approximation for pressure without any horizontal deviations.
    # Can precompute pressure for each depth level.
    ρ₀ = 1029 # from "data" text file in run directory
    g  = 9.81 # from "data"

    # depths of MITgcm levels
    pathLLC = "../inputs/GRID_LLC90/" # need to input grid earlier: see example
    fileZ = "RC" # "R" = radius, "C"= center of grid cell
    Z = read_mdsio(pathLLC,fileZ)
    Z = vec(Z)

    # standard pressures via hydrostatic balance
    Pa2dbar = 1/10000
    p = -ρ₀ .*g .* Pa2dbar .* Z # 10000 to change Pa to dbar

    # filter out zeroes.
    θz = filter(!iszero,θ[ix,iy,:])
    Sz = filter(!iszero,S[ix,iy,:])
    nz = length(θz)
    pz = p[1:nz]
    p₀ = 1000

    # ρI is in-situ density, ρR is density referenced to p0.
    (ρP,ρI,ρR) = SeaWaterDensity(θz,Sz,pz,p₀)
    σ₁ = ρR .- 1000

    extrap_low = sig1grid.<minimum(σ₁) #interpolation extends edge values to ends of grid -- will need to clear these values
    extrap_high = sig1grid.>maximum(σ₁)

    θspl = Spline1D(σ₁,θz;k=spline_order)
    θonσ = θspl(sig1grid)
    #turn extrapolated θ values into NaN
    θonσ[extrap_low] .= NaN;
    θonσ[extrap_high] .= NaN;

    Sspl = Spline1D(σ₁,Sz;k=spline_order)
    Sonσ = Sspl(sig1grid)
    Sonσ[extrap_low] .= NaN;
    Sonσ[extrap_high] .= NaN;

    pspl = Spline1D(σ₁,pz;k=spline_order)
    ponσ = pspl(sig1grid)
    ponσ[extrap_low] .= NaN;
    ponσ[extrap_high] .= NaN;

    return θonσ, Sonσ, ponσ

end
