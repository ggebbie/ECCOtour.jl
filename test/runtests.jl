using ECCOtour, MeshArrays
using Test

@testset "ECCOtour.jl" begin
    # Write your tests here.

    na = 10
    for ndupes = 1:na-1
        println(ndupes)
        a = sort(randn(na))
        b = randn(na)
        # make duplicates at end
        counter = 0
        while counter < ndupes  
            a[end-counter-1] = a[end-counter]
            counter += 1
        end
        dedup!(a,b)
        println(a)
        println(b)
        println(@test issorted(a))
    end

    path_grid="../inputs/GRID_LLC90/"
    γ = setupLLCgrid(path_grid)

    # repeated text that configures regularpoles grid
lat,lon = latlon(γ)

    # Set up Cartesian grid for interpolation.
    # Time for a structure.
    λC,λG,ϕC,ϕG,nx,ny,nyarc,nyantarc,farc,iarc,jarc,warc,fantarc,iantarc,jantarc,wantarc =
        factors4regularpoles(γ)

    # get standard levels of MITgcm
    z = depthlevels(γ)
    nz = length(z)

    # load centered longitude
    # dxc = γ.read(γ.path*"XC.data",MeshArray(γ,Float64))
    vars = Dict("XC" => γ.read(γ.path*"XC.data",MeshArray(γ,Float64)))
    dxc_regpoles = vars2regularpoles(vars,γ,nx,ny,nyarc,farc,iarc,jarc,warc,nyantarc,fantarc,iantarc,jantarc,wantarc)

    yy = 100
    for xx = 1:nx
        println(@test isapprox(dxc_regpoles["XC"][xx,yy],λC[xx], rtol=1e-6))
    end

end
