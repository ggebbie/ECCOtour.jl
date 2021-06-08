using Reemergence
using Test

@testset "Reemergence.jl" begin
    # Write your tests here.

    path_grid="../inputs/GRID_LLC90/"
    γ = setupLLCgrid(path_grid)

    lat,lon = latlon(γ)
    # Set up Cartesian grid for interpolation.
    ϕGarc,ϕCarc = latgridArctic(γ)
    ϕGantarc,ϕCantarc = latgridAntarctic(γ) 
    ϕGreg,ϕCreg = latgridRegular(γ) 
    λG = -180.0:179.0
    λC = -179.5:179.5

    ϕG = vcat(ϕGantarc,ϕGreg,ϕGarc)
    ϕC = vcat(ϕCantarc,ϕCreg,ϕCarc)

    farc,iarc,jarc,warc = prereginterp(ϕCarc,λC,γ)
    fantarc,iantarc,jantarc,wantarc = prereginterp(ϕCantarc,λC,γ)

    # Fix this
    nx = length(λC)
    ny = length(ϕC)
    nyarc = length(ϕCarc)
    nyantarc = length(ϕCantarc)

    # load centered longitude
    # dxc = γ.read(γ.path*"XC.data",MeshArray(γ,Float64))
    vars = Dict("XC" => γ.read(γ.path*"XC.data",MeshArray(γ,Float64)))
    dxc_regpoles = vars2regularpoles(vars,γ,nx,ny,nyarc,farc,iarc,jarc,warc,nyantarc,fantarc,iantarc,jantarc,wantarc)

    yy = 100
    for xx = 1:nx
        println(@test isapprox(dxc_regpoles["XC"][xx,yy],λC[xx], rtol=1e-6))
    end
    
    #format conversions
    # (γ,Γ)=GridOfOnes("CubeSphere",30,30)
    # @test isa(convert2gcmfaces(Γ["XC"]),Array)
    # @test isa(convert2array(Γ["XC"]),Array)

    # #physical oceanography
    # (ρP,ρI,ρR) = SeaWaterDensity(3.,35.5,3000.)
    # @test isapprox(ρI,1041.83267, rtol=1e-6)

    # D=collect(0.0:1.0:500.0); tmp=(1.0.-tanh.(5*(-1 .+ 2/D[end]*D)));
    # T=2.0 .+ 8.0*tmp; S=34.0 .+ 0.5*tmp;
    # (ρP,ρI,ρR) = SeaWaterDensity(T,S,D);
    # mld=MixedLayerDepth(T,S,D,"BM");
    # @test isapprox(mld,134.0)

end

end
