using Revise
using Test, ECCOtour
using MITgcmTools, MeshArrays, Statistics, Dierckx
using SigmaShift

@testset "ECCOtour.jl" begin

    pth = MeshArrays.GRID_LLC90
    γ = GridSpec("LatLonCap",pth)
    Γ = GridLoad(γ;option="full")

    @testset "SigmaShift" begin

        using GoogleDrive
        
        nf = length(γ.fSize)
        z = depthlevels(γ)
        pstdz = pressurelevels(z)
        p₀ = 1000.0; # dbar

        projectdir = dirname(Base.active_project())
        datadir = joinpath(projectdir,"data")
       
        !isdir(datadir) && mkdir(datadir)

        # download sample data set
        url = "https://docs.google.com/uc?export=download&id=1Sst5Y9AUbef1-Vk2ocBgOOiI2kudYRPx"
        filegz = google_download(url,datadir)
        println(filegz)
        run(`gunzip -f $filegz`)

        # DEFINE THE LIST OF SIGMA1 VALUES.
        sig1grid = sigma1grid()

        ## specific for state
        # the state_3d monthly-average diagnostic output
        TSroot = "state_3d_set1" # 1: θ, 2: S

        splorder = 3 # spline order

        # first filter for state_3d_set1
        filelist = searchdir(datadir,TSroot)
        # second filter for "data"
        datafile  = filter(x -> occursin("data",x),filelist)

        # Read from filelist, map to sigma-1, write to file
        fileroots = Vector{String}()
        fileroot = rstrip(datafile[1],['.','d','a','t','a'])
        push!(fileroots,fileroot)
        println(fileroots)
        varsσ = mdsio2sigma1(datadir,datadir,fileroots,γ,pstdz,sig1grid,splorder)
        @test maximum(varsσ["SALT"],NaN32) < 50.
        @test minimum(varsσ["SALT"],NaN32) ≥ 0.0
        @test maximum(varsσ["THETA"],NaN32) < 35.
        @test minimum(varsσ["THETA"],NaN32) ≥ -2.5
            
    end
    
    @testset "regularpoles" begin
        ###################################################
        # Test the mapping onto the regularpoles grid.
        pth = MeshArrays.GRID_LLC90
        γ = GridSpec("LatLonCap",pth)
        Γ = GridLoad(γ;option="full")
        nf = length(γ.fSize)

        lat,lon = latlon(γ)

        # Set up Cartesian grid for interpolation.
        # Time for a structure.
        λC,λG,ϕC,ϕG,nx,ny,nyarc,nyantarc,λarc,λantarc =
            factors4regularpoles(γ)

        # load centered longitude
        # dxc = γ.read(γ.path*"XC.data",MeshArray(γ,Float64))
        vars = Dict("XC" => γ.read(γ.path*"XC.data",MeshArray(γ,Float64)))
        dxc_regpoles = vars2regularpoles(vars,γ,nx,ny,nyarc,λarc,nyantarc,λantarc)

        yy = 100
        for xx = 1:nx
            println(@test isapprox(dxc_regpoles["XC"][xx,yy],λC[xx], rtol=1e-6))
        end

        @testset "MeshArrays.jl" begin
            ######################################
            # Test the statistics for MeshArrays.
            @test  -360.0 ≤ mean(lat,NaN) ≤ 360.0

            # can it find a NaN?
            latnan = deepcopy(lat)
            latnan[1][30,30] = NaN
            @test sum(nancount(latnan)-nancount(lat)) == 1

            # nancount work with different array size?
            lat2 = MeshArrays.gcmarray(γ,Float64,2)
            lat2[:,1] = latnan;
            lat2[:,2] = latnan;
            @test sum(nancount(lat2)) == 2*sum(nancount(latnan))

            @test maximum(latnan,NaN) ≤ 90.0

            @test minimum(lat,NaN) ≥ -90.0

            @test !isnan(std(lat,mean(lat,NaN),NaN))

            # replace NaN with zero
            @test iszero(nancount(replace!(latnan,NaN=>0.0)))
        end
    end
end
