using Revise
using Test, ECCOtour
using MITgcmTools, MeshArrays, Statistics, Dierckx
using GoogleDrive

@testset "ECCOtour.jl" begin

    pth = MeshArrays.GRID_LLC90
    γ = GridSpec("LatLonCap",pth)
    Γ = GridLoad(γ;option="full")
    nf = length(γ.fSize)
    z = depthlevels(γ)
    pstdz = pressurelevels(z)

    projectdir = dirname(Base.active_project())
    datadir = joinpath(projectdir,"data")
    srcdir = joinpath(projectdir,"src")
    
    !isdir(datadir) && mkdir(datadir)

    # download sample data sets

    # state 3d
    url = "https://docs.google.com/uc?export=download&id=1Sst5Y9AUbef1-Vk2ocBgOOiI2kudYRPx"
    filegz = google_download(url,datadir)
    cd(datadir)
    run(`tar xvzf $filegz`)

    # transport 3d # too large doesn't work due to virus scan
    #url = "https://docs.google.com/uc?export=download&id=1KKk8d_1nQFbM9xQjTelCmWTMfK3SA7U5"
    #filegz = google_download(url,datadir)

    cd(datadir)
    # workaround: use a shell script
    run(`sh $srcdir/download_google_drive.sh`)
    run(`tar xvzf trsp_3d_set1.0000000732.tar.gz`)

    ## specific for state
    # the state_3d monthly-average diagnostic output
    stateroot = "state_3d_set1" # 1: θ, 2: S

    # first filter for state_3d_set1
    filelist = searchdir(datadir,stateroot)
    # second filter for "data"
    statefile  = filter(x -> occursin("data",x),filelist)

    # Read from filelist, map to sigma-1, write to file
    fileroots = Vector{String}()
    fileroot = rstrip(statefile[1],['.','d','a','t','a'])
    push!(fileroots,fileroot)

    ## specific for transport
    transportroot = "trsp_3d_set1" # 1: θ, 2: S

    # first filter for state_3d_set1
    filelist = searchdir(datadir,transportroot)
    # second filter for "data"
    transportfile  = filter(x -> occursin("data",x),filelist)

    # Read from filelist, map to sigma-1, write to file
    fileroot = rstrip(transportfile[1],['.','d','a','t','a'])
    push!(fileroots,fileroot)

    @testset "IsopycnalSurfaces" begin
        
        p₀ = 1000.0; # dbar
        # DEFINE THE LIST OF SIGMA1 VALUES.
        sig1grid = sigma1grid()

        @testset "spline_interpolation" begin
            splorder = 3 # spline order
            varsσ = mdsio2sigma1(datadir,datadir,fileroots,γ,pstdz,sig1grid,splorder)
            for ss in eachindex(sig1grid)
                @test maximum(mask(varsσ["SALT"][:,ss],-Inf)) < 45.0
                @test minimum(mask(varsσ["SALT"][:,ss],Inf)) ≥ 0.0
                @test maximum(mask(varsσ["THETA"][:,ss],-Inf)) < 35.0
                @test minimum(mask(varsσ["THETA"][:,ss],Inf)) ≥ -3.0
                @test maximum(mask(varsσ["NVELMASS"][:,ss],-Inf)) < 3.0
                @test minimum(mask(varsσ["NVELMASS"][:,ss],Inf)) ≥ -3.0
            end

        end

        @testset "linear_interpolation" begin

            splorder = 3 # spline order
            linearinterp = true # override spline order with this optional flag
            varsσ = mdsio2sigma1(datadir,datadir,fileroots,γ,pstdz,sig1grid,splorder,linearinterp)

            for ss in eachindex(sig1grid)
                @test maximum(mask(varsσ["SALT"][:,ss],-Inf)) < 45.0
                @test minimum(mask(varsσ["SALT"][:,ss],Inf)) ≥ 0.0
                @test maximum(mask(varsσ["THETA"][:,ss],-Inf)) < 35.0
                @test minimum(mask(varsσ["THETA"][:,ss],Inf)) ≥ -3.0
                @test maximum(mask(varsσ["NVELMASS"][:,ss],-Inf)) < 3.0
                @test minimum(mask(varsσ["NVELMASS"][:,ss],Inf)) ≥ -3.0
            end
            
        end
    
        @testset "regularpoles" begin
            ####################################
            # Test the mapping onto the regularpoles grid.
            lat,lon = latlon(γ)

            # Set up Cartesian grid for interpolation.
            # Time for a structure.
            λC,λG,ϕC,ϕG,nx,ny,nyarc,nyantarc,λarc,λantarc =
                factors4regularpoles(γ)

            @testset "regularpoles 3d state" begin

                filein = fileroots[1]
                pathin = datadir

                @time varsregpoles =  mdsio2regularpoles(pathin,filein,γ,nx,ny,nyarc,λarc,nyantarc,λantarc)

                @test maximum(filter(!isnan,varsregpoles["SALT"])) < 50.
                @test minimum(filter(!isnan,varsregpoles["SALT"])) > 0.

            end

            @testset "regularpoles 3d transport" begin

                filein = fileroots[2]
                pathin = datadir

                @time varsregpoles =  mdsio2regularpoles(pathin,filein,γ,nx,ny,nyarc,λarc,nyantarc,λantarc)

                @test maximum(filter(!isnan,varsregpoles["NVELMASS"])) < 2.
                @test minimum(filter(!isnan,varsregpoles["NVELMASS"])) > -2.

            end

            @testset "regularpoles dxc" begin 
                # load centered longitude
                # dxc = γ.read(γ.path*"XC.data",MeshArray(γ,Float64))
                vars = Dict("XC" => γ.read(γ.path*"XC.data",MeshArray(γ,Float64)))
                dxc_regpoles = vars2regularpoles(vars,γ,nx,ny,nyarc,λarc,nyantarc,λantarc)

                yy = 100
                for xx = 1:nx
                    println(@test isapprox(dxc_regpoles["XC"][xx,yy],λC[xx], rtol=1e-6))
                end
            end
            
        end
        
        @testset "MeshArrays" begin
            ######################################
            # Test the statistics for MeshArrays.
            lat,lon = latlon(γ)
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

            @test maximum(mask(latnan,-Inf)) ≤ 90.0
            @test minimum(mask(lat,Inf)) ≥ -90.0

            @test !isnan(std(lat,mean(lat,NaN),NaN))

            # replace NaN with zero
            @test iszero(nancount(mask(latnan,0.0)))

            # does this run?
            faststats(lat)
            
        end
    end
end
