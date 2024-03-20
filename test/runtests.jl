using Revise
using Test
using ECCOtour
using MITgcmTools
using MeshArrays
using Statistics
using Dierckx
using NetCDF
using Downloads
#using GoogleDrive

@testset "ECCOtour.jl" begin
    
    pth = MeshArrays.GRID_LLC90
    global γ = GridSpec("LatLonCap",pth)
    Γ = GridLoad(γ;option="full")
    nf = length(γ.fSize)
    z = depthlevels(γ)
    pstdz = pressurelevels(z)

    projectdir() = ".." #dirname(Base.active_project())
    projectdir(args...) = joinpath(projectdir(), args...)
    datadir() = joinpath(projectdir(),"data")
    datadir(args...) = joinpath(datadir(), args...)
    srcdir() = joinpath(projectdir(),"src")
    srcdir(args...) = joinpath(srcdir(), args...)
    testdir() = joinpath(projectdir(),"test")
    testdir(args...) = joinpath(testdir(), args...)
    
    !ispath(datadir()) && mkdir(datadir())

    @testset "basin mask" begin
        include(testdir("test_basinmask.jl"))
    end

    cd(srcdir())
    # workaround: use a shell script
    #run(`sh $srcdir/download_google_drive.sh`)
    run(`sh download_google_drive.sh`)

    # transport > 40 MB, hits virus scanner, move to WHOI website
    url = "https://www2.whoi.edu/staff/ggebbie/wp-content/uploads/sites/146/2024/03/trsp_3d_set1.0000000732.tar_.gz"
    Downloads.download(url,"trsp_3d_set1.0000000732.tar.gz")

    run(`tar xvzf state_3d_set1.0000000732.tar.gz`)
    run(`tar xvzf trsp_3d_set1.0000000732.tar.gz`)
    
    mv("state_3d_set1.0000000732.data",datadir("state_3d_set1.0000000732.data"),force=true)
    mv("trsp_3d_set1.0000000732.data",datadir("trsp_3d_set1.0000000732.data"),force=true)
    mv("state_3d_set1.0000000732.meta",datadir("state_3d_set1.0000000732.meta"),force=true)
    mv("trsp_3d_set1.0000000732.meta",datadir("trsp_3d_set1.0000000732.meta"),force=true)

    cd(testdir())
    
    ## specific for state
    # the state_3d monthly-average diagnostic output
    stateroot = "state_3d_set1" # 1: θ, 2: S

    # first filter for state_3d_set1
    filelist = searchdir(datadir(),stateroot)
    # second filter for "data"
    statefile  = filter(x -> occursin("data",x),filelist)

    # Read from filelist, map to sigma-1, write to file
    fileroots = Vector{String}()
    fileroot = rstrip(statefile[1],['.','d','a','t','a'])
    push!(fileroots,fileroot)

    ## specific for transport
    transportroot = "trsp_3d_set1" # 1: θ, 2: S

    # first filter for state_3d_set1
    filelist = searchdir(datadir(),transportroot)
    # second filter for "data"
    transportfile  = filter(x -> occursin("data",x),filelist)

    # Read from filelist, map to sigma-1, write to file
    fileroot = rstrip(transportfile[1],['.','d','a','t','a'])
    push!(fileroots,fileroot)

    @testset "IsopycnalSurfaces" begin
        
        # DEFINE THE LIST OF SIGMA1 VALUES.
        sig1grid = sigma1grid()
        eos_mitgcm = "JMD95"

        @testset "spline_interpolation" begin
            @time varsσ = mdsio2sigma1(datadir(),datadir(),fileroots,γ,pstdz,sig1grid,splorder=3,eos=eos_mitgcm)
            for ss in eachindex(sig1grid)
                @test maximum(MeshArrays.mask(varsσ["SALT"][:,ss],-Inf)) < 45.0
                @test minimum(MeshArrays.mask(varsσ["SALT"][:,ss],Inf)) ≥ 0.0
                @test maximum(MeshArrays.mask(varsσ["THETA"][:,ss],-Inf)) < 35.0
                @test minimum(MeshArrays.mask(varsσ["THETA"][:,ss],Inf)) ≥ -3.0
                @test maximum(MeshArrays.mask(varsσ["NVELMASS"][:,ss],-Inf)) < 3.0
                @test minimum(MeshArrays.mask(varsσ["NVELMASS"][:,ss],Inf)) ≥ -3.0
            end

        end

        @testset "linear_interpolation" begin

            varsσ = mdsio2sigma1(datadir(),datadir(),fileroots,γ,pstdz,sig1grid,linearinterp=true,eos=eos_mitgcm)

            for ss in eachindex(sig1grid)
                @test maximum(MeshArrays.mask(varsσ["SALT"][:,ss],-Inf)) < 45.0
                @test minimum(MeshArrays.mask(varsσ["SALT"][:,ss],Inf)) ≥ 0.0
                @test maximum(MeshArrays.mask(varsσ["THETA"][:,ss],-Inf)) < 35.0
                @test minimum(MeshArrays.mask(varsσ["THETA"][:,ss],Inf)) ≥ -3.0
                @test maximum(MeshArrays.mask(varsσ["NVELMASS"][:,ss],-Inf)) < 3.0
                @test minimum(MeshArrays.mask(varsσ["NVELMASS"][:,ss],Inf)) ≥ -3.0
            end
            
        end
    
        @testset "regularpoles" begin
            ####################################
            # Test the mapping onto the regularpoles grid.
            lat,lon = latlon(γ)

            # Set up Cartesian grid for interpolation.
            # Time for a structure.
            rp_params = factors4regularpoles(γ)

            @testset "regularpoles 3d state" begin

                filein = fileroots[1]
                pathin = datadir()

                @time varsregpoles =  regularpoles(pathin,filein,γ,rp_params) #nx,ny,nyarc,λarc,nyantarc,λantarc)

                filesuffix = "suffix.nc"
                pathout = pathin
                filelog = srcdir("available_diagnostics.log")
                gridatts = grid_attributes()
                
                @time write(varsregpoles,
                    rp_params,
                    γ,
                    pathout,
                    filesuffix,
                    filelog,
                    gridatts)
                    
                @test maximum(filter(!isnan,varsregpoles["SALT"])) < 50.
                @test minimum(filter(!isnan,varsregpoles["SALT"])) > 0.

                # read regularpoles to see if single or double precision
                θ = ncread("../dataTHETA/THETAsuffix.nc","THETA")
                @test eltype(θ) == Float32
            end

            @testset "regularpoles 3d transport" begin

                filein = fileroots[2]
                pathin = datadir()

                @time varsregpoles = regularpoles(pathin,filein,γ,rp_params)

                @test maximum(filter(!isnan,varsregpoles["NVELMASS"])) < 2.
                @test minimum(filter(!isnan,varsregpoles["NVELMASS"])) > -2.

            end

            @testset "regularpoles dxc" begin 
                # load centered longitude
                vars = Dict("XC" => γ.read(γ.path*"XC.data",MeshArray(γ,Float64)))
                @time dxc_regpoles = ECCOtour.regularpoles(vars,γ,rp_params) 
                yy = 100
                for xx = 1:50
                    @test isapprox(dxc_regpoles["XC"][xx,yy],rp_params.λC[xx], rtol=1e-6)
                end
            end
        end
            
        @testset "get filtermatrix" begin

            # take biweekly mean. Use triangular filter with break points at:
            Δi14day = 4*14 # grid index range
            nt6hr = 4*365*27 #   length(fluxsample_point)
            i14day = 2.5:Δi14day:nt6hr+56 # goes past end of time by 14 days to be sure

            t6hr_start = 1/8 # 3Z Jan 1 1992
            Δt6hr = 6/24 # 1/4 of a day
            t6hr = range(t6hr_start,step=Δt6hr,length=nt6hr)

            # Given a full-resolution timeseries, what values at the tiepoints best represent the timeseries?
            # Determine tiepoints (in time) where fluxes are adjusted by ECCO in optimization procedure.
            Δt14day = 14 # units: days
            t14day_start = 1/2 # 12Z Jan 1 1992 by inspection of figure
            t14day = range(t14day_start,step=Δt14day,stop=t6hr[end]+14)# goes past the end by one
            nt14day = length(t14day)

            # Values are added to the t14day tiepoints, then linearly interpolated to fill gaps.
            # Careful not to store high-resolution data all at same time.
            # Make a function that takes values at the tiepoints and then makes a full-resolution timeseries.
            @time E14to6,F6to14 = get_filtermatrix(t6hr,t14day)

            daysperyear = 365.25
            fcycle = 1/(daysperyear) # units: day^{-1}
            # for removing seasonal cycle from 14-day averaged timeseries
            Ecycle,Fcycle = seasonal_matrices(fcycle,t14day)

            identity = F6to14*E14to6
            for ii in 1:size(identity,1)
                @test isapprox(identity[ii,ii],1.0,atol=0.01)
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

            @test maximum(MeshArrays.mask(latnan,-Inf)) ≤ 90.0
            @test minimum(MeshArrays.mask(lat,Inf)) ≥ -90.0

            @test !isnan(std(lat,mean(lat,NaN),NaN))

            # replace NaN with zero
            @test iszero(nancount(MeshArrays.mask(latnan,0.0)))

            # does this run?
            faststats(lat)
            
        end
    end

    @testset "landmask" begin

        msk = landmask(γ)

        # more water than land?
        @test sum(msk) < sum(iszero.(msk))

    end
    
end
