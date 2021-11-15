using Revise
using Test, ECCOtour
using MITgcmTools, MeshArrays, Statistics, Dierckx

@testset "ECCOtour.jl" begin
    # Write your tests here.

    @testset "SigmaShift.jl" begin
        #############################
        # Test dedup! function
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

        ################################
        # Idealized mapping onto sigma1
        pz = collect(0.:500.:4000.) # pressure levels
        nz = length(pz)

        θz = collect(range(20,stop=10,length=nz))
        Sz = collect(range(36,stop=35,length=nz))
        ztest = sort(rand(2:8,2))
        σ₁true = sigma1column(θz[ztest],Sz[ztest],pz[ztest])

        sig1grid = range(minimum(σ₁true),stop=maximum(σ₁true),length=20)
        splorder = 3
        σ₁=sigma1column(θz,Sz,pz)
        sgood = findall(minimum(σ₁) .<= sig1grid .<= maximum(σ₁))
        pσ = var2sigmacolumn(σ₁,pz,sig1grid[sgood],splorder)

        @test isapprox(pσ[begin],pz[ztest[begin]])
        @test isapprox(pσ[end],pz[ztest[end]])

        ################################
        # Test the mapping onto sigma1.#
        # expt = "test"
        # diagpath = pwd()
        # path_out = pwd()
        
        # γ = setupLLCgrid("grid/"))
        # nf = length(γ.fSize)

        # # get standard levels of MITgcm
        # z = depthlevels(γ)
        # pstdz = pressurelevels(z)
        # p₀ = 1000.0 ; # dbar

        # # sig1 value of interestn
        # sig1grid = 30.0;

        # TSroot = "state_3d_set1" # 1: θ, 2: S
        # splorder = 100 # spline order

        # # first filter for state_3d_set1
        # filelist = searchdir(diagpath,TSroot)

        # # second filter for "data"
        # datafilelist  = filter(x -> occursin("data",x),filelist)

        # filelist2 = searchdir(diagpath,RProot) 
        # datafilelist2  = filter(x -> occursin("data",x),filelist)

        # # make an output directory for each expteriment
        # !isdir(path_out) ? mkdir(path_out) : nothing;
        # nt = length(datafilelist)
        
        # global tt = 0
        # for datafile in datafilelist
        #     tt += 1

        #     #print timestamp
        #     year,month = timestamp_monthly_v4r4(tt)

        #     # eliminate suffix
        #     fileroot = rstrip(datafile,['.','d','a','t','a'])
        #     fileroot2 = RProot*fileroot[14:end] # a better way?
        #     fileroots = (fileroot,fileroot2)
        
        #     # Read from filelist, map to sigma-1, write to file
        #     mdsio2sigma1(diagpath,path_out,fileroots,γ,pstdz,sig1grid,splorder)
        # end
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
