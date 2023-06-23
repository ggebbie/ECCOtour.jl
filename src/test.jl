export sigma2grid, vars2sigma, mdsio2sigma, mdsio2sigma2
export densityJMD95

""" function sigma2grid()
    Standard from the Hallberg Isopycnal Model
# Arguments
- `z`: value
# Output
- `σ₂ grid`: list (vector) of σ₁ values
"""
function sigma2grid()
    σ₂grida = 24:0.05:31
    σ₂gridb = 31.02:0.02:33
    σ₂grid = vcat(σ₂grida,σ₂gridb)
    σ₂grid = σ₂grid[1:3:end]
    return σ₂grid
end

"""
    function vars2sigma(vars,p,siggrid,γ,spline_order)
    map variables onto sigma1 surfaces for gcmarrays
# Arguments
- `vars::Dict{String,MeshArrays.gcmarray{T,N,Array{T,2}}}`: dict of gcmarrays
- `p::Array{Float64,1}` : vertical profile of standard pressures
- `siggrid`: σ₁ surface values
- `γ`: grid description needed for preallocation
- `splorder`: 1-5, order of spline
- `linearinterp`: optional logical
# Output
- `varsσ::Dict{String,MeshArrays.gcmarray{T,N,Array{T,2}}}`: dict of gcmarrays of variables on sigma1 surfaces
"""
function vars2sigma(vars::Dict{String,MeshArrays.gcmarray{T,2,Matrix{T}}},pressure::Vector{Float64},siggrid::Vector{Float64},
                    p₀::Float64, γ::gcmgrid;splorder=3,linearinterp=false,eos="JMD95") where T<:AbstractFloat 

    # θ and S must exist
    !haskey(vars,"THETA") && error("θ missing")
    !haskey(vars,"SALT") && error("S missing")

    # loop over faces
    nf,nz = size(vars["THETA"])
    nσ = length(siggrid)

    # vcol = Dict with profile/column data
    # pre-allocate each key in vars
    vcol = Dict{String,Vector{T}}() # vars in a column
    varsσ = Dict{String,MeshArrays.gcmarray{T,2,Matrix{T}}}()

    for (key, value) in vars
        vcol[key] = fill(convert(T,NaN),nz)
        varsσ[key] = MeshArray(γ,T,nσ); fill!(varsσ[key],convert(T,NaN))
    end
    # allocate standard pressure by hand.
    varsσ["p"] = MeshArray(γ,T,nσ); fill!(varsσ["p"],convert(T,NaN))

    for ff = 1:nf
        nx,ny = size(vars["THETA"][ff,1])
        for xx = 1:nx
            for yy = 1:ny
                for (vcolname, vcolval) in vars
                    # vcol = Dict with profile/column data
                    vcol[vcolname] = [vcolval[ff,zz][xx,yy] for zz = 1:nz]
                end

                nw = count(notnanorzero,vcol["THETA"]) # number of wet points in column
                nwS = count(notnanorzero,vcol["SALT"])
                if nw != nwS
                    error("T,S zeroes inconsistent")
                end

                if nw > 3 #need >=n+1 points to do order-n interpolation
                    σ₀ =sigmacolumn(vcol["THETA"][1:nw],vcol["SALT"][1:nw],pressure[1:nw],p₀,eos)

                    for (vckey,vcval) in vcol
                        varσ = var2sigmacolumn(σ₀,vcval[1:nw],siggrid,splorder=splorder,linearinterp=linearinterp)
                        [varsσ[vckey][ff,ss][xx,yy] = convert(T,varσ[ss]) for ss = 1:nσ]
                    end

                    # do standard pressure by hand.
                    pσ = var2sigmacolumn(σ₁,pressure[1:nw],siggrid,splorder=splorder,linearinterp=linearinterp)
                    [varsσ["p"][ff,ss][xx,yy] = convert(T,pσ[ss]) for ss = 1:nσ]
                    
                end
            end
        end
    end
    return varsσ
end

mdsio2sigma2(pathin::String,pathout::String,fileroots::Vector{String},
             γ::gcmgrid,pstdz::Vector{Float64},siggrid::Vector{Float64}, 
             p₀::Float64, p₀_string::String;splorder=3,linearinterp=false,eos="JMD95") = mdsio2sigma(pathin,pathout,fileroots,γ,
                                                                                            pstdz,siggrid, 2000, "sigma2";
                                                                                            splorder=splorder,linearinterp=linearinterp,eos=eos) 


"""
    function mdsio2sigma(pathin::String,pathout::String,fileroots::Vector{String},γ,pstdz,siggrid;splorder=3,linearinterp=false,eos="TEOS10") 

    Take variables in a filelist, read, map to sigma1, write to file.

# Arguments
- `pathin::String`: path of input
- `pathout::String`: path of output
- `fileroots::String`: beginning of file names in pathin
- `γ`: grid description needed for preallocation
- `splorder::Integer`: optional keyword argument, default = 3, order of spline
- `linearinterp::Logical`: optional keyword argument, do linear interpolation?, default = false
- `eos::String`: optional key argument for equation of state, default = "JMD95"
"""
function mdsio2sigma(pathin::String,pathout::String,fileroots::Vector{String},γ::gcmgrid,
    pstdz::Vector{Float64},siggrid::Vector{Float64}, p₀::Float64, p₀_string::String;splorder=3,linearinterp=false,eos="JMD95") 
    # Read All Variables And Puts Them Into "Vars" Dictionary

    # ideally would be more generic and Float64 would be covered.
    vars = Dict{String,MeshArrays.gcmarray{Float32,2,Matrix{Float32}}}()
    for fileroot in fileroots
        merge!(vars,mdsio2dict(pathin,fileroot,γ))
    end

    # check for fields on the velocity (staggered) grid
    # and rotate them if necessary
    Γ = GridLoad(γ;option="full")
    rotate_velocity!(vars,Γ)

    # solve for sigma on depth levels.
    @time varsσ = vars2sigma(vars,pstdz,sig2grid,p₀, γ,splorder=splorder,linearinterp=linearinterp,eos=eos)

    fileprefix = pathout
    # use first filename to get timestamp
    filesuffix = "_on_" * p₀_string * fileroots[1][14:end]*".data"
    write_vars(varsσ,fileprefix,filesuffix)

    return varsσ
end

"""
densityJMD95(Θ,Σ,Π,Π0) from MITgcmTools.jl/PhysicalOceanography.jl SeawaterDensity, From Gael Forget

Compute potential density (ρP), in situ density (ρI), and density
referenced to PREF (Π0 in decibars) from potential temperature (Θ in °C),
salinity (Σ in psu) and pressure (Π in decibars) according to the
UNESCO / Jackett & McDougall 1994 equation of state.

Credits: code based on a Matlab implementation by B. Ferron

Reference: https://www.jodc.go.jp/info/ioc_doc/UNESCO_tech/059832eb.pdf

Check value: ρI = `1041.83267kg/m^3` for Θ=`3°Celcius`, Σ=`35psu`, Π=`3000dbar`
```
(ρP,ρI,ρR) = SeaWaterDensity(3.,35.5,3000.)
isapprox(ρI,1041.83267, rtol=1e-6)
```
"""
function densityJMD95(Θ,Σ,Π,Π0=missing)

   #square root salinity
   sqrtΣ= sqrt.(Σ)
   #compute density pure water at atm pressure
   ZR1= ((((6.536332E-9*Θ .-1.120083E-6).*Θ .+1.001685E-4).*Θ
   .-9.095290E-3).*Θ .+6.793952E-2).*Θ .+999.842594
   #seawater density atm pressure
   ZR2= (((5.3875E-9*Θ .-8.2467E-7).*Θ .+7.6438E-5).*Θ
   .-4.0899E-3).*Θ .+0.824493
   ZR3= (-1.6546E-6*Θ .+1.0227E-4).*Θ .-5.72466E-3
   ZR4= 4.8314E-4

   #potential density (referenced to the surface)
   ρP= (ZR4*Σ + ZR3.*sqrtΣ + ZR2).*Σ + ZR1

   #add the compression terms
   ZE = (-3.508914E-8*Θ .-1.248266E-8).*Θ .-2.595994E-6
   ZBW= ( 1.296821E-6*Θ .-5.782165E-9).*Θ .+1.045941E-4
   ZB = ZBW + ZE .* Σ

   ZD = -2.042967E-2
   ZC = (-7.267926E-5*Θ .+2.598241E-3).*Θ .+0.1571896
   ZAW= ((5.939910E-6*Θ .+2.512549E-3).*Θ .-0.1028859).*Θ .-4.721788
   ZA = ( ZD*sqrtΣ + ZC).*Σ + ZAW

   ZB1= (-0.1909078*Θ .+7.390729).*Θ .-55.87545
   ZA1= ((2.326469E-3*Θ .+1.553190).*Θ .-65.00517).*Θ .+1044.077
   ZKW= (((-1.361629E-4*Θ .-1.852732E-2).*Θ .-30.41638).*Θ
   .+2098.925).*Θ .+190925.6
   ZK0= (ZB1.*sqrtΣ + ZA1).*Σ + ZKW

   #in situ density
   ρI = ρP ./ (1.0 .-Π./(ZK0-Π.*(ZA-Π.*ZB)))

   #density referenced to level Π0
   if !ismissing(Π0)
      ρR = ρP ./ (1.0 .-Π0./(ZK0-Π0.*(ZA-Π0.*ZB)))
   else
      ρR = ρP
   end

   return ρP,ρI,ρR
end
