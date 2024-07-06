
"""
    cons_offset!(X::MeshArrays.gcmarray{T,2,Matrix{T}}, 
    Δ::MeshArrays.gcmarray{T,1,Matrix{T}}) where T<:AbstractFloat

Adjust the values of the 2D MeshArray `X` by adding the corresponding values from the 1D MeshArray `Δ`.

# Arguments
- `X::MeshArrays.gcmarray{T,2,Matrix{T}}`: A 2D MeshArray that will be adjusted.
- `Δ::MeshArrays.gcmarray{T,1,Matrix{T}}`: A 1D MeshArray containing the offset values to be added to `X`.

# Examples
```julia
"""
function cons_offset!(X::MeshArrays.gcmarray{T,2,Matrix{T}}, 
    Δ::MeshArrays.gcmarray{T,1,Matrix{T}}) where T<:AbstractFloat 
    for a in eachindex(X)
        X.f[a] .+= Δ.f[a[1]]
    end
end

"""
    wet_pts(Γ)

Create a mask for ocean points based on the hFacC field of a given grid Γ.

# Arguments
- `Γ`: The grid data structure containing the hFacC field.

# Output
- Returns a mask with entries equal to 1.0 for ocean points and 0.0 for non-ocean points.
"""

function wet_pts(Γ)
    ocean_mask = similar(Γ.hFacC[:,1] )
    ocean_mask[findall(Γ.hFacC[:,1] .>0.0)].=1.0
    ocean_mask[findall(Γ.hFacC[:,1] .<=0.0)].=0.0
    return Float32.(ocean_mask)
end

"""
get_msk(Γ)

Get the mask arrays for W, S, and C faces from a given grid.

Parameters:
- `Γ`: Grid object containing the hFacW, hFacS, and hFacC fields.

Returns:
- `mskW`: Mask array for W-faces.
- `mskS`: Mask array for S-faces.
- `mskC`: Mask array for C-faces.
"""

function get_msk(Γ)
    mskW = Float32.(Γ.hFacW); mskS = Float32.(Γ.hFacW); 
    mskC = Float32.(Γ.hFacW);
    @inbounds for ff in eachindex(mskC)
        mskW.f[ff][Γ.hFacW.f[ff] .> 0] .= 1.0;
        mskS.f[ff][Γ.hFacS.f[ff] .> 0] .= 1.0;
        mskC.f[ff][Γ.hFacC.f[ff] .> 0] .= 1.0;
    end
    return mskC, mskW, mskS
end

"""
    region_mask(wet_pts, λ, ϕ, latlim, lonlim)

Create a regional mask based on latitude and longitude limits.

# Arguments
- `wet_pts`: Ocean mask obtained from `wet_pts(Γ)`.
- `λ`: Latitude mesh array for each face of the grid.
- `ϕ`: Longitude mesh array for each face of the grid.
- `latlim`: Latitude limits (tuple) defining the region of interest.
- `lonlim`: Longitude limits (tuple) defining the region of interest.

# Output
- Returns a regional mask where entries are 1.0 
inside the specified latitude and longitude limits and 0.0 outside.
"""

function region_mask(wet_pts, λ, ϕ, latlim, lonlim)
    reg_msk = similar(wet_pts)
    for ff in 1:5
        in_lon = lonlim[1] .< λ[ff] .< lonlim[2]
        in_lat = latlim[1] .< ϕ[ff] .< latlim[2]
        reg_msk.f[ff] .= wet_pts.f[ff] .* in_lon .* in_lat
    end
    return reg_msk
end


"""
findlatlon(λ::MeshArray, ϕ::MeshArray, λ₀::T, ϕ₀T::)

Find the latitude and longitude coordinates closest to a given point (λ₀, ϕ₀) in a MeshArray.

Parameters:
- `λ`: Latitude mesharray.
- `ϕ`: Longitude mesharray.
- `λ₀`: Latitude of interest.
- `ϕ₀`: Longitude of interest.

Returns:
- Tuple ((face, level), (cartesian index of minimum within a face), (minimum value))
"""
function findlatlon(λ::MeshArray, ϕ::MeshArray, λ₀::T, ϕ₀::T) where T<:Real
    λ_dist = (λ .- λ₀); λ_dist = map.(x ->x^2, λ_dist);
    ϕ_dist = (ϕ .- ϕ₀); ϕ_dist = map.(x ->x^2, ϕ_dist);
    tot_dist = λ_dist .+ ϕ_dist
    return findmin(tot_dist)
end

"""
findmin(ma::MeshArray)

Find the minimum value within a MeshArray.

Parameters:
- `ma`: The MeshArray object.

Returns:
- Tuple containing ((face, level), (cartesian index of minimum within a face), (minimum value))
"""

function findmin(ma::MeshArray)
    f_list = Any[Inf]
    i_list = Any[Inf]
    val_list = Any[Inf]
    for ijk in eachindex(ma)
        temp_val, temp_idx = Base.findmin(ma[ijk])
        if val_list[1] > temp_val
            f_list[1] = ijk
            i_list[1] = temp_idx
            val_list[1] =  temp_val
        end
        
    end
    return (f_list[1], i_list[1], val_list[1])
end


"""
densityJMD95(Θ,Σ,Π,Π0) from MITgcmTools.jl/PhysicalOceanography.jl SeawaterDensity, From Gael Forget

Calculate the density of seawater using the UNESCO / Jackett & McDougall 1994 equation of state.

Parameters:
- `Θ`: Potential temperature in °C.
- `Σ`: Salinity in psu (practical salinity units).
- `Π`: Pressure in decibars.
- `Π0`: Reference pressure in decibars (optional, defaults to missing).

Returns:
- Tuple containing (ρP: potential density, ρI: in situ density, ρR: density referenced to PREF).
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

   return ρR
end



"""
get_cell_volumes(cell_area, cell_depths)

Calculate the volume of each cell in the area of interest.

Parameters:
- `cell_area`: Surface area of each cell.
- `cell_depths`: Depth of each cell.

Returns:
- `cell_volumes`: The volume of each cell
"""

function get_cell_volumes(cell_area, cell_depths)

    cell_volumes = Float32.(similar(cell_depths))

    # ff is the cartesian index (LLC90_face, depth level)
    for ff in eachindex(cell_volumes)
        f_i = ff[1]
        for k in eachindex(cell_volumes[ff])
            cell_volumes[ff][k] = cell_area[f_i][k] * cell_depths[ff][k]
        end
    end

    return cell_volumes
end


"""
get_cell_thickness(basin_mask, Δz, hFacC)

Calculate the thickness of each cell in the area of interest.

Parameters:
- `basin_mask`: Mask with entries of interest equal to 1.0 else 0.0.
- `Δz`: Possible thickness of the cells at each level.
- `hFacC`: Percentage of possible thickness of each cell.

Returns:
- `cell_depths`: Volume-weighted ocean temperature.
"""

function get_cell_thickness(basin_mask, Δz, hFacC)
    γ= basin_mask.grid
    nz = size(hFacC, 2)
    cell_depths = MeshArray(γ, Float32, nz)
    fill!(cell_depths, 0)
    for ff in eachindex(cell_depths)
        #generate the depths of each tracer cell, done in main function to avoid copying 
        tmp = basin_mask.f[ff[1]] .* (Δz[ff[2]] .* hFacC.f[ff])
        cell_depths.f[ff] .=  Float32.(tmp)
    end

    return cell_depths
end

"""
get_geothermalheating(Γ, γ)

Calculate geothermal heating in the ocean.

Parameters:
- `Γ`: Grid object containing necessary fields.
- `γ`: The grid object for geothermal heating.

Returns:
- `GTF3d`: 3D array of geothermal heating values.
"""

function get_geothermalheating(Γ, γ)
    rho = 1029; cp = 3994; 
    ocean_mask = wet_pts(Γ);
    cell_depths = get_cell_thickness(ocean_mask, Γ.DRF, Γ.hFacC); 

    bathy_mask = deepcopy(Γ.hFacC) #mask indicating bottom points
    for a in eachindex(bathy_mask)
        bathy_mask.f[a][bathy_mask.f[a] .> 0.0] .= 1.0
    end

    for ff in 1:5, k = 1:49
        bathy_mask.f[ff, k] .= bathy_mask.f[ff, k] .- bathy_mask.f[ff, k+1]
    end

    fname = "/batou/eccodrive/files/Version4/Release4/input_init/geothermalFlux.bin";
    GTF = read_bin(fname,Float32,γ)
    
    GTF3d = MeshArray(γ, Float32, 50)
    for a in eachindex(GTF3d)
        ρcpdz = cell_depths.f[a] .* (rho * cp); 
        ρcpdz[ρcpdz .== 0] .= Inf #keep empty cell empty
        GTF3d.f[a] .= (GTF.f[a[1]] .* bathy_mask.f[a]) ./ (ρcpdz) #convert to °C per second
    end

    #returns in terms of θ/s 
    return GTF3d
end

"""
calc_Wconv3D!(wFLD::MeshArrays.gcmarray{T, 2, Matrix{T}}, CONV::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:Real

Calculate convergence of W-component of velocity.

Parameters:
- `wFLD`: W-component velocity field.
- `CONV`: Convergence field for W-component.

Returns:
- Nothing. The `CONV` field is updated with the convergence values.
"""

function calc_W_conv3D!(wFLD::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    CONV::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:Real
    nz = size(wFLD, 2)
    for ff=1:5, k=1:(nz-1)
        CONV.f[ff, k] .= wFLD.f[ff, k+1] .- wFLD.f[ff, k] #in - out 
    end

    for ff=1:5
        CONV.f[ff, nz] .= -wFLD.f[ff, nz] #in - out 
    end
end

"""
calc_UV_conv3D!(uFLD::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    vFLD::MeshArrays.gcmarray{T, 2, Matrix{T}}, CONV::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:Real

Calculate convergence of U and V-components of velocity.

Parameters:
- `uFLD`: U-component velocity field.
- `vFLD`: V-component velocity field.
- `CONV`: Convergence field for U and V-components.

Returns:
- Nothing. The `CONV` field is updated with the convergence values.
"""
function calc_UV_conv3D!(uFLD::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    vFLD::MeshArrays.gcmarray{T, 2, Matrix{T}}, CONV::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:Real
        tmpU, tmpV = exch_UV_llc90(uFLD,vFLD)
        tmpU = T.(tmpU); tmpV = T.(tmpV); 
        for a in eachindex(uFLD.f)
            (s1,s2)=size(uFLD.f[a])
            @inbounds tmpU1=view(tmpU.f[a],1:s1,1:s2)
            @inbounds tmpU2=view(tmpU.f[a],2:s1+1,1:s2)
            @inbounds tmpV1=view(tmpV.f[a],1:s1,1:s2)
            @inbounds tmpV2=view(tmpV.f[a],1:s1,2:s2+1)
            @inbounds CONV.f[a] .= tmpU1-tmpU2+tmpV1-tmpV2
        end
end
    
"""
Construct the exchange grid for U and V velocities between neighboring cells.
The returned fields FLDU, FLDV can be used to compute convergences of volume, heat or salt. 
Parameters:
- `fldU`: U-component velocity field.
- `fldV`: V-component velocity field.

Returns:
- `FLDU`: U-component velocity field after exchange.
- `FLDV`: V-component velocity field after exchange.
"""
function exch_UV_llc90(fldU::MeshArrays.gcmarray{T, 2, Matrix{T}}, fldV::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:Real
    fillval = 0f0
    # Step 1
    s = size.(fldU[:, 1].f)
    nz = size(fldU, 2)
    nf = fldU.grid.nFaces
    s = vcat(s, s[3])  # Always has 5 faces in LLC90
    tp = fldU.grid.class

    FLDU = similar(fldU)
    FLDV = similar(fldV)
    (ovfW, ovfE, ovfS, ovfN, evfW, evfE, evfS, evfN) = MeshArrays.exch_cs_viewfunctions()
    for lvl = 1:nz, a = 1:nf
        FLDU.f[a, lvl] = fill(fillval, s[a][1]+1, s[a][2])
        FLDV.f[a, lvl] = fill(fillval, s[a][1], s[a][2]+1)
        FLDU.f[a, lvl][1:s[a][1], 1:s[a][2]] .= fldU.f[a, lvl]
        FLDV.f[a, lvl][1:s[a][1], 1:s[a][2]] .= fldV.f[a, lvl]

        (jW, jE, jS, jN) = MeshArrays.exch_cs_target(s[a], 1)
        (aW, aE, aS, aN, iW, iE, iS, iN) = MeshArrays.exch_cs_sources(a, s, 1)

        if (!iseven)(a)
            (aE <= nf) && (FLDU.f[a, lvl][jE[1].-1, jE[2].-1] .= ovfE(fldU.f[aE, lvl], iE[1], iE[2]))
            (aN <= nf) && (FLDV.f[a, lvl][jN[1].-1, jN[2].-1] .= ovfN(fldU.f[aN, lvl], iN[1], iN[2]))
        else
            (aE <= nf) && (FLDU.f[a, lvl][jE[1].-1, jE[2].-1] .= evfE(fldV.f[aE, lvl], iE[1], iE[2]))
            (aN <= nf) && (FLDV.f[a, lvl][jN[1].-1, jN[2].-1] .= evfN(fldV.f[aN, lvl], iN[1], iN[2]))
        end
    end
    return FLDU, FLDV
end

"""
Construct the exchange grid for U and V velocities between neighboring cells.
The returned fields FLDU, FLDV can be used to compute convergences of volume, heat or salt. 
Parameters:
- `fldU`: U-component velocity field.
- `fldV`: V-component velocity field.

Returns:
- `FLDU`: U-component velocity field after exchange.
- `FLDV`: V-component velocity field after exchange.
"""
function exch_UV_llc90(fldU::MeshArrays.gcmarray{T, 1, Matrix{T}}, fldV::MeshArrays.gcmarray{T, 1, Matrix{T}}) where T<:Real
    FLDU = similar(fldU)
    FLDV = similar(fldV)
    (ovfW, ovfE, ovfS, ovfN, evfW, evfE, evfS, evfN) = MeshArrays.exch_viewfunctions()
    nf = fldU.grid.nFaces
    fillval = 0.0

    for a in 1:5
        (s1, s2) = size(fldU.f[a])
        FLDU.f[a] = fill(fillval, s1 + 1, s2)
        FLDV.f[a] = fill(fillval, s1, s2 + 1)
        FLDU.f[a][1:s1, 1:s2] .= fldU.f[a]
        FLDV.f[a][1:s1, 1:s2] .= fldV.f[a]

        (jW, jE, jS, jN) = exch_target(s1, s2, a)
        (aW, aE, aS, aN, iW, iE, iS, iN) = exch_sources(a, s1, s2)

        if !iseven(a)
            (aE <= nf) && (FLDU.f[a][jE[1] - 1, jE[2] - 1] .= ovfE(fldU.f[aE], iE[1], iE[2]))
            (aN <= nf) && (FLDV.f[a][jN[1] - 1, jN[2] - 1] .= ovfN(fldU.f[aN], iN[1], iN[2]))
        else
            (aE <= nf) && (FLDU.f[a][jE[1] - 1, jE[2] - 1] .= evfE(fldV.f[aE], iE[1], iE[2]))
            (aN <= nf) && (FLDV.f[a][jN[1] - 1, jN[2] - 1] .= evfN(fldV.f[aN], iN[1], iN[2]))
        end
    end
    return FLDU, FLDV
end

"""
Linearly interpolates C-grid tracers to the velocity grid (uFLD, vFLD)
Parameters:
- `uFLD`: Tracer to be interpolated to U-grid
- `vFLD`: Tracer to be interpolated to V-grid
- `Γ`: Grid information.

Returns:
- `Ugrid`: Tracer interpolated to U-grid
- `Vgrid`: Tracer interpolated to V-grid
"""
function interpolate_to_lateral_faces(uFLD::MeshArrays.gcmarray{T, 1, Matrix{T}}, vFLD::MeshArrays.gcmarray{T, 1, Matrix{T}}, Γ) where T<:Real
    Ugrid = T.(similar(uFLD))
    Vgrid = T.(similar(uFLD))
    DXG = T.(Γ.DXG)
    DYG = T.(Γ.DYG)
    (DXG, DYG) = exch_UV(DXG, DYG)

    (tmpU, tmpV) = exch_UV(uFLD, vFLD)

    for a in eachindex(uFLD)
        (s1, s2) = size(uFLD.f[a])
        tmpU1 = view(tmpU.f[a], 1:s1, 1:s2)
        tmpU2 = view(tmpU.f[a], 2:s1 + 1, 1:s2)
        tmpV1 = view(tmpV.f[a], 1:s1, 1:s2)
        tmpV2 = view(tmpV.f[a], 1:s1, 2:s2 + 1)
        # Setup weights
        tmpDXG1 = view(DXG.f[a], 1:s1, 1:s2)
        tmpDXG2 = view(DXG.f[a], 2:s1 + 1, 1:s2)
        tmpDYG1 = view(DYG.f[a], 1:s1, 1:s2)
        tmpDYG2 = view(DYG.f[a], 1:s1, 2:s2 + 1)

        Ugrid.f[a] .= ((tmpU1 .* tmpDXG1) + (tmpU2 .* tmpDXG2)) ./ (tmpDXG1 .+ tmpDXG2)
        Vgrid.f[a] .= ((tmpV1 .* tmpDYG1) + (tmpV2 .* tmpDYG2)) ./ (tmpDYG1 .+ tmpDYG2)
    end

    return Ugrid, Vgrid
end

"""
Linearly interpolates C-grid tracers to the velocity grid (uFLD, vFLD)
Parameters:
- `uFLD`: Tracer to be interpolated to U-grid
- `vFLD`: Tracer to be interpolated to V-grid
- `Γ`: Grid information.

Returns:
- `Ugrid`: Tracer interpolated to U-grid
- `Vgrid`: Tracer interpolated to V-grid
"""
function interpolate_to_lateral_faces(uFLD::MeshArrays.gcmarray{T, 2, Matrix{T}}, vFLD::MeshArrays.gcmarray{T, 2, Matrix{T}}, Γ) where T<:Real
    Ugrid = T.(similar(uFLD))
    Vgrid = T.(similar(uFLD))
    DXG = T.(Γ.DXG)
    DYG = T.(Γ.DYG)
    (DXG, DYG) = exch_UV(DXG, DYG)

    (tmpU, tmpV) = exch_UV_llc90(uFLD, vFLD)

    for a in eachindex(uFLD)
        (s1, s2) = size(uFLD.f[a])
        tmpU1 = view(tmpU.f[a], 1:s1, 1:s2)
        tmpU2 = view(tmpU.f[a], 2:s1 + 1, 1:s2)
        tmpV1 = view(tmpV.f[a], 1:s1, 1:s2)
        tmpV2 = view(tmpV.f[a], 1:s1, 2:s2 + 1)
        # Setup weights
        tmpDXG1 = view(DXG.f[a[1]], 1:s1, 1:s2)
        tmpDXG2 = view(DXG.f[a[1]], 2:s1 + 1, 1:s2)
        tmpDYG1 = view(DYG.f[a[1]], 1:s1, 1:s2)
        tmpDYG2 = view(DYG.f[a[1]], 1:s1, 2:s2 + 1)

        Ugrid.f[a] .= ((tmpU1 .* tmpDXG1) + (tmpU2 .* tmpDXG2)) ./ (tmpDXG1 .+ tmpDXG2)
        Vgrid.f[a] .= ((tmpV1 .* tmpDYG1) + (tmpV2 .* tmpDYG2)) ./ (tmpDYG1 .+ tmpDYG2)
    end

    return Ugrid, Vgrid
end

nanmean(x) = mean(filter(!isnan,x))
nanmean(x,y) = mapslices(nanmean,x,dims=y)

function interpolate_to_vertical_faces!(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}, 
    ma_interp::MeshArrays.gcmarray{T, 2, Matrix{T}})  where T<:Real
    for k = 1:49 
        ma_interp[:, k+1] = (ma[:, k] .+ ma[:, k+1]) ./2 #linear interpolate to faces
    end
end