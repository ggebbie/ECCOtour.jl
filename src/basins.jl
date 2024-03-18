basinlist()=["Pacific","Atlantic","Indian","Arctic","Bering Sea",
                "South China Sea","Gulf of Mexico","Okhotsk Sea",
                "Hudson Bay","Mediterranean Sea","Java Sea","North Sea",
                "Japan Sea", "Timor Sea","East China Sea","Red Sea",
                "Gulf","Baffin Bay","GIN Seas","Barents Sea"]

"""
    function basin_mask(basin_name,γ;hemisphere=nothing,southlat=nothing,northlat=nothing,Lsmooth=nothing)

    Make a mask based on Gael Forget's definitions of ocean basins and sub-basins.
    Note: This mask contains float values. It could be more efficient with a BitArray.

# Arguments
- `basin_name::Union{String,Vector{String}}`: string options are Arctic, Atlantic, Baffin Bay, Barents Sea, Bering Sea,
East China Sea, GIN Seas, Gulf, Gulf of Mexico, Hudson Bay, indian, Japan Sea, Java Sea,
Mediterranean Sea, North Sea, Okhotsk Sea, Pacific, Red Sea, South China Sea, Timor Sea.
-`hemisphere::Symbol`: optional argument with values `:north`, `:south`, and `:both`
-'southlat::Number': optional argument specifying southerly latitude of mask
-'northlat::Number': optional argument specifying northerly latitude of mask
-`Lsmooth::Number`: smoothing lengthscale in grid points
-'southlat::Number': optional argument specifying southerly latitude of mask
-'northlat::Number': optional argument specifying northerly latitude of mask
# Output
- 'mask': space and time field of surface forcing, value of zero inside
designated lat/lon rectangle and fading to 1 outside sponge zone on each edge. This is
because this field ends up being SUBTRACTED from the total forcing
""" 
function basin_mask(basin_name)
    # open MATLAB files created elsewhere
    file = matopen(datadir("basin_grids/GRID_LLC90_"*basin_name))
    for ii in 1:5
        mask[ii] = read(file,basin_name*"_mask"*string(ii))
    end 
    close(file)
    return mask
end
function basin_mask(basin_name::String,γ)
    # this version takes one string and returns one mask.
    pth = MeshArrays.GRID_LLC90
    Γ = GridLoad(γ;option="full")
    basins=read(joinpath(pth,"v4_basin.bin"),MeshArray(γ,Float32))
    basin_list=basinlist()
    basinID=findall(basin_list.==basin_name)[1]
    basinmask=similar(basins)
    for ff in 1:5
        basinmask[ff] .= (basins[ff].==basinID) 
    end
    land2nan!(basinmask,γ)
    return basinmask
end
function basin_mask(basin_names::Vector,γ;hemisphere=nothing,southlat=nothing,northlat=nothing,Lsmooth=nothing)
    # this is the full version, take a vector of basin names, include optional arguments
    pth = MeshArrays.GRID_LLC90
    Γ = GridLoad(γ;option="full")
    basins=read(joinpath(pth,"v4_basin.bin"),MeshArray(γ,Float32))

    mask = 0 * basins # needs NaN on land
    for (ii,nn) in enumerate(basin_names)
        mask += basin_mask(nn,γ)
    end

    if !isnothing(hemisphere)
        apply_hemisphere_mask!(mask,hemisphere,γ)
    end

    if !isnothing(southlat) && !isnothing(northlat)
        apply_latitude_mask!(mask,southlat,northlat,γ)
    end
    
    if !isnothing(Lsmooth)
        mask = smooth(mask,Lsmooth,γ)
    end

    # change NaNs to zeros.
    land2zero!(mask,γ)
    return mask
end

"""
    function land2zero!(msk,γ)

# see tests for different masking from MeshArrays
that could replace this function
"""
function land2zero!(msk,γ)
    land = landmask(γ)
    for ff in eachindex(msk)
        msk[ff][land[ff]] .= zero(eltype(msk))
    end
end

"""
    function smooth(msk::MeshArrays.gcmarray,lengthscale)

    Smooth a gcmarray with a lengthscale of `X` points

    Based off Gael Forget, MeshArrays.jl
"""
function smooth(msk::MeshArrays.gcmarray,X,γ)
    Γ = GridLoad(γ;option="full")
    DXCsm=X*Γ.DXC; DYCsm=X*Γ.DYC;
    #apply smoother
    land2nan!(msk,γ)
    return msk_smooth=MeshArrays.smooth(msk,DXCsm,DYCsm,Γ);
end

"""
    function apply_hemisphere_mask!(mask,hemisphere,γ)

    overlay a hemispheric mask on a previous `mask`
    in-place function
    hemisphere options are `:north`,`:south`, and `:both`

    Note: both has not been tested with this version
"""
function apply_hemisphere_mask!(mask,hemisphere,γ)
    Γ = GridLoad(γ;option="full")
    if hemisphere == :north
        hemisphere_mask = Γ.YC .> 0.0;
    elseif hemisphere == :south #South
        hemisphere_mask = Γ.YC .< 0.0;
    elseif hemisphere == :both #both
        hemisphere_mask = Γ.YC .> 0.0 || Γ.YC .≤ 0.0; #optional argument?
    else
        error("no definition for hemisphere")
    end

    # need loop to assign/mutate mask
    for ff in 1:5
        mask[ff] .*= hemisphere_mask[ff]
    end
end

function apply_latitude_mask!(mask,southlat,northlat,γ)
    Γ = GridLoad(γ;option="full")
    if southlat < northlat
        southcondition = southlat .≤ Γ.YC;
        for ff in 1:5
            mask[ff] .*= southcondition[ff]
        end
        northcondition = Γ.YC .≤ northlat;
        for ff in 1:5
            mask[ff] .*= northcondition[ff]
        end
    else
        error("choose a southerly latitude with a value less than the northerly latitude")
    end
end
