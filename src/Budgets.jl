"""
Extracts surface wind stress components τx and τy from diagnostic files.

Parameters:
- `diagpath`: Dictionary containing paths to different experiments.
- `expname`: Name of the experiment.
- `fname`: File name of the diagnostic data.
- `γ`: Grid information.

Returns:
- `τx`: Array of τx component.
- `τy`: Array of τy component.
"""
function extract_ocnTAU(diagpath, expname, fnameτ, γ)
    flux_forcing_exps = ["seasonalclimatology", "seasonalclimatology", "climatological_tau", "clim_tau_iter0"]
    #could use 
    #  read_bin(diagpath[expname]*fnameτ,Float32,γ) and pull the last two elements
    if expname ∈ flux_forcing_exps
        @time EXF = γ.read(diagpath[expname]*fnameτ, MeshArray(γ, Float32, 10))
        τx = EXF[:, 9]
        τy = EXF[:, 10]
        return τx, τy
    else
        @time EXF = γ.read(diagpath[expname]*fnameτ, MeshArray(γ, Float32, 15))
        τx = EXF[:, 14]
        τy = EXF[:, 15]
        return τx, τy
    end
end

"""
Extracts Eulerian velocities u, v, and w from diagnostic files.

Parameters:
- `diagpath`: Dictionary containing paths to different experiments.
- `expname`: Name of the experiment.
- `fnameuvw`: File name of the diagnostic data.
- `γ`: Grid information.

Returns:
- `u`: Array of u-component of velocity.
- `v`: Array of v-component of velocity.
- `w`: Array of w-component of velocity.
"""
function extract_eulerian_velocities(diagpath, expname, fnameuvw, γ)
    UV = γ.read(diagpath[expname]*fnameuvw, MeshArray(γ, Float32, 150))
    u = UV[:, 1:50]
    v = UV[:, 51:100]
    w = UV[:, 101:150]
    return u, v, w
end

"""
Extracts Eulerian velocities u, v, and w from diagnostic files.

Parameters:
- `diagpath`: Dictionary containing paths to different experiments.
- `expname`: Name of the experiment.
- `fnameuvw`: File name of the diagnostic data.
- `γ`: Grid information.

Returns:
- `u`: Array of u-component of velocity.
- `v`: Array of v-component of velocity.
- `w`: Array of w-component of velocity.
"""
function extract_eulerian_velocities(diagpath, expname, fnameuvw, γ, level)
    UV = γ.read(diagpath[expname]*fnameuvw, MeshArray(γ, Float32, 150))
    u = UV[:, level]
    v = UV[:, 50 + level]
    w = UV[:, 100 + level]
    return u, v, w
end

"""
Extracts Eulerian and bolus velocities u, v, w, Ub, Vb, and Wb from diagnostic files.

Parameters:
- `diagpath`: Dictionary containing paths to different experiments.
- `expname`: Name of the experiment.
- `fnameuvw`: File name of the diagnostic data.
- `γ`: Grid information.
- `Γ`: Grid information.
- `mskC`: Mask data for cell centers.
- `mskW`: Mask data for west cell interfaces.
- `mskS`: Mask data for south cell interfaces.

Returns:
- `u`: Array of u-component of Eulerian velocity.
- `v`: Array of v-component of Eulerian velocity.
- `w`: Array of w-component of Eulerian velocity.
- `Ub`: Array of u-component of bolus velocity.
- `Vb`: Array of v-component of bolus velocity.
- `Wb`: Array of w-component of bolus velocity.
"""
function extract_eulerian_and_bolus_velocities(
    diagpath, expname, fnameuvw, γ, Γ,
    mskC::MeshArrays.gcmarray{T, 2, Matrix{T}},
    mskW::MeshArrays.gcmarray{T, 2, Matrix{T}},
    mskS::MeshArrays.gcmarray{T, 2, Matrix{T}}
) where T<:Real
    UV = γ.read(diagpath[expname]*fnameuvw, MeshArray(γ, Float32, 250))
    u = UV[:, 1:50]
    v = UV[:, 51:100]
    w = UV[:, 101:150]
    GM_PsiX = UV[:, 151:200]
    GM_PsiY = UV[:, 201:250]
    @time Ub, Vb, Wb = calc_bolus(GM_PsiX, GM_PsiY, Γ, mskC, mskW, mskS)
    return u, v, w, Ub, Vb, Wb
end

"""
Extracts κUθ, κVθ, Uθ, and Vθ from diagnostic files for lateral heat budget.

Parameters:
- `diagpath`: Dictionary containing paths to different experiments.
- `expname`: Name of the experiment.
- `fnameHθ`: File name of the diagnostic data.
- `γ`: Grid information.

Returns:
- `κUθ`: U-component of diffusive heat flux
- `κVθ`: V-component of diffusive heat flux
- `Uθ`:  U-component of advective heat flux
        Uses Eulerian and Bolus velocities. 
- `Vθ`:  U-component of advective heat flux
        Uses Eulerian and Bolus velocities.
"""
function extract_lateral_heatbudget(diagpath, expname, fnameHθ, γ)
    dθλ = γ.read(diagpath[expname]*fnameHθ, MeshArray(γ, Float32, 200))
    κUθ = dθλ[:, 1:50]
    κVθ = dθλ[:, 51:100]
    Uθ = dθλ[:, 101:150]
    Vθ = dθλ[:, 151:200]
    return κUθ, κVθ, Uθ, Vθ
end

"""
Extracts κzθ and wθ from diagnostic files for vertical heat budget.

Parameters:
- `diagpath`: Dictionary containing paths to different experiments.
- `expname`: Name of the experiment.
- `fnameRθ`: File name of the diagnostic data.
- `γ`: Grid information.

Returns:
- `κzθ`: Array of diffusion in the "R" direction.
- `wθ`: Array of wθ component.
        Uses Eulerian and Bolus velocities

"""
function extract_vertical_heatbudget(diagpath, expname, fnameRθ, γ)
    # Vertical convergences
    dθr = γ.read(diagpath[expname]*fnameRθ, MeshArray(γ, Float32, 150))
    wθ = dθr[:, 1:50]
    κzθ = dθr[:, 51:100] .+ dθr[:, 101:150] #implicit plus explicit
    return κzθ, wθ
end

"""
    extract_sθ(expname::String, diagpath::Dict{String, String}, γ::gcmgrid, 
               fnameS::String, fnameθ::String, 
               inv_depths::MeshArrays.gcmarray{T, 1, Matrix{T}}) where T 

Extract the potential temperature (sθ) from the given experiment.

# Arguments
- `expname::String`: The name of the experiment.
- `diagpath::Dict{String, String}`: A dictionary mapping experiment names to their respective diagnostic paths.
- `γ::gcmgrid`: The grid for the LLC90 Grid
- `fnameS::String`: The filename for the surface height (ETAN) data.
- `fnameθ::String`: The filename for the potential temperature (θ) data.
- `inv_depths::MeshArrays.gcmarray{T, 1, Matrix{T}}`: The inverse depths of the grid cells.

# Returns
- `sθ::MeshArrays.gcmarray{T, 1, Matrix{T}}`: The potential temperature at each grid cell, 
adjusted for surface height.

"""
function extract_sθ(expname::String,diagpath::Dict{String, String}, 
    γ::gcmgrid, fnameS::String, fnameθ::String, 
    inv_depths::MeshArrays.gcmarray{T, 1, Matrix{T}}) where T 
    θ = γ.read(diagpath[expname]*fnameθ,MeshArray(γ,Float32,50))
    ETAN = γ.read(diagpath[expname]*fnameS,MeshArray(γ,Float32,1))
    sθ = similar(θ)

    s1 = ETAN .* inv_depths
    s1 .+= 1
    sθ = θ .* s1
    return sθ
end

"""
Calculates the bolus velocities Ub, Vb, and Wb from GM_PsiX, GM_PsiY, and grid information.

Parameters:
- `GM_PsiX`: X-component of GM streamfunction.
- `GM_PsiY`: Y-component of GM streamfunction.
- `Γ`: Grid information.
- `mskC`: Mask data for cell centers.
- `mskW`: Mask data for west cell interfaces.
- `mskS`: Mask data for south cell interfaces.

Returns:
- `Ub`: Array of u-component of bolus velocity.
- `Vb`: Array of v-component of bolus velocity.
- `Wb`: Array of w-component of bolus velocity.
"""
function calc_bolus(
    GM_PsiX::MeshArrays.gcmarray{T, 2, Matrix{T}},
    GM_PsiY::MeshArrays.gcmarray{T, 2, Matrix{T}},
    Γ,
    mskC::MeshArrays.gcmarray{T, 2, Matrix{T}},
    mskW::MeshArrays.gcmarray{T, 2, Matrix{T}},
    mskS::MeshArrays.gcmarray{T, 2, Matrix{T}}
) where {T<:Real}
    nr = length(Γ.RC)
    γ = mskC.grid
    bolusU = MeshArray(γ, T, nr)
    bolusV = MeshArray(γ, T, nr)

    for ff=1:5, k = 1:nr-1
        bolusU.f[ff, k] .= (GM_PsiX.f[ff, k+1] .- GM_PsiX.f[ff, k]) ./ Γ.DRF[k]
        bolusV.f[ff, k] .= (GM_PsiY.f[ff, k+1] .- GM_PsiY.f[ff, k]) ./ Γ.DRF[k]
    end

    for ff=1:5
        bolusU.f[ff, nr] .= -GM_PsiX.f[ff, nr] ./ Γ.DRF[nr]
        bolusV.f[ff, nr] .= -GM_PsiY.f[ff, nr] ./ Γ.DRF[nr]
    end

    # And its vertical part
    # (seems correct, leading to 0 divergence)
    # tmpx and tmpy are the BOLUS transports
    nz_GM = size(GM_PsiX, 2)
    bolusW = MeshArray(γ, T, nz_GM)

    for a in eachindex(GM_PsiX)
        GM_PsiX.f[a] .= GM_PsiX.f[a] .* Γ.DYG.f[a[1]]
        GM_PsiY.f[a] .= GM_PsiY.f[a] .* Γ.DXG.f[a[1]]
    end

    @time calc_UV_conv3D!(GM_PsiX, GM_PsiY, bolusW)

    for a in eachindex(GM_PsiX)
        # Negative for divergence instead of convergence
        # Need to rescale bolus W
        bolusW.f[a] .= -1 .* bolusW.f[a] ./ Γ.RAC.f[a[1]]
        bolusU.f[a] .= bolusU.f[a] .* mskW.f[a]
        bolusV.f[a] .= bolusV.f[a] .* mskS.f[a]
        bolusW.f[a] .= bolusW.f[a] .* mskC.f[a]
    end

    return bolusU, bolusV, bolusW
end

function UVtoTrsp(U::MeshArrays.gcmarray{T, 2, Matrix{T}},V::MeshArrays.gcmarray{T, 2, Matrix{T}},
    G::NamedTuple, γ::gcmgrid) where {T<:Real}
    nz = size(U, 2)
    Utrsp = MeshArray(γ, T, nz)
    Vtrsp = MeshArray(γ, T, nz)

    for a in eachindex(U)
        Utrsp.f[a] .= G.DRF[a[2]]*U.f[a].*G.DYG.f[a[1]]
        Vtrsp.f[a] .= G.DRF[a[2]]*V.f[a].*G.DXG.f[a[1]]
    end
    return Utrsp, Vtrsp
end