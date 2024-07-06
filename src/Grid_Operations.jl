function sum(ma::MeshArrays.gcmarray{T, N, Matrix{T}}) where T<:AbstractFloat where N
    tmp = zeros(T, 1)
    for a in eachindex(ma)
        tmp .+= sum(T, ma.f[a])
    end

    return tmp[1]
end

function vertical_sum(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:AbstractFloat
    γ = ma.grid #grid info
    sum_temp = MeshArray(γ,T) #make a flat mesharray
    fill!(sum_temp, zero(T))
    for a in eachindex(ma)
        sum_temp.f[a[1]] .+= ma.f[a] #ijk[1] is the face
    end
    return sum_temp 
end

#must specify dims
flat_sum(x::T; dims = Inf) where T<:Matrix = vec(sum(x, dims = dims))

function lateral_sum(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:AbstractFloat
    nz = size(ma, 2)
    faces_sum = zeros(T, nz)
    for a in eachindex(ma)
        faces_sum[a[2]] += sum(ma.f[a])
    end
    return faces_sum 
end



function zonal_sum(ma::MeshArrays.gcmarray{T, 1, Matrix{T}}) where T<:Real
    temp = zeros(T, 270)

    temp[1:270] .+= flat_sum(ma.f[1], dims = 1)
    temp[1:270] .+= flat_sum(ma.f[2], dims = 1) 

    temp[1:270] .+= flat_sum(reverse(ma.f[4], dims = 1), dims = 2)
    temp[1:270] .+= flat_sum(reverse(ma.f[5], dims = 1), dims = 2) 

    return temp 
end

function zonal_sum(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:Real
    nz = size(ma, 2)
    temp = zeros(T, nz, 270)
    for lvl in 1:nz
        temp[lvl, :] .= zonal_sum(ma[:, lvl])
    end

    return temp 
end

function zonal_average(ma::MeshArrays.gcmarray{T, 1, Matrix{T}},
    weights ::MeshArrays.gcmarray{S, 1, Matrix{S}}) where {T<:Real, S<:Real}

    return zonal_sum(ma .* weights) ./ zonal_sum(weights)
end


function zonal_average(ma::MeshArrays.gcmarray{T, 2, Matrix{T}},
    weights ::MeshArrays.gcmarray{S, 2, Matrix{S}}) where {T<:Real, S<:Real}
    nz = size(ma, 2)

    numerator = zeros(T, nz, 270)
    denom = zeros(T, nz, 270)

    for lvl in 1:nz
        numerator[lvl, :] .= zonal_sum(ma[:, lvl] .* weights[:, lvl])
        denom[lvl, :] .= zonal_sum(weights[:, lvl])
    end
    return numerator ./ denom
end

function zonal_average(ma::MeshArrays.gcmarray{T, 2, Matrix{T}}) where T<:Real
    nz = size(ma, 2)

    numerator = zeros(T, nz, 270)
    denom = zeros(T, nz, 270)

    for lvl in 1:nz
        numerator[lvl, :] .= zonal_sum(ma[:, lvl])
        denom[lvl, :] .= zonal_sum((!iszero).(ma[:, lvl]))
    end
    return temp_num ./ temp_denom
end

function ma_curl(u::MeshArray,v::MeshArray,Γ::NamedTuple)

	uvcurl=similar(Γ.XC)
	fac=exchange(1.0 ./Γ.RAZ,1)
	(U,V)=exchange(u,v,1)
    (DXC,DYC)=exchange(Γ.DXC,Γ.DYC,1)
	[DXC[i].=abs.(DXC[i]) for i in eachindex(U)]
	[DYC[i].=abs.(DYC[i]) for i in eachindex(V)]
    U_DY = U .* DYC #weight by length of 
    V_DX = V .* DXC

	[U_DY.f[i][iszero.(U_DY.f[i])] .= NaN for i in eachindex(U_DY)]
	[V_DX.f[i][iszero.(V_DX.f[i])] .= NaN for i in eachindex(V_DX)]
	for i in eachindex(U)
    ucur=U_DY[i][2:end,:]
    vcur=V_DX[i][:,2:end]        
    tmpcurl=ucur[:,1:end-1]-ucur[:,2:end]
    tmpcurl=tmpcurl-(vcur[1:end-1,:]-vcur[2:end,:])
    tmpcurl=tmpcurl.*fac[i][1:end-1,1:end-1]

		##still needed:
		##- deal with corners
		##- if putCurlOnTpoints

		tmpcurl=1/4*(tmpcurl[1:end-1,2:end]+tmpcurl[1:end-1,1:end-1]+
					tmpcurl[2:end,2:end]+tmpcurl[2:end,1:end-1])

		uvcurl[i]=tmpcurl
	end
	
	return uvcurl
end
