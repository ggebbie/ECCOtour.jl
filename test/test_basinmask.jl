maskname =  ["Pacific","South China Sea","East China Sea","Okhotsk Sea","Java Sea","Japan Sea","Timor Sea"]
hemisphere = :both
Lsmooth = 5
southlat = -15
northlat = 15

msk = basin_mask(maskname,γ,southlat=southlat,northlat=northlat)
@test maximum(msk) == 1.0
@test minimum(msk) == 0.0

@test iszero(maximum(nancount(msk)))
@test sum(nancount(MeshArrays.mask(msk,NaN,0.0))) > 0

msk_smooth = basin_mask(maskname,γ,southlat=southlat,northlat=northlat,Lsmooth=Lsmooth)
maximum(nancount(msk_smooth))
# any values between 0 and 1? yes.
@test maximum(MeshArrays.mask(msk,0.0)) > maximum(MeshArrays.mask(msk_smooth,0.0)) # smooth always less than 1
