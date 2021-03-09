# interactive analysis of divergence of runs, originally 10-Feb-2021
# recreate MATLAB difference plots
# use files directly on poseidon through cifs mount

using Revise # not sure if Revise is necessary
using MITgcmTools, MeshArrays, Statistics, Reemergence

# get MITgcm / ECCOv4r4 LLC grid and depth information. Store in γ.
path_grid="../inputs/GRID_LLC90/"
γ=GridSpec("LatLonCap",path_grid)
http="https://github.com/gaelforget/GRID_LLC90"
!isdir(path_grid) ? run(`git clone $http $path_grid`) : nothing;
D=γ.read(γ.path*"Depth.data",MeshArray(γ,Float64))
tmp1=write(D); tmp2=read(tmp1,D)
show(D)

# list of experiments on poseidon
poseidonoutdir = "/poseidon/ecco/ECCOv4r4/MITgcm/exps/"
dirlist = searchdir(poseidonoutdir,"") # all files in directory
explist  = filter(x -> !occursin("README",x),dirlist) # remove README to get explist
pathlist = poseidonoutdir.*explist.*"/run/"

# old path names
# path_unadjusted = "/poseidon/ecco/ECCOv4r4/MITgcm/exps/release4_sfcunadjusted/run/"
# path_standard = "/poseidon/ecco/ECCOv4r4/MITgcm/exps/release4/run/"
# path_ff = "/poseidon/ecco/ECCOv4r4/MITgcm/exps/release4_flux_forced/run/"
# path_unadjusted_diags = "/poseidon/ecco/ECCOv4r4/MITgcm/ECCOV4/release4_sfcunadjusted/run/diags/"
# path_standard_diags = "/poseidon/ecco/ECCOv4r4/MITgcm/ECCOV4/release4/run/diags/"
# path_ff_diags = "/poseidon/ecco/ECCOv4r4/MITgcm/ECCOV4/release4_flux_forced/run/diags/"


    # what kind of averages are these?
    tsnap = 732:732:227904
    years = 1992:2017
    length(tsnap)/length(years) # 12 = monthly, 30.5 days

    pathLLC = "/home/gebbie/julia/inputs/GRID_LLC90/"
    fileZ = "RC"
    Z = read_mdsio(pathLLC,fileZ)
    Z = vec(Z)

    tt = tsnap[1]
    xstandard = read_mdsio(path_standard_diags,filename)
    nx,ny,nz,nc = size(xstandard)

# compare your choice of experiments.
# use numbers for now. can improve code late.
expcompare = (2,3)
for expno ∈ expcompare
    println(expno)

fileroot = "state_3d_set1"
filelist = searchdir(path_unadjusted_diags,fileroot) # first filter for state_3d_set1
filelist2  = filter(x -> occursin("data",x),filelist) # second filter for "data"
nt = length(filelist2)

σstd0 = NaN.*zeros(nt,nz,nc)
σmed0 = NaN.*zeros(nt,nz,nc)
σmax0 = NaN.*zeros(nt,nz,nc)
σmean0 = NaN.*zeros(nt,nz,nc)

σstdff = NaN.*zeros(nt,nz,nc)
σmedff = NaN.*zeros(nt,nz,nc)
σmaxff = NaN.*zeros(nt,nz,nc)
σmeanff = NaN.*zeros(nt,nz,nc)

tcount = 0
for fname in filelist2
    tcount += 1
    println(tcount)
    #numberzeros = 10 - length(string(tt))
    #filename = "state_3d_set1."*"0"^numberzeros*string(tt)
    #filename = filelist2[tcount]
    x0 = read_mdsio(path_unadjusted_diags,fname)
    x  = read_mdsio(path_standard_diags,fname)
    xff = read_mdsio(path_ff_diags,fname)

    # get statistics

    for zz = 1:nz
        #println(zz)
        for cc = 1:2 # cut off vertical density gradient
            tmp0 = x[:,:,zz,cc]-x0[:,:,zz,cc]
            tmpff = x[:,:,zz,cc]-xff[:,:,zz,cc]
            igood = findall(!iszero,tmp0)
            #println(length(igood))
            if length(igood) > 0
                σstd0[tcount,zz,cc] = std(tmp0[igood])
                σmed0[tcount,zz,cc] = median(abs.(tmp0[igood]))
                σmax0[tcount,zz,cc] = maximum(abs.(tmp0[igood]))
                σmean0[tcount,zz,cc] = mean(abs.(tmp0[igood]))

                σstdff[tcount,zz,cc] = std(tmpff[igood])
                σmedff[tcount,zz,cc] = median(abs.(tmpff[igood]))
                σmaxff[tcount,zz,cc] = maximum(abs.(tmpff[igood]))
                σmeanff[tcount,zz,cc] = mean(abs.(tmpff[igood]))
            end
        end
    end
end

zlev = 20
cval = 1
caseff = true
case0 = false
if cval == 1 && case0
    titlelbl = L"\theta_{129} - \theta_{ff}, "*string(floor(Int,-Z[zlev]))*"m"
    savefname = "dtheta_129vff_1992-2017_299m_nomax_12feb2021.eps"
    ylbl = L"\theta  [^{\circ}C]"
end
if cval == 2 && case0
    titlelbl = L"S_{129} - S_{ff}, "*string(floor(Int,-Z[zlev]))*"m"
    savefname = "dsalt_129vff_1992-2017_299m_nomax_12feb2021.eps"
    ylbl = L"salinity  [PSS-1978]"
end
if cval == 1 && caseff
    titlelbl = L"\theta_{129} - \theta_{ff}, "*string(floor(Int,-Z[zlev]))*"m"
    savefname = "dtheta_129vff_1992-2017_299m_nomax_12feb2021.eps"
    ylbl = L"\theta  [^{\circ}C]"
end
if cval == 2 && caseff
    titlelbl = L"S_{129} - S_{ff}, "*string(floor(Int,-Z[zlev]))*"m"
    savefname = "dsalt_129vff_1992-2017_299m_nomax_12feb2021.eps"
    ylbl = L"salinity  [PSS-1978]"
end

#figure()
if case0
    clf()
    title(titlelbl)
    plot(σstd0[:,zlev,cval],label=L"std()")
    plot(σmed0[:,zlev,cval],label=L"median(||)")
    plot(σmean0[:,zlev,cval],label=L"mean(||)")
    #plot(σmax0[:,zlev,cval],label=L"max(||)")
    xlabel("month starting in Jan. 1992")
    ylabel(ylbl)
    grid(true)
    legend()
    savefig(savefname)
end

if caseff
    #figure()
    clf()
    title(titlelbl)
    plot(σstdff[:,zlev,cval],label=L"std()")
    plot(σmedff[:,zlev,cval],label=L"median(||)")
    plot(σmeanff[:,zlev,cval],label=L"mean(||)")
    #plot(σmaxff[:,zlev,cval],label=L"max(||)")
    xlabel("month starting in Jan. 1992")
    ylabel(ylbl)
    grid(true)
    legend()
    savefig(savefname)
end

zlev = 20
#figure()
#plot(σ[:,20,1])
title(L"(\theta_{129} - \theta_{ff})/(\theta_{129} - \theta_0) "*string(floor(Int,-Z[zlev]))*"m")
plot(σstdff[:,zlev,1]./σstd0[:,zlev,1],label=L"std")
plot(σmedff[:,zlev,1]./σmed0[:,zlev,1],label=L"median")
plot(σmeanff[:,zlev,1]./σmean0[:,zlev,1],label=L"mean")
#plot(σmax[:,zlev,1],label=L"max(|\theta_{129}-\theta_{0}|)")
xlabel("month starting in Jan. 1992")
ylabel(L"\theta  [^{\circ}C]")
grid(true)
legend()
savefig("dtheta_129v0vff_1992-2017_299m_nomax_12feb2021.eps")

# make depth plot with std, median, etc.
figure()
clf()
semilogx(σ[:,1],Z,label=L"std(\theta_{129}-\theta_{0})")
semilogx(σ50[:,1],Z,label=L"median(|\theta_{129}-\theta_{0}|)")
semilogx(σmean[:,1],Z,label=L"mean(|\theta_{129}-\theta_{0}|)")
semilogx(σmax[:,1],Z,label=L"max(|\theta_{129}-\theta_{0}|)")
title("January 1992")
xlabel(L"\theta  [^{\circ}C]")
ylabel("depth [m]")
grid("true")
legend()


# plot a timeseries of T vs T and S vs S would be interesting. 

# make a timeseries of the statistics. would be cool.

# make one plan view.
