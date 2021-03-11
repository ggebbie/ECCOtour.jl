# Work in progress.
# Called from experiment_divergence.jl.
using PyPlot

zlev = 20
cval = 1

figure()
clf()
#title(titlelbl)
plot(σx[:,zlev,cval],label=L"std()")
plot(xbar[:,zlev,cval],label=L"mean(||)")
plot(xmax[:,zlev,cval],label=L"max(||)")
xlabel("month starting in Jan. 1992")
#ylabel(ylbl)
grid(true)
legend()


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
