# Work in progress.
using Revise
using MITgcmTools, MeshArrays, PyPlot, LaTeXStrings, JLD2, Reemergence

# list of experiments on poseidon
exppath = "/poseidon/ecco/ECCOv4r4/MITgcm/exps/"
runpath,diagpath = listexperiments(exppath);
# abbreviations for each experiment for labels, etc.
shortnames = expnames()

## SELECT EXPERIMENTS TO COMPARE #################################
# baseline experiment
# call 
expbase = "iter129_bulkformula"
# comparison experiment(s)
#expcompare = "iter129_fluxforced"
expcompare = "nosfcadjust"
#expcompare = "noinitadjust"

# recompute or read from file?
readfromfile = true
if readfromfile
    outpath = "../outputs/"
    outfile = outpath*"divergence_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*".jld2"
    @load outfile xbar σx xmax xmin z 
else
    include("experiment_divergence.jl")
end

# Look at both min and max to get the most extreme value.
xextreme = max.(-xmin,xmax);

# Choose the depth level and variable of interest.
zlev = 20  # there are nz=50 levels for ECCOv4r4.

for cval = 1:2 # variable 1 = theta, variable 2 = practical salinity

    ################ Linear scaling of divergence ##################
    xlbl = "month starting in Jan. 1992"
    depthlbl = string(round(Int,-z[zlev]))*"m"
    if cval == 1
        str1 = raw"$\theta_{"*shortnames[expbase]*raw"}$"
        str2 = raw"$\theta_{"*shortnames[expcompare]*raw"}$"
        titlelbl = str1*raw"$-$"*str2*","*depthlbl
        savefname = outpath*"dtheta_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_1992-2017_"*depthlbl*".eps"
        savetimefname = outpath*"dtheta_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_"*tlbl[1:3]*tlbl[5:8]*".eps"
         savelogfname = outpath*"logdtheta_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_1992-2017_"*depthlbl*".eps"
        ylbl = L"\theta [^{\circ}C]"
    elseif cval == 2
        str1 = raw"$S_{"*shortnames[expbase]*raw"}$"
        str2 = raw"$S_{"*shortnames[expcompare]*raw"}$"
        titlelbl = str1*raw"$-$"*str2*","*depthlbl
        savefname = outpath*"dsalt_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_1992-2017_"*depthlbl*".eps"
        savetimefname = outpath*"dsalt_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_"*tlbl[1:3]*tlbl[5:8]*".eps"
        savelogfname = outpath*"logdsalt_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_1992-2017_"*depthlbl*".eps"
        ylbl = L"salinity  [PSS-1978]"
    end
    zlbl  = "depth [m]"

    figure(10+cval)
    clf()
    title(titlelbl)
    plot(σx[:,zlev,cval],label=L"std()")
    plot(abs.(xbar[:,zlev,cval]),label=L"|mean()|")
    plot(xextreme[:,zlev,cval],label=L"max(||)")
    xlabel(xlbl)
    ylabel(ylbl)
    grid(true)
    legend()
    savefig(savefname)

    figure(20+cval)
    clf()
    title(titlelbl)
    semilogy(σx[:,zlev,cval],label=L"std()")
    semilogy(abs.(xbar[:,zlev,cval]),label=L"|mean()|")
    semilogy(xextreme[:,zlev,cval],label=L"max(||)")
    xlabel(xlbl)
    ylabel(ylbl)
    grid(true)
    legend()
    savefig(savelogfname)

    ## make divergence at a snapshot as a function of depth
    #  pick a time index
    tlist = (1,nt)
    for tt ∈ tlist
        tlbl = time_label(tt-1) # subtract one, months since Jan 1992
        figure(12+tt+cval)
        clf()
        semilogx(σx[tt,:,cval],z,label=L"std()")
        semilogx(xbar[tt,:,cval],z,label=L"|mean()|")
        semilogx(xextreme[tt,:,cval],z,label=L"max(||)")
        title(str1*raw"$-$"*str2*", "*tlbl)
        xlabel(ylbl)
        ylabel(zlbl)
        grid("true")
        legend()
        savefig(savetimefname)

    end
end

# Ratio of divergence 
# zlev = 20
# #figure()
# #plot(σ[:,20,1])
# title(L"(\theta_{129} - \theta_{ff})/(\theta_{129} - \theta_0) "*string(floor(Int,-Z[zlev]))*"m")
# plot(σstdff[:,zlev,1]./σstd0[:,zlev,1],label=L"std")
# plot(σmedff[:,zlev,1]./σmed0[:,zlev,1],label=L"median")
# plot(σmeanff[:,zlev,1]./σmean0[:,zlev,1],label=L"mean")
# #plot(σmax[:,zlev,1],label=L"max(|\theta_{129}-\theta_{0}|)")
# xlabel("month starting in Jan. 1992")
# ylabel(L"\theta  [^{\circ}C]")
# grid(true)
# legend()
# savefig("dtheta_129v0vff_1992-2017_299m_nomax_12feb2021.eps")


# plot a timeseries of T vs T and S vs S would be interesting. 

# make a timeseries of the statistics. would be cool.

# make one plan view.
