# Work in progress.
using Revise
using MITgcmTools, MeshArrays, PyPlot, LaTeXStrings, JLD2, Reemergence

# recompute or read from file?
readfromfile = true

# list of experiments on poseidon
exppath = "/poseidon/ecco/ECCOv4r4/MITgcm/exps/"
runpath,diagpath = listexperiments(exppath);

# abbreviations for each experiment for labels, etc.
shortnames = expnames()

## SELECT EXPERIMENTS TO COMPARE #################################
# baseline experiment
expbase = "iter129_bulkformula"

# comparison experiment(s)
#expcompare = "iter129_fluxforced"
expcompare = "nosfcadjust"
#expcompare = "noinitadjust"

outpath = "../outputs/"

if readfromfile
    outfile = outpath*"divergence_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*".jld2"
    @load outfile xbar σx xmax xmin absxbar z 
else
    include("experiment_divergence.jl")
end

# Look at both min and max to get the most extreme value.
xextreme = max.(-xmin,xmax);

# Choose the depth level and variable of interest.
zlev = 20  # there are nz=50 levels for ECCOv4r4.
zlbl  = "depth [m]"
depthlbl = string(round(Int,-z[zlev]))*"m"
for cval = 1:2 # variable 1 = theta, variable 2 = practical salinity

    ################ Linear scaling of divergence ##################
    xlbl = "month starting in Jan 1992"
    if cval == 1
        str1 = raw"$\theta_{"*shortnames[expbase]*raw"}$"
        str2 = raw"$\theta_{"*shortnames[expcompare]*raw"}$"
        titlelbl = str1*raw"$-$"*str2*","*depthlbl
        linfname = outpath*"dtheta_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_1992-2017_"*depthlbl*".eps"
        logfname = outpath*"logdtheta_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_1992-2017_"*depthlbl*".eps"
        ylbl = L"\theta"*" "*L"[^{\degree}C]"
    elseif cval == 2
        str1 = raw"$S_{"*shortnames[expbase]*raw"}$"
        str2 = raw"$S_{"*shortnames[expcompare]*raw"}$"
        titlelbl = str1*raw"$-$"*str2*","*depthlbl
        linfname = outpath*"dsalt_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_1992-2017_"*depthlbl*".eps"
        logfname = outpath*"logdsalt_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_1992-2017_"*depthlbl*".eps"
        ylbl = "salinity  [PSS-1978]"
    end

    figure(10+cval)
    clf()
    title(titlelbl)
    lev = zlev + (cval-1)*50 # 2 properties in 1 variable
    plot(σx[:,lev],label=L"std()")
    plot(absxbar[:,lev],label=L"mean(||)")
    plot(xextreme[:,lev],label=L"max(||)")
    xlabel(xlbl)
    ylabel(ylbl)
    grid(true)
    legend()
    savefig(linfname)

    figure(20+cval)
    clf()
    title(titlelbl)
    lev = zlev + (cval-1)*50
    semilogy(σx[:,lev],label=L"std()")
    semilogy(absxbar[:,lev],label=L"mean(||)")
    semilogy(xextreme[:,lev],label=L"max(||)")
    xlabel(xlbl)
    ylabel(ylbl)
    grid(true)
    legend()
    savefig(logfname)

    ## make divergence at a snapshot as a function of depth
    #  pick a time index: here, beginning and end
    tlist = (1,nt)
    for tt ∈ tlist
        tlbl = time_label(tt-1) # subtract one, months since Jan 1992
        if cval == 1
            timefname = outpath*"dtheta_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_"*tlbl[1:3]*tlbl[5:8]*".eps"
        elseif cval ==2
            
            timefname = outpath*"dsalt_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_"*tlbl[1:3]*tlbl[5:8]*".eps"
        end
        
        levs = 1+(cval-1)*50:cval*50
            
        figure(12+tt+cval)
        clf()
        semilogx(σx[tt,levs],-z,label=L"std()")
        semilogx(absxbar[tt,levs],-z,label=L"mean(||)")
        semilogx(xextreme[tt,levs],-z,label=L"max(||)")
        title(str1*raw"$-$"*str2*", "*tlbl)
        xlabel(ylbl)
        ylabel(zlbl)
        grid(true)
        legend()
        plt.gca().invert_yaxis()
        savefig(timefname)

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
