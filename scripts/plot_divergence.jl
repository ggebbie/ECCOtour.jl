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
expcompare = "iter129_fluxforced"

# recompute or read from file?
readfromfile = true
if readfromfile
    outpath = "../outputs/"
    outfile = outpath*"divergence_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*".jld2"
    @load outfile xbar σx xmax xmin z 
else
    include("experiment_divergence.jl")
end

# To improve, look at both min and max to get the most extreme value.

##### Choose the depth level and variable of interest. #########
zlev = 20
cval = 1

################ Linear scaling of divergence ##################
xlbl = "month starting in Jan. 1992"
depthlbl = string(round(Int,-z[zlev]))*"m"
if cval == 1
    str1 = raw"$\theta_{"*shortnames[expbase]*raw"}$"
    str2 = raw"$\theta_{"*shortnames[expcompare]*raw"}$"
    titlelbl = str1*raw"$-$"*str2
    savefname = outpath*"dtheta_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_1992-2017_"*depthlbl*".eps"
    ylbl = L"\theta [^{\circ}C]"
elseif cval == 2
    str1 = raw"$S_{"*shortnames[expbase]*raw"}$"
    str2 = raw"$S_{"*shortnames[expcompare]*raw"}$"
    titlelbl = str1*raw"$-$"*str2
    savefname = outpath*"dsalt_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_1992-2017_"*depthlbl*".eps"
    ylbl = L"salinity  [PSS-1978]"
end

figure(10)
clf()
title(titlelbl)
plot(σx[:,zlev,cval],label=L"std()")
plot(abs.(xbar[:,zlev,cval]),label=L"|mean()|")
plot(xmax[:,zlev,cval],label=L"max()")
xlabel(xlbl)
ylabel(ylbl)
grid(true)
legend()
savefig(savefname)

###### Log scaling of divergence ##############################
xlbl = "month starting in Jan. 1992"
depthlbl = string(round(Int,-z[zlev]))*"m"
if cval == 1
    str1 = raw"$\theta_{"*shortnames[expbase]*raw"}$"
    str2 = raw"$\theta_{"*shortnames[expcompare]*raw"}$"
    titlelbl = str1*raw"$-$"*str2
    savefname = outpath*"logdtheta_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_1992-2017_"*depthlbl*".eps"
    ylbl = L"\theta [^{\circ}C]"
elseif cval == 2
    str1 = raw"$S_{"*shortnames[expbase]*raw"}$"
    str2 = raw"$S_{"*shortnames[expcompare]*raw"}$"
    titlelbl = str1*raw"$-$"*str2
    savefname = outpath*"logdsalt_"*shortnames[expbase]*"_vs_"*shortnames[expcompare]*"_1992-2017_"*depthlbl*".eps"
    ylbl = L"salinity  [PSS-1978]"
end

figure(11)
clf()
title(titlelbl)
semilogy(σx[:,zlev,cval],label=L"std()")
semilogy(abs.(xbar[:,zlev,cval]),label=L"|mean()|")
semilogy(xmax[:,zlev,cval],label=L"max()")
xlabel(xlbl)
ylabel(ylbl)
grid(true)
legend()
savefig(savefname)

##################################
## to do: work on these other types of plots.
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

# # make depth plot with std, median, etc.
# figure()
# clf()
# semilogx(σ[:,1],Z,label=L"std(\theta_{129}-\theta_{0})")
# semilogx(σ50[:,1],Z,label=L"median(|\theta_{129}-\theta_{0}|)")
# semilogx(σmean[:,1],Z,label=L"mean(|\theta_{129}-\theta_{0}|)")
# semilogx(σmax[:,1],Z,label=L"max(|\theta_{129}-\theta_{0}|)")
# title("January 1992")
# xlabel(L"\theta  [^{\circ}C]")
# ylabel("depth [m]")
# grid("true")
# legend()


# plot a timeseries of T vs T and S vs S would be interesting. 

# make a timeseries of the statistics. would be cool.

# make one plan view.
