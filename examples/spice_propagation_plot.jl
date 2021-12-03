# Warning: setupLLCgrid no longer exists
# ggebbie, 16-Apr-2021

using Revise
using MITgcmTools, MeshArrays, Statistics, ECCOtour, PyPlot

# get MITgcm / ECCOv4r4 LLC grid and depth information. Store in γ.
path_grid="../inputs/GRID_LLC90/"
γ = setupLLCgrid(path_grid)

# list of experiments on poseidon
exppath = "/batou/ECCOv4r4/MITgcm/exps/"
runpath,diagpath = listexperiments(exppath);

# print output here
path_out = "/home/gebbie/julia/outputs/"

# abbreviations for each experiment for labels, etc.
shortnames = expnames()

## SELECT EXPERIMENTS TO ANALYZE ##############
# manually choose from available experiments listed above.
exps = ("iter129_bulkformula","nointerannual")

# to do all experiments:
# exps = keys(shortnames)
###############################################

nexps = length(exps) # number of experiments

## DEFINE THE LIST OF SIGMA1 VALUES.
sig1grida = 24:0.05:31
sig1gridb = 31.02:0.02:33
sig1grid = vcat(sig1grida,sig1gridb)
sig1grid = sig1grid[1:3:end]
nσ = length(sig1grid)

# get index for sigma1 value of interest.
sig1 = 29.9 # a choice
latsig1 = -12.0 # lat of interest

junk,isig1grid = findmin(abs.(sig1grid.-sig1))        

tecco= 1992+1/24:1/12:2018 # ecco years

ϕGreg,ϕCreg = latgridRegular(γ) 
junk,ilat = findmin(abs.(ϕCreg.-latsig1))        
λC = -179.5:179.5
# center the plot on indo-pac
cutoff = 181
λCcutoff = cat(λC[cutoff:end],λC[1:cutoff-1],dims = 1)

# make all longitudes positive
replace!(x -> x<0 ? x+=360 : x,λCcutoff)

nx = length(λC)
TSroot = "S_on_sigma1" 
cmap_seismic =get_cmap("seismic")
figure(100)
clf()
plot(λCcutoff)

# cycle through all chosen experiments
#ex = exps[1]; for manual runs
for ex in exps

    # Get list of files for salinity on sigma1
    filelist = searchdir(runpath[ex]*"sigma1/",TSroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nt = length(datafilelist)
    Ssig1slice = Array{Float32,2}(undef,nx,nt)

    # Find a latitudinal slice. i.e., 12 South. Use convert2array? Or something else? Use regular poles method.
    
    # pre-allocate
    #TS = Array{Float32, 2}(undef, nt, nz*2)

    global tt = 0
    for Sname in datafilelist
        tt += 1
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)

        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        @time S = γ.read(runpath[ex]*"sigma1/"*Sname,MeshArray(γ,Float32,nσ))

        # just take one sigma1 surface
        Ssig1 = S[:,isig1grid]

        # next slice at 12 South.
        Ssig1crop =  LLCcropC(Ssig1,γ) # regular grid by cropping        
        Ssig1slice[:,tt] = Ssig1crop[:,ilat]
    end
    #Ssig1slice = Ssig1slice[cutoff:end 1:cutoff-1,:]
    Ssig1slice = cat(Ssig1slice[cutoff:end,:],Ssig1slice[1:cutoff-1,:],dims = 1)

    figure(10)
    clf()
    iplot = 140:290; # plot indo-pac only
    lims = 35.0:.05:36.0
    contourf(λCcutoff[iplot],tecco,Ssig1slice[iplot,:]',lims,cmap=cmap_seismic)
    colorbar(label="[PSS-78]",orientation="vertical",ticks=lims)
    contour(λCcutoff[iplot],tecco,Ssig1slice[iplot,:]',lims,colors="k")


    titlelbl = ex*", "*string(latsig1)*L"\degree"*"N, "*L"\sigma_1="*string(sig1)*" kg m"*L"^{-3}"
    xlbl = "longitude "*L"[\degree E]"
    ylbl = "years"

    title(titlelbl)
    xlabel(xlbl)
    ylabel(ylbl)


    # Put output (i.e., figures) locally.
    pathoutexp = path_out*ex
    !isdir(pathoutexp) ? mkdir(pathoutexp) : nothing;
    if latsig1 > 0
        outfname = pathoutexp*"/Sonsigma1_"*shortnames[ex]*"_"*string(latsig1)*"N_"*string(sig1)*".eps"
    else
        outfname = pathoutexp*"/Sonsigma1_"*shortnames[ex]*"_"*string(abs(latsig1))*"S_"*string(sig1)*"sig1.eps"
    end
    
    savefig(outfname)
end
