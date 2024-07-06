"""
    function listexperiments(exppath)
    get a dictionary of experiments and their locations
# Arguments
- `exppath`: directory of experiments
# Output
- `runpath`: dictionary with keys=experiments, values=run paths
- `diagpath`: dictionary with keys=experiments, values=diagnostic paths
"""
function listexperiments(exppath; printavailable = false)

    # add a trailing / if needed
    exppath[end] != "/" ? exppath *= "/" : nothing
    
    explist = searchdir(exppath,"") # all files in directory
    filter!(x -> !occursin("README",x),explist) # remove README to get explist
    filter!(x -> !occursin(".out",x),explist) # remove output to get explist

    runpath = Dict(explist .=> exppath.*explist.*"/run/")
    diagpath = Dict(explist .=> exppath.*explist.*"/run/diags/")
    regpolespath = Dict(explist .=> exppath.*explist.*"/run/regularpoles/")
    
    if printavailable
        # print to screen all the available
        println("Available experiments are")
    end

    for (key,value) in runpath
        println(key)
    end
    return runpath,diagpath,regpolespath
end

"""
    function exprootdir(expt::String) 
    Root directory of the ECCOv4r4 output
"""
function exprootdir(expt::String)
    rootdir = "/batou/ECCOv4r4/exps/"*expt*"/"

    # If the experiment hasn't been copied to batou, look on poseidon.
    !isdir(rootdir) ? rootdir = "/vast/ECCOv4r4/exps/"*expt*"/" : nothing
    return rootdir
end

"""
    function exprootdir() 
    Root directory of the ECCOv4r4 output
"""
exprootdir() = "/batou/ECCOv4r4/exps"


rundir(expt::String) = exprootdir(expt)*"run/"
sig1dir(expt::String) = rundir(expt)*"sigma1/"
diagdir(expt::String) = rundir(expt)*"diags/"

vastexprootdir(expt::String) = "/vast/ECCOv4r4/exps/" * expt*"/"
vastrundir(expt::String) = vastexprootdir(expt)*"run/"
vastrundir(expt::String, runname::String) = vastexprootdir(expt)*runname * "/"

vastdiagdir(expt::String) = vastrundir(expt)*"diags/"
vastdiagdir(expt::String, runname::String) = vastrundir(expt, runname)*"diags/"

vastexprootdir() = "/vast/ECCOv4r4/exps"

update_paths_dict!(main_dict, add_dict) = [main_dict[key] = add_dict[key] for key ∈ keys(add_dict) if key ∉ keys(main_dict)]

################################################################
# get MITgcm / ECCOv4r4 LLC grid and depth information. Store in γ.
pth = MeshArrays.GRID_LLC90
γ = GridSpec("LatLonCap",pth)
Γ = GridLoad(γ;option="full")

nf = length(γ.fSize)

# get standard levels of MITgcm
p₀ = 1000.0; # dbar

Δz = read_mdsio(γ.path,"DRC")
ΔzF = vec(read_mdsio(γ.path,"DRF"))
Δz = vec(Δz)

#retrieves the array containing all basin IDs 
basins=read(joinpath(pth,"v4_basin.bin"),MeshArray(γ,Float32))

basin_list=["Pacific","Atlantic","indian","Arctic","Bering Sea",
            "South China Sea","Gulf of Mexico","Okhotsk Sea",
            "Hudson Bay","Mediterranean Sea","Java Sea","North Sea",
            "Japan Sea", "Timor Sea","East China Sea","Red Sea",
            "Gulf","Baffin Bay","GIN Seas","Barents Sea"];
            
z = depthlevels(γ);
pstdz = pressurelevels(z);

#load mask of open boundaries
mskC, mskW, mskS = get_msk(Γ)
