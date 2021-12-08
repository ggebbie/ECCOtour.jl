using Pkg, Conda

# this line doesn't seem to be hit
println("building ECCOtour package")

#Pkg.add(url="https://github.com/ggebbie/SigmaShift.jl")
#Pkg.instantiate()

if lowercase(get(ENV, "CI", "false")) == "true"    

    ENV["PYTHON"] = ""
    Pkg.build("PyCall")

    Conda.add("matplotlib",channel="conda-forge")
    Conda.add("shapely",channel="conda-forge")
    Conda.add("cartopy=0.20.0",channel="conda-forge")

end
