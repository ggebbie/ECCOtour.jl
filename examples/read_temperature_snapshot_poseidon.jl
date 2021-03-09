# Example: how to read mds format output on poseidon.whoi.edu.

path = "/poseidon/ecco/ECCOv4r4/MITgcm/exps/iter129_bulkformula/run/"
θfile = "T.0000000001.data" # it's actually potential temperature snapshot
θmds     = read_mdsio(path,θfile) # output is a 3D array

# To use a gcmarray variable type, need to have grid information in γ
θstart   = γ.read(path*θfile,MeshArray(γ,Float32))

# Do same thing with final snapshot.
θfile = "T.0000227904.data"
θend=γ.read(path*θfile,MeshArray(γ,Float32))

# operations proceed normally
Δθ = θend - θstart
maximum(ΔΘ) # etc. 
