# Example: how to read mds format output on poseidon.whoi.edu.
using MITgcmTools, MeshArrays
    
path = "/poseidon/ecco/ECCOv4r4/MITgcm/exps/iter129_bulkformula/run/"
θfile = "T.0000000001.data" # it's actually potential temperature snapshot
@time θmds     = read_mdsio(path,θfile); # output is a 3D array

# To use a gcmarray variable type, need to have grid information in γ
@time θstart   = γ.read(path*θfile,MeshArray(γ,Float32));

# Do same thing with final snapshot.
θfile = "T.0000227904.data"
θend=γ.read(path*θfile,MeshArray(γ,Float32))

# operations proceed normally
Δθ = θend - θstart
maximum(ΔΘ) # etc. 

# compare to reading monthly average fields
path = "/poseidon/ecco/ECCOv4r4/MITgcm/exps/iter129_fluxforced/run/diags/"
θfile = "state_3d_set1.0000000732.data" 
@time xmds     = read_mdsio(path,θfile); 

# To use a gcmarray variable type, need to have grid information in γ
@time xstart   = γ.read(path*θfile,MeshArray(γ,Float32));
