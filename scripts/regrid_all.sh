#!/bin/sh

julia --project=@. scripts/mdsiofiles2regularpoles.jl 1> scripts/mdsiofiles2regularpoles.out 2> scripts/mdsiofiles2regularpoles.err

julia --project=@. scripts/netcdffiles2regularpoles.jl 1> scripts/netcdffiles2regularpoles.out 2> scripts/netcdffiles2regularpoles.err

julia --project=@. scripts/regularpoles2sigma1.jl 1> scripts/regularpoles2sigma1.out 2> scripts/regularpoles2sigma1.err 
