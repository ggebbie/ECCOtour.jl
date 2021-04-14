#  First steps: 1. go into Reemergence project directory. 2. go into julia REPL package mode with `]`. 3. `activate .` 4. Backspace to return to command mode in REPL.

# make spatial plots. Use interpolation to regular grid.
# rectangular grid
longrid = -179.:2.0:179.; latgrid = -89.:2.0:89.;
f,i,j,w = prereginterp(latgrid,longrid,γ)
nx = length(longrid); ny = length(latgrid);

figure(101)
for zz = 1:nz
    βz = β[:,zz]
    βzreg = reginterp(βz,nx,ny,f,i,j,w)
    clf()
    contourf(longrid,latgrid,βzreg)

    # save it
    
end
