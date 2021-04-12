%% Examples for reading "native" MITgcm output.
%  This is "mds" format, includes a data and meta file.
%  Here, 1) read a snapshot, 2) read all timeslices, 
%  3) convert from 4D (3 spatial dimensions plus computational tiles) to 3D,
%  4) interpolate onto regular grid for plotting.

!git clone https://www.github.com/gaelforget/gcmfaces

% OPTIONAL: make your own branch
cd gcmfaces
!git checkout -b poseidon_taris
cd ..
% END OF OPTION

%add gcmfaces and MITprof directories to Matlab path:
p = genpath('gcmfaces/'); addpath(p);
%p = genpath('MITprof/'); addpath(p);

% download ECCOv4r4 grid
! wget --recursive ftp://mit.ecco-group.org/ecco_for_las/version_4/release1/nctiles_grid 
! mv mit.ecco-group.org/ecco_for_las/version_4/release1/nctiles_grid .

%load all grid variables from nctiles_grid/ into mygrid:
grid_load;

%make mygrid accessible in current workspace:
gcmfaces_global;

%display list of grid variables:
disp(mygrid);

%display one gcmfaces variable:
disp(mygrid.XC);

%fld{n}     returns the n^{th} face data (i.e., an array).
%fld(:,:,n) returns the n^{th} vertical level (i.e., a gcmfaces object).

% experiment names
expnames={'iter129_bulkformula','iter129_fluxforced','iter0_bulkformula','nointerannual','nosfcadjust','noinitadjust'}

% choose an experiment
expno = 1;

%example: MITgcm monthly-average output, available on poseidon or backup (batou)
%filename = ['/poseidon/ECCOv4r4/MITgcm/exps/',expnames{expno},'/run/diags/state_3d_set1.0000000732.data']
filename = ['/batou/ECCOv4r4/MITgcm/exps/',expnames{expno},'/run/diags/state_3d_set1.0000000732.data']

%% this is gcmfaces format with 5 faces/tiles
% according to the accompanying meta file:
% dimension 3: 1:50= theta
% dimension 3: 51:100= S
% dimension 3: 101:150= drho/dz
fld = read_bin(filename) % this is gcmfaces format with 5 faces/tiles

fldlevel = fld(:,:,10);
contourf(fldlevel{1}) % contour of one tile with eastern Atlantic

%% or convert to a standard 3D array
% convert2gcmfaces will toggle the variable type
fld3d = convert2gcmfaces(fld); 

%% How to convert this unusual 3D grid to something plottable?
% use gcmfaces/example_interp.m

%target locations:
lon=[-179.75:0.5:179.75]; lat=[-89.75:0.5:89.75];
[lat,lon] = meshgrid(lat,lon);
fld_interp=gcmfaces_interp_2d(fldlevel,lon,lat); % slow: takes ~100 s
figure(10)
contourf(fld_interp') % looks good!

%% This method reads all timeslices at once. Warning: takes a while
filename = ['/poseidon/ECCOv4r4/MITgcm/exps/',expnames{expno},'/run/diags/state_3d_set1']
fldall = rdmds2gcmfaces(filename)

%% Example: read property values on sigma-1 surfaces
% Need to check: does this work if accompanying meta files are not available?
filename = ['/batou/ECCOv4r4/MITgcm/exps/',expnames{expno},'/run/sigma1/theta_on_sigma1.0000000732.data']
ilev = 50; % pick a sigma-1 surface
fld = read_bin(filename) % this is gcmfaces format with 5 faces/tiles
fld_interp=gcmfaces_interp_2d(fld(:,:,ilev),lon,lat); % slow: takes ~100 s
figure(1)
[c,h]=contourf(fld_interp',-2:32);
clabel(c,h)
