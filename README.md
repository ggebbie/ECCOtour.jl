# ECCOtour

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ggebbie.github.io/ECCOtour.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ggebbie.github.io/ECCOtour.jl/dev)
[![Build Status](https://github.com/ggebbie/ECCOtour.jl/workflows/CI/badge.svg)](https://github.com/ggebbie/ECCOtour.jl/actions)
[![Coverage](https://codecov.io/gh/ggebbie/ECCOtour.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ggebbie/ECCOtour.jl)

Take a self-guided tour of the ECCO version 4 release 4 ocean reanalysis product. 

# Getting started

* from Emacs editor (one possible method)

Install julia-mode, julia-repl, and magit \
Skip the next 5 steps if you have already cloned the repository \
`M-x magit-clone` \
Select `u` to clone from url\
Enter ` https://github.com/ggebbie/ECCOtour.jl` as url to clone \
Select `y` in response to `remote.pushDefault' to "origin"?` \
Clone to your favorite location and rename project if necessary \
Go to any directory in the project: `C-x C-f ECCOtour.jl`\
Then activate the project and initialize a julia session: `C-c C-a`

* from the Julia REPL

`;`\
`git clone https://github.com/ggebbie/ECCOtour.jl # only do this the first time on each machine`\
or\
`git clone git@github.com:ggebbie/ECCOtour.jl.git # if you have SSH keys set up
`cd ECCOtour.jl`\
`]`\
`activate .`\
`instantiate # only do this the first time on each machine`\
To verify you are in the project environment, `]` should return `(ECCOtour) pkg>`\
Type backspace to return to command mode.

* Using an editor like Atom/Juno or Visual Studio Code, activate the environment on one of the frame panels. The default environment is @v1.x and should be changed.

# Functions
Check `export` line in src/ECCOtour.jl for list of available functions. Or check out the documentation at https://ggebbie.github.io/ECCOtour.jl/dev.

# Directory structure
- `src`: source code
- `test`: a few testing routines
- `examples`: useful code snippets
