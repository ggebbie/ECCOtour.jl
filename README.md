# ECCOonPoseidon

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ggebbie.github.io/ECCOonPoseidon.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ggebbie.github.io/ECCOonPoseidon.jl/dev)
[![Build Status](https://github.com/ggebbie/ECCOonPoseidon.jl/workflows/CI/badge.svg)](https://github.com/ggebbie/ECCOonPoseidon.jl/actions)
[![Coverage](https://codecov.io/gh/ggebbie/ECCOonPoseidon.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ggebbie/ECCOonPoseidon.jl)

Analysis of ECCO version 4 release 4 using output from Poseidon @ WHOI

# Getting started

* from Emacs editor (one possible method)

Install julia-mode, julia-repl, and magit \
Skip the next 5 steps if you have already cloned the repository \
`M-x magit-clone` \
Select `u` to clone from url\
Enter ` https://github.com/ggebbie/ECCOonPoseidon.jl` as url to clone \
Select `y` in response to `remote.pushDefault' to "origin"?` \
Clone to your favorite location and rename project if necessary \
Go to any directory in the project: `C-x C-f ECCOonPoseidon.jl`\
Then activate the project and initialize a julia session: `C-c C-a`

* from the Julia REPL

`;`\
`git clone https://github.com/ggebbie/ECCOonPoseidon.jl # only do this the first time on each machine`\
`cd ECCOonPoseidon.jl`\
`]`\
`activate .`\
`instantiate # only do this the first time on each machine`\
To verify you are in the project environment, `]` should return `(ECCOonPoseidon) pkg>`\
Type backspace to return to command mode.

* Using an editor like Atom/Juno or Visual Studio Code, activate the environment on one of the frame panels. The default environment is @v1.x and should be changed.

# Running a script (not interactively)

`cd ECCOonPoseidon.jl`\
`julia --project=@. scripts/experiment_divergence.jl`

# Functions
Check `export` line in src/ECCOonPoseidon.jl for list of available functions.

# Directory structure
- `examples`: code snippets useful for making new scripts and functions
- `test`: testing routines, currently empty
- `scripts`: production-ready scripts
