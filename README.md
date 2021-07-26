# Reemergence

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ggebbie.github.io/Reemergence.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ggebbie.github.io/Reemergence.jl/dev)
[![Build Status](https://github.com/ggebbie/Reemergence.jl/workflows/CI/badge.svg)](https://github.com/ggebbie/Reemergence.jl/actions)
[![Coverage](https://codecov.io/gh/ggebbie/Reemergence.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ggebbie/Reemergence.jl)

Tracking seawater property anomalies through the oceanic tunnel until possible reemergence

# Getting started

* from Emacs editor (one possible method)

Install julia-mode, julia-repl, and magit \
Skip the next 5 steps if you have already cloned the repository \
`M-x magit-clone` \
Select `u` to clone from url\
Enter ` https://github.com/ggebbie/Reemergence.jl` as url to clone \
Select `y` in response to `remote.pushDefault' to "origin"?` \
Clone to your favorite location and rename project if necessary \
Go to any directory in the project: `C-x C-f Reemergence.jl`\
Then activate the project and initialize a julia session: `C-c C-a`

* from the Julia REPL

`;`\
`git clone https://github.com/ggebbie/Reemergence.jl # only do this the first time on each machine`\
`cd Reemergence.jl`\
`]`\
`activate .`\
`instantiate # only do this the first time on each machine`\
To verify you are in the project environment, `]` should return `(Reemergence) pkg>`\
Type backspace to return to command mode.

* Using an editor like Atom/Juno or Visual Studio Code, activate the environment on one of the frame panels. The default environment is @v1.x and should be changed.

# Running a script (not interactively)

`cd Reemergence.jl`\
`julia --project=@. scripts/experiment_divergence.jl`

# Functions
Check `export` line in src/Reemergence.jl for list of available functions.

# Directory structure
- `examples`: code snippets useful for making new scripts and functions
- `test`: testing routines, currently empty
- `scripts`: production-ready scripts
