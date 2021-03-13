# Reemergence

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ggebbie.github.io/Reemergence.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ggebbie.github.io/Reemergence.jl/dev)
[![Build Status](https://github.com/ggebbie/Reemergence.jl/workflows/CI/badge.svg)](https://github.com/ggebbie/Reemergence.jl/actions)
[![Coverage](https://codecov.io/gh/ggebbie/Reemergence.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ggebbie/Reemergence.jl)

Tracking seawater property anomalies through the oceanic tunnel until possible reemergence

# Getting started
from the Julia REPL\
`;`\
`git clone https://github.com/ggebbie/Reemergence.jl # only do this the first time on each machine`\
`cd Reemergence.jl`\
`]`\
`activate .`\
`instantiate # only do this the first time on each machine`\
To verify you are in the project environment, `]` should return `(Reemergence) pkg>`\
Type backspace to return to command mode.

# Running a script (not interactively)

`cd Reemergence.jl`\
`julia --project=@. scripts/experiment_divergence.jl` 