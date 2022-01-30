# GISSMOReader

[![Build Status](https://github.com/RoyCCWang/GISSMOReader.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/RoyCCWang/GISSMOReader.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/RoyCCWang/GISSMOReader.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/RoyCCWang/GISSMOReader.jl)

Run `./examples/batch_graphs.jl` before `./examples/batch_graphs.jl`, and make sure the `base_dir` variable for both scripts points to the same folder.

See `./examples/batch_graphs.jl` to store the J-coupling and chemical shift values from all GISSMO entries that were hardcoded in `./src/database/GISSMO_entries.jl` (see `function getGISSMOentriesall()`). Modify that function if you wish to add new GISSMO entries than the ones already there.
