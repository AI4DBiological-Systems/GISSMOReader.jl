
# run this before graphs.jl and rest.jl variants.

import GISSMOReader

PyPlot.close("all")
fig_num = 1

Random.seed!(25)

# get the all GISSMO entries.
tmp = getGISSMOentriesall()
GISSMO_entries = extractfields(tmp, "entry")
molecule_names = extractfields(tmp, "molecule_name")

# save foldder.
#base_dir = "/home/roy/MEGAsync/data/molecules"
base_dir = "/home/roy/del/base"

molar_masses = collect( fetchentrymetadata(GISSMO_entries[i])[1] for i = 1:length( GISSMO_entries))

save_path = joinpath(base_dir, "molar_masses.jld")
JLD.save(save_path, "molecule_names", molecule_names,
                    "molar_masses", molar_masses)

# load.
load_path = save_path
dic = JLD.load(load_path)

@assert all(isfinite.(dic["molar_masses"]))
