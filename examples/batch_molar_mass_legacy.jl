
# run this before graphs.jl and rest.jl variants.

#include("../src/GISSMOReader.jl")
import GISSMOReader
import JLD

# get the all GISSMO entries.
tmp = GISSMOReader.getGISSMOentriesall()
GISSMO_entries = GISSMOReader.extractfields(tmp, "entry")
molecule_names = GISSMOReader.extractfields(tmp, "molecule_name")

# save foldder.
base_dir = "/home/roy/Documents/repo/NMRData/input/GISSMO_data"

molar_masses = collect( GISSMOReader.fetchentrymetadata(GISSMO_entries[i])[1] for i = 1:length( GISSMO_entries))

save_path = joinpath(base_dir, "molar_masses.jld")
JLD.save(save_path, "molecule_names", molecule_names,
                    "molar_masses", molar_masses)

# load.
load_path = save_path
dic = JLD.load(load_path)

@assert all(isfinite.(dic["molar_masses"]))
