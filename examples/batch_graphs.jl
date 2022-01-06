# run batch_molar_mass.jl first
# re-run this script if there are network issues.
#   The GISSMO server might be blocking you for downloading so many entries.
#   This script will skip downloading entries that were already downloaded.


import GISSMOReader
#import DelimitedFiles

#include("../src/coupling/parse2.jl")



entries = GISSMOReader.getGISSMOentriesall()
#entries = getGISSMOentries(["L-Lysine";])

# storage for the experiment 1D 1H files.
#base_dir = "/home/roy/MEGAsync/data/molecules"
base_dir = "/home/roy/del/base"

# storage for the extract J-coupling and chemical shift values, stored in a JLD file per GISSMO entry.
save_dir = "/home/roy/MEGAsync/inputs/NMR/molecules/"


δ_lb = 0.1
δ_ub = 0.1
unique_cs_tol = 1e-6


GISSMOReader.downloadGISSMOentries(entries,
base_dir,
save_dir;
δ_lb = δ_lb,
δ_ub = δ_ub,
unique_cs_tol = unique_cs_tol)

# julia batch_graphs.jl > /home/roy/MEGAsync/inputs/NMR/molecules/output.txt
