#### new scheme.

#STR_file_root = "/home/roy/MEGAsync/data/NMR/GISSMO_str_files"

include("../src/GISSMOReader.jl")
import .GISSMOReader

import JSON
import InfoZIP

base_URL = "http://gissmo.nmrfam.wisc.edu/"
unique_cs_tol = 1e-6
zero_tol_sigdigits = 6

#JLD_save_folder = "/home/roy/Documents/repo/NMRData/input/compounds"
coupling_info_folder = "/home/roy/Documents/repo/NMRData/intput/coupling_info"
edata_save_folder = "/home/roy/MEGAsync/data/NMR/GISSMO_edata"

JSON_folder = "/home/roy/Documents/repo/GISSMOReader.jl/GISSMO_info"
JSON_filename = "GISSMO_entries_metadata_Jan_2022.json"
JSON_load_path = joinpath(JSON_folder, JSON_filename)

db_dict = JSON.parse(open(JSON_load_path))


# ###### single entry.
#
# #entry_key = "http://gissmo.nmrfam.wisc.edu/entry/bmse000861/simulation_1"
# entry_key = "http://gissmo.nmrfam.wisc.edu/entry/Maybridge_Ro3_Fragment_10_A08/simulation_1"
# dict = db_dict[entry_key]
#
# molecule_name = dict["molecule_name"]
# eData_URL = "$(base_URL)$(dict["GISSMO experiment URL"])"
#
# simulation_label = split(entry_key, "/")[end]
# entry_label = split(entry_key, "/")[end-1]
#
# coupling_info_folder = "/home/roy/Documents/repo/NMRData/intput/coupling_info"
# edata_save_folder = "/home/roy/MEGAsync/data/NMR/GISSMO_edata"
#
# GISSMOReader.downloadGISSMOedata(molecule_name,
#     simulation_label,
#     entry_label,
#     coupling_info_folder,
#     edata_save_folder,
#     eData_URL;
#     unique_cs_tol = 1e-6,
#     zero_tol_sigdigits = 6)
#
#
# # test loading.
# load_path = "/home/roy/Documents/repo/NMRData/coupling_info/Maybridge_Ro3_Fragment_10_A08_simulation_1.json"
# H_IDs, H_css, J_IDs, J_vals = GISSMOReader.loadcouplinginfojson(load_path)
#
# @assert 1==2

############# all entries in the json.

function getedatakeys(JSON_load_path)

    db_dict = JSON.parse(open(JSON_load_path))

    entry_keys = Vector{String}(undef, 0)
    for (entry_key, dict) in db_dict
        push!(entry_keys, entry_key)
    end

    return entry_keys
end

function batchdownloadGISSMOedata(JSON_load_path,
    entry_keys::Vector{String},
    start_index::Int,
    coupling_info_folder,
    edata_save_folder;
    unique_cs_tol = 1e-6,
    zero_tol_sigdigits = 6,
    keep_zip_file = false)

    db_dict = JSON.parse(open(JSON_load_path))

    for i = start_index:length(entry_keys)

        println("Starting entry_key $(i)")

        entry_key = entry_keys[i]
        dict = db_dict[entry_key]

        molecule_name = dict["molecule_name"]
        eData_URL = "$(base_URL)$(dict["GISSMO experiment URL"])"

        simulation_label = split(entry_key, "/")[end]
        entry_label = split(entry_key, "/")[end-1]



        GISSMOReader.downloadGISSMOedata(molecule_name,
            simulation_label,
            entry_label,
            coupling_info_folder,
            edata_save_folder,
            eData_URL;
            unique_cs_tol = unique_cs_tol,
            zero_tol_sigdigits = zero_tol_sigdigits,
            keep_zip_file = keep_zip_file)
    end
end

key_list = getedatakeys(JSON_load_path)

st_ind = 1
batchdownloadGISSMOedata(JSON_load_path,
    key_list,
    st_ind,
    coupling_info_folder,
    edata_save_folder;
    unique_cs_tol = 1e-6,
    zero_tol_sigdigits = 6,
    keep_zip_file = false)
