include("../src/GISSMOReader.jl")
import .GISSMOReader

import Unicode
import JSON
import JSON3

output_folder = "/home/roy/Documents/repo/NMRData/input/compound_mapping"
mkpath(output_folder)
save_path = joinpath(output_folder, "GISSMO_names.json")

JSON_folder = "/home/roy/Documents/repo/GISSMOReader.jl/GISSMO_info"
JSON_filename = "GISSMO_entries_metadata_Jan_2022.json"
JSON_load_path = joinpath(JSON_folder, JSON_filename)



function createGISSMOcompoundnamesdict(JSON_load_path)
    db_dict = JSON.parse(open(JSON_load_path))

    dict_out = Dict()
    for (key, entry_dict) in db_dict

        compound_name = entry_dict["molecule_name"]

        U = split(key, "/")
        entry_name = U[end-1]
        simulation_name = U[end]

        key_out = "$(compound_name) $(entry_name) $(simulation_name)"
        file_name = "$(entry_name)_$(simulation_name).json"
        notes = key

        dict_out[key_out] = Dict( "file name" => file_name, "notes" => notes)
    end

    return dict_out
end

dict_out = createGISSMOcompoundnamesdict(JSON_load_path)

stringdata = JSON.json(dict_out)

open(save_path, "w") do f
    JSON3.pretty(f, stringdata)
    println(f)
end

#### common compounds.



key_names = ["DSS";
"L-Tyrosine";
"L-Leucine";
"Glycine";
"L-Isoleucine";
"L-Valine";
"L-Tryptophan";
"L-Threonine";
"L-Serine";
"L-Phenylalanine";
"L-Methionine";
"L-Glutamine";
"L-(+) Lactic acid";
"D-(+)-Glucose";
"Caffeine";
"L-Histidine";
"L-Proline";
"Choline";
"L-Asparagine";
"L-Aspartic acid";
"L-Lysine";
"L-Alanine";
"L-Arginine";
"Folate";
"L-Cysteine";
"Putrescine";
"Succinic acid";
"Purine";
"L-Glutamic acid";
"Ethanol";
"Epinephrine";
"Norepinephrine";
"3,4-Dihydroxy-L-phenylalanine";
"Dopamine";
"5-Hydroxy-L-tryptophan";
"Serotonin";
"L-Kynurenine";
"Tryptamine";
"Creatine";
"Creatinine";
"Myo-Inositol";
"Glycerol";
"Gamma-Aminobutyric acid";
"3-Hydroxybutyrate";
"HEPES";
"Acetylcholine";
"DL-Alanine";
"Methanol";
"Triethanolamine";
"DL-Serine";
"Diethanolamine";
"N-Alpha-Acetyl-L-Lysine";
"Phenylacetylglycine";
"2-Phenylethanol"]

search_names = ["DSS";
"L-Tyrosine";
"L-Leucine";
"Glycine";
"L-Isoleucine";
"L-Valine";
"L-Tryptophan";
"L-Threonine";
"L-Serine";
"L-Phenylalanine";
"L-Methionine";
"L-Glutamine";
"Lactic acid";
"Glucose";
"Caffeine";
"L-Histidine";
"L-Proline";
"Choline";
"L-Asparagine";
"L-Aspartic acid";
"L-Lysine";
"L-Alanine";
"L-Arginine";
"Folate";
"L-Cysteine";
"Putrescine";
"Succinic acid";
"Purine";
"L-Glutamic acid";
"Ethanol";
"Epinephrine";
"Norepinephrine";
"Dihydroxy-L-phenylalanine";
"Dopamine";
"5-Hydroxy-L-tryptophan";
"Serotonin";
"L-Kynurenine";
"Tryptamine";
"Creatine";
"Creatinine";
"Myo-Inositol";
"Glycerol";
"Gamma-Aminobutyric acid";
"3-Hydroxybutyrate";
"HEPES";
"Acetylcholine";
"DL-Alanine";
"Methanol";
"Triethanolamine";
"DL-Serine";
"Diethanolamine";
"Acetyl-L-lysine";
"Phenylacetylglycine";
"2-Phenylethanol"]


function createGISSMOcompoundnamesdict(JSON_load_path, search_names::Vector{String})
    db_dict = JSON.parse(open(JSON_load_path))

    dict_out = Dict()
    for (key, entry_dict) in db_dict

        compound_name = entry_dict["molecule_name"]

        U = split(key, "/")
        entry_name = U[end-1]
        simulation_name = U[end]

        key_out = "$(compound_name) $(entry_name) $(simulation_name)"
        file_name = "$(entry_name)_$(simulation_name).json"
        notes = key

        compound_name0 = Unicode.normalize(compound_name, casefold = true)

        if any(occursin(Unicode.normalize(search_names[i], casefold = true), compound_name0) for i = 1:length(search_names))

            dict_out[key_out] = Dict( "file name" => file_name, "notes" => notes)
        end
    end

    return dict_out
end

dict_out_subset = createGISSMOcompoundnamesdict(JSON_load_path, search_names)


stringdata = JSON.json(dict_out_subset)

save_path = joinpath(output_folder, "subset_names.json")
open(save_path, "w") do f
    JSON3.pretty(f, stringdata)
    println(f)
end
