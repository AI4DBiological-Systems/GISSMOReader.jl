
function constructBMRBlist(dict)
    N = length(dict)
    entries = Vector{String}(undef, N)
    molecule_names = Vector{String}(undef, N)

    i = 0
    for (key, value) in dict
        #
        #println("typeof(key) = ", typeof(key))
        #println("typeof(value[1]) = ", typeof(value[1]))

        i += 1
        entries[i] = String(key)
        molecule_names[i] = value[1]
    end

    return entries, molecule_names
end



### fetch the current list of molecules in the BMRB metabolomics database.
function getchBMRBmetabolomicslist(folder_path::String = joinpath(homedir(),"MEGAsync/data/NMR/BMRB"))

    a = "http://api.bmrb.io/current/search/get_all_values_for_tag/Chem_comp.Name?database=metabolomics"
    HTTP_result = HTTP.request("GET", a)
    dict = JSON3.read(HTTP_result.body)
    BMRB_list, list_names = constructBMRBlist(dict)

    # store queried results.
    isdir(folder_path) || mkdir(folder_path)

    save_path = joinpath(folder_path, "list.jld")
    JLD.save(save_path, "BMRB_list", BMRB_list,
                        "list_names", list_names)

    return BMRB_list, list_names
end




function downloadBMRBdata(target_names::Vector{String},
                BMRB_list::Vector{String},
                list_names::Vector{String},
                GISSMO_dir::String,
                BMRB_zip_dir::String,
                BMRB_data_dir::String;
                required_experiment_substring::String = "1D 1H",
                delay_time = 0.2)

    #
    N_molecules = length(target_names)

    entry_set = Vector{Vector{String}}(undef, N_molecules)
    query_success_flags = falses(N_molecules)
    entry_success_flags_set = Vector{BitArray{1}}(undef, N_molecules)

    for i = 1:N_molecules

        display_string = "Processing $(i)/$(N_molecules), $(target_names[i])"
        println(display_string)

        entry_set[i], query_success_flags[i],
            entry_success_flags_set[i] = downloadBMRBdata(target_names[i],
                BMRB_list,
                list_names,
                GISSMO_dir,
                BMRB_zip_dir,
                BMRB_data_dir;
                required_experiment_substring = required_experiment_substring,
                delay_time = delay_time)

        println()
        println()
    end

    return entry_set, query_success_flags, entry_success_flags_set
end

# for one target molecule.
# delay_time in seconds. To respect REST API rate limit.
function downloadBMRBdata(target_name::String,
                BMRB_list::Vector{String},
                list_names::Vector{String},
                GISSMO_dir::String,
                BMRB_zip_dir::String,
                BMRB_data_dir::String;
                required_experiment_substring::String = "1D 1H",
                delay_time = 0.2)


    ### search for the entries that contain the target name.
    inds, target_exist_flag = searchmoleculelist(list_names, target_name)
    search_result_entries = BMRB_list[inds]

    N = length(search_result_entries)
    entry_success_flags = falses(N)

    for n = 1:N

        entry = search_result_entries[n]
        println("Processing entry ", entry)

        # get metadata, so we can figure out which experiments to keep.
        URLs, all_metadata_JSON, inds_1D1H = queryBMRBexperimentdataURL(entry; target_string = required_experiment_substring)
        metadata_JSON = all_metadata_JSON[inds_1D1H]

        # download, extract, keep only 1D 1H entries. unzip_dir not used.
        unzip_dir = downloadandextractBMRBentry(entry, BMRB_zip_dir, BMRB_data_dir; delay = delay_time)
        entry_dir, experiment_paths, data_move_success_flag = keeponly1D1Hexperiments(entry, BMRB_zip_dir, BMRB_data_dir, unzip_dir, URLs)

        if !data_move_success_flag

            # report bad entry, skip.
            entry_success_flags[n] = false

        else
            # save raw metadata to storage.
            save_path = joinpath(entry_dir, "metadata.json")
            open(save_path, "w") do f
                JSON3.pretty(f, JSON3.write(metadata_JSON))
                println(f)
            end

            # cache the 1D 1H experiment directories for this entry.
            save_path = joinpath(entry_dir, "experiment_dirs.jld")
            JLD.save(save_path, "experiment_dirs", experiment_paths,
                                "molecule_name", target_name)

            ### store the sample metadata for each experiment.
            N_experiments = length(metadata_JSON)
            @assert length(experiment_paths) == N_experiments

            for j = 1:N_experiments

                # concentration.
                X = metadata_JSON[j][:Sample_component]
                solute_mM, ref_mM, ref_molecule_name,
                    solubility_gL_solute, molar_mass_solute,
                    solubility_gL_ref, molar_mass_ref = convertBMRBexperimentconcentrationinfo(X, target_name, GISSMO_dir)

                # conditions.
                pH = metadata_JSON[j][:Sample_condition_variable][:ph]
                temperature = metadata_JSON[j][:Sample_condition_variable][:temperature]

                # store.
                # debug.

                @task begin;
                    save_path = joinpath(experiment_paths[j], "metadata.jld")
                    JLD.save(save_path, "solute_mM", solute_mM,
                                    "ref_mM", ref_mM,
                                    "ref_molecule_name", ref_molecule_name,
                                    "solubility_gL_solute", solubility_gL_solute,
                                    "molar_mass_solute", molar_mass_solute,
                                    "solubility_gL_ref", solubility_gL_ref,
                                    "molar_mass_ref", molar_mass_ref,
                                    "pH", pH,
                                    "temperature", temperature)
                    # debug.
                    println("save_path= ", save_path)
                    

                    save_path = joinpath(experiment_paths[j], "metadata.json")
                    open(save_path, "w") do f
                        JSON3.pretty(f, JSON3.write(metadata_JSON[j]))
                        println(f)
                    end
                end
                sleep(0.1)
            end

            entry_success_flags[n] = true
        end

    end

    return search_result_entries, target_exist_flag, entry_success_flags
end

# base_dir is the parent folder where the BMRB entries can be found.
function loadBMRBdata(entries::Vector{String},
                        base_dir::String)

    #
    full_paths = Vector{String}(undef, 0)

    #BMRB_entries = readdir(base_path)

    # loop over each specified entry.
    for i = 1:length(entries)
        entry = entries[i]

        # see if the entry exists on our local storage.
        entry_path = joinpath(base_dir, entry)
        if isdir(entry_path)

            experiment_tags = readdir(entry_path)

            # loop over the experiment for the current entry.
            entry_paths = collect( joinpath(entry_path, experiment_tags[k]) for k = 1:length(experiment_tags) )

            # push.
            push!(full_paths, entry_paths...)
        end
    end

    # remove paths to files.
    keep_flags = collect( !isfile(full_paths[i]) for i = 1:length(full_paths) )
    full_paths = full_paths[keep_flags]

    return full_paths
end

# X = metadata_JSON[i][:Sample_component]
function convertBMRBexperimentconcentrationinfo(X,
            molecule_name::String, GISSMO_dir::String)

    ### concentration
    ref_mM = NaN
    solute_mM = NaN

    ref_ind = findfirst(xx->(occursin("Reference", xx[:Type]) || occursin("reference", xx[:Type])), X)
    if length(ref_ind) < 1
        println("Warning, cannot find reference entry in provided sample composition metadata. Assigning NaN as reference concentration.")
        ref_mM = NaN
    end

    solute_ind = findfirst(xx->(occursin("Solute", xx[:Type]) || occursin("solute", xx[:Type])), X)
    if length(solute_ind) < 1
        println("Warning, cannot find solute entry in provided sample composition metadata. Assigning NaN as solute concentration.")
        solute_mM = NaN
    end

    A = X[solute_ind]

    solubility_gL_solute = solubilityinwater(molecule_name)
    molar_mass_solute = getmolarmass(GISSMO_dir, molecule_name)

    solute_mM = convertstringtomM(A[:Concentration_val],
                A[:Concentration_val_units], solubility_gL_solute, molar_mass_solute)

    #
    B = X[ref_ind]
    ref_molecule_name = B[:Mol_common_name]

    solubility_gL_ref = solubilityinwater(ref_molecule_name)

    ref_entry = getGISSMOentry(ref_molecule_name)
    molar_mass_ref = getmolarmass(GISSMO_dir, ref_molecule_name)

    ref_mM = convertstringtomM(B[:Concentration_val],
                B[:Concentration_val_units], solubility_gL_ref, molar_mass_ref)

    return solute_mM, ref_mM, ref_molecule_name,
        solubility_gL_solute, molar_mass_solute,
        solubility_gL_ref, molar_mass_ref
end

# download to zip dir and unzip to data dir. delay in seconds.
function downloadandextractBMRBentry(entry, BMRB_zip_dir, BMRB_data_dir; delay = 0.2)

    # download entry to zip directory..
    URL = "https://bmrb.io/metabolomics/mol_summary/zip_entry_directory.php?id=$(entry)"
    dest = joinpath(BMRB_zip_dir, "$(entry).zip")


    t = @task begin; isfile(dest) || download(URL, dest); println("done download ", entry); end
    schedule(t); wait(t)
    sleep(delay) # respect REST API query rate limit.

    # unzip to dataset directory.
    unzip_dir = joinpath(BMRB_data_dir, "$(entry)_unzip")
    t = @task begin; isdir(unzip_dir) || mkdir(unzip_dir); end
    schedule(t); wait(t)

    t = @task begin; InfoZIP.unzip(dest, unzip_dir); println("done unzip ", entry); end
    schedule(t); wait(t)

    return unzip_dir
end

# keep only 1D1H experiments from the unzipped directory.
function keeponly1D1Hexperiments(entry,
            BMRB_zip_dir, BMRB_data_dir, unzip_dir, URLs)

    ### delete directories that aren't 1D 1H.
    # move each of the experiment paths.
    entry_dir = joinpath(BMRB_data_dir, "$(entry)")
    t = @task begin; ispath(entry_dir) || mkpath(entry_dir); end
    schedule(t); wait(t)


    experiment_tmp_paths = collect( getlocal1D1Hexperimentpath(URLs[k], entry , unzip_dir) for k = 1:length(URLs) )
    problem_flags = collect( !ispath(experiment_tmp_paths[i]) for i = 1:length(experiment_tmp_paths) )

    experiment_paths = collect( getlocal1D1Hexperimentpath(URLs[k], entry , entry_dir) for k = 1:length(URLs) )



    if any(problem_flags)

        # zip archive and metadata URL not the same.
        # report bad entry. Remove the unzip folder without moving.
        println("Path issue for ", entry, " does not agree with metadata URL. Skip.")

        println("experiment_tmp_paths ", experiment_tmp_paths)
        println("experiment_paths ", experiment_paths)
        println("URLs ", URLs)
        println("problem_flags ", problem_flags)

        #rm(unzip_dir; force = true, recursive = true)

        return entry_dir, experiment_paths, false
    end

    t = @task begin; mvfolders(experiment_tmp_paths, experiment_paths); println("done moving ", entry); end
    schedule(t); wait(t)

    rm(unzip_dir; force = true, recursive = true)

    return entry_dir, experiment_paths, true
end

######### work-around for experiment ID mismatch with download URL issue.

# if we've X[k][:attribute] has target_string as a sub-string, then return k.
# return a vector of all k's that match.
function searchJSON3arrayformatch(X, attribute::Symbol, target_string::String)
    N = length(X)

    inds = Vector{Int}(undef, 0)

    for i = 1:length(X)
        if occursin(target_string, X[i][attribute])
            push!(inds, i)
        end
    end

    return inds
end

# assumes mkpath(dest) will work.
function mvfolders(experiment_tmp_paths, experiment_paths)

    for k = 1:length(experiment_paths)

        src = experiment_tmp_paths[k]
        dest = experiment_paths[k]

        #println("src ", src)
        #println("dest ", dest)

        t = @task begin; ispath(dest) || mkpath(dest); end
        schedule(t); wait(t)

        mv(src, dest; force = true)
    end
end

# get path to the 1D 1H experiment reference in the URL.
function getlocal1D1Hexperimentpath(URL::String,
            entry::String, entry_dir::String)::String

    # keep only the folder names after the entry in the URL.
    q = split(URL, '/')
    l = findfirst(xx->xx==entry, q)

    # assemble local path.
    experiment_path = entry_dir
    for i = l+1:length(q)
        experiment_path = joinpath(experiment_path, q[i])

        # if !ispath(experiment_path)
        #     # error.
        #     return ""
        # end
    end

    return experiment_path
end


# search for 1D 1H for the current entry.
function queryBMRBexperimentdataURL(entry::String;
                        target_string::String = "1D 1H")
    a = "http://api.bmrb.io/v2/entry/$(entry)/experiments"

    HTTP_result = HTTP.request("GET", a)
    experiments = JSON3.read(HTTP_result.body)

    inds = searchJSON3arrayformatch(experiments, :Name, target_string)

    URLs = Vector{String}(undef, 0)
    for k in inds
        tmp = experiments[k][:Experiment_file]
        inds2 = searchJSON3arrayformatch(tmp, :description, "experiment directory")

        for j in inds2
            push!(URLs, tmp[j][:url])
        end
    end

    return URLs, experiments, inds
end
