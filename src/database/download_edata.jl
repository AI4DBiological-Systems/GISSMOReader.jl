function downloadGISSMOedata(molecule_name::String,
    simulation_label,
    entry_label,
    coupling_info_folder::String,
    edata_save_folder::String,
    eData_URL::String;
    unique_cs_tol::T = 1e-6,
    zero_tol_sigdigits = 6,
    keep_zip_file = false) where T <: Real

    # Download zip file.
    entry_save_folder = joinpath(edata_save_folder, "$(entry_label)_$(simulation_label)")
    mkpath(entry_save_folder)

    file_name = "nmredata.zip"
    dest = joinpath(entry_save_folder, file_name)
    t = @task begin; isfile(dest) || download(eData_URL, dest); end
    schedule(t); wait(t)

    #q = filter(x->endswith(x, ".zip"), readdir(entry_save_folder))

    # unzip.
    InfoZIP.unzip(dest, entry_save_folder)

    # remove zip file.
    if !keep_zip_file
        rm(dest)
    end

    # get SDF filename.
    tmp = filter(isdir, readdir(entry_save_folder; join=true))
    @assert length(tmp) == 1 # should be only one folder, which is the unzipped one.
    folder_path = tmp[1]

    tmp = filter(x->endswith(x, ".sdf"), readdir(folder_path))
    @assert length(tmp) == 1 # should be exactly 1 SDF file.
    SDF_name = tmp[1]
    SDF_path = joinpath(folder_path, SDF_name)

    file_strings = readlines(SDF_path)

    #### find the section on chemical shifts.
    cs_string = "> <NMREDATA_ASSIGNMENT>"

    H_IDs, H_css = parsecssfromNMReDATA(file_strings)

    ##### find the section for J-coupling. J_string = "> <NMREDATA_J>"
    J_IDs, J_vals = GISSMOReader.parseJsfromNMReDATA(file_strings)

    # store.
    mkpath(coupling_info_folder)
    save_path = joinpath(coupling_info_folder, "$(entry_label)_$(simulation_label).json")

    cs_dict = cstodict(H_IDs, H_css)
    J_dict = Jcouplingtodict(J_IDs, J_vals)
    # dict1 = Dict("H_IDs" => H_IDs, "H_css" => H_css,
    #             "J_IDs" => J_IDs, "J_vals" => J_vals)
    dict1 = Dict("chemical shift" => cs_dict, "J-coupling" => J_dict)
    stringdata = JSON.json(dict1)

    open(save_path, "w") do f
        JSON3.pretty(f, stringdata)
        println(f)
    end

    return save_path
end

function Jcouplingtodict(J_IDs::Vector{Tuple{Int,Int}}, J_vals::Vector{T}) where T

    return collect( Dict("ID1" => J_IDs[i][1], "ID2" => J_IDs[i][2], "value" => J_vals[i]) for i = 1:length(J_vals) )
end

function cstodict(H_IDs::Vector{Int}, H_css::Vector{T}) where T

    return collect( Dict("ID" => H_IDs[i], "value" => H_css[i]) for i = 1:length(H_css) )
end

# come back to this later.
function tempfunc()
    #J_string = "> <NMREDATA_J>"
    J_IDs, J_vals = parseJsfromNMReDATA(file_strings)

    ### attempt on auto-grouping.

    J_IDs_sys, J_vals_sys, H_IDs_sys = partitionJcouplings(J_IDs, J_vals)

    css_sys = getcssys(H_IDs_sys, H_IDs, H_css)

    𝐽_IDs_sys = getJsys(J_IDs_sys, J_vals_sys, H_IDs_sys)

    H_singlets, cs_singlets,
        cs_singlets_compact = createsingletsysems(H_IDs, H_IDs_sys,
                                    H_css;
                                    zero_tol = unique_cs_tol)
    #
    ### remove spin groups that only have one unique chem shift, but have J-coupling between the atoms.
    # Each of such groups are are really one singlet group.
    removeisolatedcs!(css_sys,
                    cs_singlets,
                    cs_singlets_compact,
                    H_singlets,
                    H_IDs_sys,
                    J_IDs_sys,
                    J_vals_sys,
                    𝐽_IDs_sys;
                    zero_tol_sigdigits = zero_tol_sigdigits)
    #
    println("Entry and simulation: ", "$(entry_label)_$(simulation_label)")
    println("cs_singlets: ", cs_singlets)
    println("cs_singlets_compact: ", cs_singlets_compact)
    println("H_singlets: ", H_singlets)

    println("css_sys: ", css_sys)
    println("J_vals_sys: ", J_vals_sys)
    println("𝐽_IDs_sys: ", 𝐽_IDs_sys)
    println()
    #### end code for moving spin group to singlet group.

    #
    cs_LUT, p_cs_sys = constructLUTcss(css_sys; zero_tol = unique_cs_tol)


    cs_len_sys = getcslengthfromLUT(cs_LUT)

    # ### cs and J-values upper and lower bounds.
    # perturblbfunc = xx->additiveperturbation(xx, δ_lb)
    # perturbubfunc = xx->additiveperturbation(xx, δ_ub)
    #
    # #
    # J_lb_sys = getJperturbed(J_vals_sys, perturblbfunc)
    # J_ub_sys = getJperturbed(J_vals_sys, perturbubfunc)
    #
    # cs_lb_sys = getJperturbed(css_sys, perturblbfunc)
    # cs_ub_sys = getJperturbed(css_sys, perturbubfunc)

    # store.
    save_path = joinpath(JLD_save_folder, "$(entry_label)_$(simulation_label).jld")

    JLD.save(save_path,
        "css_sys", css_sys,
        "J_IDs_sys", 𝐽_IDs_sys,
        "J_vals_sys", J_vals_sys,
        "cs_singlets", cs_singlets,
        "cs_LUT", cs_LUT)
end
