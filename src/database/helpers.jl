
#### query.
# # switch to exact search.
# function searchmoleculelist(list_names::Vector{String},
#                 target_strings::Vector{String})
#
#     M = length(target_strings)
#     N = length(list_names)
#
#     inds = collect(1:N)
#
#     target_exist_flags = falses(M)
#
#
#     keep_flags = falses(N)
#     for i = 1:N
#         db_string = Unicode.normalize(list_names[i], casefold=true)
#
#         for m = 1:M
#
#             search_string = Unicode.normalize(target_strings[m], casefold=true)
#
#             if db_string == search_string
#                 keep_flags[i] = true
#                 target_exist_flags[m] = true
#             end
#         end
#     end
#
#     return inds[keep_flags], target_exist_flags
# end

# single entry.
function searchmoleculelist(list_names::Vector{String},
                target_string::String)

    N = length(list_names)

    inds = collect(1:N)

    target_exist_flag = false

    search_string = Unicode.normalize(target_string, casefold=true)

    keep_flags = falses(N)
    for i = 1:N

        db_string = Unicode.normalize(list_names[i], casefold=true)

        if search_string == db_string
            keep_flags[i] = true
            target_exist_flag = true
        end

    end

    return inds[keep_flags], target_exist_flag
end

########## molar mass look up.

function getmolarmasssolvent(name)::Float64
    #

    if name == "D2O"
        return 20.028 #g/mol, D2O.
    end

    return NaN
end

# in grams/mol.
function getmolarmass(name)::Float64
    #

    if name == "DSS"
        return 218.32177
    end

    if name == "TSP"
        return 163.94
    end

    if name == "TMSP"
        return 146.26
    end

    if name == "TMS"
        return 88.22
    end

    return NaN
end



function loadGISSMOmolecule(base_path, file_name)

   load_path = joinpath(base_path, "$(file_name).jld")

   dic = JLD.load(load_path)

   css_sys = dic["css_sys"]
   J_IDs_sys = dic["J_IDs_sys"]
   J_vals_sys = dic["J_vals_sys"]
   cs_singlets = dic["cs_singlets"]

   J_lb_sys = dic["J_lb_sys"]
   J_ub_sys = dic["J_ub_sys"]
   cs_lb_sys = dic["cs_lb_sys"]
   cs_ub_sys = dic["cs_ub_sys"]
   cs_LUT = dic["cs_LUT"]

   ref_mM = dic["ref_mM"]
   solute_mM = dic["solute_mM"]
   solubility_gL_solute = dic["solubility_gL_solute"]
   solubility_gL_ref = dic["solubility_gL_ref"]
   molar_mass_ref = dic["molar_mass_ref"]
   molar_mass_solute = dic["molar_mass_solute"]
   pH = dic["pH"]
   temperature = dic["temperature"]
   ref_molecule_name = dic["ref_molecule_name"]

   return css_sys, J_IDs_sys, J_vals_sys, cs_singlets,
            J_lb_sys, J_ub_sys, cs_lb_sys, cs_ub_sys, cs_LUT,
            ref_mM, solute_mM, solubility_gL_solute,
            solubility_gL_ref, molar_mass_ref,
            molar_mass_solute, pH, temperature, ref_molecule_name
end
