
# assume 25 Celcius. From PubChem.
# returns solubility in grams / litre.
function solubilityinwater(molecule_name)::Float64
    #

    if Unicode.normalize(molecule_name, casefold = true) == Unicode.normalize("L-Tyrosine", casefold=true)
        # https://pubchem.ncbi.nlm.nih.gov/compound/Tyrosine#section=Melting-Point
        #return 0.453

        # https://www.sigmaaldrich.com/content/dam/sigma-aldrich/docs/Sigma-Aldrich/Product_Information_Sheet/t3754pis.pdf
        return 0.45

    elseif Unicode.normalize(molecule_name, casefold = true) == Unicode.normalize("L-Tryptophan", casefold=true)
        # https://pubchem.ncbi.nlm.nih.gov/compound/Tryptophan#section=Solubility
        return 11.4

    elseif Unicode.normalize(molecule_name, casefold = true) == Unicode.normalize("L-Glutamic acid", casefold=true)
        # https://pubchem.ncbi.nlm.nih.gov/compound/Glutamic-acid#section=Solubility
        return 8570 * 1e-3

    elseif Unicode.normalize(molecule_name, casefold = true) == Unicode.normalize("Caffeine", casefold=true)
        # https://pubchem.ncbi.nlm.nih.gov/compound/Caffeine#section=Solubility
        return (2.16*10+4) * 1e-3
    end

    # TODO
    # 30: Epinephrine
    # solute_concentration = NaN
    # reference_concentration = 0.5
    #
    # 31: (-)-Norepinephrine
    # solute_concentration = NaN
    # reference_concentration = 0.5
    #
    # 32: 3,4-Dihydroxy-L-phenylalanine
    # solute_concentration = NaN
    # reference_concentration = 0.5


    return NaN
end

########### fetch metadata from GISSMO or BMRB.



# look up from JLD.
function getmolarmass(base_dir, target_name::String)
    load_path = joinpath(base_dir, "molar_masses.jld")
    dic = JLD.load(load_path)

    molecule_names = dic["molecule_names"]
    molar_masses = dic["molar_masses"]

    ind = findfirst(xx->xx==target_name, molecule_names)

    if typeof(ind) == Nothing
        return NaN
    end

    return molar_masses[ind]
end

function fetchentrymetadata(entry::String)

    # get molecular weight (natural abundance).
    b = "http://api.bmrb.io/v2/entry/$(entry)"
    HTTP_result = HTTP.request("GET", b)
    json_data = JSON3.read(HTTP_result.body)
    frame_data = json_data[Symbol(entry)][:saveframes]

    ind = findfirst(xx->(xx[:category] == "chem_comp"), frame_data)
    ind2 = findfirst(xx->(xx[1] == "Formula_weight"), frame_data[ind][:tags])
    molar_mass_solute = parse(Float64, frame_data[ind][:tags][ind2][2]) # grams per mol.

    return molar_mass_solute, json_data
end

# was called parseconcentrationfromNMReDATA
function convertsampleconcentration(compound_names,
                                    concentration_value_strings,
                                    concentration_units,
                                    compound_types,
                                    solute_concentration_string,
                                    solute_concentration_unit,
                                    molar_mass_solute::Float64,
                                    solubility_gL_solute::Float64)

    ### case if solute is saturated.
    solute_mM = convertstringtomM(solute_concentration_string,
                                solute_concentration_unit,
                                solubility_gL_solute,
                                molar_mass_solute)

    ### find reference.
    ref_mM = NaN
    ref_name = ""
    molar_mass_ref = NaN
    solubility_gL_ref = NaN

    ind = findfirst(xx->(occursin("Reference", xx) || occursin("reference", xx)), compound_types)
    if length(ind) < 1
        println("Warning, cannot find reference entry in provided sample composition metadata. Assigning NaN as reference concentration.")
        ref_mM = NaN

    else

        # convert reference.
        ref_name = compound_names[ind]
        molar_mass_ref = getmolarmass(ref_name)
        solubility_gL_ref = solubilityinwater(ref_name)

        ref_mM = convertstringtomM(concentration_value_strings[ind],
                                        concentration_units[ind],
                                        solubility_gL_ref,
                                        molar_mass_ref)
    end

    return solute_mM, ref_mM, molar_mass_ref, solubility_gL_ref, ref_name
end

function convertstringtomM( value_string,
                            unit_string,
                            solubility_gL::T,
                            molar_mass::T) where T

    #
    x_mM = NaN

    ret_val = tryparse(T, value_string)

    if typeof(ret_val) != T
        # fall back method: use saturation value.
        x_mM = convertsaturationtomM(solubility_gL, molar_mass)
    else
        x_mM = ret_val
    end

    # if saturated, set unit to a known case of convertconcentrationtomM().
    if occursin("saturated", value_string) || occursin("Saturated", value_string)
        unit_string = "mM"
    end

    ### solute conversion.
    x_mM = convertconcentrationtomM(x_mM, unit_string, molar_mass)

    return x_mM
end

function convertsaturationtomM(solubility_gL::T, molar_mass::T)::T where T
    # # old
    # solubility_g_100mL = solubility_gL/10
    # x_M = solubility_g_100mL/molar_mass # molar concentration.
    # x_mM = x_M*1000 # millimolar concentration.

    x_M = solubility_gL/molar_mass # molar concentration.
    x_mM = x_M*1000 # millimolar concentration.


    return x_mM
end

function convertconcentrationtomM(x::T, unit, molar_mass::T)::T where T <: Real

    if unit == "uM"
        return x * 1e-3

    elseif unit == "%"
        # see http://abacus.bates.edu/~ganderso/biology/resources/dilutions.html#percenttomolarity

        x_M = x*10/molar_mass # molar concentration.
        return x_M*1000 # millimolar concentration.

    elseif unit == "mM"

        return x
    end

    return NaN # unknown case.
end

# assumes temperature is always in Kelvins.
function parseconditionsfromNMReDATA(file_strings::Vector{String})
    # pH.
    header_string = "> <NMREDATA_PH>"

    ind = findfirst(xx->xx==header_string, file_strings)
    @assert typeof(ind) == Int # must not be empty.

    pH = NaN
    val = tryparse(Float64, file_strings[ind+1])
    if typeof(val) == Float64
        pH = val
    end

    # temperature.
    header_string = "> <NMREDATA_TEMPERATURE>"

    ind = findfirst(xx->xx==header_string, file_strings)
    @assert typeof(ind) == Int # must not be empty.

    s = file_strings[ind+1]
    tokens = split(s, ' ')

    temperature = NaN
    val = tryparse(Float64, tokens[1])
    if typeof(val) == Float64
        temperature = val
    end

    return pH, temperature
end


function parsecompositionfromNMReDATA(file_strings::Vector{String})

    ### find the solute concentration.
    header_string = "> <NMREDATA_CONCENTRATION>"

    ind = findfirst(xx->xx==header_string, file_strings)
    @assert typeof(ind) == Int # must not be empty.

    # split and parse.
    s = file_strings[ind+1]
    tokens = split(s, ' ')

    ### solute concentration. Possibly numerical (string) value or "saturated".
    # current strategy: when in doubt, use 100%.
    solute_concentration_unit = tokens[2]
    solute_concentration_string = tokens[1]


    ### find the mixture composotion info
    header_string = "> <NMREDATA_SOLVENT>"

    ind = findfirst(xx->xx==header_string, file_strings)
    @assert typeof(ind) == Int # must not be empty.

    # split and parse.
    s = file_strings[ind+1]
    tokens = split(s, '/')
    tmp = tokens[1:end-1]

    tokens = split(tokens[end], ' ')

    compound_names = [tmp; tokens[1]]

    concentration_values = split(tokens[2], ':')

    concentration_units = split(tokens[3], ':')
    compound_types = split(tokens[4], ':')
    @assert length(compound_names) == length(concentration_values) == length(concentration_units) == length(compound_types)



    return compound_names, concentration_values, concentration_units, compound_types, solute_concentration_string, solute_concentration_unit
end


### convert concentration values to millimolar.
function converttomM(   compound_types,
                        concentration_units,
                        concentration_values,
                        molar_mass_solute::T,
                        molar_mass_solvent::T,
                        molar_mass_ref::T) where T

    concentration_values_mM = collect( parse(T, concentration_values[i]) for i = 1:length(concentration_values))

    for i = 1:length(compound_types)
        cv = concentration_values_mM[i]

        molar_mass = molar_mass_solute
        if compound_types[i] == "Solvent"
            molar_mass = molar_mass_solvent

        elseif compound_types[i] == "Reference"
            molar_mass = molar_mass_ref
        end

        if compound_types[i] == "Reference" || compound_types[i] == "Solvent"

            if concentration_units[i] == "uM"
                concentration_values_mM[i] = cv/1000
            elseif concentration_units[i] == "%"
                # see http://abacus.bates.edu/~ganderso/biology/resources/dilutions.html#percenttomolarity

                concentration_values_mM[i] = cv*10/molar_mass # molar concentration.
                concentration_values_mM[i] = concentration_values_mM[i]*1000 # millimolar concentration.
            end
        else
            concentration_values_mM[i] = NaN
        end
    end

    return concentration_values_mM
end
