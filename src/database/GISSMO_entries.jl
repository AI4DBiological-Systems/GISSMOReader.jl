
#### used in batch_graphs.jl.

mutable struct GISSMOEntryType{}
    #
    entry::String
    molecule_name::String
    eDATA_URL::String # includes FID data, cs, J-coupling values.
    STAR_URL::String # contains concentration info.

    κ_range
    # each row corresponds to a spin_group.
    # each consecutive pairs of columns correspond to the start
    #   and end of an interval for an unique κ.
    # a negative number means the entire spin group has one κ.
    # see cs_config.txt to help visualize κ_range.

    ι_range::Vector{Tuple{Float64,Float64}}
    #N_protons::Int
end
# TODO update graphs and batch_graphs to pull concentration data from STAR_URL.
#  near Sample_component.Sample_ID

function getGISSMOentriesall()::Vector{GISSMOEntryType}

    entries = Vector{GISSMOEntryType}(undef, 0)

    tmp = GISSMOEntryType("bmse000795",
        "DSS",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000795/simulation_1/bmse000795_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000795/simulation_1/bmse000795-simulation_1.str",
        ((-1.2),(-1.2)),
        [(-0.2, 3.1);] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000051",
        "L-Tyrosine",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000051/simulation_1/bmse000051_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000051/simulation_1/bmse000051-simulation_1.str",
        ((-1.2),(-1.2)),
        [(2.5, 3.6); (6.4, 7.2)] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000042",
        "L-Leucine",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000042/simulation_1/bmse000042_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000042/simulation_1/bmse000042-simulation_1.str",
        ((0.7, 1.3, 2.5, 3.9),),
        [(2.5, 3.6); (6.4, 7.2)] )
    push!(entries, tmp)

    # tmp = GISSMOEntryType("bmse000089",
    #     "Glycine",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000089/simulation_1/bmse000089_simulation_1_nmredata.zip",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000089/simulation_1/bmse000089-simulation_1.str")
    # push!(entries, tmp)
    tmp = GISSMOEntryType("bmse000977",
        "Glycine",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000977/simulation_1/bmse000977_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000977/simulation_1/bmse000977-simulation_1.str",
        ((-1.2),),
        [(3.3, 3.7);] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000041",
        "L-Isoleucine",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000041/simulation_1/bmse000041_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000041/simulation_1/bmse000041-simulation_1.str",
        ((0.8, 1.65, 3.0, 3.8),),
        [(0.8, 3.8);] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000860",
        "L-Valine",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000860/simulation_1/bmse000860_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000860/simulation_1/bmse000860-simulation_1.str",
        ((0.7, 1.5, 2.8, 3.7),),
        [(0.7, 3.7);] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000050",
        "L-Tryptophan",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000050/simulation_1/bmse000050_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000050/simulation_1/bmse000050-simulation_1.str",
        ((-1.2),(-1.2),(-1.2)),
        [(3.0, 4.2); (6.9, 8.0)] )
    push!(entries, tmp)

    # tmp = GISSMOEntryType("bmse000049",
    #     "L-Threonine",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000049/simulation_1/bmse000049_simulation_1_nmredata.zip",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000049/simulation_1/bmse000049-simulation_1.str",
    #     [1.0 2.5 3.8 3.9 4.5 ])
    # push!(entries, tmp)
    tmp = GISSMOEntryType("bmse000859",
        "L-Threonine",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000859/simulation_1/bmse000859_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000859/simulation_1/bmse000859-simulation_1.str",
        ((0.9, 2.2, 4.4),),
        [(1.1, 1.5); (3.4, 4.4)] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000867",
        "L-Serine",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000867/simulation_1/bmse000867_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000867/simulation_1/bmse000867-simulation_1.str",
        ((3.7, 3.88, 4.02),),
        [(3.7, 4.02);] )
    push!(entries, tmp)
    # [2/3; 1.0] # κ.

    tmp = GISSMOEntryType("bmse000900",
        "L-Phenylalanine",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000900/simulation_1/bmse000900_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000900/simulation_1/bmse000900-simulation_1.str",
        ((3.0, 3.5, 4.2), (7.2, 7.6)),
        [(3.0, 4.1); (7.1, 7.6)] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000044",
        "L-Methionine",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000044/simulation_1/bmse000044_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000044/simulation_1/bmse000044-simulation_1.str",
        ((1.9, 3.66, 4.0), (-1.2)),
        [(1.9, 4.0);] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000038",
        "L-Glutamine",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000038/simulation_1/bmse000038_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000038/simulation_1/bmse000038-simulation_1.str",
        ((-1.2),),
        [(1.82, 2.7); (3.6, 3.9)] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000208",
        "L-(+) Lactic acid",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000208/simulation_1/bmse000208_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000208/simulation_1/bmse000208-simulation_1.str",
        ((1.0, 2.8, 4.4),),
        [(1.23, 1.4);] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000015",
        "D-(+)-Glucose",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000015/simulation_1/bmse000015_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000015/simulation_1/bmse000015-simulation_1.str",
        ((3.0, 3.3, 4.0, 5.0), (3.0, 4.5, 5.5)),
        [(3.15, 3.96); (4.57, 4.7); (5.16, 5.28)] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000206",
        "Caffeine",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000206/simulation_1/bmse000206_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000206/simulation_1/bmse000206-simulation_1.str",
        ((-1.2), (-1.2), (-1.2), (-1.2)),
        [(3.1, 3.55); (3.73, 4.0); (7.75, 8.0)] )
    push!(entries, tmp)

    # tmp = GISSMOEntryType("bmse000039",
    #     "L-Histidine",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000039/simulation_1/bmse000039_simulation_1_nmredata.zip",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000039/simulation_1/bmse000039-simulation_1.str",
    #     ((-1.2), (-1.2)),
    #     [(3.0, 3.35); (3.8, 4.17); (6.83, 8.2)] )
    # push!(entries, tmp)
    tmp = GISSMOEntryType("bmse000976",
        "L-Histidine",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000976/simulation_1/bmse000976_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000976/simulation_1/bmse000976-simulation_1.str",
        #((-1.2), (-1.2), (-1.2)),
        ((3.0, 3.5, 4.2), (6.9, 7.4), (7.8, 8.3)),
        [(3.0, 3.35); (3.8, 4.17); (6.83, 8.2)] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000947",
        "L-Proline",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000947/simulation_1/bmse000947_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000947/simulation_1/bmse000947-simulation_1.str",
        ((1.5, 3.7, 4.6),),
        [(1.85, 2.15); (2.25, 2.42); (3.24, 3.48); (4.05, 4.18) ] )
    push!(entries, tmp)

    # tmp = GISSMOEntryType("bmse000285", # 400 MHz.
    #     "Choline",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000285/simulation_1/bmse000285_simulation_1_nmredata.zip",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000285/simulation_1/bmse000285-simulation_1.str")
    # push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000953", # 600 MHz.
        "Choline",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000953/simulation_1/bmse000953_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000953/simulation_1/bmse000953-simulation_1.str",
        ((3.2, 3.8, 4.1), (-1.2)),
        [(2.97, 3.67); (3.96, 4.16)] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000912",
        "L-Asparagine",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000912/simulation_1/bmse000912_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000912/simulation_1/bmse000912-simulation_1.str",
        ((2.74, 3.5, 4.1),),
        [(2.68, 3.15); (3.82, 4.16)] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000875",
        "L-Aspartic acid",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000875/simulation_1/bmse000875_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000875/simulation_1/bmse000875-simulation_1.str",
        ((2.4, 3.4, 4.4),),
        [(2.5, 3.0); (3.68, 4.06)] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000914",
        "L-Lysine",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000914/simulation_1/bmse000914_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000914/simulation_1/bmse000914-simulation_1.str",
        ((-1.2),),
        [(1.3, 2.0); (2.95, 3.09); (3.66, 3.83)] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000028",
        "L-Alanine",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000028/simulation_1/bmse000028_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000028/simulation_1/bmse000028-simulation_1.str",
        ((-1.2),),
        [(1.2, 1.6); (3.57, 3.99)] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000899",
        "L-Arginine",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000899/simulation_1/bmse000899_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000899/simulation_1/bmse000899-simulation_1.str",
        ((-1.2),),
        [(1.52, 2.05); (3.12, 3.38); (3.65, 3.87)] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000299",
        "Folate",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000299/simulation_1/bmse000299_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000299/simulation_1/bmse000299-simulation_1.str",
        ((-1.2),(-1.2),(-1.2),(-1.2),(-1.2)),
        [(1.82, 2.62); (3.62, 4.4); (5.98, 6.46); (7.2, 7.52); (8.25, 8.52)] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000975",
        "L-Cysteine",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000975/simulation_1/bmse000975_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000975/simulation_1/bmse000975-simulation_1.str",
        ((2.9, 3.5, 4.1),),
        [(2.95, 3.12); (3.86, 4.04)] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000109",
        "Putrescine",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000109/simulation_1/bmse000109_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000109/simulation_1/bmse000109-simulation_1.str",
        ((-1.2),),
        [(1.56, 1.95); (2.85, 3.24)] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000183",
        "Succinic acid",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000183/simulation_1/bmse000183_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000183/simulation_1/bmse000183-simulation_1.str",
        ((-1.2),),
        [(2.22, 2.54);] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000454",
        "Purine",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000454/simulation_1/bmse000454_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000454/simulation_1/bmse000454-simulation_1.str",
        ((-1.2),(-1.2)),
        [(8.41, 8.61); (8.7, 9.09)] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000913",
        "L-Glutamic acid",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000913/simulation_1/bmse000913_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000913/simulation_1/bmse000913-simulation_1.str",
        ((1.9, 3.0, 4.13),),
        [(1.93, 2.54); (3.61, 3.84)] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000976",
        "Ethanol",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000297/simulation_1/bmse000297_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000297/simulation_1/bmse000297-simulation_1.str",
        #((-1.2),),
        ((1.0, 2.5, 4.0),),
        [(1.0, 1.5); (3.0, 4.0)] )
    push!(entries, tmp)

    # tmp = GISSMOEntryType("bmse000316",
    #     "Epinephrine",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000316/simulation_1/bmse000316_simulation_1_nmredata.zip",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000316/simulation_1/bmse000316-simulation_1.str",
    #     ((3.0, 4.0, 5.2), (-1.2), (-1.2)) )
    # push!(entries, tmp)
    #
    # tmp = GISSMOEntryType("bmse000404",
    #     "(-)-Norepinephrine",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000404/simulation_1/bmse000404_simulation_1_nmredata.zip",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000404/simulation_1/bmse000404-simulation_1.str",
    #     ((-1.2), (2.9, 4.0, 5.1)) )
    # push!(entries, tmp)
    #
    # # L-DOPA
    # tmp = GISSMOEntryType("bmse000322",
    #     "3,4-Dihydroxy-L-phenylalanine",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000322/simulation_1/bmse000322_simulation_1_nmredata.zip",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000322/simulation_1/bmse000322-simulation_1.str",
    #     ((-1.2),(-1.2)) )
    # push!(entries, tmp)

    # tmp = GISSMOEntryType("bmse000933",
    #     "Dopamine",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000933/simulation_1/bmse000933_simulation_1_nmredata.zip",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000933/simulation_1/bmse000933-simulation_1.str")
    # push!(entries, tmp)
    #
    # tmp = GISSMOEntryType("bmse000457",
    #     "5-Hydroxy-L-tryptophan",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000457/simulation_1/bmse000457_simulation_1_nmredata.zip",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000457/simulation_1/bmse000457-simulation_1.str")
    # push!(entries, tmp)
    #
    # tmp = GISSMOEntryType("bmse000757",
    #     "Serotonin",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000757/simulation_1/bmse000757_simulation_1_nmredata.zip",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000757/simulation_1/bmse000757-simulation_1.str")
    # push!(entries, tmp)
    #
    # tmp = GISSMOEntryType("bmse000172",
    #     "L-Kynurenine",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000172/simulation_1/bmse000172_simulation_1_nmredata.zip",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000172/simulation_1/bmse000172-simulation_1.str")
    # push!(entries, tmp)
    #
    # tmp = GISSMOEntryType("bmse000207",
    #     "Tryptamine",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000207/simulation_1/bmse000207_simulation_1_nmredata.zip",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000207/simulation_1/bmse000207-simulation_1.str")
    # push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000950",
        "Creatine",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000950/simulation_1/bmse000950_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000950/simulation_1/bmse000950-simulation_1.str",
        ((-1.2),(-1.2)),
        [(2.96, 3.09); (3.85, 3.98)] )
    push!(entries, tmp)

    tmp = GISSMOEntryType("bmse000155",
        "Creatinine",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000155/simulation_1/bmse000155_simulation_1_nmredata.zip",
        "http://gissmo.nmrfam.wisc.edu/entry/bmse000155/simulation_1/bmse000155-simulation_1.str",
        ((-1.2),(-1.2)),
        [(2.95, 3.1); (3.97, 4.1)] )
    push!(entries, tmp)

    # # uses TMS!
    # tmp = GISSMOEntryType("bmse000923",
    #     "L-Thyroxine",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000923/simulation_1/bmse000923_simulation_1_nmredata.zip",
    #     "http://gissmo.nmrfam.wisc.edu/entry/bmse000923/simulation_1/bmse000923-simulation_1.str")
    # push!(entries, tmp)

    return entries
end

function getGISSMOentriesDMEM()::Vector{GISSMOEntryType}

    molecule_names = Vector{String}(undef, 0)
    push!(molecule_names, "L-Tyrosine")
    push!(molecule_names, "L-Leucine")
    push!(molecule_names, "Glycine")
    push!(molecule_names, "L-Isoleucine")
    push!(molecule_names, "L-Valine")
    push!(molecule_names, "L-Tryptophan")
    push!(molecule_names, "L-Threonine")
    push!(molecule_names, "L-Serine")
    push!(molecule_names, "L-Phenylalanine")
    push!(molecule_names, "L-Methionine")
    push!(molecule_names, "L-Glutamine")
    push!(molecule_names, "D-(+)-Glucose")
    push!(molecule_names, "L-Histidine")
    push!(molecule_names, "Choline")
    push!(molecule_names, "L-Lysine")
    push!(molecule_names, "L-Arginine")

    return getGISSMOentries(molecule_names)
end

function extractfields(X::Vector{GISSMOEntryType}, z::String)::Vector{String}

    return collect( getproperty(X[i], Symbol(z)) for i = 1:length(X))
end
# example usage:
# DMEM_entries = getGISSMOentriesDMEM()
# A = extractfields(DMEM_entries, "molecule_name")
# B = extractfields(DMEM_entries, "entry")
# display([A B])

function exportallGISSMOnamestofile(save_path)
    entries = getGISSMOentriesall()
    #molecule_names = collect( entries[i].molecule_name for i = 1:length(entries))
    molecule_names = extractfields(entries, "molecule_name")

    sort!(molecule_names)

    DelimitedFiles.writedlm(save_path, molecule_names)
end

function getGISSMOentry(target_name::String)::GISSMOEntryType

    entries = getGISSMOentriesall()
    #molecule_names = collect( entries[i].molecule_name for i = 1:length(entries))
    molecule_names = extractfields(entries, "molecule_name")

    inds, _ = searchmoleculelist(molecule_names, target_name)

    k = 0
    if length(inds) > 1
        k = inds[1]
        message_string = "Warning, $(target_name) found multiple entries. Using the first one: entry $(entries[k].entry)."
        println(message_string)

    elseif isempty(inds)
        message_string = "Warning, $(target_name) did not find any entries."
        println(message_string)

        return GISSMOEntryType("", "", "", "")

    else
        k = inds[1]
    end

    return entries[k]
end


function getGISSMOentries(target_names::Vector{String})::Vector{GISSMOEntryType}

    entries = collect( getGISSMOentry(target_names[i]) for i = 1:length(target_names) )

    return entries
end



##### for batch_graphs.jl.

function downloadGISSMOentries(entries::Vector{GISSMOEntryType},
                            base_dir::String,
                            save_dir::String;
                            δ_lb::T = 0.1,
                            δ_ub::T = 0.1,
                            unique_cs_tol::T = 1e-6,
                            delay_time = 0.2,
                            default_cs_delta::T = 0.1) where T <: Real
    #
    N = length(entries)

    out_string = Vector{String}(undef, 0)

    for i = 1:N

        print_string = "$(i): $(entries[i].molecule_name)"
        println(print_string)

        sleep(delay_time) # replace with retry later.
        css_sys, cs_singlets_compact = downloadGISSMOentry(entries[i],
            base_dir, save_dir;
                δ_lb = δ_lb,
                δ_ub = δ_ub,
                unique_cs_tol = unique_cs_tol)

        println()

        addentrydefaultconfigtobuffer!(out_string,
                entries[i].molecule_name,
                css_sys,
                cs_singlets_compact;
                default_cs_delta = default_cs_delta)
    end

    # special cases: add D2O.
    entry = "0"
    molecule_name = "D2O"
    singlet_ppm = 4.8
    css_sys, cs_singlets_compact = specialsingletentry(entry,
        molecule_name,
        base_dir,
        save_dir,
        singlet_ppm;
            δ_lb = δ_lb,
            δ_ub = δ_ub,
            unique_cs_tol = unique_cs_tol)

    # add to buffer.
    addentrydefaultconfigtobuffer!(out_string,
            molecule_name,
            css_sys,
            cs_singlets_compact;
            default_cs_delta = default_cs_delta)

    # write buffer to file.
    file_path = joinpath(save_dir, "default_cs_config.txt")
    DelimitedFiles.writedlm(file_path, out_string)

    return nothing
end

# zero_tol_sigdigits and unique_cs_tol are used to determine unique cs with tolerance.
# base_dir is where the GISSMO entry will be downloaded to.
# save_dir is where the cs and J-value JLD files will be saved.
function downloadGISSMOentry(X::GISSMOEntryType,
                            base_dir::String,
                            save_dir::String;
                            δ_lb::T = 0.1,
                            δ_ub::T = 0.1,
                            unique_cs_tol::T = 1e-6,
                            zero_tol_sigdigits = 6) where T <: Real

    # parse.
    molecule_name = X.molecule_name
    entry = X.entry
    eDATA_URL = X.eDATA_URL
    STAR_URL = X.STAR_URL

    # storage.
    record_name = "$(molecule_name)_$(entry)"
    data_dir = joinpath(base_dir, record_name)
    t = @task begin; isdir(data_dir) || mkdir(data_dir); end
    schedule(t); wait(t)

    # Download the STAR file.
    file_name = "$(record_name).str"
    dest = joinpath(data_dir, file_name)
    t = @task begin; isfile(dest) || download(STAR_URL, dest); end
    schedule(t); wait(t)

    # Download zip file.
    file_name = "$(record_name).zip"
    dest = joinpath(data_dir, file_name)
    t = @task begin; isfile(dest) || download(eDATA_URL, dest); end
    schedule(t); wait(t)

    #q = filter(x->endswith(x, ".zip"), readdir(data_dir))

    # unzip.
    InfoZIP.unzip(dest, data_dir)

    # remove zip file.
    #rm(dest)

    # get SDF filename.
    tmp = filter(isdir, readdir(data_dir; join=true))
    @assert length(tmp) == 1 # should be only one folder, which is the unzipped one.
    folder_path = tmp[1]

    tmp = filter(x->endswith(x, ".sdf"), readdir(folder_path))
    @assert length(tmp) == 1 # should be exactly 1 SDF file.
    SDF_name = tmp[1]
    SDF_path = joinpath(folder_path, SDF_name)

    file_strings = readlines(SDF_path)



    #### get absolute concentrations.
    solubility_gL_solute = solubilityinwater(molecule_name)
    molar_mass_solute = getmolarmass(base_dir, molecule_name)

    compound_names, concentration_values, concentration_units,
    compound_types, solute_concentration_string,
        solute_concentration_unit = parsecompositionfromNMReDATA(file_strings)

    #
    solute_mM, ref_mM, molar_mass_ref, solubility_gL_ref,
    ref_molecule_name = convertsampleconcentration(compound_names, concentration_values, concentration_units,
                            compound_types,
                            solute_concentration_string,
                            solute_concentration_unit,
                            molar_mass_solute, solubility_gL_solute)

    # assemble output. Units in mM<.
    println("solute_concentration = ", solute_mM)
    println("reference_concentration = ", ref_mM)

    #### find pH info.
    pH, temperature = parseconditionsfromNMReDATA(readlines(SDF_path))

    #### find the section on chemical shifts.
    cs_string = "> <NMREDATA_ASSIGNMENT>"

    H_IDs, H_css = parsecssfromNMReDATA(file_strings)


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
    println("Entry: ", entry)
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

    ### cs and J-values upper and lower bounds.
    perturblbfunc = xx->additiveperturbation(xx, δ_lb)
    perturbubfunc = xx->additiveperturbation(xx, δ_ub)

    #
    J_lb_sys = getJperturbed(J_vals_sys, perturblbfunc)
    J_ub_sys = getJperturbed(J_vals_sys, perturbubfunc)

    cs_lb_sys = getJperturbed(css_sys, perturblbfunc)
    cs_ub_sys = getJperturbed(css_sys, perturbubfunc)

    # store.
    #save_path = "$(save_dir)/$(record_name).jld"
    save_path = "$(save_dir)$(record_name).jld"

    JLD.save(save_path,
        "css_sys", css_sys,
        "J_IDs_sys", 𝐽_IDs_sys,
        "J_vals_sys", J_vals_sys,
        "cs_singlets", cs_singlets,

        "J_lb_sys", J_lb_sys,
        "J_ub_sys", J_ub_sys,
        "cs_lb_sys", cs_lb_sys,
        "cs_ub_sys", cs_ub_sys,
        "cs_LUT", cs_LUT,

        "ref_mM", ref_mM,
        "solute_mM", solute_mM,
        "solubility_gL_solute", solubility_gL_solute,
        "solubility_gL_ref", solubility_gL_ref,
        "molar_mass_ref", molar_mass_ref,
        "molar_mass_solute", molar_mass_solute,
        "pH", pH,
        "temperature", temperature,
        "ref_molecule_name", ref_molecule_name)

    return css_sys, cs_singlets_compact
end

function specialsingletentry(  entry,
                    molecule_name,
                    base_dir::String,
                    save_dir::String,
                    singlet_ppm;
                    unique_cs_tol = 1e-6,
                    δ_lb = 0.1,
                    δ_ub = 0.1)
    #
    H_IDs = Vector{Int}(undef, 1)
    H_IDs[1] = 1

    H_css = Vector{Float64}(undef, 1)
    H_css[1] = singlet_ppm # we're setting 4.8 ppm for D2O.

    J_IDs = Vector{Tuple{Int,Int}}(undef, 0)

    J_vals = Vector{Float64}(undef, 0)


    J_IDs_sys, J_vals_sys, H_IDs_sys = partitionJcouplings(J_IDs, J_vals)



    css_sys = getcssys(H_IDs_sys, H_IDs, H_css)

    𝐽_IDs_sys = getJsys(J_IDs_sys, J_vals_sys, H_IDs_sys)


    H_singlets, cs_singlets,
        cs_singlets_compact = createsingletsysems(H_IDs, H_IDs_sys,
                                    H_css;
                                    zero_tol = unique_cs_tol)

    cs_LUT, p_cs_sys = constructLUTcss(css_sys; zero_tol = unique_cs_tol)


    cs_len_sys = getcslengthfromLUT(cs_LUT)


    perturblbfunc = xx->additiveperturbation(xx, δ_lb)
    perturbubfunc = xx->additiveperturbation(xx, δ_ub)

    J_lb_sys = getJperturbed(J_vals_sys, perturblbfunc)
    J_ub_sys = getJperturbed(J_vals_sys, perturbubfunc)

    cs_lb_sys = getJperturbed(css_sys, perturblbfunc)
    cs_ub_sys = getJperturbed(css_sys, perturbubfunc)

    record_name = "$(molecule_name)_$(entry)"
    save_path = "$(save_dir)$(record_name).jld"

    JLD.save(save_path,
        "css_sys", css_sys,
        "J_IDs_sys", 𝐽_IDs_sys,
        "J_vals_sys", J_vals_sys,
        "cs_singlets", cs_singlets,

        "J_lb_sys", J_lb_sys,
        "J_ub_sys", J_ub_sys,
        "cs_lb_sys", cs_lb_sys,
        "cs_ub_sys", cs_ub_sys,
        "cs_LUT", cs_LUT,

        "ref_mM", NaN,
        "solute_mM", NaN,
        "solubility_gL_solute", NaN,
        "solubility_gL_ref", NaN,
        "molar_mass_ref", NaN,
        "molar_mass_solute", NaN,
        "pH", NaN,
        "temperature", NaN,
        "ref_molecule_name", NaN)

    return css_sys, cs_singlets_compact
end

function addentrydefaultconfigtobuffer!(out_string::Vector{String},
                                        molecule_name::String,
                                        css_sys::Vector{Vector{T}},
                                        cs_singlets_compact::Vector{T};
                                        default_cs_delta = 0.1) where T <: Real
    #
    # save heading to template.
    push!(out_string, molecule_name)

    # save the default cs configuration info.
    N_spin_sys = length(css_sys)
    for i = 1:N_spin_sys
        push!(out_string, "Group $(i)")

        for k = 1:length(css_sys[i])
            line_string = "cs $(css_sys[i][k]): $(default_cs_delta), $(k)"
            push!(out_string, line_string)
        end

        push!(out_string, "end group")
    end

    for i = 1:length(cs_singlets_compact)
        push!(out_string, "Group $(i+N_spin_sys)")
        line_string = "cs $(cs_singlets_compact[i]): $(default_cs_delta), 1"
        push!(out_string, line_string)
        push!(out_string, "end group")
    end

    push!(out_string, "end compound")
    push!(out_string, "")

    return nothing
end



function getGISSMOexperimentpaths(record_names, record_entries, dir_list)

    exp_paths = Vector{String}(undef, length(record_names))
    for k = 1:length(record_names)

        inds = findall(xx->occursin("$(record_names[k])_$(record_entries[k])", xx), dir_list)
        @assert length(inds) > 0

        current_folder_path = dir_list[inds[1]]
        next_folder = filter(isdir, readdir(current_folder_path; join = true))[1]

        while !isempty(next_folder) && !occursin("pdata", next_folder)

            current_folder_path = next_folder
            next_folder = filter(isdir, readdir(current_folder_path; join = true))[1]
        end

        exp_paths[k] = current_folder_path
    end

    return exp_paths
end
