#### new scheme.

#STR_file_root = "/home/roy/MEGAsync/data/NMR/GISSMO_str_files"

include("../src/GISSMOReader.jl")
import .GISSMOReader

import JSON3
import JSON
import InfoZIP

base_URL = "http://gissmo.nmrfam.wisc.edu/"
unique_cs_tol = 1e-6
zero_tol_sigdigits = 6

JLD_save_folder = "/home/roy/Documents/repo/NMRData/input/compounds"
edata_save_folder = "/home/roy/MEGAsync/data/NMR/GISSMO_edata"

JSON_folder = "/home/roy/Documents/repo/GISSMOReader.jl/GISSMO_info"
JSON_filename = "GISSMO_entries_metadata_Jan_2022.json"
JSON_load_path = joinpath(JSON_folder, JSON_filename)

db_dict = JSON.parse(open(JSON_load_path))
# TODO here. debug getcssys(). on key index 85.

###### single entry.

# IDs are the labelling used in GISSMO.

# entry_key = "http://gissmo.nmrfam.wisc.edu/entry/bmse000042/simulation_1" # L-Leucine.
# entry_key = "http://gissmo.nmrfam.wisc.edu/entry/bmse000860/simulation_1" # L-Valine.
# entry_key = "http://gissmo.nmrfam.wisc.edu/entry/bmse000041/simulation_1" # L-isoleucine.
# entry_key = "http://gissmo.nmrfam.wisc.edu/entry/bmse000297/simulation_1" # Ethanol
# entry_key = "http://gissmo.nmrfam.wisc.edu/entry/bmse000949/simulation_1" # D-Carnitine # problem.

#entry_key = "http://gissmo.nmrfam.wisc.edu/entry/bmse000113/simulation_1" # single resonance, multiple J-coupling.
entry_key = "http://gissmo.nmrfam.wisc.edu/entry/Maybridge_Ro3_Fragment_12_G10/simulation_1"

#entry_key = "http://gissmo.nmrfam.wisc.edu/entry/bmse000795/simulation_1" # DSS

dict = db_dict[entry_key]

molecule_name = dict["molecule_name"]
eData_URL = "$(base_URL)$(dict["GISSMO experiment URL"])"

simulation_label = split(entry_key, "/")[end]
entry_label = split(entry_key, "/")[end-1]

#JLD_save_folder = "/home/roy/Documents/repo/NMRData/input/compounds"
edata_save_folder = "/home/roy/MEGAsync/data/NMR/GISSMO_edata"
coupling_info_folder = "/home/roy/Documents/repo/GISSMOReader.jl/coupling_info"

# GISSMOReader.downloadGISSMOedata(molecule_name,
#     simulation_label,
#     entry_label,
#     JLD_save_folder,
#     edata_save_folder,
#     eData_URL;
#     unique_cs_tol = 1e-6,
#     zero_tol_sigdigits = 6)

#mkpath(JLD_save_folder)
mkpath(coupling_info_folder)

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
#rm(dest)

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

H_IDs, H_css_GISSMO = GISSMOReader.parsecssfromNMReDATA(file_strings)


#J_string = "> <NMREDATA_J>"
J_IDs, J_vals = GISSMOReader.parseJsfromNMReDATA(file_strings)


#######

# do this for ethanol first.

load_path = "/home/roy/MEGAsync/inputs/NMR/debug/ethanol_eq.json"
a1, a2, a3, a4 = GISSMOReader.loadcouplinginfojson(load_path)

@assert 1==4

import NMRSpectraSimulator

Id = NMRSpectraSimulator.getsingleId()
Ix = NMRSpectraSimulator.getsingleIx()
Iy_no_im = NMRSpectraSimulator.getsingleIynoim()
Iz = NMRSpectraSimulator.getsingleIz()


# pair-wise terms.

push!(J_IDs, (24,23))
push!(J_vals, -14.0)

N = length(H_IDs)

search_list = [22; 23; 24; 25]

inds = findall(xx->(any(xx[1] .== search_list) || any(xx[2] .== search_list)), J_IDs)
J_IDs2 = J_IDs[inds]
J_vals2 = J_vals[inds]
#fill!(J_vals2, 1.23)

J_IDs2_flat = unique([collect( J_IDs2[i][1] for i = 1:length(J_IDs2) ); collect( J_IDs2[i][2] for i = 1:length(J_IDs2) )])
H_IDs2 = filter( xx->any(xx .== J_IDs2_flat), H_IDs)

dict_H_IDs2_to_g = Dict(H_IDs2 .=> collect(1:length(H_IDs2)))
J_inds2 = collect( (dict_H_IDs2_to_g[J_IDs2[i][1]], dict_H_IDs2_to_g[J_IDs2[i][2]]) for i = 1:length(J_IDs2) )

# pair-wise terms.
H1s = collect( NMRSpectraSimulator.getmultiIjIk(Id, Ix, Iy_no_im, Iz, J_inds2[i][1], J_inds2[i][2], length(H_IDs2), 2*π*J_vals2[i]) for i = 1:length(J_inds2) )
H1 = sum(H1s)

using LinearAlgebra
s1,Q1 = eigen(H1)

#####
push!(J_IDs, (24,23))
inds = findall(xx->!isapprox(-14,xx), J_vals)
J_IDs = J_IDs[inds]
J_vals = J_vals[inds]

N = length(H_IDs)

search_list = [22; 23; 24; 25]

inds = findall(xx->(any(xx[1] .== search_list) || any(xx[2] .== search_list)), J_IDs)
J_IDs2 = J_IDs[inds]
J_vals2 = J_vals[inds]
#fill!(J_vals2, 1.23)

J_IDs2_flat = unique([collect( J_IDs2[i][1] for i = 1:length(J_IDs2) ); collect( J_IDs2[i][2] for i = 1:length(J_IDs2) )])
H_IDs2 = filter( xx->any(xx .== J_IDs2_flat), H_IDs)

dict_H_IDs2_to_g = Dict(H_IDs2 .=> collect(1:length(H_IDs2)))
J_inds2 = collect( (dict_H_IDs2_to_g[J_IDs2[i][1]], dict_H_IDs2_to_g[J_IDs2[i][2]]) for i = 1:length(J_IDs2) )

# pair-wise terms.
H2s = collect( NMRSpectraSimulator.getmultiIjIk(Id, Ix, Iy_no_im, Iz, J_inds2[i][1], J_inds2[i][2], length(H_IDs2), 2*π*J_vals2[i]) for i = 1:length(J_inds2) )
H2 = sum(H2s)

using LinearAlgebra
s2,Q2 = eigen(H2)

# # single terms.
# ωIzs = collect( ω0[j] .* Iz for j = 1:N )
# H2 = Kronecker.kroneckersum(ωIzs...)


@assert 1==2



H1 = zeros(T, 2^N, 2^N)
for j = 1:N
    for k = j+1:N
        #H += (2*π*J[j,k]) .* getmultiIjIk(Id, Ix, Iy_no_im, Iz, j, k, N)
        H1 += getmultiIjIk(Id, Ix, Iy_no_im, Iz, j, k, N, 2*π*J[j,k])

    end
end

# single terms.
ωIzs = collect( ω0[j] .* Iz for j = 1:N )
H2 = Kronecker.kroneckersum(ωIzs...)

@assert 1==43


# chemical equivalence: see if chemical shifts are the same.
H_css = round.(H_css_GISSMO, sigdigits = zero_tol_sigdigits)

A = collect(zip(H_IDs, H_css))
dict_ID_to_cs = Dict(H_IDs .=> H_css)


itr = pairs(H_css)
ua = Dict{valtype(itr),keytype(itr)}()
for (idx, val) in itr
    ua[val] = idx
end

A2 = unique(xx -> round(xx[2], sigdigits = zero_tol_sigdigits), A)

function getuniquecsIDs(A2, dict_ID_to_cs, H_IDs)

    H_IDs_g = collect(1:length(H_IDs))
    dict_H_IDs_to_g = Dict(H_IDs .=> H_IDs_g)

    A2_IDs = Vector{Vector{Int}}(undef, length(A2))
    A2_IDs_g = Vector{Vector{Int}}(undef, length(A2))
    for i = 1:length(A2)

        keep_flags = collect( isapprox( dict_ID_to_cs[H_IDs[k]], A2[i][2]) for k = 1:length(H_IDs) )
        #inds_A2[i] = collect(1:length(H_IDs))[keep_flags]
        A2_IDs[i] = H_IDs[keep_flags]
        A2_IDs_g[i] = H_IDs_g[keep_flags]
    end

    return A2_IDs, A2_IDs_g
end

A2_IDs, A2_IDs_g = getuniquecsIDs(A2, dict_ID_to_cs, H_IDs)


# re-index.
import Graphs

function constructspinsystemsgraphforcompound(H_IDs, J_IDs)
    g = Graphs.SimpleGraph(length(H_IDs))

    # re-number the spin atom labels to start at 1.
    J_IDs_reind = Vector{Tuple{Int,Int}}(undef, 0)
    for k = 1:length(J_IDs)

        i = findfirst(xx->xx==J_IDs[k][1], H_IDs)
        j = findfirst(xx->xx==J_IDs[k][2], H_IDs)

        push!(J_IDs_reind, (i,j))
        Graphs.add_edge!(g, i, j)
    end

    return g, J_IDs_reind
end

dict_H_IDs_to_g = Dict(H_IDs .=> collect(1:length(H_IDs)))
dict_g_to_H_IDs = Dict( collect(1:length(H_IDs)) .=> H_IDs)
g, J_IDs_g = constructspinsystemsgraphforcompound(H_IDs, J_IDs)
Q_g = Graphs.connected_components(g)

function mapQtoHID(Q_g::Vector{Vector{Int}}, H_IDs)

    dict_g_to_H_IDs = Dict( collect(1:length(H_IDs)) .=> H_IDs)

    Q = Vector{Vector{Int}}(undef, length(Q_g))
    for i = 1:length(Q_g)

        Q[i] = Vector{Int}(undef, length(Q_g[i]))
        for k = 1:length(Q_g[i])
            Q[i][k] = dict_g_to_H_IDs[Q_g[i][k]]
        end
    end

    return Q
end
Q = mapQtoHID(Q_g, H_IDs)

# given a subgraph, see if the vertices (spins) are magnetically equivalent.

function getJfromdict(i::Int, j::Int, dict::Dict{Tuple{Int,Int},T} ) where T
    #
    if haskey(dict, (i,j))
        return dict[(i,j)]
    end

    if haskey(dict, (j,i))
        return dict[(j,i)]
    end

    return zero(T)
end

# magnetic equivalence. For every pair of chemically equivalent spins, see if their J-coupling value to any other another spin is identifical.
h = Graphs.SimpleGraph(length(H_IDs)) # edge means the pair of nodes are magnetically equivalent.

# For each pair of singletons, see if they are chemically equivalent.
# for each

dict_J_ID_to_val = Dict(J_IDs .=> J_vals)


i = 2
chem_eq_IDs = A2_IDs[i]

# Updates h by placing edges between nodes that have magnetic equivalence.
"""
Given a list of spin IDs `chem_eq_IDs` and a graph `h` where the nodes are IDs (re-labeled to start at 1), place edges in graph `h` between the IDs specified by `chem_eq_IDs` that have magnetic equivalence.

The list of spin IDs should have the same chemical shift.
"""
function processmagneticeq!(h,
    g,
    chem_eq_IDs::Vector{Int},
    dict_H_IDs_to_g,
    dict_g_to_H_IDs,
    dict_ID_to_cs,
    dict_J_ID_to_val)

    for k0 = 1:length(chem_eq_IDs)
        ID_0 = chem_eq_IDs[k0]
        node_0 = dict_H_IDs_to_g[ID_0]

        neighbour_nodes_0 = Graphs.neighbors(g, node_0)
        neighbour_IDs_0 = collect( dict_g_to_H_IDs[neighbour_nodes_0[l]] for l = 1:length(neighbour_nodes_0) )

        cs_0 = dict_ID_to_cs[ID_0]

        for k1 = k0+1:length(chem_eq_IDs)

            ID_1 = chem_eq_IDs[k1]
            node_1 = dict_H_IDs_to_g[ID_1]

            neighbour_nodes_1 = Graphs.neighbors(g, node_1)
            neighbour_IDs_1 = collect( dict_g_to_H_IDs[neighbour_nodes_1[l]] for l = 1:length(neighbour_nodes_1) )

            # println("neighbour_IDs_0 = ", neighbour_IDs_0)
            # println("neighbour_IDs_1 = ", neighbour_IDs_1)
            #
            # println("ID_0 = ", ID_0)
            # println("ID_1 = ", ID_1)

            # get the IDs that are neighbours to both ID_0 and ID_1.
            #J_test_IDs = intersect(neighbour_IDs_0, neighbour_IDs_1)
            J_test_IDs = union(neighbour_IDs_0, neighbour_IDs_1)

            # See if ID_1 and ID_0 are magnetically equivalent using J_test_IDs.
            is_magnetically_equivalent = true
            for j = 1:length(J_test_IDs)
                ID_2 = J_test_IDs[j]

                if (ID_2 != ID_0) && (ID_2 != ID_1)
                    J_20 = getJfromdict(ID_2, ID_0, dict_J_ID_to_val)
                    J_21 = getJfromdict(ID_2, ID_1, dict_J_ID_to_val)

                    # println("ID_2 = ", ID_2)
                    # println("J_20 = ", J_20)
                    # println("J_21 = ", J_21)

                    if !isapprox(J_20, J_21)
                        is_magnetically_equivalent = false
                    end
                end
            end
            #println()

            # store the equivalence in h.
            if is_magnetically_equivalent
                Graphs.add_edge!(h, node_0, node_1)
            end
        end
    end

    return nothing
end


# # apply.
# processmagneticeq!(h, g, A2_IDs[4], dict_H_IDs_to_g,
#     dict_g_to_H_IDs, dict_ID_to_cs, dict_J_ID_to_val)
#
# # inspect manually to see if magnetic equivalence is recorded correctly. Consult GISSMO website. for spin matrix to assist with this task.
# connected_g = Graphs.connected_components(h)
# connected_IDs = collect( collect( dict_g_to_H_IDs[connected_g[i][j]] for j = 1:length(connected_g[i])) for i = 1:length(connected_g) )
#
# display(connected_IDs)
#
# @assert 55==1

function processmagneticeq!(h,
    g,
    A2_IDs::Vector{Vector{Int}},
    dict_H_IDs_to_g,
    dict_g_to_H_IDs,
    dict_ID_to_cs,
    dict_J_ID_to_val)

    for i = 1:length(A2_IDs)

        processmagneticeq!(h, g, A2_IDs[i], dict_H_IDs_to_g,
            dict_g_to_H_IDs, dict_ID_to_cs, dict_J_ID_to_val)
    end

    return nothing
end

# before.
connected_g = Graphs.connected_components(h)
connected_IDs = collect( collect( dict_g_to_H_IDs[connected_g[i][j]] for j = 1:length(connected_g[i])) for i = 1:length(connected_g) )
display(connected_IDs)

# apply.
processmagneticeq!(h, g, A2_IDs, dict_H_IDs_to_g,
    dict_g_to_H_IDs, dict_ID_to_cs, dict_J_ID_to_val)

# after.
# inspect manually to see if magnetic equivalence is recorded correctly. Consult GISSMO website. for spin matrix to assist with this task.
connected_g = Graphs.connected_components(h)
connected_IDs = collect( collect( dict_g_to_H_IDs[connected_g[i][j]] for j = 1:length(connected_g[i])) for i = 1:length(connected_g) )

display(connected_IDs)

@assert 55==4

function removemagneticequivalence(H_IDs, H_css_in::Vector{T}, J_IDs, J_vals;
    zero_tol_sigdigits = 6) where T

    H_css = round.(H_css_in, sigdigits = zero_tol_sigdigits)

    A = collect(zip(H_IDs, H_css))
    #dict_ID_to_cs = Dict(H_IDs .=> H_css)
    mapping_H_IDs_to_ind_A = Dict(H_IDs .=> collect(1:length(H_IDs)))


    # These are the atoms to keep, together with their unique chemical shift values.
    A2 = unique(xx -> round(xx[2], sigdigits = zero_tol_sigdigits), A)

    J_IDs_out = Vector{Tuple{Int,Int}}(undef, 0)
    for k = 1:length(J_IDs)
        id1, id2 = J_IDs[k]

        ind1 = mapping_H_IDs_to_ind_A[id1]
        cs_val1, ind1_A2 = findmin(xx->abs(xx[2]-A[ind1][2]), A2)

        ind2 = mapping_H_IDs_to_ind_A[id2]
        cs_val2, ind2_A2 = findmin(xx->abs(xx[2]-A[ind2][2]), A2)

        if !isapprox(cs_val1, cs_val2)
            push!(J_IDs_out, A[ind1][1], A[ind2][1])
        end

    end

    return A2, J_IDs_out
end

A2, J_IDs_out = removeidenticalcs(H_IDs, H_css, J_IDs)

### attempt on auto-grouping.

J_IDs_sys, J_vals_sys, H_IDs_sys = GISSMOReader.partitionJcouplings(J_IDs, J_vals)

css_sys = GISSMOReader.getcssys(H_IDs_sys, H_IDs, H_css)

𝐽_IDs_sys = GISSMOReader.getJsys(J_IDs_sys, J_vals_sys, H_IDs_sys)

H_singlets, cs_singlets,
    cs_singlets_compact = GISSMOReader.createsingletsysems(H_IDs, H_IDs_sys,
                                H_css;
                                zero_tol = unique_cs_tol)
#
### remove spin groups that only have one unique chem shift, but have J-coupling between the atoms.
# Each of such groups are are really one singlet group.
GISSMOReader.removeisolatedcs!(css_sys,
                cs_singlets,
                cs_singlets_compact,
                H_singlets,
                H_IDs_sys,
                J_IDs_sys,
                J_vals_sys,
                𝐽_IDs_sys;
                zero_tol_sigdigits = zero_tol_sigdigits)
#



@assert 1==2
JLD.save(save_path,
    "css_sys", css_sys,
    "J_IDs_sys", 𝐽_IDs_sys,
    "J_vals_sys", J_vals_sys,
    "cs_singlets", cs_singlets,
    "cs_LUT", cs_LUT)
