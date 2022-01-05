
function parsemetadatafromNMReDATA(file_strings::Vector{String})

    header_string = "> <NMREDATA_SOLVENT>"

    H_IDs = Vector{Int}(undef, 0)
    H_css = Vector{Float64}(undef, 0)

    record_flag = false

    for i = 1:length(file_strings)

        s = file_strings[i]

        if record_flag

            #tokens = split("_My input.string", (' ','_',',','.'))
            tokens = split(s, ',')
            #println("tokens = ", tokens)

            if length(tokens) >= 2
                tmp1 = tryparse(Int, tokens[1])
                tmp2 = tryparse(Float64, tokens[2])

                if typeof(tmp1) == Int && typeof(tmp2) == Float64
                    push!(H_IDs, tmp1)
                    push!(H_css, tmp2)
                else
                    record_flag = false
                end
            else
                record_flag = false
            end

        end


        if s == header_string
            record_flag = true
        end


    end

    return H_IDs, H_css
end

function parsecssfromNMReDATA(file_strings::Vector{String})

    header_string = "> <NMREDATA_ASSIGNMENT>"

    H_IDs = Vector{Int}(undef, 0)
    H_css = Vector{Float64}(undef, 0)

    record_flag = false

    for i = 1:length(file_strings)

        s = file_strings[i]

        if record_flag

            #tokens = split("_My input.string", (' ','_',',','.'))
            tokens = split(s, ',')
            #println("tokens = ", tokens)

            if length(tokens) >= 2
                tmp1 = tryparse(Int, tokens[1])
                tmp2 = tryparse(Float64, tokens[2])

                if typeof(tmp1) == Int && typeof(tmp2) == Float64
                    push!(H_IDs, tmp1)
                    push!(H_css, tmp2)
                else
                    record_flag = false
                end
            else
                record_flag = false
            end

        end


        if s == header_string
            record_flag = true
        end


    end

    return H_IDs, H_css
end


function parseJsfromNMReDATA(file_strings::Vector{String})

    header_string = "> <NMREDATA_J>"

    J_IDs = Vector{Tuple{Int,Int}}(undef, 0)
    J_vals = Vector{Float64}(undef, 0)

    record_flag = false

    for i = 1:length(file_strings)

        s = file_strings[i]

        if record_flag

            tokens = split(s, (' ',','))
            # println("tokens = ", tokens)
            # tokens = SubString{String}["26", "", "27", "", "-12.000000", "\\"]

            if length(tokens) >= 5
                tmp1 = tryparse(Int, tokens[1])
                tmp2 = tryparse(Int, tokens[3])
                tmp3 = tryparse(Float64, tokens[5])

                # println("tmp1 = ", tmp1)
                # println("tmp2 = ", tmp2)
                # println("tmp3 = ", tmp3)

                if typeof(tmp1) == Int && typeof(tmp2) == Int &&
                        typeof(tmp3) == Float64

                    push!(J_IDs, (tmp1,tmp2))
                    push!(J_vals, tmp3)
                else
                    record_flag = false
                end
            else
                record_flag = false
            end

        end


        if s == header_string
            record_flag = true
        end


    end

    return J_IDs, J_vals
end


function partitionJcouplings(  J_IDs::Vector{Tuple{Int,Int}},
                                J_vals::Vector{T}) where T <: Real

    #
    @assert length(J_IDs) == length(J_vals)

    H_IDs_sys = Vector{Vector{Int}}(undef, 0)
    J_IDs_sys = Vector{Vector{Tuple{Int,Int}}}(undef, 0)
    J_vals_sys = Vector{Vector{T}}(undef, 0)

    for i = 1:length(J_IDs)
        a, b = J_IDs[i]

        # fetch fresh H_IDs_sys.
        H_IDs_sys = collect(finduniqueIDs(J_IDs_sys[m]) for m = 1:length(J_IDs_sys))

        #println("H_IDs_sys = ", H_IDs_sys)

        presence_flag = falses(length(H_IDs_sys))

        # search through existing cliques
        for m = 1:length(H_IDs_sys)

            # search within the current clique.
            for l = 1:length(H_IDs_sys[m])
                presence_flag[m] = presence_flag[m] | (H_IDs_sys[m][l] == a) | (H_IDs_sys[m][l] == b)
            end
        end

        # update
        if any(presence_flag)

            js = findall(presence_flag)

            # add to the first existing subsystem.
            push!(J_IDs_sys[js[1]], J_IDs[i])
            push!(J_vals_sys[js[1]], J_vals[i])

            if length(js) > 1
                # merge all existing subsystems that share the hydrogens in this ID.
                consolidatearray!(J_IDs_sys, js)
                consolidatearray!(J_vals_sys, js)
            end

        else
            # start new subsystem to the end of J_IDs_sys.
            push!(J_IDs_sys, Vector{Tuple{Int,Int}}(undef,1))
            J_IDs_sys[end][1] = J_IDs[i]

            push!(J_vals_sys, Vector{T}(undef,1))
            J_vals_sys[end][1] = J_vals[i]
        end

        ## force update on H_IDs_sys. Need this for for 2-spin systems.
        # fetch fresh H_IDs_sys.
        H_IDs_sys = collect(finduniqueIDs(J_IDs_sys[m]) for m = 1:length(J_IDs_sys))

    end

    return J_IDs_sys, J_vals_sys, H_IDs_sys
end


function finduniqueIDs(J_IDs::Vector{Tuple{Int,Int}})

    out = Vector{Int}(undef, 0)
    for i = 1:length(J_IDs)

        push!(out, J_IDs[i][1])
        push!(out, J_IDs[i][2])
    end

    return unique(out)
end

function consolidatearray!(X::Vector{Vector{T}}, inds::Vector{Int})::Nothing where T
    @assert length(inds) > 1
    j = inds[1]
    @assert 1 <= j <= length(X)

    for i = 2:length(inds)
        @assert 1 <= inds[i] <= length(X)

        push!(X[j], X[i]...)
    end

    deleteat!(X, inds[2:end])

    return nothing
end

function getcssys(  H_IDs_sys::Vector{Vector{Int}},
                    H_IDs::Vector{Int},
                    H_css::Vector{T}) where T <: Real

    @assert length(H_IDs) == length(H_css)

    N_subsystems = length(H_IDs_sys)

    # for sanity-checking that H_IDs_sys have unique finest-grain elements.
    processed_flags = falses(length(H_IDs))

    # output.
    css_sys = Vector{Vector{T}}(undef, N_subsystems)

    for m = 1:length(H_IDs_sys)

        L = length(H_IDs_sys[m])

        css_sys[m] = Vector{T}(undef, L)

        for l = 1:L

            t = H_IDs_sys[m][l]

            k = findfirst(xx->xx==t, H_IDs)
            @assert typeof(k) != Nothing
            @assert processed_flags[k] == false

            # store.
            css_sys[m][l] = H_css[k]

            # book keep.
            processed_flags[k] = true
        end
    end

    return css_sys
end

function getJsys(  J_IDs_sys::Vector{Vector{Tuple{Int,Int}}},
                    J_vals_sys::Vector{Vector{T}},
                    H_IDs_sys::Vector{Vector{Int}}) where T <: Real

    N_subsystems = length(H_IDs_sys)
    @assert length(J_IDs_sys) == length(J_vals_sys) == N_subsystems

    # for sanity-checking that H_IDs_sys have unique finest-grain elements.
    #processed_flags = falses(length(H_IDs))

    # output.
    𝐽_IDs_sys = Vector{Vector{Tuple{Int,Int}}}(undef, N_subsystems)

    for m = 1:length(J_IDs_sys)

        # parse & set up.
        Hs = H_IDs_sys[m]
        Js = J_IDs_sys[m]

        N = length(Js)

        # allocate.
        𝐽_IDs_sys[m] = Vector{Tuple{Int,Int}}(undef, N)

        for i = 1:N

            a, b = Js[i]

            k = findfirst(xx->xx==a, Hs)
            @assert typeof(k) != Nothing

            z = findfirst(xx->xx==b, Hs)
            @assert typeof(z) != Nothing


            # store.
            𝐽_IDs_sys[m][i] = (k,z)
        end
    end

    return 𝐽_IDs_sys
end

function createsingletsysems(   H_IDs,
                                H_IDs_multiplets_sys,
                                H_css::Vector{T};
                                zero_tol = 1e-6) where T
    #
    M = length(H_IDs_multiplets_sys)

    H_set = Set(H_IDs)

    for m = 1:M
        A = Set(H_IDs_multiplets_sys[m])

        H_set = setdiff(H_set, A)
    end

    # remainder goes into a new set each, depending on their chemical shift.
    R = collect(H_set)
    R_cs = Vector{T}(undef, length(R))

    for i = 1:length(R)

        k = findfirst(xx->xx==R[i], H_IDs)
        @assert typeof(k) != Nothing

        R_cs[i] = H_css[k]
    end


    # consolidate groups with similar cs (values within zero_tol).
    cs_singlets = Vector{Vector{T}}(undef, 0)
    H_singlets = Vector{Vector{Int}}(undef, 0)

    while length(R_cs) > 0
        t = R_cs[1]

        inds = findall(xx->isnumericallyclose(xx, t, zero_tol), R_cs)

        push!(cs_singlets, R_cs[inds])
        push!(H_singlets, R[inds])

        deleteat!(R_cs, inds)
    end

    # compact representation.
    cs_singlets_compact = collect( cs_singlets[m][1] for m = 1:length(cs_singlets) )

    return H_singlets, cs_singlets, cs_singlets_compact
end

function constructLUTcss(   css_sys::Vector{Vector{T}};
                            zero_tol = 1e-6) where T <: Real
    #
    N_subsystems = length(css_sys)

    cs_LUT = Vector{Vector{Vector{Int}}}(undef, N_subsystems)
    p_cs = Vector{Vector{T}}(undef, N_subsystems)

    for m = 1:N_subsystems

        # set up.
        cs = css_sys[m]
        N = length(css_sys[m])
        available_flags = trues(N)

        # allocate.
        cs_LUT[m] = Vector{Vector{Int}}(undef,0)
        p_cs[m] = Vector{T}(undef,0)

        for i = 1:N
            if available_flags[i]

                # search.
                t = cs[i]

                inds = findall(xx->isnumericallyclose(xx, t, zero_tol), cs)

                # store.
                push!(cs_LUT[m], inds)
                push!(p_cs[m], t)

                # book keep.
                for l in inds
                    available_flags[l] = false
                end
            end
        end
    end

    return cs_LUT, p_cs
end

### remove spin groups that have only one chem shift, but non-zero J-coupling.
# Those are really singlets.
function removeisolatedcs!(css_sys::Vector{Vector{T}},
                            cs_singlets::Vector{Vector{T}},
                            cs_singlets_compact::Vector{T},
                            H_singlets,
                            H_IDs_sys::Vector{Vector{Int}},
                            J_IDs_sys::Vector{Vector{Tuple{Int,Int}}},
                            J_vals_sys::Vector{Vector{T}},
                            𝐽_IDs_sys::Vector{Vector{Tuple{Int,Int}}};
                            zero_tol_sigdigits::Int = 6,
                            zero_tol = 1e-10) where T

    # add approximate spin groups as singlet groups.
    keep_flags = trues(length(css_sys))
    for i = 1:length(css_sys)
        #
        tmp = unique(x -> round(x, sigdigits = zero_tol_sigdigits), css_sys[i])
        if length(tmp) == 1
            keep_flags[i] = false

            new_cs_group = tmp .* ones(T, length(css_sys[i]))
            push!(cs_singlets, new_cs_group)

            push!(cs_singlets_compact, tmp[1])

            push!(H_singlets, H_IDs_sys[i])
        end
    end

    # remove those sping groups from spin system-related quantities.
    delete_flags = collect( !keep_flags[i] for i = 1:length(keep_flags) )

    deleteat!(H_IDs_sys, delete_flags)
    deleteat!(css_sys, delete_flags)

    deleteat!(J_IDs_sys, delete_flags)
    deleteat!(J_vals_sys, delete_flags)
    deleteat!(𝐽_IDs_sys, delete_flags)

    # merge singlet groups if neccessary.
    unique_css = unique(x -> round(x, sigdigits = zero_tol_sigdigits), cs_singlets_compact)
    if length(unique_css) != length(cs_singlets_compact)
        # we can merge at least two singlet groups.

        # take care of the cs_singlets.
        for i = 1:length(unique_css)

            x = unique_css[i]

            merge_flags = collect( isapprox(cs_singlets_compact[i], x; atol = zero_tol) for i = 1:length(cs_singlets_compact))
            merge_inds = collect(1:length(cs_singlets_compact))[merge_flags]
            @assert length(merge_inds) > 0

            # add to the keep_ind.
            keep_ind = merge_inds[1]
            for j = 2:length(merge_inds)

                k = merge_inds[j]
                push!(cs_singlets[keep_ind], cs_singlets[k]...)
                push!(H_singlets[keep_ind], H_singlets[k]...)
            end

            if length(merge_inds) > 1
                # delete the merge_inds that isn't the keep_ind (which is the first entry).
                delete_flags = [false; merge_flags[2:end]]
                deleteat!(cs_singlets, delete_flags)
                deleteat!(H_singlets, delete_flags)
            end
        end

        # update compact version.
        resize!(cs_singlets_compact, length(cs_singlets))
        for i = 1:length(cs_singlets)
            cs_singlets_compact[i] = cs_singlets[i][1]
        end
    end

    return nothing
end

function mergefirst(x::Vector{T}, a::T; tol::T)::Int where T
    for i = 1:length(x)
        if isapprox(x[i], a, atol = tol)
            return i
        end
    end

    return 0
end

#### perturb.

function additiveperturbation(x::T, δ::T)::T where T
    return x + δ
end

# function multiplicativeperturbation(x::T, δ::T)::T where T
#     return x*δ
# end

function getJperturbed(  J_vals_sys::Vector{Vector{T}},
                        perturbfunc::Function)::Vector{Vector{T}} where T <: Real

    #
    N_subsystems = length(J_vals_sys)

    J_p_sys = Vector{Vector{T}}(undef, N_subsystems)

    for m = 1:length(J_vals_sys)
        N = length(J_vals_sys[m])

        J_p_sys[m] = Vector{T}(undef, N)

        for i = 1:N
            J_p_sys[m][i] = perturbfunc(J_vals_sys[m][i])
        end
    end

    return J_p_sys
end
