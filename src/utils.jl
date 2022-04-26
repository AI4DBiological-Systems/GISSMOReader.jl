function isnumericallyclose(x::T, y::T, tol = eps(T)*2) where T
    if abs(x-y) < tol
        return true
    end

    return false
end

function getcslengthfromLUT(cs_LUT_sys::Vector{Vector{Vector{Int}}})

    N_subsystems = length(cs_LUT_sys)

    cs_len_sys = Vector{Int}(undef, N_subsystems)

    for m = 1:N_subsystems
        cs_LUT = cs_LUT_sys[m]

        cs_len_sys[m] = maximum(maximum(cs_LUT[i]) for i = 1:length(cs_LUT))
    end

    return cs_len_sys
end

function loadcouplinginfojson(save_path::String)::Tuple{Vector{Int}, Vector{Float64}, Vector{Tuple{Int,Int}}, Vector{Float64}}

    J = JSON.parsefile(save_path)

    H_IDs = Vector{Int}(undef, 0)
    H_css = Vector{Float64}(undef, 0)
    for dict in J["chemical shift"]
        id = convert(Int, dict["ID"])
        value = convert(Float64, dict["value"])
        push!(H_IDs, id)
        push!(H_css, value)
    end

    J_IDs = Vector{Tuple{Int,Int}}(undef, 0)
    J_vals = Vector{Float64}(undef, 0)
    for dict in J["J-coupling"]
        id1 = convert(Int, dict["ID1"])
        id2 = convert(Int, dict["ID2"])
        value = convert(Float64, dict["value"])
        push!(J_IDs, (id1, id2))
        push!(J_vals, value)
    end

    return H_IDs, H_css, J_IDs, J_vals
end
