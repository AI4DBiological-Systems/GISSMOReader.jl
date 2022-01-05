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