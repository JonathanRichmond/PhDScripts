"""
Visibility utility functions

Author: Jonathan Richmond
C: 7/16/25
"""

using MBD

export checkExclusion, getLunarExclusionAngles, visibility

function checkExclusion(r_O::Matrix{Float64}, r_T::Matrix{Float64}, primaryPos::Vector{Float64}, sigma::Vector{Float64})
    nObs::Int64 = size(r_O, 2)
    nTargs::Int64 = size(r_T, 2)
    r_OP::Matrix{Float64} = primaryPos .- r_O
    norm_OP::Matrix{Float64} = sqrt.(sum(r_OP .^ 2, dims = 1))
    r_OT::Array{Float64, 3} = reshape(r_T, 3, 1, nTargs) .- reshape(r_O, 3, nObs, 1)
    dots::Array{Float64, 3} = sum(r_OT .* reshape(r_OP, 3, nObs, 1), dims = 1)
    norm_OT::Array{Float64, 3} = sqrt.(sum(r_OT .^ 2, dims = 1))
    kappa::Array{Float64, 3} = acos.(clamp.(dots ./ (norm_OT .* reshape(norm_OP, 1, nObs, 1)), -1.0, 1.0))

    return dropdims(abs.(kappa) .<= reshape(sigma, 1, nObs, 1), dims = 1)
end

function getLunarExclusionAngles(r_O::Matrix{Float64}, MoonPos::Vector{Float64}, SunPos::Vector{Float64}, R_M::Float64, lstar::Float64, sigmaM::Float64)
    r_OM::Matrix{Float64} = MoonPos .- r_O
    norm_OM::Matrix{Float64} = sqrt.(sum(r_OM .^ 2, dims = 1))
    r_OS::Matrix{Float64} = SunPos .- r_O
    norm_OS::Matrix{Float64} = sqrt.(sum(r_OS .^ 2, dims = 1))
    dots::Matrix{Float64} = sum(r_OM .* r_OS, dims = 1)
    Chi::Vector{Float64} = vec(acos.(clamp.(dots ./ (norm_OM .*norm_OS), -1.0, 1.0)))
    MoonOffset::Vector{Float64} = atan.(R_M ./ (lstar .* vec(norm_OM)))

    return sigmaM/pi .* Chi .+ MoonOffset
end

function visibility(r_O::Matrix{Float64}, r_T::Matrix{Float64}, EarthPos::Vector{Float64}, MoonPos::Vector{Float64}, SunPos::Vector{Float64}, EarthAngle::Float64, MoonAngle::Float64, SunAngle::Float64, r_M::Float64, lstar::Float64, objR::Float64, Cd::Float64, a4::Float64)
    nObs::Int64 = size(r_O, 2)
    nTargs::Int64 = size(r_T, 2)
    # sigmaM::Vector{Float64} = getLunarExclusionAngles(r_O, MoonPos, SunPos, r_M*lstar, lstar, MopnAngle)
    isEarthExcluded::Matrix{Bool} = checkExclusion(r_O, r_T, EarthPos, fill(EarthAngle, nObs))
    # isMoonExcluded::Matrix{Bool} = checkExclusion(r_O, r_T, MoonPos, sigmaM)
    isMoonExcluded::Matrix{Bool} = checkExclusion(r_O, r_T, MoonPos, fill(MoonAngle, nObs))
    isSunExcluded::Matrix{Bool} = checkExclusion(r_O, r_T, SunPos, fill(SunAngle, nObs))
    isExcluded::Matrix{Bool} = isEarthExcluded .| isMoonExcluded .| isSunExcluded
    r_T_rs::Array{Float64, 3} = reshape(r_T, 3, 1, nTargs)
    r_TO::Array{Float64, 3} = reshape(r_O, 3, nObs, 1) .- r_T_rs
    norm_TO::Array{Float64, 3} = sqrt.(sum(r_TO .^ 2, dims = 1))
    r_TS::Array{Float64, 3} = reshape(SunPos, 3, 1, 1) .- r_T_rs
    norm_TS::Array{Float64, 3} = sqrt.(sum(r_TS .^ 2, dims = 1))
    dots::Array{Float64, 3} = sum(r_TO .* r_TS, dims = 1)
    alpha::Array{Float64, 3} = acos.(clamp.(dots ./ (norm_TO .* norm_TS), -1.0, 1.0))
    brightness::Matrix{Float64} = dropdims(-26.832 .- 2.5 .* log10.(objR^2*Cd/(3*pi^2*a4^2) .* norm_TS .^ 2 .* (sin.(alpha) .+ (pi .- alpha) .* cos.(alpha)) ./ norm_TO .^ 2), dims = 1)
    brightness[isExcluded] .= 100.0

    return brightness
end
