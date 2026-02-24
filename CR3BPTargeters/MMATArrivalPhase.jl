"""
Arrival phase targeter for CR3BP MMATs

Author: Jonathan LeFevre Richmond
C: 1/30/26
U: 2/23/26
"""

using MBD, LinearAlgebra
using ..EphemeridesLoader: FrameTransformations

export MMATArrivalPhaseTargeter, MMATEnv
export centralDifference_ext, centralDifference_int, correct_ext, correct_int
export evaluateConstraints_ext, evaluateConstraints_int

struct MMATEnv
    EMDynamicsModel::MBD.CR3BPDynamicsModel
    SDynamicsModel::MBD.KDynamicsModel
    SEDynamicsModel::MBD.CR3BPDynamicsModel
    SArrDynamicsModel::MBD.CR3BPDynamicsModel

    Earth::MBD.BodyData
    arrBody::MBD.BodyData
    Moon::MBD.BodyData
    Sun::MBD.BodyData

    charValues::NamedTuple

    momentumPropagator::MBD.Propagator
    orbitDepartureEvent::DifferentialEquations.ContinuousCallback
    propagator::MBD.Propagator
    P1DistanceEvent::DifferentialEquations.ContinuousCallback
    P2DistanceEvent::DifferentialEquations.ContinuousCallback

    EarthSoI::Float64
    EMMomentumDiff::Float64
    arrSoI::Float64
    MoonSoI::Float64
    SArrMomentumDiff::Float64

    days::Vector{Float64}
    epochs::Vector{String}
    frame::FrameTransformations.FrameSystem
    initialEpoch::String
    initialEpochTime::Float64

    oe_E::Vector{Float64}
    oe_arrBody::Vector{Float64}
end

"""
    MMATArrivalPhaseTargeter(dynamicsModel, kappa)

MMAT arrival phase targeter object

# Arguments
- `dynamicsModel::KDynamicsModel`: Keplerian dynamics model object
- `kappa::Float64`: Central differencing step
"""
struct MMATArrivalPhaseTargeter
    dynamicsModel::MBD.KDynamicsModel                                   # Keplerian dynamics model object
    step::Float64                                                       # Central differencing step

    function MMATArrivalPhaseTargeter(dynamicsModel::MBD.KDynamicsModel, kappa::Float64)
        this = new(dynamicsModel, kappa)

        return this
    end
end

"""
    centralDifference_ext(targeter, env, X, oe_bridge_peri, q_arr_SoI)

Return numerical Jacobian for exterior transfer

# Arguments
- `targeter::MMATArrivalPhaseTargeter`: MMAT arrival phase targeter object
- `env::MMATEnv`: MMAT environment object
- `X::Vector{Float64}`: Free variable vector
- `oe_bridge_peri::Float64`: Bridge arc orbital elements at periapsis
- `q_arr_SoI::Vector{Float64}`: Arrival SoI state in CR3BP rotating frame [ndim]
"""
function centralDifference_ext(targeter::MMATArrivalPhaseTargeter, env::MMATEnv, X::Vector{Float64}, oe_bridge_peri::Vector{Float64}, q_arr_SoI::Vector{Float64})
    DF::Matrix{Float64} = Matrix{Float64}(undef, 3, 3)
    for j::Int64 in 1:length(X)
        X_pos::Vector{Float64} = copy(X)
        X_neg::Vector{Float64} = copy(X)
        X_pos[j] += targeter.step
        X_neg[j] -= targeter.step
        F_pos::Vector{Float64} = evaluateConstraints_ext(targeter, env, X_pos, oe_bridge_peri, q_arr_SoI)
        F_neg::Vector{Float64} = evaluateConstraints_ext(targeter, env, X_neg, oe_bridge_peri, q_arr_SoI)
        DF[:,j] = (F_pos-F_neg)./(2*targeter.step)
    end

    return DF
end

"""
    centralDifference_int(targeter, env, X, oe_dep_SoI, q_arr_SoI, a_bridge, e_bridge)

Return numerical Jacobian for interior transfer

# Arguments
- `targeter::MMATArrivalPhaseTargeter`: MMAT arrival phase targeter object
- `env::MMATEnv`: MMAT environment object
- `X::Vector{Float64}`: Free variable vector
- `oe_dep_SoI::Vector{Float64}`: Departure arc orbital elements at SoI
- `q_arr_SoI::Vector{Float64}`: Arrival SoI state in CR3BP rotating frame [ndim]
- `a_bridge::Float64`: Bridge arc semimajor axis [km]
- `e_bridge::Float64`: Bridge arc eccentricity
"""
function centralDifference_int(targeter::MMATArrivalPhaseTargeter, env::MMATEnv, X::Vector{Float64}, oe_dep_SoI::Vector{Float64}, q_arr_SoI::Vector{Float64}, a_bridge::Float64, e_bridge::Float64)
    DF::Matrix{Float64} = Matrix{Float64}(undef, 3, 3)
    for j::Int64 in 1:length(X)
        X_pos::Vector{Float64} = copy(X)
        X_neg::Vector{Float64} = copy(X)
        X_pos[j] += targeter.step
        X_neg[j] -= targeter.step
        F_pos::Vector{Float64} = evaluateConstraints_int(targeter, env, X_pos, oe_dep_SoI, q_arr_SoI, a_bridge, e_bridge)
        F_neg::Vector{Float64} = evaluateConstraints_int(targeter, env, X_neg, oe_dep_SoI, q_arr_SoI, a_bridge, e_bridge)
        DF[:,j] = (F_pos-F_neg)./(2*targeter.step)
    end

    return DF
end

"""
    correct_ext(targeter, env, X, oe_bridge_peri, q_arr_SoI)

Return corrected multiple shooter problem for exterior transfer

# Arguments
- `targeter::MMATArrivalPhaseTargeter`: MMAT arrival phase targeter object
- `env::MMATEnv`: MMAT environment object
- `X::Vector{Float64}`: Free variable vector ([theta_bridge_int, theta_M_SoI, theta_arr_int])
- `oe_bridge_peri::Vector{Float64}`: Bridge arc orbital elements at periapsis
- `q_arr_SoI::Vector{Float64}`: Arrival SoI state in CR3BP rotating frame [ndim]
"""
function correct_ext(targeter::MMATArrivalPhaseTargeter, env::MMATEnv, X::Vector{Float64}, oe_bridge_peri::Vector{Float64}, q_arr_SoI::Vector{Float64})
    F::Vector{Float64} = -evaluateConstraints_ext(targeter, env, X, oe_bridge_peri, q_arr_SoI)
    iter::Int64 = 1
    while (LinearAlgebra.norm(F) >= 1E-5) && (iter <= 50)
        DF::Matrix{Float64} = centralDifference_ext(targeter, env, X, oe_bridge_peri, q_arr_SoI)
        solver = LinearAlgebra.qr(DF, LinearAlgebra.ColumnNorm())
        dX::Vector{Float64} = solver\F
        X += dX
        F = -evaluateConstraints_ext(targeter, env, X, oe_bridge_peri, q_arr_SoI)
        iter += 1
    end

    (LinearAlgebra.norm(F) >= 1E-5) ? throw(ErrorException("Corrections algorithm could not converge")) : (return X)
end

"""
    correct_int(targeter, env, X, oe_dep_SoI, q_arr_SoI, a_bridge, e_bridge)

Return corrected multiple shooter problem for interior transfer

# Arguments
- `targeter::MMATArrivalPhaseTargeter`: MMAT arrival phase targeter object
- `env::MMATEnv`: MMAT environment object
- `X::Vector{Float64}`: Free variable vector ([theta_bridge_int, theta_M_SoI, theta_arr_int])
- `oe_dep_SoI::Vector{Float64}`: Departure arc orbital elements at SoI
- `q_arr_SoI::Vector{Float64}`: Arrival SoI state in CR3BP rotating frame [ndim]
- `a_bridge::Float64`: Bridge arc semimajor axis [km]
- `e_bridge::Float64`: Bridge arc eccentricity
"""
function correct_int(targeter::MMATArrivalPhaseTargeter, env::MMATEnv, X::Vector{Float64}, oe_dep_SoI::Vector{Float64}, q_arr_SoI::Vector{Float64}, a_bridge::Float64, e_bridge::Float64)
    F::Vector{Float64} = -evaluateConstraints_int(targeter, env, X, oe_dep_SoI, q_arr_SoI, a_bridge, e_bridge)
    iter::Int64 = 1
    while (LinearAlgebra.norm(F) >= 1E-5) && (iter <= 50)
        DF::Matrix{Float64} = centralDifference_int(targeter, env, X, oe_dep_SoI, q_arr_SoI, a_bridge, e_bridge)
        solver = LinearAlgebra.qr(DF, LinearAlgebra.ColumnNorm())
        dX::Vector{Float64} = solver\F
        X += dX
        F = -evaluateConstraints_int(targeter, env, X, oe_dep_SoI, q_arr_SoI, a_bridge, e_bridge)
        iter += 1
    end

    (LinearAlgebra.norm(F) >= 1E-5) ? throw(ErrorException("Corrections algorithm could not converge")) : (return X)
end

"""
    evaluateConstraints_ext(targeter, env, X, oe_bridge_peri, q_arr_SoI)

Return evaluated constraint vector for exterior transfer

# Arguments
- `targeter::MMATArrivalPhaseTargeter`: MMAT arrival phase targeter object
- `env::MMATEnv`: MMAT environment object
- `X::Vector{Float64}`: Free variable vector
- `oe_bridge_peri::Vector{Float64}`: Bridge arc orbital elements at periapsis
- `q_arr_SoI::Vector{Float64}`: Arrival SoI state in CR3BP rotating frame [ndim]
"""
function evaluateConstraints_ext(targeter::MMATArrivalPhaseTargeter, env::MMATEnv, X::Vector{Float64}, oe_bridge_peri::Vector{Float64}, q_arr_SoI::Vector{Float64})
    Q_bridge_int::Vector{Float64} = getCartesianState(targeter.dynamicsModel, append!(oe_bridge_peri[1:5], X[1]))
    q_arr_SoI_SI::Vector{Float64} = rotatingToSunEclipJ2000(env.SArrDynamicsModel, env.frame, env.initialEpochTime, [q_arr_SoI], [X[2]])[1]
    @views Q_arr_SoI_SI = similar(q_arr_SoI_SI)
    Q_arr_SoI_SI[1:3] .= q_arr_SoI_SI[1:3].*env.charValues.SArr.lstar
    Q_arr_SoI_SI[4:6] .= q_arr_SoI_SI[4:6].*env.charValues.SArr.lstar./env.charValues.SArr.tstar
    oe_arr_SoI::Vector{Float64} = getOrbitalElements(targeter.dynamicsModel, Q_arr_SoI_SI)
    Q_arr_int::Vector{Float64} = getCartesianState(targeter.dynamicsModel, append!(oe_arr_SoI[1:5], X[3]))

    return Q_arr_int[1:3]-Q_bridge_int[1:3]
end

"""
    evaluateConstraints_int(targeter, env, X, oe_dep_SoI, q_arr_SoI, a_bridge, e_bridge)

Return evaluated constraint vector for interior transfer

# Arguments
- `targeter::MMATArrivalPhaseTargeter`: MMAT arrival phase targeter object
- `env::MMATEnv`: MMAT environment object
- `X::Vector{Float64}`: Free variable vector
- `oe_dep_SoI::Vector{Float64}`: Departure arc orbital elements at SoI
- `q_arr_SoI::Vector{Float64}`: Arrival SoI state in CR3BP rotating frame [ndim]
- `a_bridge::Float64`: Bridge arc semimajor axis [km]
- `e_bridge::Float64`: Bridge arc eccentricity
"""
function evaluateConstraints_int(targeter::MMATArrivalPhaseTargeter, env::MMATEnv, X::Vector{Float64}, oe_dep_SoI::Vector{Float64}, q_arr_SoI::Vector{Float64}, a_bridge::Float64, e_bridge::Float64)
    Q_dep_int::Vector{Float64} = getCartesianState(targeter.dynamicsModel, append!(oe_dep_SoI[1:5], X[1]))
    q_arr_SoI_SI::Vector{Float64} = rotatingToSunEclipJ2000(env.SArrDynamicsModel, env.frame, env.initialEpochTime, [q_arr_SoI], [X[2]])[1]
    @views Q_arr_SoI_SI = similar(q_arr_SoI_SI)
    Q_arr_SoI_SI[1:3] .= q_arr_SoI_SI[1:3].*env.charValues.SArr.lstar
    Q_arr_SoI_SI[4:6] .= q_arr_SoI_SI[4:6].*env.charValues.SArr.lstar./env.charValues.SArr.tstar
    oe_arr_SoI::Vector{Float64} = getOrbitalElements(targeter.dynamicsModel, Q_arr_SoI_SI)
    oe_bridge_peri::Vector{Float64} = [a_bridge, e_bridge, oe_arr_SoI[3], oe_arr_SoI[4], oe_arr_SoI[5], 0.0]
    Q_bridge_int::Vector{Float64} = getCartesianState(targeter.dynamicsModel, append!(oe_bridge_peri[1:5], X[3]))

    return Q_bridge_int[1:3]-Q_dep_int[1:3]
end
