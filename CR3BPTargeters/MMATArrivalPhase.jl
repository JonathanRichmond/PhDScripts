"""
Arrival phase targeter for CR3BP MMATs

Author: Jonathan LeFevre Richmond
C: 1/30/26
U: 2/7/26
"""

using MBD, FrameTransformations, LinearAlgebra

export MMATArrivalPhaseTargeter, MMATEnv
export centralDifference, correct, evaluateConstraints

struct MMATEnv
    EMDynamicsModel::MBD.CR3BPDynamicsModel
    SDynamicsModel::MBD.KDynamicsModel
    SEDynamicsModel::MBD.CR3BPDynamicsModel
    SMDynamicsModel::MBD.CR3BPDynamicsModel

    Earth::MBD.BodyData
    Mars::MBD.BodyData
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
    MarsSoI::Float64
    MoonSoI::Float64
    SMMomentumDiff::Float64

    days::Vector{Float64}
    epochs::Vector{String}
    frame::FrameTransformations.FrameSystem
    initialEpoch::String
    initialEpochTime::Float64

    oe_E::Vector{Float64}
    oe_M::Vector{Float64}
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
    centralDifference(targeter, env, X, oe_bridge_peri, q_arr_SoI)

Return numerical Jacobian

# Arguments
- `targeter::MMATArrivalPhaseTargeter`: MMAT arrival phase targeter object
- `env::MMATEnv`: MMAT environment object
- `X::Vector{Float64}`: Free variable vector
- `oe_bridge_peri::Float64`: Bridge arc orbital elements at periapsis
- `q_arr_SoI::Vector{Float64}`: Arrival SoI state in CR3BP rotating frame [ndim]
"""
function centralDifference(targeter::MMATArrivalPhaseTargeter, env::MMATEnv, X::Vector{Float64}, oe_bridge_peri::Vector{Float64}, q_arr_SoI::Vector{Float64})
    DF::Matrix{Float64} = Matrix{Float64}(undef, 3, 3)
    for j::Int64 in 1:length(X)
        X_pos::Vector{Float64} = copy(X)
        X_neg::Vector{Float64} = copy(X)
        X_pos[j] += targeter.step
        X_neg[j] -= targeter.step
        F_pos::Vector{Float64} = evaluateConstraints(targeter, env, X_pos, oe_bridge_peri, q_arr_SoI)
        F_neg::Vector{Float64} = evaluateConstraints(targeter, env, X_neg, oe_bridge_peri, q_arr_SoI)
        DF[:,j] = (F_pos-F_neg)./(2*targeter.step)
    end

    return DF
end

"""
    correct(targeter, env, X, oe_bridge_peri, q_arr_SoI)

Return corrected multiple shooter problem

# Arguments
- `targeter::MMATArrivalPhaseTargeter`: MMAT arrival phase targeter object
- `env::MMATEnv`: MMAT environment object
- `X::Vector{Float64}`: Free variable vector ([theta_bridge_int, theta_M_SoI, theta_arr_int])
- `oe_bridge_peri::Vector{Float64}`: Bridge arc orbital elements at periapsis
- `q_arr_SoI::Vector{Float64}`: Arrival SoI state in CR3BP rotating frame [ndim]
"""
function correct(targeter::MMATArrivalPhaseTargeter, env::MMATEnv, X::Vector{Float64}, oe_bridge_peri::Vector{Float64}, q_arr_SoI::Vector{Float64})
    F::Vector{Float64} = -1.0 .* evaluateConstraints(targeter, env, X, oe_bridge_peri, q_arr_SoI)
    iter::Int64 = 1
    while (LinearAlgebra.norm(F) >= 1E-5) && (iter <= 50)
        DF::Matrix{Float64} = centralDifference(targeter, env, X, oe_bridge_peri, q_arr_SoI)
        solver = LinearAlgebra.qr(DF, LinearAlgebra.ColumnNorm())
        dX::Vector{Float64} = solver\F
        alpha::Float64 = 1.0#(LinearAlgebra.norm(F) > 5000) ? 0.1 : 1.0
        X += alpha.*dX
        F = -1.0 .* evaluateConstraints(targeter, env, X, oe_bridge_peri, q_arr_SoI)
        iter += 1
    end

    (LinearAlgebra.norm(F) >= 1E-5) ? throw(ErrorException("Corrections algorithm could not converge")) : (return X)
end

"""
    evaluateConstraints(targeter, env, X, oe_bridge_peri, q_arr_SoI)

Return evaluated constraint vector

# Arguments
- `targeter::MMATArrivalPhaseTargeter`: MMAT arrival phase targeter object
- `env::MMATEnv`: MMAT environment object
- `X::Vector{Float64}`: Free variable vector
- `oe_bridge_peri::Vector{Float64}`: Bridge arc orbital elements at periapsis
- `q_arr_SoI::Vector{Float64}`: Arrival SoI state in CR3BP rotating frame [ndim]
"""
function evaluateConstraints(targeter::MMATArrivalPhaseTargeter, env::MMATEnv, X::Vector{Float64}, oe_bridge_peri::Vector{Float64}, q_arr_SoI::Vector{Float64})
    Q_bridge_int::Vector{Float64} = getCartesianState(targeter.dynamicsModel, append!(oe_bridge_peri[1:5], X[1]))
    q_arr_SoI_SI::Vector{Float64} = rotatingToSunEclipJ2000(env.SMDynamicsModel, env.frame, env.initialEpochTime, [q_arr_SoI], [X[2]])[1]
    @views Q_arr_SoI_SI = similar(q_arr_SoI_SI)
    Q_arr_SoI_SI[1:3] .= q_arr_SoI_SI[1:3].*env.charValues.SM.lstar
    Q_arr_SoI_SI[4:6] .= q_arr_SoI_SI[4:6].*env.charValues.SM.lstar./env.charValues.SM.tstar
    oe_arr_SoI::Vector{Float64} = getOrbitalElements(targeter.dynamicsModel, Q_arr_SoI_SI)
    oe_arr_int::Vector{Float64} = append!(oe_arr_SoI[1:5], X[3])
    Q_arr_int::Vector{Float64} = getCartesianState(targeter.dynamicsModel, oe_arr_int)

    return Q_arr_int[1:3]-Q_bridge_int[1:3]
end
