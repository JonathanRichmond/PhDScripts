"""
Period perpendicular crossing targeter for BCR4BP planar orbits

Author: Jonathan Richmond
C: 4/9/25
U: 6/5/25
"""

using MBD, DifferentialEquations, StaticArrays

export PlanarPerpP12Targeter
export correct, getMonodromy, getPeriod, propagateState#, getIndividualPeriodicOrbit, interpOrbit

"""
    PlanarPerpP12Targeter(dynamicsModel)

BCR4BP P1-P2 planar perpendicular crossing period targeter object

# Arguments
- `dynamicsModel::BCR4BP12DynamicsModel`: BCR4BP P1-P2 dynamics model object
"""
struct PlanarPerpP12Targeter <: MBD.AbstractTargeter
    dynamicsModel::MBD.BCR4BP12DynamicsModel                            # BCR4BP P1-P2 dynamics model object

    function PlanarPerpP12Targeter(dynamicsModel::MBD.BCR4BP12DynamicsModel)
        this = new(dynamicsModel)

        return this
    end
end

"""
    correct(targeter, q0, targetP, synRevs; tol, JTol)

Return corrected BCR4BP P1-P2 multiple shooter problem object

# Arguments
- `targeter::PlanarPerpP12Targeter`: BCR4BP P1-P2 planar perpendicular crossing period targeter object
- `q0::Vector{Float64}`: Initial state guess [ndim]
- `targetP::Float64`: Target period [ndim]
- `synRevs::Int64`: Number of synodic revolutions
- `tol::Float64`: Convergence tolerance (default = 1E-11)
- `JTol::Float64`: Jacobian accuracy tolerance (default = 2E-3)
"""
function correct(targeter::PlanarPerpP12Targeter, q0::Vector{Float64}, targetP::Float64, synRevs::Int64; tol::Float64 = 1E-11, JTol::Float64 = 2E-3)
    halfPeriod::Float64 = targetP/2
    qPCGuess::Vector{Float64} = propagateState(targeter, q0, [0, halfPeriod])
    originNode = MBD.BCR4BP12Node(0.0, q0, targeter.dynamicsModel)
    originNode.state.name = "Initial State"
    setFreeVariableMask!(originNode.state, [true, false, false, false, true, false, false])
    terminalNode = MBD.BCR4BP12Node(halfPeriod, qPCGuess, targeter.dynamicsModel)
    terminalNode.state.name = "Target State"
    segment = MBD.BCR4BP12Segment(halfPeriod, originNode, terminalNode)
    problem = MBD.BCR4BP12MultipleShooterProblem()
    addSegment!(problem, segment)
    addConstraint!(problem, MBD.BCR4BP12ContinuityConstraint(segment))
    addConstraint!(problem, MBD.BCR4BP12StateConstraint(terminalNode, [2, 4, 7], [0, 0, q0[7]-synRevs*pi]))
    checkJacobian(problem; relTol = JTol)
    shooter = MBD.BCR4BP12MultipleShooter(tol)
    # shooter.printProgress = true
    solution::MBD.BCR4BP12MultipleShooterProblem = MBD.solve!(shooter, problem)

    return solution
end

"""
    getMonodromy(targeter, solution)

Return orbit monodromy matrix

# Arguments
- `targeter::PlanarPerpP12Targeter`: BCR4BP P1-P2 planar perpendicular crossing period targeter object
- `solution::BCR4BP12MultipleShooterProblem`: Solved BCR4BP P1-P2 multiple shooter problem object
"""
function getMonodromy(targeter::PlanarPerpP12Targeter, solution::MBD.BCR4BP12MultipleShooterProblem)
    propagator = MBD.Propagator(equationType = MBD.STM)
    n_simple::Int64 = getStateSize(targeter.dynamicsModel, MBD.SIMPLE)
    n_STM::Int64 = getStateSize(targeter.dynamicsModel, MBD.STM)
    renormalizeEvent::DifferentialEquations.DiscreteCallback = DifferentialEquations.PeriodicCallback(MBD.renormalize!, pi/10)
    Rs::Vector{Matrix{Float64}} = []
    orbit::MBD.BCR4BP12Arc = propagateWithPeriodicEvent(propagator, renormalizeEvent, appendExtraInitialConditions(targeter.dynamicsModel, solution.nodes[1].state.data, MBD.STM), [0, getPeriod(targeter, solution)], targeter.dynamicsModel, [targeter.dynamicsModel, Rs])
    endState::StaticArrays.SVector{n_STM, Float64} = StaticArrays.SVector{n_STM, Float64}(getStateByIndex(orbit, -1))
    M::Matrix{Float64} = reshape(endState[n_simple+1:n_STM], (n_simple,n_simple))
    for R::Matrix{Float64} in reverse(Rs)
        M *= R
    end

    return M
end

"""
    getPeriod(targeter, solution)

Return orbit period

# Arguments
- `targeter::PlanarPerpP12Targeter`: BCR4BP P1-P2 planar perpendicular crossing period targeter object
- `solution::BCR4BP12MultipleShooterProblem`: Solved BCR4BP P1-P2 multiple shooter problem object
"""
function getPeriod(targeter::PlanarPerpP12Targeter, solution::MBD.BCR4BP12MultipleShooterProblem)
    return 2*solution.segments[1].TOF.data[1]
end

"""
    propagateState(targeter, q_simple, tSpan)

Return propagated state

# Arguments
- `targeter::PlanarPerpP12Targeter`: BCR4BP P1-P2 planar perpendicular crossing period targeter object
- `q_simple::Vector{Float64}`: Simple state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
"""
function propagateState(targeter::PlanarPerpP12Targeter, q_simple::Vector{Float64}, tSpan::Vector{Float64})
    propagator = MBD.Propagator()
    arc::MBD.BCR4BP12Arc = propagate(propagator, q_simple, tSpan, targeter.dynamicsModel)

    return getStateByIndex(arc, -1)
end
