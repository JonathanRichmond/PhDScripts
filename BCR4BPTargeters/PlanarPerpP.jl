"""
Period perpendicular crossing targeter for BCR4BP planar orbits

Author: Jonathan Richmond
C: 4/9/25
"""

using MBD

export PlanarPerpP12Targeter
export correct, propagateState, getPeriod#, getIndividualPeriodicOrbit, getMonodromy, interpOrbit

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
    correct(targeter, q0, targetP; tol)

Return corrected BCR4BP P1-P2 multiple shooter problem object

# Arguments
- `targeter::PlanarPerpP12Targeter`: BCR4BP P1-P2 planar perpendicular crossing period targeter object
- `q0::Vector{Float64}`: Initial state guess [ndim]
- `targetP::Float64`: Target period [ndim]
- `tol::Float64`: Convergence tolerance (default = 1E-11)
"""
function correct(targeter::PlanarPerpP12Targeter, q0::Vector{Float64}, targetP::Float64, tol::Float64 = 1E-11)
    halfPeriod::Float64 = targetP/2
    qPCGuess::Vector{Float64} = propagateState(targeter, q0, [0, halfPeriod])
    originNode = MBD.BCR4BP12Node(0.0, q0, targeter.dynamicsModel)
    originNode.state.name = "Initial State"
    setFreeVariableMask!(originNode.state, [true, false, false, false, true, false, false])
    terminalNode = MBD.BCR4BP12Node(halfPeriod, qPCGuess, targeter.dynamicsModel)
    terminalNode.state.name = "Target State"
    segment = MBD.BCR4BP12Segment(halfPeriod, originNode, terminalNode)
    setFreeVariableMask!(segment.TOF, [false])
    problem = MBD.BCR4BP12MultipleShooterProblem()
    addSegment!(problem, segment)
    addConstraint!(problem, MBD.BCR4BP12ContinuityConstraint(segment))
    addConstraint!(problem, MBD.BCR4BP12StateConstraint(terminalNode, [2, 4], [0.0, 0.0]))
    shooter = MBD.BCR4BP12MultipleShooter(tol)
    # shooter.printProgress = true
    solution::MBD.BCR4BP12MultipleShooterProblem = MBD.solve!(shooter, problem)

    return solution
end

"""
    getPeriod(targeter, solution)

Return orbit period

# Arguments
- `targeter::PlanarPerpP12Targeter`: BCR4BP P1-P2 planar perpendicular crossing period targeter object
- `solution::BR4BP12MultipleShooterProblem`: Solved BCR4BP P1-P2 multiple shooter problem object
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
