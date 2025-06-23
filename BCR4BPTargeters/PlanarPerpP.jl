"""
Period perpendicular crossing targeter for BCR4BP planar orbits

Author: Jonathan Richmond
C: 4/9/25
U: 6/23/25
"""

using MBD, DifferentialEquations, Logging, StaticArrays

export PlanarPerpP12Targeter
export correct, getMonodromy, getPeriod, getResonantOrbit, propagateState

"""
    PlanarPerpP12Targeter(dynamicsModel)

BCR4BP P1-P2 planar perpendicular crossing period targeter object

# Arguments
- `dynamicsModel::BCR4BP12DynamicsModel`: BCR4BP P1-P2 dynamics model object
"""
struct PlanarPerpP12Targeter
    dynamicsModel::MBD.BCR4BP12DynamicsModel                            # BCR4BP P1-P2 dynamics model object

    function PlanarPerpP12Targeter(dynamicsModel::MBD.BCR4BP12DynamicsModel)
        this = new(dynamicsModel)

        return this
    end
end

"""
    correct(targeter, q0, targetP; tol, JTol)

Return corrected BCR4BP P1-P2 multiple shooter problem object

# Arguments
- `targeter::PlanarPerpP12Targeter`: BCR4BP P1-P2 planar perpendicular crossing period targeter object
- `q0::Vector{Float64}`: Initial state guess [ndim]
- `targetP::Float64`: Target period [ndim]
- `tol::Float64`: Convergence tolerance (default = 1E-11)
- `JTol::Float64`: Jacobian accuracy tolerance (default = 2E-3)
"""
function correct(targeter::PlanarPerpP12Targeter, q0::Vector{Float64}, targetP::Float64; tol::Float64 = 1E-11, JTol::Float64 = 2E-3)
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
    getResonantOrbit(targeter, initialOrbit, theta40, p, q, boundingBoxJumpCheck; Deltaeps, tol, refTol, JTol)

Return synodic-resonant periodic orbit

# Arguments
- `targeter::PlanarPerpP12Targeter`: BCR4BP P1-P2 planar perpendicular crossing period targeter object
- `initialOrbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit initial guess
- `theta40::Float64`: Initial P4 angle [ndim]
- `p::Int64`: Orbit revolutions
- `q::Int64`: Synodic revolutions
- `boundingBoxJumpCheck::BoundingBoxJumpCheck`: Bounding box jump check object
- `Deltaeps::Float64`: Homotopy parameter step (Sun mass)
- `tol::Float64`: Convergence tolerance (default = 1E-11)
- `refTol::Float64`: Refined solution convergence tolerance (default = 1E-11)
- `JTol::Float64`: Jacobian accuracy tolerance (default = 2E-3)
"""
function getResonantOrbit(targeter::PlanarPerpP12Targeter, initialOrbit::MBD.CR3BPPeriodicOrbit, theta40::Float64, p::Int64, q::Int64, boundingBoxJumpCheck::MBD.BoundingBoxJumpCheck; Deltaeps = 0.001, tol = 1E-11, refTol = 1E-11, JTol = 2E-3)
    systemData0 = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
    systemData0.P4Mass = 0.0
    dynamicsModel0 = MBD.BCR4BP12DynamicsModel(systemData0)
    targeter0 = PlanarPerpP12Targeter(dynamicsModel0)
    systemData1 = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
    systemData1.P4Mass = 0.001*targeter.dynamicsModel.systemData.P4Mass
    dynamicsModel1 = MBD.BCR4BP12DynamicsModel(systemData1)
    targeter1 = PlanarPerpP12Targeter(dynamicsModel1)
    initialStateGuess::Vector{Float64} = push!(copy(initialOrbit.initialCondition), theta40)
    targetP::Float64 = q*getSynodicPeriod(targeter.dynamicsModel)
    solution0::MBD.BCR4BP12MultipleShooterProblem = correct(targeter0, initialStateGuess, targetP, tol = tol, JTol = JTol)
    Logging.@info "Converged Homotopy Orbit 0:\n\tIC:\t$(solution0.nodes[1].state.data[1:7])\n\tP:\t$(getPeriod(targeter0, solution0))\n"
    solution1::MBD.BCR4BP12MultipleShooterProblem = correct(targeter1, copy(solution0.nodes[1].state.data[1:7]), targetP, tol = tol, JTol = JTol)
    Logging.@info "Converged Homotopy Orbit 1:\n\tIC:\t$(solution1.nodes[1].state.data[1:7])\n\tP:\t$(getPeriod(targeter1, solution1))\n"
    continuationEngine = MBD.P4MassContinuationEngine(solution0, solution1, Deltaeps, 0.1, tol = tol, JTol = JTol)
    addJumpCheck!(continuationEngine, boundingBoxJumpCheck)
    homotopyEndCheck = MBD.HomotopyEndCheck(1.0)
    addEndCheck!(continuationEngine, homotopyEndCheck)
    # continuationEngine.printProgress = true
    solutions::MBD.BCR4BP12ContinuationFamily = doContinuation!(continuationEngine, solution0, solution1)
    Logging.@info "Converged Last Homotopy Orbit:\n\tIC:\t$(solutions.nodes[end][1].state.data[1:7])\n\tP:\t$(2*solutions.segments[end][1].TOF.data[1])\n"
    solution::MBD.BCR4BP12MultipleShooterProblem = correct(targeter, solutions.nodes[end][1].state.data[1:7], targetP, tol = refTol, JTol = JTol)
    orbit = MBD.BCR4BP12PeriodicOrbit(targeter.dynamicsModel, solution.nodes[1].state.data[1:7], getPeriod(targeter, solution), getMonodromy(targeter, solution))
    println("Converged $p:$q BCR4BP Orbit:\n\tIC:\t$(orbit.initialCondition)\n\tP:\t$(orbit.period)\n")

    return orbit
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
