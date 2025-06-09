"""
Period perpendicular crossing targeter for BCR4BP planar orbits

Author: Jonathan Richmond
C: 4/9/25
U: 6/9/25
"""

using MBD, DifferentialEquations, Logging, StaticArrays

export PlanarPerpP12Targeter
export correct, getMonodromy, getPeriod, getResonantOrbit, homotopy, propagateState

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
    getResonantOrbit(targeter, initialOrbit, p, q; Deltaeps, tol, refTol, JTol)

Return synodic-resonant periodic orbit

# Arguments
- `targeter::PlanarPerpP12Targeter`: BCR4BP P1-P2 planar perpendicular crossing period targeter object
- `initialOrbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit initial guess
- `p::Int64`: Orbit revolutions
- `q::Int64`: Synodic revolutions
- `Deltaeps::Float64`: Homotopy parameter step (Sun mass)
- `tol::Float64`: Convergence tolerance (default = 1E-10)
- `refTol::Float64`: Refined solution convergence tolerance (default = 1E-11)
- `JTol::Float64`: Jacobian accuracy tolerance (default = 2E-3)
"""
function getResonantOrbit(targeter::PlanarPerpP12Targeter, initialOrbit::MBD.CR3BPPeriodicOrbit, p::Int64, q::Int64; Deltaeps = 0.001, tol = 1E-11, refTol = 1E-11, JTol = 2E-3)
    initialStateGuess::Vector{Float64} = push!(copy(initialOrbit.initialCondition), 0.0)
    solution::MBD.BCR4BP12MultipleShooterProblem = homotopy(targeter, initialStateGuess, p, q, Deltaeps = Deltaeps, tol = tol, JTol = JTol)
    refinedSolution::MBD.BCR4BP12MultipleShooterProblem = correct(targeter, solution.nodes[1].state.data[1:7], q*getSynodicPeriod(targeter.dynamicsModel), q, tol = refTol, JTol = JTol)
    println("Converged $p:$q BCR4BP Orbit:\n\tIC:\t$(refinedSolution.nodes[1].state.data[1:7])\n\tP:\t$(getPeriod(targeter, refinedSolution))\n")
    orbit = MBD.BCR4BP12PeriodicOrbit(dynamicsModel, refinedSolution.nodes[1].state.data[1:7], getPeriod(targeter, refinedSolution), getMonodromy(targeter, refinedSolution))

    return orbit
end

"""
    homotopy(targeter, CR3BPq0, p, q; Deltaeps, tol, JTol)

Return continued BCR4BP P1-P2 solution

# Arguments
- `targeter::PlanarPerpP12Targeter`: BCR4BP P1-P2 planar perpendicular crossing period targeter object
- `CR3BPq0::Vector{Float64}`: CR3BP initial condition [ndim]
- `p::Int64`: Orbit revolutions
- `q::Int64`: Synodic revolutions
- `Deltaeps::Float64`: Homotopy parameter step (Sun mass)
- `tol::Float64`: Convergence tolerance (default = 1E-11)
- `JTol::Float64`: Jacobian accuracy tolerance (default = 2E-3)
"""
function homotopy(targeter::PlanarPerpP12Targeter, CR3BPq0::Vector{Float64}, p::Int64, q::Int64; Deltaeps::Float64 = 0.001, tol::Float64 = 1E-11, JTol::Float64 = 2E-3)
    initialStateGuess::Vector{Float64} = copy(CR3BPq0)
    homoSystemData::MBD.BCR4BPSystemData = MBD.shallowClone(targeter.dynamicsModel.systemData)
    homoDynamicsModel = MBD.BCR4BP12DynamicsModel(homoSystemData)
    homoTargeter = PlanarPerpP12Targeter(homoDynamicsModel)
    SunGravParam::Float64 = copy(homoSystemData.primaryData[3].gravParam)
    homoSystemData.primaryData[3].gravParam = 0.0
    targetP::Float64 = q*getSynodicPeriod(homoDynamicsModel)
    solution::MBD.BCR4BP12MultipleShooterProblem = correct(homoTargeter, CR3BPq0, targetP, q; tol = tol, JTol = JTol)
    println("\nConverged $p:$q CR3BP Orbit:\n\tIC:\t$(solution.nodes[1].state.data[1:7])\n\tP:\t$(getPeriod(homoTargeter, solution))\n")
    
    epsilon = Deltaeps
    while epsilon < 1.0
        homoSystemData.primaryData[3].gravParam = epsilon*SunGravParam
        targetP = q*getSynodicPeriod(homoDynamicsModel)
        initialStateGuess = solution.nodes[1].state.data[1:7]
        solution = correct(homoTargeter, initialStateGuess, targetP, q; tol = tol, JTol = JTol)
        Logging.@info "Converged $p:$q Homotopy Orbit (epsilon = $epsilon):\n\tIC:\t$(solution.nodes[1].state.data[1:7])\n\tP:\t$(getPeriod(homoTargeter, solution))\n"
        epsilon += Deltaeps
    end

    return solution
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
