"""
Period state continuity targeter for BCR4BP spatial orbits

Author: Jonathan Richmond
C: 6/25/25
U: 7/9/25
"""

using MBD, DifferentialEquations, Logging, StaticArrays

export SpatialContP12Targeter
export correct, getMonodromy, getPeriod, getResonantOrbit, propagateState

"""
    SpatialContP12Targeter(dynamicsModel)

BCR4BP P1-P2 spatial state continuity period targeter object

# Arguments
- `dynamicsModel::BCR4BP12DynamicsModel`: BCR4BP P1-P2 dynamics model object
"""
struct SpatialContP12Targeter
    dynamicsModel::MBD.BCR4BP12DynamicsModel                            # BCR4BP P1-P2 dynamics model object

    function SpatialContP12Targeter(dynamicsModel::MBD.BCR4BP12DynamicsModel)
        this = new(dynamicsModel)

        return this
    end
end

"""
    correct(targeter, qVector, tVector, numNodes; tol, JTol)

Return corrected BCR4BP P1-P2 multiple shooter problem object

# Arguments
- `targeter::SpatialContP12Targeter`: BCR4BP P1-P2 spatial state continuity period targeter object
- `qVector::Vector{Vector{Float64}}`: Node state guesses [ndim]
- `tVector::Vector{Float64}`: Node time guesses [ndim]
- `numNodes::Int64`: Number of nodes
- `tol::Float64`: Convergence tolerance (default = 1E-11)
- `JTol::Float64`: Jacobian accuracy tolerance (default = 2E-3)
"""
function correct(targeter::SpatialContP12Targeter, qVector::Vector{Vector{Float64}}, tVector::Vector{Float64}, numNodes::Int64; tol::Float64 = 1E-11, JTol::Float64 = 2E-3)
    nodes::Vector{MBD.BCR4BP12Node} = []
    for n = 1:numNodes
        node = MBD.BCR4BP12Node(tVector[n], qVector[n], targeter.dynamicsModel)
        node.state.name = "Node "*string(n)*" State"
        node.epoch.name = "Node "*string(n)*" Epoch"
        push!(nodes, node)
    end    
    setFreeVariableMask!(nodes[1].state, [true, false, true, false, true, false, false])
    problem = MBD.BCR4BP12MultipleShooterProblem()
    segments::Vector{MBD.BCR4BP12Segment} = []
    for s = 1:numNodes-1
        segment = MBD.BCR4BP12Segment(tVector[s+1]-tVector[s], nodes[s], nodes[s+1])
        segment.TOF.name = "Segment "*string(s)*" TOF"
        setFreeVariableMask!(segment.TOF, [false])
        push!(segments, segment)
    end
    map(s -> addSegment!(problem, s), segments)
    map(s -> addConstraint!(problem, MBD.BCR4BP12ContinuityConstraint(s)), segments)
    addConstraint!(problem, MBD.BCR4BP12StateMatchConstraint(nodes[1].state, nodes[end].state, [1, 2, 3, 4, 6]))
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
- `targeter::SpatialContP12Targeter`: BCR4BP P1-P2 spatial state continuity period targeter object
- `solution::BCR4BP12MultipleShooterProblem`: Solved BCR4BP P1-P2 multiple shooter problem object
"""
function getMonodromy(targeter::SpatialContP12Targeter, solution::MBD.BCR4BP12MultipleShooterProblem)
    propagator = MBD.Propagator(equationType = MBD.STM)
    nStates::Int64 = getStateSize(targeter.dynamicsModel, MBD.STM)
    nSegs::Int64 = length(solution.segments)
    renormalizeEvent::DifferentialEquations.DiscreteCallback = DifferentialEquations.PeriodicCallback(MBD.renormalize!, 0.1)
    STMs::Vector{Matrix{Float64}} = []
    for s::Int64 = 1:nSegs
        Rs::Vector{Matrix{Float64}} = []
        propSegment::MBD.BCR4BP12Arc = propagateWithPeriodicEvent(propagator, renormalizeEvent, appendExtraInitialConditions(targeter.dynamicsModel, solution.segments[s].originNode.state.data, MBD.STM), [0, solution.segments[s].TOF.data[1]], targeter.dynamicsModel, [targeter.dynamicsModel, Rs])
        endState::StaticArrays.SVector{nStates, Float64} = StaticArrays.SVector{nStates, Float64}(getStateByIndex(propSegment, -1))
        STM::Matrix{Float64} = reshape(endState[8:56], (7,7))
        for R::Matrix{Float64} in reverse(Rs)
            STM *= R
        end
        push!(STMs, STM)
    end
    M::Matrix{Float64} = Matrix(STMs[end])
    for Phi::Matrix{Float64} in reverse(STMs[1:end-1])
        M *= Phi
    end

    return M
end

"""
    getPeriod(targeter, segments)

Return orbit period

# Arguments
- `targeter::SpatialContP12Targeter`: BCR4BP P1-P2 spatial state continuity period targeter object
- `segments::Vector{BCR4BP12Segment}`: BCR4BP P1-P2 segment objects
"""
function getPeriod(targeter::SpatialContP12Targeter, segments::Vector{MBD.BCR4BP12Segment})
    TOF::Float64 = 0.0
    for s::Int64 = 1:length(segments)
        TOF += segments[s].TOF.data[1]
    end

    return TOF
end

"""
    getPeriod(targeter, solution)

Return orbit period

# Arguments
- `targeter::SpatialContP12Targeter`: BCR4BP P1-P2 spatial state continuity period targeter object
- `solution::BCR4BP12MultipleShooterProblem`: Solved BCR4BP P1-P2 multiple shooter problem object
"""
function getPeriod(targeter::SpatialContP12Targeter, solution::MBD.BCR4BP12MultipleShooterProblem)
    return getPeriod(targeter, solution.segments)
end

"""
    getResonantOrbit(targeter, initialOrbit, numSegs, theta40, p, q, boundingBoxJumpCheck; targeteps, Deltaeps, tol, refTol, JTol)

Return synodic-resonant periodic orbit

# Arguments
- `targeter::SpatialContP12Targeter`: BCR4BP P1-P2 spatial state continuity period targeter object
- `initialOrbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit initial guess
- `numSegs::Int64`: NUmber of segments
- `theta40::Float64`: Initial P4 angle [ndim]
- `p::Int64`: Orbit revolutions
- `q::Int64`: Synodic revolutions
- `boundingBoxJumpCheck::BoundingBoxJumpCheck`: Bounding box jump check object
- `targeteps::Float64`: Target homotopy parameter (default = 1.0)
- `Deltaeps::Float64`: Homotopy parameter step (Sun mass, default = 0.001)
- `tol::Float64`: Convergence tolerance (default = 1E-11)
- `refTol::Float64`: Refined solution convergence tolerance (default = 1E-11)
- `JTol::Float64`: Jacobian accuracy tolerance (default = 2E-3)
"""
function getResonantOrbit(targeter::SpatialContP12Targeter, initialOrbit::MBD.CR3BPPeriodicOrbit, numSegs::Int64, theta40::Float64, p::Int64, q::Int64, boundingBoxJumpCheck::MBD.BoundingBoxJumpCheck; targeteps::Float64 = 1.0, Deltaeps = 0.001, tol = 1E-11, refTol = 1E-11, JTol = 2E-3)
    systemData0 = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
    systemData0.P4Mass = 0.0
    dynamicsModel0 = MBD.BCR4BP12DynamicsModel(systemData0)
    targeter0 = SpatialContP12Targeter(dynamicsModel0)
    propagator = MBD.Propagator()
    initialStateGuess::Vector{Float64} = push!(copy(initialOrbit.initialCondition), theta40)
    targetTime::Float64 = p*initialOrbit.period/2
    orbitArc::MBD.BCR4BP12Arc = propagate(propagator, initialStateGuess, [0, initialOrbit.period], dynamicsModel0)
    guessArc = MBD.BCR4BP12Arc(dynamicsModel0)
    currentTime::Float64 = 0.0
    push!(guessArc.states, getStateByIndex(orbitArc, 1))
    push!(guessArc.times, getTimeByIndex(orbitArc, 1))
    while currentTime+initialOrbit.period <= targetTime
        append!(guessArc.states, orbitArc.states[2:end])
        append!(guessArc.times, orbitArc.times[2:end].+currentTime)
        currentTime = getTimeByIndex(guessArc, -1)
    end
    if (targetTime-currentTime) > 1E-10
        partialArc::MBD.BCR4BP12Arc = propagate(propagator, initialStateGuess, [0, targetTime-currentTime], dynamicsModel0)
        append!(guessArc.states, partialArc.states[2:end])
        append!(guessArc.times, partialArc.times[2:end].+currentTime)
    end
    numStates::Int64 = getStateCount(guessArc)
    numNodes::Int64 = numSegs+1
    indices::Vector{Int64} = round.(Int64, range(1, numStates, numNodes))
    stateGuesses::Vector{Vector{Float64}} = [getStateByIndex(guessArc, indices[i]) for i = eachindex(indices)]
    timeGuesses::Vector{Float64} = [getTimeByIndex(guessArc, indices[i]) for i = eachindex(indices)]
    solution0::MBD.BCR4BP12MultipleShooterProblem = correct(targeter0, stateGuesses, timeGuesses, numNodes, tol = tol, JTol = JTol)
    Logging.@info "Converged Homotopy Orbit 0:\n\tIC:\t$(solution0.nodes[1].state.data[1:7])\n\tP:\t$(getPeriod(targeter0, solution0))\n"
    systemData1 = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
    systemData1.P4Mass = 0.001*targeter.dynamicsModel.systemData.P4Mass
    dynamicsModel1 = MBD.BCR4BP12DynamicsModel(systemData1)
    targeter1 = SpatialContP12Targeter(dynamicsModel1)
    solution1::MBD.BCR4BP12MultipleShooterProblem = correct(targeter1, [solution0.nodes[n].state.data[1:7] for n = 1:numNodes], [solution0.nodes[n].epoch.data[1] for n = 1:numNodes], numNodes, tol = tol, JTol = JTol)
    Logging.@info "Converged Homotopy Orbit 1:\n\tIC:\t$(solution1.nodes[1].state.data[1:7])\n\tP:\t$(getPeriod(targeter1, solution1))\n"
    continuationEngine = MBD.P4MassContinuationEngine(solution0, solution1, Deltaeps, 0.1, tol = tol, JTol = JTol)
    addJumpCheck!(continuationEngine, boundingBoxJumpCheck)
    homotopyEndCheck = MBD.HomotopyEndCheck(targeteps)
    addEndCheck!(continuationEngine, homotopyEndCheck)
    # continuationEngine.printProgress = true
    solutions::MBD.BCR4BP12ContinuationFamily = doContinuation!(continuationEngine, solution0, solution1)
    (continuationEngine.dataInProgress.converging == false) && throw(ErrorException("Homotopy could not converge: Failed at $(continuationEngine.dataInProgress.previousSolution.nodes[1].dynamicsModel.systemData.P4Mass/targeter.dynamicsModel.systemData.primaryData[3].mass) with IC = $(continuationEngine.dataInProgress.previousSolution.nodes[1].state.data[1:7])"))
    Logging.@info "Converged Last Homotopy Orbit:\n\tIC:\t$(solutions.nodes[end][1].state.data[1:7])\n\tP:\t$(getPeriod(targeter, solutions.segments[end]))\n"
    systemDataEnd = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
    systemDataEnd.P4Mass = targeteps*targeter.dynamicsModel.systemData.P4Mass
    dynamicsModelEnd = MBD.BCR4BP12DynamicsModel(systemDataEnd)
    targeterEnd = SpatialContP12Targeter(dynamicsModelEnd)
    solution::MBD.BCR4BP12MultipleShooterProblem = correct(targeterEnd, [solutions.nodes[end][n].state.data[1:7] for n = 1:numNodes], [solutions.nodes[end][n].epoch.data[1] for n = 1:numNodes], numNodes, tol = refTol, JTol = JTol)
    orbit = MBD.BCR4BP12PeriodicOrbit(targeterEnd.dynamicsModel, [solution.nodes[n].state.data[1:7] for n = 1:numNodes], [solution.nodes[n].epoch.data[1] for n = 1:numNodes], getPeriod(targeterEnd, solution), getMonodromy(targeterEnd, solution))
    println("Converged $p:$q BCR4BP Orbit:\n\tIC:\t$(orbit.initialCondition)\n\tP:\t$(orbit.period)\n\tStab.:\t$(getStabilityIndex(orbit))\n")

    return orbit
end

"""
    propagateState(targeter, q_simple, tSpan)

Return propagated state

# Arguments
- `targeter::SpatialContP12Targeter`: BCR4BP P1-P2 spatial state continuity period targeter object
- `q_simple::Vector{Float64}`: Simple state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
"""
function propagateState(targeter::SpatialContP12Targeter, q_simple::Vector{Float64}, tSpan::Vector{Float64})
    propagator = MBD.Propagator()
    arc::MBD.BCR4BP12Arc = propagate(propagator, q_simple, tSpan, targeter.dynamicsModel)

    return getStateByIndex(arc, -1)
end
