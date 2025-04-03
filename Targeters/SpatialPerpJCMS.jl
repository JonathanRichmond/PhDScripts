"""
Jacobi constant perpendicular crossing multiple shooter for spatial orbits

Author: Jonathan Richmond
C: 2/26/25
U: 4/2/25
"""

using MBD, CSV, DataFrames, LinearAlgebra, StaticArrays

export SpatialPerpJCMSTargeter
export correct, doContinuation!, getPeriod, propagateState, tryConverging!#, getIndividualPeriodicOrbit, getMonodromy, interpOrbit

"""
    SpatialPerpJCMSTargeter(dynamicsModel)

CR3BP spatial Jacobi constant perpendicular crossing multiple shooter object

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
"""
struct SpatialPerpJCMSTargeter <: MBD.AbstractTargeter
    dynamicsModel::MBD.CR3BPDynamicsModel                               # CR3BP dynamics model object

    function SpatialPerpJCMSTargeter(dynamicsModel::MBD.CR3BPDynamicsModel)
        this = new(dynamicsModel)

        return this
    end
end

"""
    correct(targeter, q0, tSpan, numSegs, targetJC; tol)

Return corrected CR3BP multiple shooter problem object

# Arguments
- `targeter::SpatialPerpJCMSTargeter`: CR3BP spatial Jacobi constant perpendicular crossing multiple shooter object
- `q0::Vector{Float64}`: Initial state guess [ndim]
- `tSpan::Vector{Float64}`: Time span guess [ndim]
- `numSegs::Int64`: Number of segments
- `targetJC::Float64`: Target Jacobi constant
- `tol::Float64`: Convergence tolerance (default = 1E-11)
"""
function correct(targeter::SpatialPerpJCMSTargeter, q0::Vector{Float64}, tSpan::Vector{Float64}, numSegs::Int64, targetJC::Float64, tol::Float64 = 1E-11)
    halfPeriodGuess::Float64 = (tSpan[2]-tSpan[1])/2
    tPCGuess::Float64 = tSpan[1]+halfPeriodGuess
    propagator = MBD.Propagator()
    arc::MBD.CR3BPArc = propagate(propagator, q0, [tSpan[1], tPCGuess], targeter.dynamicsModel)
    numStates::Int64 = getStateCount(arc)
    numNodes::Int64 = numSegs/2+1
    indices::Vector{Int64} = round.(Int64, range(1, numStates, numNodes))
    stateGuesses::Vector{Vector{Float64}} = [getStateByIndex(arc, indices[i]) for i = eachindex(indices)]
    timeGuesses::Vector{Float64} = [getTimeByIndex(arc, indices[i]) for i = eachindex(indices)]
    nodes::Vector{MBD.CR3BPNode} = []
    for n = 1:numNodes
        node = MBD.CR3BPNode(timeGuesses[n], stateGuesses[n], dynamicsModel)
        node.state.name = "Node "*string(n)*" State"
        node.epoch.name = "Node "*string(n)*" Epoch"
        push!(nodes, node)
    end    
    setFreeVariableMask!(nodes[1].state, [true, false, true, false, true, false])
    problem = MBD.CR3BPMultipleShooterProblem()
    segments::Vector{MBD.CR3BPSegment} = []
    for s = 1:numNodes-1
        segment = MBD.CR3BPSegment(timeGuesses[s+1]-timeGuesses[s], nodes[s], nodes[s+1])
        segment.TOF.name = "Segment "*string(s)*" TOF"
        push!(segments, segment)
    end
    map(s -> setFreeVariableMask!(s.TOF, [false]), segments[1:end-1])
    map(s -> addSegment!(problem, s), segments)
    map(s -> addConstraint!(problem, MBD.CR3BPContinuityConstraint(s)), segments)
    addConstraint!(problem, MBD.CR3BPStateConstraint(nodes[end], [2, 4, 6], [0.0, 0.0, 0.0]))
    addConstraint!(problem, MBD.JacobiConstraint(nodes[1], targetJC))
    checkJacobian(problem)
    shooter = MBD.CR3BPMultipleShooter(tol)
    # shooter.printProgress = true
    shooter.maxIterations = 50
    solution::MBD.CR3BPMultipleShooterProblem = MBD.solve!(shooter, problem)

    return solution
end

"""
    doContinuation!(targeter, multipleShooterContinuationEngine, initialGuess1, initialGuess2, numSegs; tol)

Return family of solutions

# Arguments
- `targeter::SpatialPerpJCMSTargeter`: CR3BP spatial Jacobi constant perpendicular crossing multiple shooter object
- `multipleShooterContinuationEngine::CR3BPMultipleShooterContinuationEngine`: CR3BP multiple shooter continuation engine object
- `initialGuess1::CR3BPMultipleShooterProblem`: First member of family
- `initialGuess2::CR3BPMultipleShooterProblem`: Second member of family
- `numSegs::Int64`: Number of segments
- `tol::Float64`: Convergence tolerance (default = 1E-11)
"""
function doContinuation!(targeter::SpatialPerpJCMSTargeter, multipleShooterContinuationEngine::MBD.CR3BPMultipleShooterContinuationEngine, initialGuess1::MBD.CR3BPMultipleShooterProblem, initialGuess2::MBD.CR3BPMultipleShooterProblem, numSegs::Int64, tol::Float64 = 1E-11)
    isempty(multipleShooterContinuationEngine.endChecks) && throw(ErrorException("Cannot do continuation without at least one end check"))
    resetEngine!(multipleShooterContinuationEngine, initialGuess1, initialGuess2)
    multipleShooterContinuationEngine.corrector.printProgress = multipleShooterContinuationEngine.printProgress
    multipleShooterContinuationEngine.printProgress && println("Converging initial guesses...")
    multipleShooterContinuationEngine.dataInProgress.twoPreviousSolution = convergeInitialSolution(multipleShooterContinuationEngine, initialGuess1)
    multipleShooterContinuationEngine.dataInProgress.previousSolution = convergeInitialSolution(multipleShooterContinuationEngine, initialGuess2)
    multipleShooterContinuationEngine.dataInProgress.numIterations = multipleShooterContinuationEngine.corrector.recentIterationCount
    push!(multipleShooterContinuationEngine.dataInProgress.family.nodes, [MBD.shallowClone(multipleShooterContinuationEngine.dataInProgress.twoPreviousSolution.nodes[n]) for n = 1:length(multipleShooterContinuationEngine.dataInProgress.twoPreviousSolution.nodes)], [MBD.shallowClone(multipleShooterContinuationEngine.dataInProgress.previousSolution.nodes[n]) for n = 1:length(multipleShooterContinuationEngine.dataInProgress.previousSolution.nodes)])
    push!(multipleShooterContinuationEngine.dataInProgress.family.segments, [MBD.shallowClone(multipleShooterContinuationEngine.dataInProgress.twoPreviousSolution.segments[s]) for s = 1:length(multipleShooterContinuationEngine.dataInProgress.twoPreviousSolution.segments)], [MBD.shallowClone(multipleShooterContinuationEngine.dataInProgress.previousSolution.segments[s]) for s = 1:length(multipleShooterContinuationEngine.dataInProgress.previousSolution.segments)])
    multipleShooterContinuationEngine.dataInProgress.initialGuess = initialGuess2
    multipleShooterContinuationEngine.dataInProgress.converging = true
    multipleShooterContinuationEngine.dataInProgress.forceEndContinuation = false
    multipleShooterContinuationEngine.dataInProgress.currentStepSize = multipleShooterContinuationEngine.stepSizeGenerator.initialStepSize
    while (!endContinuation(multipleShooterContinuationEngine, multipleShooterContinuationEngine.dataInProgress) && !multipleShooterContinuationEngine.dataInProgress.forceEndContinuation)
        multipleShooterContinuationEngine.printProgress && println("\nConverging family member $(getNumSteps(multipleShooterContinuationEngine.dataInProgress)+1)...")
        tryConverging!(targeter, multipleShooterContinuationEngine, numSegs, tol)
        while (!multipleShooterContinuationEngine.dataInProgress.converging && !multipleShooterContinuationEngine.dataInProgress.forceEndContinuation)
            tryConverging!(targeter, multipleShooterContinuationEngine, numSegs, tol)
        end
        if (multipleShooterContinuationEngine.storeIntermediateMembers && multipleShooterContinuationEngine.dataInProgress.converging)
            println("\tConverged Initial State: $(multipleShooterContinuationEngine.dataInProgress.previousSolution.nodes[1].state.data[1:6])")
            push!(multipleShooterContinuationEngine.dataInProgress.family.nodes, [MBD.shallowClone(multipleShooterContinuationEngine.dataInProgress.previousSolution.nodes[n]) for n = 1:length(multipleShooterContinuationEngine.dataInProgress.previousSolution.nodes)])
            push!(multipleShooterContinuationEngine.dataInProgress.family.segments, [MBD.shallowClone(multipleShooterContinuationEngine.dataInProgress.previousSolution.segments[s]) for s = 1:length(multipleShooterContinuationEngine.dataInProgress.previousSolution.segments)])
        end
    end
    if (!multipleShooterContinuationEngine.dataInProgress.converging && (getNumSteps(multipleShooterContinuationEngine.dataInProgress) == 2))
        throw(ErrorException("Could not converge any solutions beyond initial guess"))
    end
    if (!multipleShooterContinuationEngine.storeIntermediateMembers && (getNumSteps(multipleShooterContinuationEngine.dataInProgress) > 2))
        push!(multipleShooterContinuationEngine.dataInProgress.family.nodes, [MBD.shallowClone(multipleShooterContinuationEngine.dataInProgress.previousSolution.nodes[n]) for n = 1:length(multipleShooterContinuationEngine.dataInProgress.previousSolution.nodes)])
        push!(multipleShooterContinuationEngine.dataInProgress.family.segments, [MBD.shallowClone(multipleShooterContinuationEngine.dataInProgress.previousSolution.segments[s]) for s = 1:length(multipleShooterContinuationEngine.dataInProgress.previousSolution.segments)])
    end

    return multipleShooterContinuationEngine.dataInProgress.family
end

"""
    getPeriod(targeter, segments)

Return orbit period

# Arguments
- `targeter::SpatialPerpJCMSTargeter`: CR3BP spatial Jacobi constant perpendicular crossing multiple shooter object
- `segments::Vector{CR3BPSegment}`: CR3BP segment objects
"""
function getPeriod(targeter::SpatialPerpJCMSTargeter, segments::Vector{MBD.CR3BPSegment})
    TOF::Float64 = 0.0
    for s = 1:length(segments)
        TOF += segments[s].TOF.data[1]
    end

    return 2*TOF
end

"""
    getPeriod(targeter, solution)

Return orbit period

# Arguments
- `targeter::SpatialPerpJCMSTargeter`: CR3BP spatial Jacobi constant perpendicular crossing multiple shooter object
- `solution::CR3BPMultipleShooterProblem`: Solved CR3BP multiple shooter problem object
"""
function getPeriod(targeter::SpatialPerpJCMSTargeter, solution::MBD.CR3BPMultipleShooterProblem)
    TOF::Float64 = 0.0
    for s = 1:length(solution.segments)
        TOF += solution.segments[s].TOF.data[1]
    end

    return 2*TOF
end

"""
    propagateState(targeter, q_simple, tSpan)

Return propagated state

# Arguments
- `targeter::SpatialPerpJCMSTargeter`: CR3BP spatial Jacobi constant perpendicular crossing multiple shooter object
- `q_simple::Vector{Float64}`: Simple state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
"""
function propagateState(targeter::SpatialPerpJCMSTargeter, q_simple::Vector{Float64}, tSpan::Vector{Float64})
    propagator = MBD.Propagator()
    arc::MBD.CR3BPArc = propagate(propagator, q_simple, tSpan, targeter.dynamicsModel)

    return getStateByIndex(arc, -1)
end

"""
    tryConverging!(targeter, multipleShooterContinuationEngine, numSegs; tol)

Return updated CR3BP multiple shooter continuation engine object

# Arguments
- `targeter::SpatialPerpJCMSTargeter`: CR3BP spatial Jacobi constant perpendicular crossing multiple shooter object
- `multipleShooterContinuationEngine::CR3BPMultipleShooterContinuationEngine`: CR3BP multiple shooter continuation engine object
- `numSegs::Int64`: Number of segments
- `tol::Float64`: Convergence tolerance (default = 1E-11)
"""
function tryConverging!(targeter::SpatialPerpJCMSTargeter, multipleShooterContinuationEngine::MBD.CR3BPMultipleShooterContinuationEngine, numSegs::Int64, tol::Float64 = 1E-11)
    stateStep::Vector{Float64} = computeStateStep(multipleShooterContinuationEngine, multipleShooterContinuationEngine.dataInProgress)
    timeStep::Float64 = computeTimeStep(multipleShooterContinuationEngine, multipleShooterContinuationEngine.dataInProgress)
    updateStepSize!(multipleShooterContinuationEngine.stepSizeGenerator, multipleShooterContinuationEngine.dataInProgress)
    multipleShooterContinuationEngine.printProgress && println("\tCurrent step size: $(multipleShooterContinuationEngine.dataInProgress.currentStepSize)")
    try
        twoPreviousConvergedSolution::MBD.CR3BPMultipleShooterProblem = multipleShooterContinuationEngine.dataInProgress.twoPreviousSolution
        previousConvergedSolution::MBD.CR3BPMultipleShooterProblem = multipleShooterContinuationEngine.dataInProgress.previousSolution
        targJC::Float64 = getJacobiConstant(targeter.dynamicsModel, previousConvergedSolution.nodes[1].state.data[1:6])+multipleShooterContinuationEngine.dataInProgress.currentStepSize
        println("\tTarget Jacobi Constant: $targJC")
        multipleShooterContinuationEngine.dataInProgress.twoPreviousSolution = MBD.deepClone(multipleShooterContinuationEngine.dataInProgress.previousSolution)
        multipleShooterContinuationEngine.dataInProgress.previousSolution = correct(targeter, previousConvergedSolution.nodes[1].state.data+stateStep.*multipleShooterContinuationEngine.dataInProgress.currentStepSize, [0, getPeriod(targeter, previousConvergedSolution)+timeStep*multipleShooterContinuationEngine.dataInProgress.currentStepSize], numSegs, targJC, tol)
        multipleShooterContinuationEngine.dataInProgress.converging = true
        for jumpCheck::MBD.AbstractContinuationJumpCheck in multipleShooterContinuationEngine.jumpChecks
            if typeof(jumpCheck) == MBD.BoundingBoxJumpCheck
                for (index::MBD.Variable, value::Int16) in multipleShooterContinuationEngine.dataInProgress.previousSolution.freeVariableIndexMap
                    if index.name == jumpCheck.paramName
                        addBounds!(jumpCheck, multipleShooterContinuationEngine.dataInProgress.previousSolution, index, jumpCheck.paramBounds)
                        multipleShooterContinuationEngine.dataInProgress.converging = isFamilyMember(jumpCheck, multipleShooterContinuationEngine.dataInProgress)
                        !multipleShooterContinuationEngine.dataInProgress.converging && println("\tSolution jumped")
                        removeBounds!(jumpCheck, multipleShooterContinuationEngine.dataInProgress.previousSolution, index)
                    end
                end
            end
        end
        if multipleShooterContinuationEngine.dataInProgress.converging
            multipleShooterContinuationEngine.dataInProgress.numIterations = multipleShooterContinuationEngine.corrector.recentIterationCount
            multipleShooterContinuationEngine.dataInProgress.nextGuess = MBD.deepClone(multipleShooterContinuationEngine.dataInProgress.previousSolution)
        else
            multipleShooterContinuationEngine.dataInProgress.twoPreviousSolution = twoPreviousConvergedSolution
            multipleShooterContinuationEngine.dataInProgress.previousSolution = previousConvergedSolution
        end
    catch err
        multipleShooterContinuationEngine.dataInProgress.converging = false
        println("\tFailed to converge")
        # @error exception = (err, catch_backtrace())
    end
end
