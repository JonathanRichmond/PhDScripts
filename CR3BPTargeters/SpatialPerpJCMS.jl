"""
Jacobi constant perpendicular crossing multiple shooter for CR3BP spatial orbits

Author: Jonathan Richmond
C: 2/26/25
U: 7/16/25
"""

using MBD, CSV, DataFrames, DifferentialEquations, LinearAlgebra, StaticArrays, Statistics

export SpatialPerpJCMSTargeter
export correct, getIndividualPeriodicOrbit, getMonodromy, getPeriod, interpOrbit, propagateState

"""
    SpatialPerpJCMSTargeter(dynamicsModel)

CR3BP spatial Jacobi constant perpendicular crossing multiple shooter object

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
"""
struct SpatialPerpJCMSTargeter
    dynamicsModel::MBD.CR3BPDynamicsModel                               # CR3BP dynamics model object

    function SpatialPerpJCMSTargeter(dynamicsModel::MBD.CR3BPDynamicsModel)
        this = new(dynamicsModel)

        return this
    end
end

"""
    correct(targeter, qVector, tVector, numNodes, targetJC; tol, JTol)

Return corrected CR3BP multiple shooter problem object

# Arguments
- `targeter::SpatialPerpJCMSTargeter`: CR3BP spatial Jacobi constant perpendicular crossing multiple shooter object
- `qVector::Vector{Vector{Float64}}`: Node state guesses [ndim]
- `tVector::Vector{Float64}`: Node time guesses [ndim]
- `numNodes::Int64`: Number of nodes
- `targetJC::Float64`: Target Jacobi constant
- `tol::Float64`: Convergence tolerance (default = 1E-11)
- `JTol::Float64`: Jacobian accuracy tolerance (default = 2E-3)
"""
function correct(targeter::SpatialPerpJCMSTargeter, qVector::Vector{Vector{Float64}}, tVector::Vector{Float64}, numNodes::Int64, targetJC::Float64; tol::Float64 = 1E-11, JTol::Float64 = 2E-3)
    nodes::Vector{MBD.CR3BPNode} = []
    for n = 1:numNodes
        node = MBD.CR3BPNode(tVector[n], qVector[n], targeter.dynamicsModel)
        node.state.name = "Node "*string(n)*" State"
        node.epoch.name = "Node "*string(n)*" Epoch"
        push!(nodes, node)
    end    
    setFreeVariableMask!(nodes[1].state, [true, false, true, false, true, false])
    problem = MBD.CR3BPMultipleShooterProblem()
    segments::Vector{MBD.CR3BPSegment} = []
    for s = 1:numNodes-1
        segment = MBD.CR3BPSegment(tVector[s+1]-tVector[s], nodes[s], nodes[s+1])
        segment.TOF.name = "Segment "*string(s)*" TOF"
        push!(segments, segment)
    end
    map(s -> addSegment!(problem, s), segments)
    map(s -> addConstraint!(problem, MBD.CR3BPContinuityConstraint(s)), segments)
    addConstraint!(problem, MBD.CR3BPStateConstraint(nodes[end], [2, 4, 6], [0.0, 0.0, 0.0]))
    addConstraint!(problem, MBD.JacobiConstraint(nodes[1], targetJC))
    checkJacobian(problem, relTol = JTol)
    shooter = MBD.CR3BPMultipleShooter(tol)
    # shooter.printProgress = true
    solution::MBD.CR3BPMultipleShooterProblem = MBD.solve!(shooter, problem)
    map(s -> updateTerminalNodeEpoch!(s), solution.segments)

    return solution
end

"""
    getIndividualPeriodicOrbit(targeter, family, orbit)

Return periodic orbit object

# Arguments
- `targeter::SpatialPerpJCMSTargeter`: CR3BP spatial perpendicular crossing Jacobi constant multiple shooter object
- `family::CR3BPContinuationFamily`: CR3BP continuation family object
- `orbit::Int64`: Orbit identifier
"""
function getIndividualPeriodicOrbit(targeter::SpatialPerpJCMSTargeter, family::MBD.CR3BPContinuationFamily, orbit::Int64)
    numNodes::Int64 = length(family.nodes[orbit])
    period::Float64 = getPeriod(targeter, family.segments[orbit])
    orbitStates::Vector{Vector{Float64}} = append!([family.nodes[orbit][n].state.data[1:6] for n = 1:numNodes-1], [family.nodes[orbit][n].state.data[1:6].*[1, -1, 1, -1, 1, -1] for n = numNodes:-1:1])
    orbitTimes::Vector{Float64} = append!([family.nodes[orbit][n].epoch.data[1] for n = 1:numNodes-1], [period-family.nodes[orbit][n].epoch.data[1] for n = numNodes:-1:1]) 

    return MBD.CR3BPMSPeriodicOrbit(targeter.dynamicsModel, orbitStates, orbitTimes, period, getMonodromy(targeter, family.segments[orbit]))
end

"""
    getMonodromy(targeter, segments)

Return orbit monodromy matrix

# Arguments
- `targeter::SpatialPerpJCMSTargeter`: CR3BP spatial Jacobi constant perpendicular crossing multiple shooter object
- `segments::Vector{CR3BPSegment}`: CR3BP segment objects
"""
function getMonodromy(targeter::SpatialPerpJCMSTargeter, segments::Vector{MBD.CR3BPSegment})
    propagator = MBD.Propagator(equationType = MBD.STM)
    nStates::Int64 = getStateSize(targeter.dynamicsModel, MBD.STM)
    nSegs::Int64 = length(segments)
    renormalizeEvent::DifferentialEquations.DiscreteCallback = DifferentialEquations.PeriodicCallback(MBD.renormalize!, 0.1)
    STMs::Vector{Matrix{Float64}} = []
    for s::Int64 = 1:nSegs
        Rs::Vector{Matrix{Float64}} = []
        propSegment::MBD.CR3BPArc = propagateWithPeriodicEvent(propagator, renormalizeEvent, appendExtraInitialConditions(targeter.dynamicsModel, segments[s].originNode.state.data, MBD.STM), [0,segments[s].TOF.data[1]], targeter.dynamicsModel, [targeter.dynamicsModel, Rs])
        endState::StaticArrays.SVector{nStates, Float64} = StaticArrays.SVector{nStates, Float64}(getStateByIndex(propSegment, -1))
        STM::Matrix{Float64} = reshape(endState[7:42], (6,6))
        for R::Matrix{Float64} in reverse(Rs)
            STM *= R
        end
        push!(STMs, STM)
    end
    for s::Int64 = nSegs:-1:1
        Rs::Vector{Matrix{Float64}} = []
        propSegment::MBD.CR3BPArc = propagateWithPeriodicEvent(propagator, renormalizeEvent, appendExtraInitialConditions(targeter.dynamicsModel, segments[s].terminalNode.state.data.*[1, -1, 1, -1, 1, -1], MBD.STM), [0, segments[s].TOF.data[1]], targeter.dynamicsModel, [targeter.dynamicsModel, Rs])
        endState::StaticArrays.SVector{nStates, Float64} = StaticArrays.SVector{nStates, Float64}(getStateByIndex(propSegment, -1))
        STM::Matrix{Float64} = reshape(endState[7:42], (6,6))
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
    getMonodromy(targeter, solution)

Return orbit monodromy matrix

# Arguments
- `targeter::SpatialPerpJCMSTargeter`: CR3BP spatial Jacobi constant perpendicular crossing multiple shooter object
- `solution::CR3BPMultipleShooterProblem`: Solved CR3BP multiple shooter problem object
"""
function getMonodromy(targeter::SpatialPerpJCMSTargeter, solution::MBD.CR3BPMultipleShooterProblem)
    return getMonodromy(targeter, solution.segments)
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
    for s::Int64 = 1:length(segments)
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
    return getPeriod(targeter, solution.segments)
end

"""
    interpOrbit(targeter, fileName, paramName, paramValue, numNodes; choiceIndex, printProgress, tol, JTol)

Return interpolated periodic orbit object via bisection

# Arguments
- `targeter::SpatialPerpJCMSTargeter`: CR3BP spatial Jacobi constant perpendicular crossing multiple shooter object
- `fileName::String`: Family data CSV file
- `paramName::String`: Desired parameter name
- `paramValue::Float64`: Desired parameter value
- `numNodes::Int64`: Number of nodes
- `choiceIndex::Int64`: Desired orbit option index (default = 1)
- `printProgress::Bool`: Print progress? (default = false)
- `tol::Float64`: Convergence tolerance (default = 1E-11)
- `JTol::Float64`: Jacobian accuracy tolerance (default = 2E-3)
"""
function interpOrbit(targeter::SpatialPerpJCMSTargeter, fileName::String, paramName::String, paramValue::Float64, numNodes::Int64; choiceIndex::Int64 = 1, printProgress::Bool = false, tol::Float64 = 1E-11, JTol::Float64 = 2E-3)
    familyData::DataFrames.DataFrame = DataFrames.DataFrame(CSV.File(fileName))
    nMem::Int16 = Int16(size(familyData, 1))
    !any(occursin.(paramName, names(familyData))) && throw(ErrorException("Parameter not supported"))
    increasing::Bool = ((familyData[2,paramName]-familyData[1,paramName]) > 0) ? true : false
    switchIndices::Vector{Int16} = [1]
    for d::Int16 in Int16(3):nMem
        if (familyData[d,paramName]-familyData[d-1,paramName] > 0) != increasing
            increasing = !increasing
            push!(switchIndices, d-1)
        end
    end
    push!(switchIndices, nMem)
    solutionIndices::Vector{Int16} = []
    for s::Int16 in Int16(2):Int16(length(switchIndices))
        rangeData::DataFrames.DataFrame = familyData[switchIndices[s-1]:switchIndices[s],:]
        closestIndex::Int16 = Int16(argmin(abs.(rangeData[:,paramName].-paramValue)))
        (abs(rangeData[closestIndex,paramName]-paramValue) <= 1E-2) && push!(solutionIndices, closestIndex+switchIndices[s-1]-1)
    end
    isempty(solutionIndices) && throw(ErrorException("Family does not contain member with desired value"))
    nSol::Int16 = Int16(length(solutionIndices))
    index::Int16 = solutionIndices[1]
    if nSol > Int16(1)
        println("$nSol orbit options exist; choosing selected index...")
        index = solutionIndices[choiceIndex]
    end
    propagator = MBD.Propagator()
    if abs(familyData[index,paramName]-paramValue) <= 1E-8
        println("Orbit already exists in database!")
        orbitData::DataFrames.DataFrameRow = familyData[index,:]
        initialCondition::Vector{Float64} = [orbitData[p] for p in ["x", "y", "z", "xdot", "ydot", "zdot"]]
        period::Float64 = orbitData["Period"]
        JC::Float64 = orbitData["JC"]
        guessArc::MBD.CR3BPArc = propagate(propagator, initialCondition, [0, period/2], targeter.dynamicsModel)
        numStates::Int64 = getStateCount(guessArc)
        indices::Vector{Int64} = round.(Int64, range(1, numStates, numNodes))
        stateGuesses::Vector{Vector{Float64}} = [getStateByIndex(guessArc, indices[i]) for i = eachindex(indices)]
        timeGuesses::Vector{Float64} = [getTimeByIndex(guessArc, indices[i]) for i = eachindex(indices)]
        solution::MBD.CR3BPMultipleShooterProblem = correct(targeter, stateGuesses, timeGuesses, numNodes, JC; tol, JTol)
        period = getPeriod(targeter, solution)
        orbitStates::Vector{Vector{Float64}} = append!([solution.nodes[n].state.data[1:6] for n = 1:numNodes-1], [solution.nodes[n].state.data[1:6].*[1, -1, 1, -1, 1, -1] for n = numNodes:-1:1])
        orbitTimes::Vector{Float64} = append!([solution.nodes[n].epoch.data[1] for n = 1:numNodes-1], [period-solution.nodes[n].epoch.data[1] for n = numNodes:-1:1]) 
        orbit = MBD.CR3BPMSPeriodicOrbit(targeter.dynamicsModel, orbitStates, orbitTimes, period, getMonodromy(targeter, solution))
    else
        println("Orbit does not already exist in database; attempting bisection...")
        lowerBound::Int16 = index-Int16(1)
        upperBound::Int16 = index+Int16(1)
        (index == Int16(1)) && (lowerBound = 1)
        (index == nMem) && (upperBound = nMem)
        boundsData::DataFrames.DataFrame = DataFrames.sort(familyData[lowerBound:upperBound,:], paramName)
        lowerData::DataFrames.DataFrameRow = boundsData[1,:]
        upperData::DataFrames.DataFrameRow = boundsData[size(boundsData, 1),:]
        midInitialCondition::Vector{Float64} = [lowerData[p]+0.5*(upperData[p]-lowerData[p]) for p = ["x", "y", "z", "xdot", "ydot", "zdot"]]
        midPeriod::Float64 = lowerData["Period"]+0.5*(upperData["Period"]-lowerData["Period"])
        midJC::Float64 = lowerData["JC"]+0.5*(upperData["JC"]-lowerData["JC"])
        currentError::Float64 = abs(lowerData[paramName]-paramValue)
        iter::Int16 = Int16(1)
        while (currentError > 1E-8) && (iter <= 20)
            guessArc = propagate(propagator, midInitialCondition, [0, midPeriod/2], targeter.dynamicsModel)
            numStates = getStateCount(guessArc)
            indices = round.(Int64, range(1, numStates, numNodes))
            stateGuesses = [getStateByIndex(guessArc, indices[i]) for i = eachindex(indices)]
            timeGuesses = [getTimeByIndex(guessArc, indices[i]) for i = eachindex(indices)]
            solution = correct(targeter, stateGuesses, timeGuesses, numNodes, midJC; tol, JTol)
            newInitialCondition::Vector{Float64} = solution.nodes[1].state.data[1:6]
            newPeriod::Float64 = getPeriod(targeter, solution)
            newJC::Float64 = getJacobiConstant(targeter.dynamicsModel, newInitialCondition)
            newMonodromy::Matrix{Float64} = getMonodromy(targeter, solution)
            newEigenvalues::Vector{Complex{Float64}} = LinearAlgebra.eigvals(newMonodromy)
            newStability::Float64 = LinearAlgebra.norm(newEigenvalues, Inf)
            midData::DataFrames.DataFrameRow = DataFrames.DataFrame("x" => newInitialCondition[1], "y" => newInitialCondition[2], "z" => newInitialCondition[3], "xdot" => newInitialCondition[4], "ydot" => newInitialCondition[5], "zdot" => newInitialCondition[6], "Period" => newPeriod, "JC" => newJC, "Stability Index" => newStability)[1,:]
            currentValue::Float64 = midData[paramName]
            currentError = abs(currentValue-paramValue)
            printProgress && println("Current parameter value: "*paramName*" = $currentValue")
            (sign(currentValue-paramValue) == sign(lowerData[paramName]-paramValue)) ? (lowerData = midData) : (upperData = midData)
            midInitialCondition = [lowerData[p]+0.5*(upperData[p]-lowerData[p]) for p = ["x", "y", "z", "xdot", "ydot", "zdot"]]
            midPeriod = lowerData["Period"]+0.5*(upperData["Period"]-lowerData["Period"])
            midJC = lowerData["JC"]+0.5*(upperData["JC"]-lowerData["JC"])
            iter += Int16(1)
        end
        (iter > 50) && throw(ErrorException("Bisection failed"))
        period = getPeriod(targeter, solution)
        orbitStates = append!([solution.nodes[n].state.data[1:6] for n = 1:numNodes-1], [solution.nodes[n].state.data[1:6].*[1, -1, 1, -1, 1, -1] for n = numNodes:-1:1])
        orbitTimes = append!([solution.nodes[n].epoch.data[1] for n = 1:numNodes-1], [period-solution.nodes[n].epoch.data[1] for n = numNodes:-1:1]) 
        orbit = MBD.CR3BPMSPeriodicOrbit(targeter.dynamicsModel, orbitStates, orbitTimes, period, getMonodromy(targeter, solution))
    end

    return orbit
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
