"""
Initial x-position targeter for CR3BP planar lunar free returns

Author: Jonathan Richmond
C: 7/1/25
"""

using MBD, CSV, DataFrames, DifferentialEquations, LinearAlgebra, StaticArrays

export PlanarPerpJCTargeter
export correct, getIndividualSolution, getTOF, interpSolution, propagateState

"""
    PlanarLFRXTargeter(dynamicsModel)

CR3BP planar lunar free return x-position targeter object

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
"""
struct PlanarLFRXTargeter
    dynamicsModel::MBD.CR3BPDynamicsModel                               # CR3BP dynamics model object

    function PlanarLFRXTargeter(dynamicsModel::MBD.CR3BPDynamicsModel)
        this = new(dynamicsModel)

        return this
    end
end

"""
    correct(targeter, qVector, tVector, numNodes, targetx0, targeth, targetgamma; tol, JTol)

Return corrected CR3BP multiple shooter problem object

# Arguments
- `targeter::PlanarLFRXTargeter`: CR3BP planar lunar free return initial x-position targeter object
- `qVector::Vector{Vector{Float64}}`: Node state guesses [ndim]
- `tVector::Vector{Float64}`: Node time guesses [ndim]
- `numNodes::Int64`: Number of nodes
- `targetx0::Float64`: Target initial x-position
- `targeth::Float64`: Target altitude [ndim]
- `targetgamma::Float64`: Target flight path angle [rad]
- `tol::Float64`: Convergence tolerance (default = 1E-11)
- `JTol::Float64`: Jacobian accuracy tolerance (default = 2E-3)
"""
function correct(targeter::PlanarLFRXTargeter, qVector::Vector{Vector{Float64}}, tVector::Vector{Float64}, numNodes::Int64, targetx0::Float64, targeth::Float64, targetgamma::Float64; tol::Float64 = 1E-11, JTol::Float64 = 2E-3)
    nodes::Vector{MBD.CR3BPNode} = []
    for n = 1:numNodes
        node = MBD.CR3BPNode(tVector[n], qVector[n], targeter.dynamicsModel)
        node.state.name = "Node "*string(n)*" State"
        node.epoch.name = "Node "*string(n)*" Epoch"
        push!(nodes, node)
    end    
    setFreeVariableMask!(nodes[1].state, [true, false, false, false, true, false])
    problem = MBD.CR3BPMultipleShooterProblem()
    segments::Vector{MBD.CR3BPSegment} = []
    for s = 1:numNodes-1
        segment = MBD.CR3BPSegment(tVector[s+1]-tVector[s], nodes[s], nodes[s+1])
        segment.TOF.name = "Segment "*string(s)*" TOF"
        push!(segments, segment)
    end
    map(s -> addSegment!(problem, s), segments)
    map(s -> addConstraint!(problem, MBD.CR3BPContinuityConstraint(s)), segments)
    addConstraint!(problem, MBD.CR3BPStateConstraint(nodes[1], [1], [targetx0]))
    addConstraint!(problem, MBD.CR3BPAltitudeConstraint(nodes[end], 1, targeth))
    addConstraint!(problem, MBD.CR3BPFlightPathAngleConstraint(nodes[end], 1, targetgamma))
    checkJacobian(problem, relTol = JTol)
    shooter = MBD.CR3BPMultipleShooter(tol)
    # shooter.printProgress = true
    solution::MBD.CR3BPMultipleShooterProblem = MBD.solve!(shooter, problem)
    map(s -> updateTerminalNodeEpoch!(s), solution.segments)

    return solution
end

"""
    getIndividualSolution(targeter, family, solution)

Return solution object

# Arguments
- `targeter::PlanarLFRXTargeter`: CR3BP planar lunar free return initial x-position targeter object
- `family::CR3BPContinuationFamily`: CR3BP continuation family object
- `solution::Int64`: Solution identifier
"""
function getIndividualSolution(targeter::PlanarLFRXTargeter, family::MBD.CR3BPContinuationFamily, solution::Int64)
    problem = MBD.CR3BPMultipleShooterProblem()
    numNodes::Int64 = length(family.nodes[solution])
    problem.nodes = [MBD.shallowClone(family.nodes[solution][n]) for n = 1:numNodes]
    problem.segments = [MBD.shallowClone(family.segments[solution][s]) for s = 1:numNodes-1]

    return problem
end

"""
    getTOF(targeter, segments)

Return trajectory time-of-flight

# Arguments
- `targeter::PlanarLFRXTargeter`: CR3BP planar lunar free return intial x-position targeter object
- `segments::Vector{CR3BPSegment}`: CR3BP segment objects
"""
function getTOF(targeter::PlanarLFRXTargeter, segments::Vector{MBD.CR3BPSegment})
    TOF::Float64 = 0.0
    for s::Int64 = 1:length(segments)
        TOF -= segments[s].TOF.data[1]
    end

    return TOF
end

"""
    getTOF(targeter, solution)

Return trajectory time-of-flight

# Arguments
- `targeter::PlanarLFRXTargeter`: CR3BP planar lunar free return initial x-position targeter object
- `solution::CR3BPMultipleShooterProblem`: Solved CR3BP multiple shooter problem object
"""
function getTOF(targeter::PlanarLFRXTargeter, solution::MBD.CR3BPMultipleShooterProblem)
    return getTOF(targeter, solution.segments)
end

"""
    interpSolution(targeter, fileName, paramName, paramValue, numNodes, targeth, targetgamma; choiceIndex, printProgress, tol, JTol)

Return interpolated trajectory solution object via bisection

# Arguments
- `targeter::PlanarLFRXTargeter`: CR3BP planar lunar free return initial x-position targeter object
- `fileName::String`: Family data CSV file
- `paramName::String`: Desired parameter name
- `paramValue::Float64`: Desired parameter value
- `numNodes::Int64`: Number of nodes
- `targeth::Float64`: Target altitude [ndim]
- `targetgamma::Float64`: Target flight path angle [rad]
- `choiceIndex::Int64`: Desired orbit option index (default = 1)
- `printProgress::Bool`: Print progress? (default = false)
- `tol::Float64`: Convergence tolerance (default = 1E-11)
- `JTol::Float64`: Jacobian accuracy tolerance (default = 2E-3)
"""
function interpSolution(targeter::PlanarLFRXTargeter, fileName::String, paramName::String, paramValue::Float64, numNodes::Int64, targeth::Float64, targetgamma::Float64; choiceIndex::Int64 = 1, printProgress::Bool = false, tol::Float64 = 1E-11, JTol::Float64 = 2E-3)
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
        println("$nSol transfer options exist; choosing selected index...")
        index = solutionIndices[choiceIndex]
    end
    propagator = MBD.Propagator()
    if abs(familyData[index,paramName]-paramValue) <= 1E-8
        println("Transfer already exists in database!")
        orbitData::DataFrames.DataFrameRow = familyData[index,:]
        initialCondition::Vector{Float64} = [orbitData[p] for p in ["x", "y", "z", "xdot", "ydot", "zdot"]]
        TOF::Float64 = orbitData["TOF"]
        guessArc::MBD.CR3BPArc = propagate(propagator, initialCondition, [0, TOF], targeter.dynamicsModel)
        numStates::Int64 = getStateCount(guessArc)
        indices::Vector{Int64} = round.(Int64, range(1, numStates, numNodes))
        stateGuesses::Vector{Vector{Float64}} = [getStateByIndex(guessArc, indices[i]) for i = eachindex(indices)]
        timeGuesses::Vector{Float64} = [getTimeByIndex(guessArc, indices[i]) for i = eachindex(indices)]
        solution::MBD.CR3BPMultipleShooterProblem = correct(targeter, stateGuesses, timeGuesses, numNodes, initialCondition[1], targeth, targetgamma; tol, JTol)
    else
        println("Transfer does not already exist in database; attempting bisection...")
        lowerBound::Int16 = index-Int16(1)
        upperBound::Int16 = index+Int16(1)
        (index == Int16(1)) && (lowerBound = 1)
        (index == nMem) && (upperBound = nMem)
        boundsData::DataFrames.DataFrame = DataFrames.sort(familyData[lowerBound:upperBound,:], paramName)
        lowerData::DataFrames.DataFrameRow = boundsData[1,:]
        upperData::DataFrames.DataFrameRow = boundsData[size(boundsData, 1),:]
        midInitialCondition::Vector{Float64} = [lowerData[p]+0.5*(upperData[p]-lowerData[p]) for p = ["x", "y", "z", "xdot", "ydot", "zdot"]]
        midTOF::Float64 = lowerData["TOF"]+0.5*(upperData["TOF"]-lowerData["TOF"])
        currentError::Float64 = abs(lowerData[paramName]-paramValue)
        iter::Int16 = Int16(1)
        while (currentError > 1E-8) && (iter <= 20)
            guessArc = propagate(propagator, midInitialCondition, [0, midTOF], targeter.dynamicsModel)
            numStates = getStateCount(guessArc)
            indices = round.(Int64, range(1, numStates, numNodes))
            stateGuesses = [getStateByIndex(guessArc, indices[i]) for i = eachindex(indices)]
            timeGuesses = [getTimeByIndex(guessArc, indices[i]) for i = eachindex(indices)]
            solution = correct(targeter, stateGuesses, timeGuesses, numNodes, midInitialCondition[1], targeth, targetgamma; tol, JTol)
            newInitialCondition::Vector{Float64} = solution.nodes[1].state.data[1:6]
            newTOF::Float64 = -getTOF(targeter, solution)
            newJC::Float64 = getJacobiConstant(targeter.dynamicsModel, newInitialCondition)
            midData::DataFrames.DataFrameRow = DataFrames.DataFrame("x" => newInitialCondition[1], "y" => newInitialCondition[2], "z" => newInitialCondition[3], "xdot" => newInitialCondition[4], "ydot" => newInitialCondition[5], "zdot" => newInitialCondition[6], "TOF" => newTOF, "JC" => newJC)[1,:]
            currentValue::Float64 = midData[paramName]
            currentError = abs(currentValue-paramValue)
            printProgress && println("Current parameter value: "*paramName*" = $currentValue")
            (sign(currentValue-paramValue) == sign(lowerData[paramName]-paramValue)) ? (lowerData = midData) : (upperData = midData)
            midInitialCondition = [lowerData[p]+0.5*(upperData[p]-lowerData[p]) for p = ["x", "y", "z", "xdot", "ydot", "zdot"]]
            midTOF = lowerData["TOF"]+0.5*(upperData["TOF"]-lowerData["TOF"])
            iter += Int16(1)
        end
        (iter > 50) && throw(ErrorException("Bisection failed"))
    end

    return solution
end

"""
    propagateState(targeter, q_simple, tSpan)

Return propagated state

# Arguments
- `targeter::PlanarLFRXTargeter`: CR3BP planar lunar free return initial x-position targeter object
- `q_simple::Vector{Float64}`: Simple state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
"""
function propagateState(targeter::PlanarLFRXTargeter, q_simple::Vector{Float64}, tSpan::Vector{Float64})
    propagator = MBD.Propagator()
    arc::MBD.CR3BPArc = propagate(propagator, q_simple, tSpan, targeter.dynamicsModel)

    return getStateByIndex(arc, -1)
end
