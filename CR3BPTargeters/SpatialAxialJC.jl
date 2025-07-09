"""
Jacobi constant axial crossing targeter for CR3BP spatial orbits

Author: Jonathan Richmond
C: 7/8/25
"""

using MBD, CSV, DataFrames, DifferentialEquations, LinearAlgebra, StaticArrays

export SpatialAxialJCTargeter
export correct, getIndividualPeriodicOrbit, getMonodromy, getPeriod, interpOrbit, propagateState

"""
    SpatialAxialJCTargeter(dynamicsModel)

CR3BP spatial axial crossing Jacobi constant targeter object

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
"""
struct SpatialAxialJCTargeter
    dynamicsModel::MBD.CR3BPDynamicsModel                               # CR3BP dynamics model object

    function SpatialAxialJCTargeter(dynamicsModel::MBD.CR3BPDynamicsModel)
        this = new(dynamicsModel)

        return this
    end
end

"""
    correct(targeter, q0, tSpan, targetJC; tol, JTol)

Return corrected CR3BP multiple shooter problem object

# Arguments
- `targeter::SpatialAxialJCTargeter`: CR3BP spatial axial crossing Jacobi constant targeter object
- `q0::Vector{Float64}`: Initial state guess [ndim]
- `tSpan::Vector{Float64}`: Time span guess [ndim]
- `targetJC::Float64`: Target Jacobi constant
- `tol::Float64`: Convergence tolerance (default = 1E-11)
- `JTol::Float64`: Jacobian accuracy tolerance (default = 2E-3)
"""
function correct(targeter::SpatialAxialJCTargeter, q0::Vector{Float64}, tSpan::Vector{Float64}, targetJC::Float64; tol::Float64 = 1E-11, JTol::Float64 = 2E-3)
    halfPeriodGuess::Float64 = (tSpan[2]-tSpan[1])/2
    tPCGuess::Float64 = tSpan[1]+halfPeriodGuess
    qPCGuess::Vector{Float64} = propagateState(targeter, q0, [tSpan[1], tPCGuess])
    originNode = MBD.CR3BPNode(tSpan[1], q0, targeter.dynamicsModel)
    originNode.state.name = "Initial State"
    setFreeVariableMask!(originNode.state, [true, false, false, false, true, true])
    terminalNode = MBD.CR3BPNode(tPCGuess, qPCGuess, targeter.dynamicsModel)
    terminalNode.state.name = "Target State"
    segment = MBD.CR3BPSegment(halfPeriodGuess, originNode, terminalNode)
    problem = MBD.CR3BPMultipleShooterProblem()
    addSegment!(problem, segment)
    addConstraint!(problem, MBD.CR3BPContinuityConstraint(segment))
    addConstraint!(problem, MBD.JacobiConstraint(originNode, targetJC))
    addConstraint!(problem, MBD.CR3BPStateConstraint(terminalNode, [2, 3, 4], [0.0, 0.0, 0.0]))
    checkJacobian(problem; relTol = JTol)
    shooter = MBD.CR3BPMultipleShooter(tol)
    # shooter.printProgress = true
    solution::MBD.CR3BPMultipleShooterProblem = MBD.solve!(shooter, problem)
    updateTerminalNodeEpoch!(segment)

    return solution
end

"""
    getIndividualPeriodicOrbit(targeter, family, orbit)

Return periodic orbit object

# Arguments
- `targeter::SpatialAxialJCTargeter`: CR3BP spatial axial crossing Jacobi constant targeter object
- `family::CR3BPContinuationFamily`: CR3BP continuation family object
- `orbit::Int64`: Orbit identifier
"""
function getIndividualPeriodicOrbit(targeter::SpatialAxialJCTargeter, family::MBD.CR3BPContinuationFamily, orbit::Int64)
    period::Float64 = 2*family.segments[orbit][1].TOF.data[1]
    propagator = MBD.Propagator(equationType = MBD.STM)
    nStates::Int64 = getStateSize(targeter.dynamicsModel, MBD.STM)
    halfOrbit::MBD.CR3BPArc = propagate(propagator, appendExtraInitialConditions(targeter.dynamicsModel, family.nodes[orbit][1].state.data, MBD.STM), [0, family.segments[orbit][1].TOF.data[1]], targeter.dynamicsModel)
    PCState::StaticArrays.SVector{nStates, Float64} = StaticArrays.SVector{nStates, Float64}(getStateByIndex(halfOrbit, -1))
    PCSTM::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}(reshape(PCState[7:42], (6,6)))
    G::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}(1.0, 0, 0, 0, 0, 0, 0, -1.0, 0, 0, 0, 0, 0, 0, -1.0, 0, 0, 0, 0, 0, 0, -1.0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 1.0)
    Omega::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}(0, -1.0, 0, 1.0, 0, 0, 0, 0, 0)
    A::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([zeros(Float64, (3,3)) -1 .*LinearAlgebra.I; LinearAlgebra.I -2 .*Omega])
    B::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([-2 .*Omega LinearAlgebra.I; -1 .*LinearAlgebra.I zeros(Float64, (3,3))])
    monodromy::Matrix{Float64} = G*A*(PCSTM')*B*G*PCSTM

    return MBD.CR3BPPeriodicOrbit(targeter.dynamicsModel, getStateByIndex(halfOrbit, 1)[1:6], period, monodromy)
end

"""
    getMonodromy(targeter, solution)

Return orbit monodromy matrix

# Arguments
- `targeter::SpatialAxialJCTargeter`: CR3BP spatial axial crossing Jacobi constant targeter object
- `solution::CR3BPMultipleShooterProblem`: Solved CR3BP multiple shooter problem object
"""
function getMonodromy(targeter::SpatialAxialJCTargeter, solution::MBD.CR3BPMultipleShooterProblem)
    propagator = MBD.Propagator(equationType = MBD.STM)
    n_simple::Int64 = getStateSize(targeter.dynamicsModel, MBD.SIMPLE)
    n_STM::Int64 = getStateSize(targeter.dynamicsModel, MBD.STM)
    renormalizeEvent::DifferentialEquations.DiscreteCallback = DifferentialEquations.PeriodicCallback(MBD.renormalize!, pi/10)
    Rs::Vector{Matrix{Float64}} = []
    halfOrbit::MBD.CR3BPArc = propagateWithPeriodicEvent(propagator, renormalizeEvent, appendExtraInitialConditions(targeter.dynamicsModel, solution.nodes[1].state.data, MBD.STM), [0, solution.segments[1].TOF.data[1]], targeter.dynamicsModel, [targeter.dynamicsModel, Rs])
    PCState::StaticArrays.SVector{n_STM, Float64} = StaticArrays.SVector{n_STM, Float64}(getStateByIndex(halfOrbit, -1))
    PCSTM::Matrix{Float64} = reshape(PCState[n_simple+1:n_STM], (n_simple,n_simple))
    for R::Matrix{Float64} in reverse(Rs)
        PCSTM *= R
    end
    G::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}(1.0, 0, 0, 0, 0, 0, 0, -1.0, 0, 0, 0, 0, 0, 0, -1.0, 0, 0, 0, 0, 0, 0, -1.0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 1.0)
    Omega::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}(0, -1.0, 0, 1.0, 0, 0, 0, 0, 0)
    A::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([zeros(Float64, (3,3)) -1 .*LinearAlgebra.I; LinearAlgebra.I -2 .*Omega])
    B::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([-2 .*Omega LinearAlgebra.I; -1 .*LinearAlgebra.I zeros(Float64, (3,3))])

    return G*A*(PCSTM')*B*G*PCSTM
end

"""
    getPeriod(targeter, solution)

Return orbit period

# Arguments
- `targeter::SpatialAxialJCTargeter`: CR3BP spatial axial crossing Jacobi constant targeter object
- `solution::CR3BPMultipleShooterProblem`: Solved CR3BP multiple shooter problem object
"""
function getPeriod(targeter::SpatialAxialJCTargeter, solution::MBD.CR3BPMultipleShooterProblem)
    return 2*solution.segments[1].TOF.data[1]
end

"""
    interpOrbit(targeter, fileName, paramName, paramValue; choiceIndex, printProgress, tol, JTol)

Return interpolated periodic orbit object via bisection

# Arguments
- `targeter::SpatialAxialJCTargeter`: CR3BP spatial axial crossing Jacobi constant targeter object
- `fileName::String`: Family data CSV file
- `paramName::String`: Desired parameter name
- `paramValue::Float64`: Desired parameter value
- `choiceIndex::Int64`: Desired orbit option index (default = 1)
- `printProgress::Bool`: Print progress? (default = false)
- `tol::Float64`: Convergence tolerance (default = 1E-11)
- `JTol::Float64`: Jacobian accuracy tolerance (default = 2E-3)
"""
function interpOrbit(targeter::SpatialAxialJCTargeter, fileName::String, paramName::String, paramValue::Float64; choiceIndex::Int64 = 1, printProgress::Bool = false, tol::Float64 = 1E-11, JTol::Float64 = 2E-3)
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
    propagator = MBD.Propagator(equationType = MBD.STM)
    nStates::Int64 = getStateSize(targeter.dynamicsModel, MBD.STM)
    G::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}(1.0, 0, 0, 0, 0, 0, 0, -1.0, 0, 0, 0, 0, 0, 0, -1.0, 0, 0, 0, 0, 0, 0, -1.0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 1.0)
    Omega::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}(0, -1.0, 0, 1.0, 0, 0, 0, 0, 0)
    A::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([zeros(Float64, (3,3)) -1 .*LinearAlgebra.I; LinearAlgebra.I -2 .*Omega])
    B::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([-2 .*Omega LinearAlgebra.I; -1 .*LinearAlgebra.I zeros(Float64, (3,3))])
    if abs(familyData[index,paramName]-paramValue) <= 1E-8
        println("Orbit already exists in database!")
        orbitData::DataFrames.DataFrameRow = familyData[index,:]
        initialCondition::Vector{Float64} = [orbitData[p] for p in ["x", "y", "z", "xdot", "ydot", "zdot"]]
        period::Float64 = orbitData["Period"]
        halfOrbit::MBD.CR3BPArc = propagate(propagator, appendExtraInitialConditions(targeter.dynamicsModel, initialCondition, MBD.STM), [0, period/2], targeter.dynamicsModel)
        PCState::StaticArrays.SVector{nStates, Float64} = StaticArrays.SVector{nStates, Float64}(getStateByIndex(halfOrbit, -1))
        PCSTM::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}(reshape(PCState[7:42], (6,6)))
        monodromy::Matrix{Float64} = G*A*(PCSTM')*B*G*PCSTM
        orbit = MBD.CR3BPPeriodicOrbit(targeter.dynamicsModel, initialCondition, period, monodromy)
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
            midSolution::MBD.CR3BPMultipleShooterProblem = correct(targeter, midInitialCondition, [0, midPeriod], midJC; tol, JTol)
            newInitialCondition::Vector{Float64} = midSolution.nodes[1].state.data[1:6]
            newPeriod::Float64 = getPeriod(targeter, midSolution)
            newJC::Float64 = getJacobiConstant(targeter.dynamicsModel, newInitialCondition)
            newMonodromy::Matrix{Float64} = getMonodromy(targeter, midSolution)
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
        halfOrbit = propagate(propagator, appendExtraInitialConditions(targeter.dynamicsModel, midInitialCondition, MBD.STM), [0, midPeriod/2], targeter.dynamicsModel)
        PCState = StaticArrays.SVector{nStates, Float64}(getStateByIndex(halfOrbit, -1))
        PCSTM = StaticArrays.SMatrix{6, 6, Float64}(reshape(PCState[7:42], (6,6)))
        monodromy = G*A*(PCSTM')*B*G*PCSTM
        orbit = MBD.CR3BPPeriodicOrbit(targeter.dynamicsModel, midInitialCondition, midPeriod, monodromy)
    end

    return orbit
end

"""
    propagateState(targeter, q_simple, tSpan)

Return propagated state

# Arguments
- `targeter::SpatialAxialJCTargeter`: CR3BP spatial axial crossing Jacobi constant targeter object
- `q_simple::Vector{Float64}`: Simple state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
"""
function propagateState(targeter::SpatialAxialJCTargeter, q_simple::Vector{Float64}, tSpan::Vector{Float64})
    propagator = MBD.Propagator()
    arc::MBD.CR3BPArc = propagate(propagator, q_simple, tSpan, targeter.dynamicsModel)

    return getStateByIndex(arc, -1)
end
