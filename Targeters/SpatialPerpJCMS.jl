"""
Jacobi constant perpendicular crossing multiple shooter for spatial orbits

Author: Jonathan Richmond
C: 2/26/25
U: 3/6/25
"""

using MBD, CSV, DataFrames, LinearAlgebra, StaticArrays

export SpatialPerpJCMSTargeter
export correct, propagateState, getPeriod#, getIndividualPeriodicOrbit, getMonodromy, interpOrbit

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
    nodes::Vector{MBD.CR3BPNode} = [MBD.CR3BPNode(timeGuesses[n], stateGuesses[n], dynamicsModel) for n = 1:numNodes]
    setFreeVariableMask!(nodes[1].state, [true, false, true, false, true, false])
    problem = MBD.CR3BPMultipleShooterProblem()
    segments::Vector{MBD.CR3BPSegment} = [MBD.CR3BPSegment(timeGuesses[n+1]-timeGuesses[n], nodes[n], nodes[n+1]) for n = 1:(numNodes-1)]
    map(s -> setFreeVariableMask!(s.TOF, [false]), segments[1:(end-1)])
    map(s -> addSegment!(problem, s), segments)
    map(s -> addConstraint!(problem, MBD.CR3BPContinuityConstraint(s)), segments)
    addConstraint!(problem, MBD.CR3BPStateConstraint(nodes[end], [2, 4, 6], [0.0, 0.0, 0.0]))
    addConstraint!(problem, MBD.JacobiConstraint(nodes[1], targetJC))
    # println("Free variables: $(getNumFreeVariables!(problem))")
    # println("Constraints: $(getNumConstraints(problem))")
    checkJacobian(problem)
    shooter = MBD.CR3BPMultipleShooter(tol)
    # shooter.printProgress = true
    shooter.maxIterations = 50
    solution::MBD.CR3BPMultipleShooterProblem = MBD.solve!(shooter, problem)

    return solution
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
