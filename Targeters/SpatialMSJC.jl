"""
Jacobi constant multiple shooter for spatial orbits

Author: Jonathan Richmond
C: 2/26/25
"""

using MBD, CSV, DataFrames, LinearAlgebra, StaticArrays

export SpatialMSJCTargeter
export correct, propagateState, getPeriod#, getIndividualPeriodicOrbit, getMonodromy, interpOrbit

"""
    SpatialMSJCTargeter(dynamicsModel)

CR3BP spatial Jacobi constant multiple shooter object

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
"""
struct SpatialMSJCTargeter <: MBD.AbstractTargeter
    dynamicsModel::MBD.CR3BPDynamicsModel                               # CR3BP dynamics model object

    function SpatialMSJCTargeter(dynamicsModel::MBD.CR3BPDynamicsModel)
        this = new(dynamicsModel)

        return this
    end
end

"""
    correct(targeter, q0, tSpan, numNodes, targetJC; tol)

Return corrected CR3BP multiple shooter problem object

# Arguments
- `targeter::SpatialMSTargeter`: CR3BP spatial Jacobi constant multiple shooter object
- `q0::Vector{Float64}`: Initial state guess [ndim]
- `tSpan::Vector{Float64}`: Time span guess [ndim]
- `numNodes::Int64`: Number of nodes
- `targetJC::Float64`: Target Jacobi constant
- `tol::Float64`: Convergence tolerance (default = 1E-11)
"""
function correct(targeter::SpatialMSJCTargeter, q0::Vector{Float64}, tSpan::Vector{Float64}, numNodes::Int64, targetJC::Float64, tol::Float64 = 1E-11)
    timeGuesses::Vector{Float64} = collect(LinRange(tSpan[1], tSpan[2], numNodes+1))
    stateGuesses::Vector{Vector{Float64}} = [q0]
    [push!(stateGuesses, propagateState(targeter, q0, [0, timeGuesses[n]])) for n = 2:(numNodes+1)]
    nodes::Vector{MBD.CR3BPNode} = [MBD.CR3BPNode(timeGuesses[n], stateGuesses[n], dynamicsModel) for n = 1:(numNodes+1)]
    problem = MBD.CR3BPMultipleShooterProblem()
    segments::Vector{MBD.CR3BPSegment} = [MBD.CR3BPSegment(timeGuesses[n+1]-timeGuesses[n], nodes[n], nodes[n+1]) for n = 1:numNodes]
    map(s -> addSegment!(problem, s), segments)
    map(s -> addConstraint!(problem, MBD.CR3BPContinuityConstraint(s)), segments)
    addConstraint!(problem, MBD.StateMatchConstraint(nodes[1].state, nodes[end].state, [1, 2, 3, 5, 6]))
    addConstraint!(problem, MBD.JacobiConstraint(nodes[1], targetJC))
    println(getNumFreeVariables!(problem))
    println(getNumConstraints(problem))
    checkJacobian(problem)
    shooter = MBD.CR3BPMultipleShooter(tol)
    shooter.printProgress = true
    solution::MBD.CR3BPMultipleShooterProblem = MBD.solve!(shooter, problem)

    return solution
end

"""
    getPeriod(targeter, solution)

Return orbit period

# Arguments
- `targeter::SpatialMSJCTargeter`: CR3BP spatial Jacobi constant multiple shooter object
- `solution::CR3BPMultipleShooterProblem`: Solved CR3BP multiple shooter problem object
"""
function getPeriod(targeter::SpatialMSJCTargeter, solution::MBD.CR3BPMultipleShooterProblem)
    TOF::Float64 = 0.0
    for s = 1:length(solution.segments)
        TOF += solution.segments[s].TOF.data[1]
    end

    return TOF
end

"""
    propagateState(targeter, q_simple, tSpan)

Return propagated state

# Arguments
- `targeter::SpatialMSJCTargeter`: CR3BP spatial Jacobi constant multiple shooter object
- `q_simple::Vector{Float64}`: Simple state vector [ndim]
- `tSpan::Vector{Float64}`: Time span [ndim]
"""
function propagateState(targeter::SpatialMSJCTargeter, q_simple::Vector{Float64}, tSpan::Vector{Float64})
    propagator = MBD.Propagator()
    arc::MBD.CR3BPArc = propagate(propagator, q_simple, tSpan, targeter.dynamicsModel)

    return getStateByIndex(arc, -1)
end
