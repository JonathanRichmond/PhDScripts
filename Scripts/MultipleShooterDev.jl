"""
Script for multiple shooter code development

Author: Jonathan Richmond
C: 2/26/25
U: 4/2/25
"""
module MSDev
println()

using MBD, MATLAB

include("../Targeters/SpatialPerpJCMS.jl")
include("../Utilities/Export.jl")

systemData = MBD.CR3BPSystemData("Earth", "Moon")
dynamicsModel = MBD.CR3BPDynamicsModel(systemData)

Earth = systemData.primaryData[1]
Moon = systemData.primaryData[2]
L1::Vector{Float64} = getEquilibriumPoint(dynamicsModel, 1)
L2::Vector{Float64} = getEquilibriumPoint(dynamicsModel, 2)
L3::Vector{Float64} = getEquilibriumPoint(dynamicsModel, 3)
L4::Vector{Float64} = getEquilibriumPoint(dynamicsModel, 4)
L5::Vector{Float64} = getEquilibriumPoint(dynamicsModel, 5)

targeter = SpatialPerpJCMSTargeter(dynamicsModel)

initialStateGuess::Vector{Float64} = [0.850522, 0, -0.177117, 0, 0.255761, 0] # [0.881995, 0, -0.24036, 0, 0.144858, 0]
tSpanGuess::Vector{Float64} = [0, 18.0561] # [0, 20.3756]
targetJC::Float64 = 3.0098 # 2.9833
numSegs::Int64 = 28
solution1::MBD.CR3BPMultipleShooterProblem = correct(targeter, initialStateGuess, tSpanGuess, numSegs, targetJC, 1E-10)
println("Converged Orbit 1:\n\tState:$(solution1.nodes[1].state.data[1:6])\n\tPeriod: $(getPeriod(targeter, solution1))\n\tJC: $(getJacobiConstant(dynamicsModel, solution1.nodes[1].state.data[1:6]))")

solution2::MBD.CR3BPMultipleShooterProblem = correct(targeter, initialStateGuess, tSpanGuess, numSegs, targetJC-1E-5, 1E-10)
println("\nConverged Orbit 2:\n\tState:$(solution2.nodes[1].state.data[1:6])\n\tPeriod: $(getPeriod(targeter, solution2))\n\tJC: $(getJacobiConstant(dynamicsModel, solution2.nodes[1].state.data[1:6]))")

engine = MBD.CR3BPMultipleShooterContinuationEngine(solution1, solution2, "Jacobi Constant", 1, -5E-6, -1E-4)
q0JumpCheck = MBD.BoundingBoxJumpCheck("Node 1 State", [0.85 0.91; -0.5 -0.17; 0.08 0.26])
addJumpCheck!(engine, q0JumpCheck)
numStepsEndCheck = MBD.NumberStepsContinuationEndCheck(2500)
addEndCheck!(engine, numStepsEndCheck)

println()
solutions::MBD.CR3BPContinuationFamily = doContinuation!(targeter, engine, solution1, solution2, numSegs, 1E-10)
println("\nLast Converged Orbit:\n\tState:$(solutions.nodes[end][1].state.data[1:6])\n\tPeriod: $(getPeriod(targeter, solutions.segments[end]))\n\tJC: $(getJacobiConstant(dynamicsModel, solutions.nodes[end][1].state.data[1:6]))")

exportSolution::Int64 = length(solutions.nodes)
propagator = MBD.Propagator()
mf = MATLAB.MatFile("Output/CR3BPTraj.mat", "w")
for s::Int64 = 1:numSegs/2
    elapsedT::Float64 = 0.0
    if s > 1
        for j::Int64 = 2:s
            elapsedT += solutions.segments[exportSolution][s-1].TOF.data[1]
        end
    end
    arc::MBD.CR3BPArc = propagate(propagator, solutions.nodes[exportSolution][s].state.data[1:6], [elapsedT, elapsedT+solutions.segments[exportSolution][s].TOF.data[1]], dynamicsModel)
    nStates::Int64 = getStateCount(arc)
    x::Vector{Float64} = zeros(Float64, nStates)
    y::Vector{Float64} = zeros(Float64, nStates)
    z::Vector{Float64} = zeros(Float64, nStates)
    xdot::Vector{Float64} = zeros(Float64, nStates)
    ydot::Vector{Float64} = zeros(Float64, nStates)
    zdot::Vector{Float64} = zeros(Float64, nStates)
    t::Vector{Float64} = zeros(Float64, nStates)
    for j::Int64 in 1:nStates
        state::Vector{Float64} = getStateByIndex(arc, j)
        x[j] = state[1]
        y[j] = state[2]
        z[j] = state[3]
        xdot[j] = state[4]
        ydot[j] = state[5]
        zdot[j] = state[6]
        t[j] = getTimeByIndex(arc, j)
    end
    exportCR3BPTrajectory(x, y, z, xdot, ydot, zdot, t, mf, Symbol("trajCR3BP"*string(s)))
end
for s::Int64 = 1:numSegs/2
    elapsedT::Float64 = 0.0
    for j = 1:length(solutions.segments[exportSolution])
        elapsedT += solutions.segments[exportSolution][j].TOF.data[1]
    end
    if s > 1
        for j::Int64 = 2:s
            elapsedT += solutions.segments[exportSolution][end+2-s].TOF.data[1]
        end
    end
    arc::MBD.CR3BPArc = propagate(propagator, solutions.nodes[exportSolution][end-s].state.data[1:6].*[1, -1, 1, -1, 1, -1], [elapsedT+solutions.segments[exportSolution][end+1-s].TOF.data[1], elapsedT], dynamicsModel)
    nStates::Int64 = getStateCount(arc)
    x::Vector{Float64} = zeros(Float64, nStates)
    y::Vector{Float64} = zeros(Float64, nStates)
    z::Vector{Float64} = zeros(Float64, nStates)
    xdot::Vector{Float64} = zeros(Float64, nStates)
    ydot::Vector{Float64} = zeros(Float64, nStates)
    zdot::Vector{Float64} = zeros(Float64, nStates)
    t::Vector{Float64} = zeros(Float64, nStates)
    for j::Int64 in 1:nStates
        state::Vector{Float64} = getStateByIndex(arc, j)
        x[j] = state[1]
        y[j] = state[2]
        z[j] = state[3]
        xdot[j] = state[4]
        ydot[j] = state[5]
        zdot[j] = state[6]
        t[j] = getTimeByIndex(arc, j)
    end
    exportCR3BPTrajectory(x, y, z, xdot, ydot, zdot, t, mf, Symbol("trajCR3BP"*string(s+numSegs)))
end
MATLAB.close(mf)

println()
end
