"""
Script for Earth-Moon CR3BP spatial cycler orbit family
    *Initial condition provided by Gomez*

Author: Jonathan Richmond
C: 2/26/25
U: 3/4/25
"""
module EMSCycler
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

initialStateGuess::Vector{Float64} = [0.8327, 0, 0.08783, 0, 0.1718, 0]
tSpanGuess::Vector{Float64} = [0, 6.7911]
targetJC::Float64 = 3.12618
numNodes::Int64 = 6
solution1::MBD.CR3BPMultipleShooterProblem = correct(targeter, initialStateGuess, tSpanGuess, numNodes, targetJC)
println("Converged Orbit 1:\n\tState:$(solution1.nodes[1].state.data[1:6])\n\tPeriod: $(getPeriod(targeter, solution1))\n\tJC: $(getJacobiConstant(dynamicsModel, solution1.nodes[1].state.data[1:6]))")

# solution2::MBD.CR3BPMultipleShooterProblem = correct(targeter, initialStateGuess+[0, 0, 0.0001, 0, 0, 0], tSpanGuess, 2, targetJC-1E-5, 1E-10)
# println("\nConverged Orbit 2:\n\tState:$(solution2.nodes[1].state.data[1:6])\n\tPeriod: $(getPeriod(targeter, solution2))\n\tJC: $(getJacobiConstant(dynamicsModel, solution2.nodes[1].state.data[1:6]))")

exportSolution::MBD.CR3BPMultipleShooterProblem = solution1
propagator = MBD.Propagator()
mf = MATLAB.MatFile("Output/CR3BPTraj.mat", "w")
for s::Int64 = 1:numNodes/2
    elapsedT::Float64 = 0.0
    if s > 1
        for j::Int64 = 2:s
            elapsedT += exportSolution.segments[s-1].TOF.data[1]
        end
    end
    arc::MBD.CR3BPArc = propagate(propagator, exportSolution.nodes[s].state.data[1:6], [elapsedT, elapsedT+exportSolution.segments[s].TOF.data[1]], dynamicsModel)
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
for s::Int64 = 1:numNodes/2
    elapsedT::Float64 = 0.0
    for j = 1:length(exportSolution.segments)
        elapsedT += exportSolution.segments[j].TOF.data[1]
    end
    if s > 1
        for j::Int64 = 2:s
            elapsedT += exportSolution.segments[end+2-s].TOF.data[1]
        end
    end
    arc::MBD.CR3BPArc = propagate(propagator, exportSolution.nodes[end-s].state.data[1:6].*[1, -1, 1, -1, 1, -1], [elapsedT+exportSolution.segments[end+1-s].TOF.data[1], elapsedT], dynamicsModel)
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
    exportCR3BPTrajectory(x, y, z, xdot, ydot, zdot, t, mf, Symbol("trajCR3BP"*string(s+numNodes)))
end
MATLAB.close(mf)

println()
end
