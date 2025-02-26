"""
Script for Earth-Moon CR3BP spatial cycler orbit family
    *Initial condition provided by Gomez*

Author: Jonathan Richmond
C: 2/26/25
"""
module EMSCycler
println()

using MBD, MATLAB

include("../Targeters/SpatialMSJC.jl")
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

targeter = SpatialMSJCTargeter(dynamicsModel)

initialStateGuess::Vector{Float64} = [0.8190590215962604, 0, 0, 0, 0.19936048867848094, 0]
tSpanGuess::Vector{Float64} = [0, 6.526254920424554]
targetJC::Float64 = getJacobiConstant(dynamicsModel, initialStateGuess)
solution1::MBD.CR3BPMultipleShooterProblem = correct(targeter, initialStateGuess, tSpanGuess, 2, targetJC, 1E-10)
println("Converged Orbit 1:\n\tState:$(solution1.nodes[1].state.data[1:6])\n\tPeriod: $(getPeriod(targeter, solution1))\n\tJC: $(getJacobiConstant(dynamicsModel, solution1.nodes[1].state.data[1:6]))")

# solution2::MBD.CR3BPMultipleShooterProblem = correct(targeter, initialStateGuess+[0, 0, 0.0001, 0, 0, 0], tSpanGuess, 2, targetJC-1E-5, 1E-10)
# println("\nConverged Orbit 2:\n\tState:$(solution2.nodes[1].state.data[1:6])\n\tPeriod: $(getPeriod(targeter, solution2))\n\tJC: $(getJacobiConstant(dynamicsModel, solution2.nodes[1].state.data[1:6]))")

propagator = MBD.Propagator()
arc::MBD.CR3BPArc = propagate(propagator, solution1.nodes[1].state.data[1:6], [0, getPeriod(targeter, solution1)], dynamicsModel)
nStates::Int64 = getStateCount(arc)
x::Vector{Float64} = zeros(Float64, nStates)
y::Vector{Float64} = zeros(Float64, nStates)
z::Vector{Float64} = zeros(Float64, nStates)
xdot::Vector{Float64} = zeros(Float64, nStates)
ydot::Vector{Float64} = zeros(Float64, nStates)
zdot::Vector{Float64} = zeros(Float64, nStates)
t::Vector{Float64} = zeros(Float64, nStates)
for s::Int64 in 1:nStates
    state::Vector{Float64} = getStateByIndex(arc, s)
    x[s] = state[1]
    y[s] = state[2]
    z[s] = state[3]
    xdot[s] = state[4]
    ydot[s] = state[5]
    zdot[s] = state[6]
    t[s] = getTimeByIndex(arc, s)
end

mf = MATLAB.MatFile("Output/CR3BPTraj.mat", "w")
exportCR3BPTrajectory(x, y, z, xdot, ydot, zdot, t, mf, :trajCR3BP)
MATLAB.close(mf)

println()
end
