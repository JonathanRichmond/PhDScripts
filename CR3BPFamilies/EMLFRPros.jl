"""
Script for Earth-Moon CR3BP lunar free return prograde orbit family

Author: Jonathan Richmond
C: 7/1/25
"""
module EMLFRPro
println()

using MBD, GLMakie, Logging

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../CR3BPTargeters/PlanarLFRX.jl")
include("../Utilities/Export.jl")
include("../Utilities/Plot.jl")

systemData = MBD.CR3BPSystemData("Earth", "Moon")
dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
Earth::MBD.BodyData, Moon::MBD.BodyData = systemData.primaryData[1], systemData.primaryData[2]
LagrangePoints::Vector{Vector{Float64}} = [getEquilibriumPoint(dynamicsModel, l) for l = 1:5]

propagator = MBD.Propagator()
targeter = PlanarLFRXTargeter(dynamicsModel)

h_E::Float64 = 185
h_M::Float64 = 30000
gamma::Float64 = 0.0
initialStateGuess::Vector{Float64} = getPrimaryState(dynamicsModel, 2)+[(Moon.bodyRadius+h_M)/getCharLength(systemData), 0, 0, 0, -1.0, 0]
tSpanGuess::Vector{Float64} = [0, -5*24*3600/getCharTime(systemData)]

numSegs::Int64 = 2
guessArc::MBD.CR3BPArc = propagate(propagator, initialStateGuess, tSpanGuess, dynamicsModel)
numStates::Int64 = getStateCount(guessArc)
numNodes::Int64 = numSegs+1
indices::Vector{Int64} = round.(Int64, range(1, numStates, numNodes))
stateGuesses::Vector{Vector{Float64}} = [getStateByIndex(guessArc, indices[i]) for i = eachindex(indices)]
timeGuesses::Vector{Float64} = [getTimeByIndex(guessArc, indices[i]) for i = eachindex(indices)]
solution1::MBD.CR3BPMultipleShooterProblem = correct(targeter, stateGuesses, timeGuesses, numNodes, initialStateGuess[1], h_E/getCharLength(systemData), gamma)
println("Converged Transfer 1:\n\tIC:\t$(solution1.nodes[1].state.data[1:6])\n\tTOF:\t$(getTOF(targeter, solution1))\n\tJC:\t$(getJacobiConstant(dynamicsModel, solution1.nodes[1].state.data[1:6]))\n")

solution2::MBD.CR3BPMultipleShooterProblem = correct(targeter, stateGuesses, timeGuesses, numNodes, initialStateGuess[1]-1E-5, h_E/getCharLength(systemData), gamma)
println("Converged Transfer 2:\n\tIC:\t$(solution2.nodes[1].state.data[1:6])\n\tTOF:\t$(getTOF(targeter, solution2))\n\tJC:\t$(getJacobiConstant(dynamicsModel, solution2.nodes[1].state.data[1:6]))\n")

println("Continuing transfers...")
continuationEngine = MBD.CR3BPNaturalParameterContinuationEngine(solution1, solution2, "Node 1 State", 1, -1E-5, -1E-2)
ydot0JumpCheck = MBD.BoundingBoxJumpCheck("Node 1 State", [NaN NaN; -2.5 -1.0])
addJumpCheck!(continuationEngine, ydot0JumpCheck)
MoonEndCheck = MBD.CR3BPPrimarySurfaceContinuationEndCheck(dynamicsModel, 2)
addEndCheck!(continuationEngine, MoonEndCheck)
solutions::MBD.CR3BPContinuationFamily = doContinuation!(continuationEngine, solution1, solution2)
lastTransfer::MBD.CR3BPMultipleShooterProblem = getIndividualSolution(targeter, solutions, getNumMembers(solutions))
println("Last Converged Transfer:\n\tIC:\t$(lastTransfer.nodes[1].state.data[1:6])\n\tTOF:\t$(getTOF(targeter, lastTransfer))\n\tJC:\t$(getJacobiConstant(dynamicsModel, lastTransfer.nodes[1].state.data[1:6]))\n")

# println("\nExporting family data...")
# fullExportCR3BPFamily(solutions, "FamilyData/CR3BPEMLFRPros.mat", "FamilyData/CR3BPEMLFRPros.csv")

# println("\nTesting interpolation...")
# testTransfer::MBD.CR3BPMultipleShooterProblem = interpSolution(targeter, "FamilyData/CR3BPEMLFRPros.csv", "x", 1.07, numNodes, h_E/getCharLength(systemData), gamma)
# println("Test Transfer:\n\tIC:\t$(testTransfer.nodes[1].state.data[1:6])\n\tP:\t$(getTOF(targeter, testTransfer))\n\tJC:\t$(getJacobiConstant(dynamicsModel, testTransfer.nodes[1].state.data[1:6]))\n")

# println("Plotting transfer...")
# plotTransfer::Int64 = getNumMembers(solutions)
# orbitArc::MBD.CR3BPArc = propagate(propagator, solutions.nodes[plotTransfer][1].state.data[1:6], [0, -getTOF(targeter, solutions.segments[plotTransfer])], dynamicsModel)
# xData::Vector{Float64}, yData::Vector{Float64} = Vector{Float64}(undef, getStateCount(orbitArc)), Vector{Float64}(undef, getStateCount(orbitArc))
# for s::Int64 in 1:getStateCount(orbitArc)
#     xData[s], yData[s] = orbitArc.states[s][1], orbitArc.states[s][2]
# end
# (figure, axis) = set2DPlotParameters(L"Earth-Moon Planar LFR ($JC=%$(round(getJacobiConstant(dynamicsModel, solutions.nodes[plotTransfer][1].state.data[1:6]); digits = 4))$)", L"$x$ [ndim]", L"$y$ [ndim]")
# GLMakie.lines!(axis, xData, yData, color = :white, label = L"\mathrm{Spatial\ Cycler\ Orbit}")
# GLMakie.scatter!(axis, getPrimaryState(dynamicsModel, 1)[1], getPrimaryState(dynamicsModel, 1)[2], color = :blue, markerspace = :data, markersize = 2*Earth.bodyRadius/getCharLength(systemData), label = L"\mathrm{Earth}" => (; markersize = 20))
# GLMakie.scatter!(axis, getPrimaryState(dynamicsModel, 2)[1], getPrimaryState(dynamicsModel, 2)[2], color = :gray, markerspace = :data, markersize = 2*Moon.bodyRadius/getCharLength(systemData), label = L"\mathrm{Moon}" => (; markersize = 20))
# GLMakie.Legend(figure[1,2], axis)

println()
end
