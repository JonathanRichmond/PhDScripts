"""
Script for Earth-Moon CR3BP spatial lunar free return prograde orbit family

Author: Jonathan Richmond
C: 7/21/25
U: 7/23/25
"""
module EMLFRProSpat
println()

using MBD, GLMakie, Logging, MATLAB

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../CR3BPTargeters/SpatialLFRZ.jl")
include("../Utilities/Export.jl")
include("../Utilities/Plot.jl")

systemData = MBD.CR3BPSystemData("Earth", "Moon")
dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
Earth::MBD.BodyData, Moon::MBD.BodyData = systemData.primaryData[1], systemData.primaryData[2]
LagrangePoints::Vector{Vector{Float64}} = [getEquilibriumPoint(dynamicsModel, l) for l = 1:5]

propagator = MBD.Propagator()
targeter = SpatialLFRZTargeter(dynamicsModel)

h_E::Float64 = 185
h_M::Float64 = 16000
gamma::Float64 = 0.0
initialStateGuess::Vector{Float64} = getPrimaryState(dynamicsModel, 2)+[(Moon.bodyRadius+h_M)/getCharLength(systemData), 0, 0, 0, -1.13, 0]
tSpanGuess::Vector{Float64} = [0, -4.35*24*3600/getCharTime(systemData)]

numSegs::Int64 = 2
guessArc::MBD.CR3BPArc = propagate(propagator, initialStateGuess, tSpanGuess, dynamicsModel)
numStates::Int64 = getStateCount(guessArc)
numNodes::Int64 = numSegs+1
indices::Vector{Int64} = round.(Int64, range(1, numStates, numNodes))
stateGuesses::Vector{Vector{Float64}} = [getStateByIndex(guessArc, indices[i]) for i = eachindex(indices)]
timeGuesses::Vector{Float64} = [getTimeByIndex(guessArc, indices[i]) for i = eachindex(indices)]
solution1::MBD.CR3BPMultipleShooterProblem = correct(targeter, stateGuesses, timeGuesses, numNodes, initialStateGuess[6], h_M/getCharLength(systemData), h_E/getCharLength(systemData), gamma)
println("Converged Transfer 1:\n\tIC:\t$(solution1.nodes[1].state.data[1:6])\n\tTOF:\t$(getTOF(targeter, solution1))\n\tJC:\t$(getJacobiConstant(dynamicsModel, solution1.nodes[1].state.data[1:6]))\n")

solution2::MBD.CR3BPMultipleShooterProblem = correct(targeter, stateGuesses, timeGuesses, numNodes, initialStateGuess[6]+1E-5, h_M/getCharLength(systemData), h_E/getCharLength(systemData), gamma)
println("Converged Transfer 2:\n\tIC:\t$(solution2.nodes[1].state.data[1:6])\n\tTOF:\t$(getTOF(targeter, solution2))\n\tJC:\t$(getJacobiConstant(dynamicsModel, solution2.nodes[1].state.data[1:6]))\n")

println("Continuing transfers...")
continuationEngine = MBD.CR3BPNaturalParameterContinuationEngine(solution1, solution2, "Node 1 State", 2, 1E-5, 1E-2)
ydot0JumpCheck = MBD.BoundingBoxJumpCheck("Node 1 State", [NaN NaN; NaN NaN; -1.5 -1.0])
addJumpCheck!(continuationEngine, ydot0JumpCheck)
z0EndCheck = MBD.BoundingBoxContinuationEndCheck("Node 1 State", [NaN NaN; 0 sin(pi/3)*(Moon.bodyRadius+h_M)/getCharLength(systemData); NaN NaN])
addEndCheck!(continuationEngine, z0EndCheck)
solutions::MBD.CR3BPContinuationFamily = doContinuation!(continuationEngine, solution1, solution2)
lastTransfer::MBD.CR3BPMultipleShooterProblem = getIndividualSolution(targeter, solutions, getNumMembers(solutions))
println("Last Converged Transfer:\n\tIC:\t$(lastTransfer.nodes[1].state.data[1:6])\n\tTOF:\t$(getTOF(targeter, lastTransfer))\n\tJC:\t$(getJacobiConstant(dynamicsModel, lastTransfer.nodes[1].state.data[1:6]))\n")

# println("\nExporting family data...")
# fullExportCR3BPFamily(solutions, "FamilyData/CR3BPEMLFRProsSpatial.mat", "FamilyData/CR3BPEMLFRProsSpatial.csv", :LFRPros)

# println("\nTesting interpolation...")
# testTransfer::MBD.CR3BPMultipleShooterProblem = interpSolution(targeter, "FamilyData/CR3BPEMLFRProsSpatial.csv", "z", sin(pi/6)*(Moon.bodyRadius+h_M)/getCharLength(systemData), numNodes, h_M/getCharLength(systemData), h_E/getCharLength(systemData), gamma)
# println("Test Transfer:\n\tIC:\t$(testTransfer.nodes[1].state.data[1:6])\n\tTOF:\t$(getTOF(targeter, testTransfer))\n\tJC:\t$(getJacobiConstant(dynamicsModel, testTransfer.nodes[1].state.data[1:6]))\n")

# println("Plotting transfer...")
# plotTransfer::Int64 = getNumMembers(solutions)
# orbitArc::MBD.CR3BPArc = propagate(propagator, solutions.nodes[plotTransfer][1].state.data[1:6], [0, -getTOF(targeter, solutions.segments[plotTransfer])], dynamicsModel)
# xData::Vector{Float64}, yData::Vector{Float64}, zData::Vector{Float64} = Vector{Float64}(undef, getStateCount(orbitArc)), Vector{Float64}(undef, getStateCount(orbitArc)), Vector{Float64}(undef, getStateCount(orbitArc))
# for s::Int64 in 1:getStateCount(orbitArc)
#     xData[s], yData[s], zData[s] = orbitArc.states[s][1], orbitArc.states[s][2], orbitArc.states[s][3]
# end
# (figure, axis) = set3DPlotParameters(L"Earth-Moon Spatial LFR ($JC=%$(round(getJacobiConstant(dynamicsModel, solutions.nodes[plotTransfer][1].state.data[1:6]); digits = 4))$)", L"$x$ [ndim]", L"$y$ [ndim]", L"$z$ [ndim]")
# GLMakie.lines!(axis, xData, yData, zData, color = :white, label = L"\mathrm{Spatial\ Lunar\ Free\ Return}")
# GLMakie.scatter!(axis, getPrimaryState(dynamicsModel, 1)[1], getPrimaryState(dynamicsModel, 1)[2], getPrimaryState(dynamicsModel, 1)[3], color = :blue, markerspace = :data, markersize = 2*Earth.bodyRadius/getCharLength(systemData), label = L"\mathrm{Earth}" => (; markersize = 20))
# GLMakie.scatter!(axis, getPrimaryState(dynamicsModel, 2)[1], getPrimaryState(dynamicsModel, 2)[2], getPrimaryState(dynamicsModel, 2)[3], color = :gray, markerspace = :data, markersize = 2*Moon.bodyRadius/getCharLength(systemData), label = L"\mathrm{Moon}" => (; markersize = 20))
# GLMakie.Legend(figure[1,2], axis)

# mf = MATLAB.MatFile("Output/CR3BPTraj.mat", "w")
# exportCR3BPTrajectory(testTransfer, mf, :CR3BPTraj)
# MATLAB.close(mf)

println()
end
