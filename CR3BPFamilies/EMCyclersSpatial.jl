"""
Script for Earth-Moon CR3BP spatial cycler orbit family

Author: Jonathan Richmond
C: 7/2/25
U: 7/16/25
"""
module EMCyclerSpat
println()

using MBD, GLMakie, Logging

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../CR3BPTargeters/SpatialPerpJCMS.jl")
include("../Utilities/Export.jl")
include("../Utilities/Plot.jl")

systemData = MBD.CR3BPSystemData("Earth", "Moon")
dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
Earth::MBD.BodyData, Moon::MBD.BodyData = systemData.primaryData[1], systemData.primaryData[2]
LagrangePoints::Vector{Vector{Float64}} = [getEquilibriumPoint(dynamicsModel, l) for l = 1:5]

propagator = MBD.Propagator()
targeter = SpatialPerpJCMSTargeter(dynamicsModel)

numSegs::Int64 = 4
initialStateGuess::Vector{Float64} = [0.81923488, 0, -0.00995508, 0, 0.19902041, 0]
tSpanGuess::Vector{Float64} = [0, 6.5282364/2]
guessArc::MBD.CR3BPArc = propagate(propagator, initialStateGuess, tSpanGuess, dynamicsModel)
numStates::Int64 = getStateCount(guessArc)
numNodes::Int64 = numSegs/2+1
indices::Vector{Int64} = round.(Int64, range(1, numStates, numNodes))
stateGuesses::Vector{Vector{Float64}} = [getStateByIndex(guessArc, indices[i]) for i = eachindex(indices)]
timeGuesses::Vector{Float64} = [getTimeByIndex(guessArc, indices[i]) for i = eachindex(indices)]
targetJC::Float64 = getJacobiConstant(dynamicsModel, initialStateGuess)
solution1::MBD.CR3BPMultipleShooterProblem = correct(targeter, stateGuesses, timeGuesses, numNodes, targetJC)
println("Converged Orbit 1:\n\tIC:\t$(solution1.nodes[1].state.data[1:6])\n\tP:\t$(getPeriod(targeter, solution1))\n\tJC:\t$(getJacobiConstant(dynamicsModel, solution1.nodes[1].state.data[1:6]))\n")

solution2::MBD.CR3BPMultipleShooterProblem = correct(targeter, stateGuesses, timeGuesses, numNodes, targetJC-1E-4)
println("Converged Orbit 2:\n\tIC:\t$(solution2.nodes[1].state.data[1:6])\n\tP:\t$(getPeriod(targeter, solution2))\n\tJC:\t$(getJacobiConstant(dynamicsModel, solution2.nodes[1].state.data[1:6]))\n")

println("Continuing orbits...")
continuationEngine = MBD.JacobiConstantContinuationEngine(solution1, solution2, -1E-4, -1E-2)
z0JumpCheck = MBD.BoundingBoxJumpCheck("Node 1 State", [NaN NaN; -0.2 0; NaN NaN])
addJumpCheck!(continuationEngine, z0JumpCheck)
numStepsEndCheck = MBD.NumberStepsContinuationEndCheck(1000)
addEndCheck!(continuationEngine, numStepsEndCheck)
family = MBD.CR3BPMSOrbitFamily(dynamicsModel)
solutions::MBD.CR3BPContinuationFamily = doContinuation!(continuationEngine, solution1, solution2)
lastOrbit::MBD.CR3BPMSPeriodicOrbit = getIndividualPeriodicOrbit(targeter, solutions, getNumMembers(solutions))
println("Last Converged Orbit:\n\tIC:\t$(lastOrbit.initialCondition)\n\tP:\t$(lastOrbit.period)\n\tJC:\t$(getJacobiConstant(lastOrbit))\n")

for s::Int64 in 1:getNumMembers(solutions)
    orbit::MBD.CR3BPMSPeriodicOrbit = getIndividualPeriodicOrbit(targeter, solutions, s)
    push!(family.initialConditions, orbit.initialCondition)
    push!(family.periods, orbit.period)
    push!(family.monodromies, orbit.monodromy)
    push!(family.nodeStates, orbit.nodeStates)
    push!(family.nodeEpochs, orbit.nodeEpochs)
end
eigenSort!(family)

# println("\nExporting family data...")
# fullExportCR3BPFamily(family, "FamilyData/CR3BPEMCyclersSpatial.mat", "FamilyData/CR3BPEMCyclersSpatial.csv", :CyclersSpatial)

# println("\nTesting interpolation...")
# testOrbit::MBD.CR3BPMSPeriodicOrbit = interpOrbit(targeter, "FamilyData/CR3BPEMCyclersSpatial.csv", "JC", 3.1, numNodes)
# println("Test Orbit:\n\tIC:\t$(testOrbit.initialCondition)\n\tP:\t$(testOrbit.period)\n\tJC:\t$(getJacobiConstant(testOrbit))\n")

# println("Plotting orbit...")
# plotOrbit::Int64 = getNumMembers(family)
# orbitArc::MBD.CR3BPArc = propagate(propagator, family.initialConditions[plotOrbit], [0, family.periods[plotOrbit]], dynamicsModel)
# xData::Vector{Float64}, yData::Vector{Float64}, zData::Vector{Float64} = Vector{Float64}(undef, getStateCount(orbitArc)), Vector{Float64}(undef, getStateCount(orbitArc)), Vector{Float64}(undef, getStateCount(orbitArc))
# for s::Int64 in 1:getStateCount(orbitArc)
#     xData[s], yData[s], zData[s] = orbitArc.states[s][1], orbitArc.states[s][2], orbitArc.states[s][3]
# end
# (figure, axis) = set3DPlotParameters(L"Earth-Moon Spatial Cycler ($JC=%$(round(getJacobiConstant(dynamicsModel, family.initialConditions[plotOrbit]); digits = 4))$)", L"$x$ [ndim]", L"$y$ [ndim]", L"$z$ [ndim]")
# GLMakie.lines!(axis, xData, yData, zData, color = :white, label = L"\mathrm{Spatial\ Cycler\ Orbit}")
# GLMakie.scatter!(axis, LagrangePoints[1][1], LagrangePoints[1][2], LagrangePoints[1][3], color = :red, markersize = 5, label = L"$L_{1}$" => (; markersize = 20))
# GLMakie.scatter!(axis, LagrangePoints[2][1], LagrangePoints[2][2], LagrangePoints[2][3], color = :orange, markersize = 5, label = L"$L_{2}$" => (; markersize = 20))
# GLMakie.scatter!(axis, getPrimaryState(dynamicsModel, 2)[1], getPrimaryState(dynamicsModel, 2)[2], getPrimaryState(dynamicsModel, 2)[3], color = :gray, markerspace = :data, markersize = 2*Moon.bodyRadius/getCharLength(systemData), label = L"\mathrm{Moon}" => (; markersize = 20))
# GLMakie.Legend(figure[1,2], axis)

println()
end
