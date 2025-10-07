"""
Script for Sun-B1 CR3BP L1 Lyapunov orbit family

Author: Jonathan Richmond
C: 10/6/25
"""
module SB1L1Lyap
println()

using MBD, GLMakie, Logging

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../CR3BPTargeters/PlanarPerpJC.jl")
include("../Utilities/Export.jl")
include("../Utilities/Plot.jl")

systemData = MBD.CR3BPSystemData("Sun", "Earth_Barycenter")
dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
Sun::MBD.BodyData, B1::MBD.BodyData = systemData.primaryData[1], systemData.primaryData[2]
Earth = MBD.BodyData("Earth")
LagrangePoints::Vector{Vector{Float64}} = [getEquilibriumPoint(dynamicsModel, l) for l = 1:5]

propagator = MBD.Propagator()
targeter = PlanarPerpJCTargeter(dynamicsModel)

(initialStateGuess::Vector{Float64}, tSpanGuess::Vector{Float64}) = getLinearVariation(dynamicsModel, 2, LagrangePoints[1], [0.0001, 0, 0])
targetJC::Float64 = getJacobiConstant(dynamicsModel, initialStateGuess)
solution1::MBD.CR3BPMultipleShooterProblem = correct(targeter, initialStateGuess, tSpanGuess, targetJC)
println("Converged Orbit 1:\n\tIC:\t$(solution1.nodes[1].state.data[1:6])\n\tP:\t$(getPeriod(targeter, solution1))\n\tJC:\t$(getJacobiConstant(dynamicsModel, solution1.nodes[1].state.data[1:6]))\n")

solution2::MBD.CR3BPMultipleShooterProblem = correct(targeter, initialStateGuess, tSpanGuess, targetJC-1E-5)
println("Converged Orbit 2:\n\tIC:\t$(solution2.nodes[1].state.data[1:6])\n\tP:\t$(getPeriod(targeter, solution2))\n\tJC:\t$(getJacobiConstant(dynamicsModel, solution2.nodes[1].state.data[1:6]))\n")

println("Continuing orbits...")
continuationEngine = MBD.JacobiConstantContinuationEngine(solution1, solution2, -1E-5, -1E-2; boundRadiusScaleFactor = 5000.0)
ydot0JumpCheck = MBD.BoundingBoxJumpCheck("Initial State", [NaN NaN; -2.5 0])
addJumpCheck!(continuationEngine, ydot0JumpCheck)
numEndCheck = MBD.NumberStepsContinuationEndCheck(100)
addEndCheck!(continuationEngine, numEndCheck)
family = MBD.CR3BPOrbitFamily(dynamicsModel)
solutions::MBD.CR3BPContinuationFamily = doContinuation!(continuationEngine, solution1, solution2)
lastOrbit::MBD.CR3BPPeriodicOrbit = getIndividualPeriodicOrbit(targeter, solutions, getNumMembers(solutions))
println("Last Converged Orbit:\n\tIC:\t$(lastOrbit.initialCondition)\n\tP:\t$(lastOrbit.period)\n\tJC:\t$(getJacobiConstant(lastOrbit))\n")

for s::Int64 in 1:getNumMembers(solutions)
    orbit::MBD.CR3BPPeriodicOrbit = getIndividualPeriodicOrbit(targeter, solutions, s)
    push!(family.initialConditions, orbit.initialCondition)
    push!(family.periods, orbit.period)
    push!(family.monodromies, orbit.monodromy)
end
eigenSort!(family)

println("\nExporting family data...")
fullExportCR3BPFamily(family, "FamilyData/CR3BPSB1L1Lyapunovs.mat", "FamilyData/CR3BPSB1L1Lyapunovs.csv", :L1Lyapunovs)

println("\nTesting interpolation...")
testOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(targeter, "FamilyData/CR3BPSB1L1Lyapunovs.csv", "JC", 3.0)
println("Test Orbit:\n\tIC:\t$(testOrbit.initialCondition)\n\tP:\t$(testOrbit.period)\n\tJC:\t$(getJacobiConstant(testOrbit))\n")

println("Plotting orbit...")
plotOrbit::Int64 = getNumMembers(family)
orbitArc::MBD.CR3BPArc = propagate(propagator, family.initialConditions[plotOrbit], [0, family.periods[plotOrbit]], dynamicsModel)
xData::Vector{Float64}, yData::Vector{Float64} = Vector{Float64}(undef, getStateCount(orbitArc)), Vector{Float64}(undef, getStateCount(orbitArc))
for s::Int64 in 1:getStateCount(orbitArc)
    xData[s], yData[s] = orbitArc.states[s][1], orbitArc.states[s][2]
end
(figure, axis) = set2DPlotParameters(L"Sun-$B_{1}$ $L_{1}$ Lyapunov ($JC=%$(round(getJacobiConstant(dynamicsModel, family.initialConditions[plotOrbit]); digits = 4))$)", L"$x$ [ndim]", L"$y$ [ndim]")
GLMakie.lines!(axis, xData, yData, color = :white, label = L"$L_{2}$ Lyapunov Orbit")
GLMakie.scatter!(axis, LagrangePoints[1][1], LagrangePoints[1][2], color = :red, markersize = 5, label = L"$L_{1}$" => (; markersize = 20))
GLMakie.scatter!(axis, LagrangePoints[2][1], LagrangePoints[2][2], color = :orange, markersize = 5, label = L"$L_{2}$" => (; markersize = 20))
GLMakie.scatter!(axis, getPrimaryState(dynamicsModel, 2)[1], getPrimaryState(dynamicsModel, 2)[2], color = :blue, markerspace = :data, markersize = 2*Earth.bodyRadius/getCharLength(systemData), label = L"$B_{1}$" => (; markersize = 20))
GLMakie.Legend(figure[1,2], axis)

println()
end