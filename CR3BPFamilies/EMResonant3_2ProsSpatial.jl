"""
Script for Earth-Moon CR3BP 3:2 prograde spatial resonant orbit family
    ***NEED TO CLASSIFY***

Author: Jonathan Richmond
C: 7/1/25
U: 7/2/25
"""
module EMRes32ProSpat
println()

using MBD, GLMakie, Logging

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../CR3BPTargeters/SpatialPerpVy.jl")
include("../Utilities/Export.jl")
include("../Utilities/Plot.jl")

systemData = MBD.CR3BPSystemData("Earth", "Moon")
dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
Earth::MBD.BodyData, Moon::MBD.BodyData = systemData.primaryData[1], systemData.primaryData[2]
LagrangePoints::Vector{Vector{Float64}} = [getEquilibriumPoint(dynamicsModel, l) for l = 1:5]

propagator = MBD.Propagator()
targeter = SpatialPerpVyTargeter(dynamicsModel)

initialStateGuess::Vector{Float64} = [0.732397, 0, -0.147068, 0, 0.523851, 0]
tSpanGuess::Vector{Float64} = [0, 11.3186]
solution1::MBD.CR3BPMultipleShooterProblem = correct(targeter, initialStateGuess, tSpanGuess, initialStateGuess[5])
println("Converged Orbit 1:\n\tIC:\t$(solution1.nodes[1].state.data[1:6])\n\tP:\t$(getPeriod(targeter, solution1))\n\tJC:\t$(getJacobiConstant(dynamicsModel, solution1.nodes[1].state.data[1:6]))\n")

solution2::MBD.CR3BPMultipleShooterProblem = correct(targeter, initialStateGuess, tSpanGuess, initialStateGuess[5]-1E-4)
println("Converged Orbit 2:\n\tIC:\t$(solution2.nodes[1].state.data[1:6])\n\tP:\t$(getPeriod(targeter, solution2))\n\tJC:\t$(getJacobiConstant(dynamicsModel, solution2.nodes[1].state.data[1:6]))\n")

println("Continuing orbits...")
continuationEngine = MBD.CR3BPNaturalParameterContinuationEngine(solution1, solution2, "Initial State", 3, -1E-4, -1E-2)
ydot0JumpCheck = MBD.BoundingBoxJumpCheck("Initial State", [NaN NaN; -0.5 0; NaN NaN])
addJumpCheck!(continuationEngine, ydot0JumpCheck)
planeEndCheck = MBD.BoundingBoxContinuationEndCheck("Initial State", [NaN NaN; -0.5 0; NaN NaN])
addEndCheck!(continuationEngine, planeEndCheck)
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

# println("\nExporting family data...")
# fullExportCR3BPFamily(family, "FamilyData/CR3BPEMResonant3_2ProsSpatial.mat", "FamilyData/CR3BPEMResonant3_2ProsSpatial.csv")

# println("\nTesting interpolation...")
# testOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(targeter, "FamilyData/CR3BPEMResonant3_2ProsSpatial.csv", "JC", 3.0)
# println("Test Orbit:\n\tIC:\t$(testOrbit.initialCondition)\n\tP:\t$(testOrbit.period)\n\tJC:\t$(getJacobiConstant(testOrbit))\n")

# println("Plotting orbit...")
# plotOrbit::Int64 = getNumMembers(family)
# orbitArc::MBD.CR3BPArc = propagate(propagator, family.initialConditions[plotOrbit], [0, family.periods[plotOrbit]], dynamicsModel)
# xData::Vector{Float64}, yData::Vector{Float64} = Vector{Float64}(undef, getStateCount(orbitArc)), Vector{Float64}(undef, getStateCount(orbitArc))
# for s::Int64 in 1:getStateCount(orbitArc)
#     xData[s], yData[s] = orbitArc.states[s][1], orbitArc.states[s][2]
# end
# (figure, axis) = set3DPlotParameters(L"Earth-Moon $3:2$ Prograde Spatial Resonant Orbit ($JC=%$(round(getJacobiConstant(dynamicsModel, family.initialConditions[plotOrbit]); digits = 4))$)", L"$x$ [ndim]", L"$y$ [ndim]", L"$z$ [ndim]")
# GLMakie.lines!(axis, xData, yData, color = :white, label = L"$3:2$ Prograde Spatial Resonant Orbit")
# GLMakie.scatter!(axis, LagrangePoints[1][1], LagrangePoints[1][2], LagrangePoints[1][3], color = :red, markersize = 5, label = L"$L_{1}$" => (; markersize = 20))
# GLMakie.scatter!(axis, LagrangePoints[2][1], LagrangePoints[2][2], LagrangePoints[2][3], color = :orange, markersize = 5, label = L"$L_{2}$" => (; markersize = 20))
# GLMakie.scatter!(axis, LagrangePoints[3][1], LagrangePoints[3][2], LagrangePoints[3][3], color = :yellow, markersize = 5, label = L"$L_{3}$" => (; markersize = 20))
# GLMakie.scatter!(axis, LagrangePoints[4][1], LagrangePoints[4][2], LagrangePoints[4][3], color = :green, markersize = 5, label = L"$L_{4}$" => (; markersize = 20))
# GLMakie.scatter!(axis, LagrangePoints[5][1], LagrangePoints[5][2], LagrangePoints[5][3], color = :purple, markersize = 5, label = L"$L_{5}$" => (; markersize = 20))
# GLMakie.scatter!(axis, getPrimaryState(dynamicsModel, 1)[1], getPrimaryState(dynamicsModel, 1)[2], getPrimaryState(dynamicsModel, 1)[3], color = :blue, markerspace = :data, markersize = 2*Earth.bodyRadius/getCharLength(systemData), label = L"\mathrm{Earth}" => (; markersize = 20))
# GLMakie.scatter!(axis, getPrimaryState(dynamicsModel, 2)[1], getPrimaryState(dynamicsModel, 2)[2], getPrimaryState(dynamicsModel, 2)[3], color = :gray, markerspace = :data, markersize = 2*Moon.bodyRadius/getCharLength(systemData), label = L"\mathrm{Moon}" => (; markersize = 20))
# GLMakie.Legend(figure[1,2], axis)

println()
end
