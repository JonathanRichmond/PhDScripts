"""
Script for Sun-Mars CR3BP L2 southern halo orbit family

Author: Jonathan LeFevre Richmond
C: 1/24/26
"""
# module SML2Halo
println("Running SML2Halos.jl...\n")

using MBD, GLMakie, Logging

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../CR3BPTargeters/SpatialPerpVy.jl")
include("../Utilities/Export.jl")
include("../Utilities/Plot.jl")

systemData = MBD.CR3BPSystemData("Sun", "Mars")
dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
Sun::MBD.BodyData, Mars::MBD.BodyData = systemData.primaryData[1], systemData.primaryData[2]
LagrangePoints::Vector{Vector{Float64}} = [getEquilibriumPoint(dynamicsModel, l) for l = 1:5]

propagator = MBD.Propagator()
targeter = SpatialPerpVyTargeter(dynamicsModel)

initialStateGuess::Vector{Float64} = [1.0039864420872824, 0, 0, 0, 0.004615276462168598, 0]
tSpanGuess::Vector{Float64} = [0, 3.0914797887407874]
solution1::MBD.CR3BPMultipleShooterProblem = correct(targeter, initialStateGuess, tSpanGuess, initialStateGuess[5])
println("Converged Orbit 1:\n\tIC:\t$(solution1.nodes[1].state.data[1:6])\n\tP:\t$(getPeriod(targeter, solution1))\n\tJC:\t$(getJacobiConstant(dynamicsModel, solution1.nodes[1].state.data[1:6]))\n")

solution2::MBD.CR3BPMultipleShooterProblem = correct(targeter, initialStateGuess+[0, 0, 0.0002, 0, 0, 0], tSpanGuess, initialStateGuess[5]+1E-5)
println("Converged Orbit 2:\n\tIC:\t$(solution2.nodes[1].state.data[1:6])\n\tP:\t$(getPeriod(targeter, solution2))\n\tJC:\t$(getJacobiConstant(dynamicsModel, solution2.nodes[1].state.data[1:6]))\n")

println("Continuing orbits...")
continuationEngine = MBD.CR3BPNaturalParameterContinuationEngine(solution1, solution2, "Initial State", 3, 1E-5, 1E-2)
z0JumpCheck = MBD.BoundingBoxJumpCheck("Initial State", [NaN NaN; 0 0.005; NaN NaN])
addJumpCheck!(continuationEngine, z0JumpCheck)
MarsEndCheck = MBD.CR3BPPrimarySurfaceContinuationEndCheck(dynamicsModel, 2)
addEndCheck!(continuationEngine, MarsEndCheck)
numEndCheck = MBD.NumberStepsContinuationEndCheck(3000)
addEndCheck!(continuationEngine, numEndCheck)
family = MBD.CR3BPOrbitFamily(dynamicsModel)
solutions::MBD.CR3BPContinuationFamily = doContinuation!(continuationEngine, solution1, solution2)
lastOrbit::MBD.CR3BPPeriodicOrbit = getIndividualPeriodicOrbit(targeter, solutions, getNumMembers(solutions))
println("\nLast Converged Orbit:\n\tIC:\t$(lastOrbit.initialCondition)\n\tP:\t$(lastOrbit.period)\n\tJC:\t$(getJacobiConstant(lastOrbit))\n")

for s::Int64 in 1:getNumMembers(solutions)
    orbit::MBD.CR3BPPeriodicOrbit = getIndividualPeriodicOrbit(targeter, solutions, s)
    push!(family.initialConditions, orbit.initialCondition)
    push!(family.periods, orbit.period)
    push!(family.monodromies, orbit.monodromy)
end
eigenSort!(family)

# println("\nExporting family data...")
# fullExportCR3BPFamily(family, "FamilyData/CR3BPSML2Halos.mat", "FamilyData/CR3BPSML2Halos.csv", :L2Halos)

# println("\nTesting interpolation...")
# testOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(targeter, "FamilyData/CR3BPSML2Halos.csv", "JC", 3.0001)
# println("Test Orbit:\n\tIC:\t$(testOrbit.initialCondition)\n\tP:\t$(testOrbit.period)\n\tJC:\t$(getJacobiConstant(testOrbit))\n")

# println("Plotting orbit...")
# plotOrbit::Int64 = 200
# orbitArc::MBD.CR3BPArc = propagate(propagator, family.initialConditions[plotOrbit], [0, family.periods[plotOrbit]], dynamicsModel)
# xData::Vector{Float64}, yData::Vector{Float64}, zData::Vector{Float64} = Vector{Float64}(undef, getStateCount(orbitArc)), Vector{Float64}(undef, getStateCount(orbitArc)), Vector{Float64}(undef, getStateCount(orbitArc))
# for s::Int64 in 1:getStateCount(orbitArc)
#     xData[s], yData[s], zData[s] = orbitArc.states[s][1], orbitArc.states[s][2], orbitArc.states[s][3]
# end
# (figure, axis) = set3DPlotParameters(L"Sun-Mars $L_{2}$ Halo ($JC=%$(round(getJacobiConstant(dynamicsModel, family.initialConditions[plotOrbit]); digits = 6))$)", L"$x$ [ndim]", L"$y$ [ndim]", L"$z$ [ndim]")
# GLMakie.lines!(axis, xData, yData, zData, color = :white, label = L"$L_{2}$ Halo Orbit")
# GLMakie.scatter!(axis, LagrangePoints[1][1], LagrangePoints[1][2], LagrangePoints[1][3], color = :orange, markersize = 5, label = L"$L_{1}$" => (; markersize = 20))
# GLMakie.scatter!(axis, LagrangePoints[2][1], LagrangePoints[2][2], LagrangePoints[2][3], color = :green, markersize = 5, label = L"$L_{2}$" => (; markersize = 20))
# GLMakie.scatter!(axis, getPrimaryState(dynamicsModel, 2)[1], getPrimaryState(dynamicsModel, 2)[2], getPrimaryState(dynamicsModel, 2)[3], color = :red, markerspace = :data, markersize = 2*Mars.bodyRadius/getCharLength(systemData), label = L"\mathrm{Mars}" => (; markersize = 20))
# GLMakie.Legend(figure[1,2], axis)

println()
# end
