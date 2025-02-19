"""
Script for Earth-Moon CR3BP L1 Lyapunov orbit family

Author: Jonathan Richmond
C: 2/4/25
U: 2/17/25
"""
module EML1Lyap
println()

using MBD, MATLAB

include("../Targeters/PlanarPerpJC.jl")
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

targeter = PlanarPerpJCTargeter(dynamicsModel)

(initialStateGuess::Vector{Float64}, tSpanGuess::Vector{Float64}) = getLinearVariation(dynamicsModel, 1, L1, [0.005, 0, 0])
targetJC::Float64 = getJacobiConstant(dynamicsModel, initialStateGuess)
solution1::MBD.CR3BPMultipleShooterProblem = correct(targeter, initialStateGuess, tSpanGuess, targetJC)
println("Converged Orbit 1:\n\tState:$(solution1.nodes[1].state.data[1:6])\n\tPeriod: $(getPeriod(targeter, solution1))\n\tJC: $(getJacobiConstant(dynamicsModel, solution1.nodes[1].state.data[1:6]))")

solution2::MBD.CR3BPMultipleShooterProblem = correct(targeter, initialStateGuess, tSpanGuess, targetJC-1E-4)
println("\nConverged Orbit 2:\n\tState:$(solution2.nodes[1].state.data[1:6])\n\tPeriod: $(getPeriod(targeter, solution2))\n\tJC: $(getJacobiConstant(dynamicsModel, solution2.nodes[1].state.data[1:6]))")

engine = MBD.JacobiConstantContinuationEngine(solution1, solution2, -1E-4, -1E-2)
ydot0JumpCheck = MBD.BoundingBoxJumpCheck("Initial State", [NaN NaN; -2.5 0])
addJumpCheck!(engine, ydot0JumpCheck)
MoonEndCheck = MBD.BoundingBoxContinuationEndCheck("Initial State", [L1[1] getPrimaryPosition(dynamicsModel, 2)[1]-Moon.bodyRadius/getCharLength(systemData); NaN NaN])
addEndCheck!(engine, MoonEndCheck)

println()
family = MBD.CR3BPOrbitFamily(dynamicsModel)
solutions::MBD.CR3BPContinuationFamily = doContinuation!(engine, solution1, solution2)
lastOrbit::MBD.CR3BPPeriodicOrbit = getIndividualPeriodicOrbit(targeter, solutions, getNumMembers(solutions))
println("\nLast Converged Orbit:\n\tState:$(lastOrbit.initialCondition)\n\tPeriod: $(lastOrbit.period)\n\tJC: $(getJacobiConstant(lastOrbit))")

println()
for s::Int64 in 1:getNumMembers(solutions)
    orbit::MBD.CR3BPPeriodicOrbit = getIndividualPeriodicOrbit(targeter, solutions, s)
    push!(family.initialConditions, orbit.initialCondition)
    push!(family.periods, orbit.period)
    push!(family.monodromies, orbit.monodromy)
end
eigenSort!(family)

# exportData(family, "FamilyData/EML1Lyapunovs.csv")

println("\nTesting interpolation...")
testOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(targeter, "FamilyData/EML1Lyapunovs.csv", "JC", 3.0)
println("\nTest Orbit:\n\tState:$(testOrbit.initialCondition)\n\tPeriod: $(testOrbit.period)\n\tJC: $(getJacobiConstant(testOrbit))")

mf = MATLAB.MatFile("Output/CR3BPTraj.mat", "w")
exportCR3BPOrbit(testOrbit, dynamicsModel, mf)
MATLAB.close(mf)

println()
end
