"""
Script for BCR4BP Earth-Moon planar orbits

Author: Jonathan Richmond
C: 4/9/25
"""
module EMPlanar
println()

using MBD, MATLAB

include("../BCR4BPTargeters/PlanarPerpP.jl")
include("../Utilities/Export.jl")

systemData = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
dynamicsModel = MBD.BCR4BP12DynamicsModel(systemData)
CR3BPSystemData = MBD.CR3BPSystemData("Earth", "Moon")
CR3BPDynamicsModel = MBD.CR3BPDynamicsModel(CR3BPSystemData)

Earth = systemData.primaryData[1]
Moon = systemData.primaryData[2]

targeter = PlanarPerpP12Targeter(dynamicsModel)

initialStateGuess::Vector{Float64} = append!(getEquilibriumPoint(CR3BPDynamicsModel, 1), [0.0, 0.0, 0.0, 0.0])
targetP::Float64 = getSynodicPeriod(dynamicsModel)
solution1::MBD.BCR4BP12MultipleShooterProblem = correct(targeter, initialStateGuess, targetP)
# println("Converged Orbit 1:\n\tState:$(solution1.nodes[1].state.data[1:6])\n\tPeriod: $(getPeriod(targeter, solution1))\n\tJC: $(getJacobiConstant(dynamicsModel, solution1.nodes[1].state.data[1:6]))")

println()
end
