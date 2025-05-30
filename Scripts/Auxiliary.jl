"""
Auxiliary script for code development

Author: Jonathan Richmond
U: 5/28/25
"""
module Aux
println()

using MBD

include("../CR3BPTargeters/PlanarPerpJC.jl")

systemData = MBD.CR3BPSystemData("Earth", "Moon")
dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
L1::Vector{Float64} = getEquilibriumPoint(dynamicsModel, 1)

targeter = PlanarPerpJCTargeter(dynamicsModel)
propagator = MBD.Propagator(equationType = MBD.STM)

(initialStateGuess::Vector{Float64}, tSpanGuess::Vector{Float64}) = getLinearVariation(dynamicsModel, 1, L1, [0.005, 0, 0])
targetJC::Float64 = getJacobiConstant(dynamicsModel, initialStateGuess)
solution1::MBD.CR3BPMultipleShooterProblem = correct(targeter, initialStateGuess, tSpanGuess, targetJC)
M::Matrix{Float64} = getMonodromy(targeter, solution1)
println(M)

println()
end
