"""
Auxiliary script for code development

Author: Jonathan Richmond
U: 5/28/25
"""
module Aux
println()

using MBD, DifferentialEquations, LinearAlgebra

include("../CR3BPTargeters/PlanarPerpJC.jl")

function reconstructSTM(Q::Matrix{Float64}, rs::Vector{Matrix{Float64}})
    totalR = LinearAlgebra.I
    for R::Matrix{Float64} in reverse(Rs)
        totalR *= R
    end

    return Q*totalR
end

systemData = MBD.CR3BPSystemData("Earth", "Moon")
dynamicsModel = MBD.CR3BPDynamicsModel(systemData)

targeter = PlanarPerpJCTargeter(dynamicsModel)
propagator = MBD.Propagator(equationType = MBD.STM)

L1::Vector{Float64} = getEquilibriumPoint(dynamicsModel, 1)
(initialStateGuess::Vector{Float64}, tSpanGuess::Vector{Float64}) = getLinearVariation(dynamicsModel, 1, L1, [0.005, 0, 0])
renormalizeEvent = DifferentialEquations.PeriodicCallback(renormalize!, tSpanGuess[2]/10)
Rs::Vector{Matrix{Float64}} = []
arc::MBD.CR3BPArc = propagate(propagator, appendExtraInitialConditions(dynamicsModel, initialStateGuess, MBD.STM), tSpanGuess, dynamicsModel)
Phi::Matrix{Float64} = reshape(getStateByIndex(arc, -1)[7:42], (6, 6))
arcRecon::MBD.CR3BPArc = propagateWithPeriodicEvent(propagator, renormalizeEvent, appendExtraInitialConditions(dynamicsModel, initialStateGuess, MBD.STM), tSpanGuess, dynamicsModel, [dynamicsModel, Rs])
PhiRecon::Matrix{Float64} = reconstructSTM(reshape(getStateByIndex(arcRecon, -1)[7:42], (6, 6)), Rs)
println(Phi)
println(PhiRecon)

println()
end
