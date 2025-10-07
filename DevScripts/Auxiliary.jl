"""
Auxiliary script for code development

Author: Jonathan Richmond
U: 10/1/25
"""
module Aux
println()

using MBD, Logging

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../BCR4BPTargeters/PlanarPerpP.jl")
include("../CR3BPTargeters/PlanarPerpJC.jl")

EMSystemData = MBD.CR3BPSystemData("Earth", "Moon")
EMSSystemData = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
EMDynamicsModel = MBD.CR3BPDynamicsModel(EMSystemData)
EMSDynamicsModel = MBD.BCR4BP12DynamicsModel(EMSSystemData)

CR3BPTargeter = PlanarPerpJCTargeter(EMDynamicsModel)
BCR4BPTargeter = PlanarPerpP12Targeter(EMSDynamicsModel)

familyFile::String = "FamilyData/CR3BPEML2Lyapunovs.csv"
p::Int64, q::Int64 = 1, 1
compOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", getSynodicPeriod(EMSDynamicsModel)*q/p; choiceIndex = 1)
println("Converged $p:$q CR3BP Orbit:\n\tIC:\t$(compOrbit.initialCondition)\n\tP:\t$(compOrbit.period)\n\tStab.:\t$(getStabilityIndex(compOrbit))\n\tJC:\t$(getJacobiConstant(compOrbit))\n")
q0JumpCheck = MBD.BoundingBoxJumpCheck("Node 1 State", [0.9 1.0; 1.0 2.0])
orbit::MBD.BCR4BP12PeriodicOrbit = getResonantOrbit(BCR4BPTargeter, compOrbit, 4*q, 0.0, p, q, q0JumpCheck)

println()
end
