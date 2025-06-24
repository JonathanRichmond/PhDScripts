"""
Script for comparing CR3BP and BCR4BP invariant and pseudo-manifolds

Author: Jonathan Richmond
C: 6/10/25
U: 6/24/25
"""
module ManComp
println()

using MBD, Logging, MATLAB

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../BCR4BPTargeters/SpatialPerpP.jl")
include("../CR3BPTargeters/SpatialPerpVy.jl")
include("../Utilities/Export.jl")

systemData = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
dynamicsModel = MBD.BCR4BP12DynamicsModel(systemData)
CR3BPSystemData = MBD.CR3BPSystemData("Earth", "Moon")
CR3BPDynamicsModel = MBD.CR3BPDynamicsModel(CR3BPSystemData)
Earth::MBD.BodyData, Moon::MBD.BodyData = systemData.primaryData[1], systemData.primaryData[2]

propagator = MBD.Propagator()
targeter = SpatialPerpP12Targeter(dynamicsModel)
CR3BPTargeter = SpatialPerpVyTargeter(CR3BPDynamicsModel)

familyFile::String = "FamilyData/CR3BPEML2Halos.csv"
p::Int64, q::Int64 = 2, 1
compOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", getSynodicPeriod(dynamicsModel)*q/p; choiceIndex = 1)
println("Converged $p:$q CR3BP Orbit:\n\tIC:\t$(compOrbit.initialCondition)\n\tP:\t$(compOrbit.period)\n\tStability:\t$(getStabilityIndex(compOrbit))\n")
orbitArc::MBD.CR3BPArc = propagate(propagator, compOrbit.initialCondition, [0, compOrbit.period/2], CR3BPDynamicsModel)
halfState::Vector{Float64} = getStateByIndex(orbitArc, -1)
compOrbit = MBD.CR3BPPeriodicOrbit(CR3BPDynamicsModel, halfState, compOrbit.period, Matrix{Float64}(compOrbit.monodromy))
q0JumpCheck = MBD.BoundingBoxJumpCheck("Initial State", [0.9 1.3; -0.1 0; -0.5 0])
orbit::MBD.BCR4BP12PeriodicOrbit = getResonantOrbit(targeter, compOrbit, 0.0, p, q, q0JumpCheck)

propTime::Float64 = pi*2
BCR4BPPosManifold::MBD.BCR4BP12Manifold = getManifoldByArclength(orbit, "Unstable", "Positive", 25/get12CharLength(systemData), 100*p)
BCR4BPPosManifold.TOF = propTime
BCR4BPPosManifoldArcs::Vector{MBD.BCR4BP12ManifoldArc} = stopCrashes(BCR4BPPosManifold)
BCR4BPNegManifold::MBD.BCR4BP12Manifold = getManifoldByArclength(orbit, "Unstable", "Negative", 25/get12CharLength(systemData), 100*p)
BCR4BPNegManifold.TOF = propTime
BCR4BPNegManifoldArcs::Vector{MBD.BCR4BP12ManifoldArc} = stopCrashes(BCR4BPNegManifold)
BCR4BPManifoldArcs::Vector{MBD.BCR4BP12ManifoldArc} = append!(BCR4BPPosManifoldArcs, BCR4BPNegManifoldArcs)
CR3BPPosManifold::MBD.CR3BPManifold = getManifoldByArclength(compOrbit, "Unstable", "Positive", 25/getCharLength(CR3BPSystemData), 100)
CR3BPPosManifold.TOF = propTime
CR3BPPosManifoldArcs::Vector{MBD.CR3BPManifoldArc} = stopCrashes(CR3BPPosManifold)
CR3BPNegManifold::MBD.CR3BPManifold = getManifoldByArclength(compOrbit, "Unstable", "Negative", 25/getCharLength(CR3BPSystemData), 100)
CR3BPNegManifold.TOF = propTime
CR3BPNegManifoldArcs::Vector{MBD.CR3BPManifoldArc} = stopCrashes(CR3BPNegManifold)
CR3BPManifoldArcs::Vector{MBD.CR3BPManifoldArc} = append!(CR3BPPosManifoldArcs, CR3BPNegManifoldArcs)
pseudoManifold0::MBD.BCR4BPPseudoManifold = getPseudoManifoldByArclength(compOrbit, dynamicsModel, 0.0, 50)
pseudoManifold0.TOF = pi*1.5
theta40::Vector{Float64} = collect(range(0, 2*pi, 21))
pseudoManifoldArcs::Vector{MBD.BCR4BPPseudoManifoldArc} = []
for t::Float64 in theta40[1:end-1]
    pseudoManifold::MBD.BCR4BPPseudoManifold = getPseudoManifoldByArclength(compOrbit, dynamicsModel, t, 50)
    pseudoManifold.TOF = pseudoManifold0.TOF
    pseudoManifoldArcsSub::Vector{MBD.BCR4BPPseudoManifoldArc} = stopCrashes(pseudoManifold)
    append!(pseudoManifoldArcs, pseudoManifoldArcsSub)
end

mf = MATLAB.MatFile("Output/ManifoldComparison.mat", "w")
exportCR3BPOrbit(compOrbit, mf, :CR3BPCompOrbit)
exportBCR4BP12Orbit(orbit, mf, :BCR4BPOrbit)
exportCR3BPManifold(CR3BPPosManifold, CR3BPManifoldArcs, mf, :CR3BPManifold)
exportBCR4BP12Manifold(BCR4BPPosManifold, BCR4BPManifoldArcs, mf, :BCR4BPManifold)
exportPseudoManifold(pseudoManifold0, pseudoManifoldArcs, mf, :pseudoManifold)
MATLAB.close(mf)

println()
end
