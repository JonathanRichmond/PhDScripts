"""
Script for comparing CR3BP and BCR4BP invariant and pseudo-manifolds

Author: Jonathan Richmond
C: 6/10/25
U: 6/18/25
"""
module ManComp
println()

using MBD, Logging, MATLAB

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../BCR4BPTargeters/PlanarPerpP.jl")
include("../CR3BPTargeters/PlanarPerpJC.jl")
include("../Utilities/Export.jl")

systemData = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
dynamicsModel = MBD.BCR4BP12DynamicsModel(systemData)
CR3BPSystemData = MBD.CR3BPSystemData("Earth", "Moon")
CR3BPDynamicsModel = MBD.CR3BPDynamicsModel(CR3BPSystemData)
Earth::MBD.BodyData, Moon::MBD.BodyData = systemData.primaryData[1], systemData.primaryData[2]

propagator = MBD.Propagator()
targeter = PlanarPerpP12Targeter(dynamicsModel)
CR3BPTargeter = PlanarPerpJCTargeter(CR3BPDynamicsModel)

CR3BPP::Float64, targetP::Float64 = 6.2840049277323855, getSynodicPeriod(dynamicsModel)

familyFile::String = "FamilyData/CR3BPEML2Lyapunovs.csv"
p1::Int64, q1::Int64 = 2, 1
guessOrbit1::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", CR3BPP*q1/p1)
orbit1::MBD.BCR4BP12PeriodicOrbit = getResonantOrbit(targeter, guessOrbit1, 0.0, p1, q1, Deltaeps = 0.01, tol = 1E-11, refTol = 1E-11, JTol = 2E-3)
compOrbit1::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", targetP*q1/p1)

propTime::Float64 = pi*5
BCR4BPPosManifold::MBD.BCR4BP12Manifold = getManifoldByArclength(orbit1, "Unstable", "Positive", 25/get12CharLength(systemData), 100*p1)
BCR4BPPosManifold.TOF = propTime
BCR4BPPosManifoldArcs::Vector{MBD.BCR4BP12ManifoldArc} = stopCrashes(BCR4BPPosManifold)
BCR4BPNegManifold::MBD.BCR4BP12Manifold = getManifoldByArclength(orbit1, "Unstable", "Negative", 25/get12CharLength(systemData), 100*p1)
BCR4BPNegManifold.TOF = propTime
BCR4BPNegManifoldArcs::Vector{MBD.BCR4BP12ManifoldArc} = stopCrashes(BCR4BPNegManifold)
BCR4BPManifoldArcs::Vector{MBD.BCR4BP12ManifoldArc} = append!(BCR4BPPosManifoldArcs, BCR4BPNegManifoldArcs)
CR3BPPosManifold::MBD.CR3BPManifold = getManifoldByArclength(compOrbit1, "Unstable", "Positive", 25/getCharLength(CR3BPSystemData), 100)
CR3BPPosManifold.TOF = propTime
CR3BPPosManifoldArcs::Vector{MBD.CR3BPManifoldArc} = stopCrashes(CR3BPPosManifold)
CR3BPNegManifold::MBD.CR3BPManifold = getManifoldByArclength(compOrbit1, "Unstable", "Negative", 25/getCharLength(CR3BPSystemData), 100)
CR3BPNegManifold.TOF = propTime
CR3BPNegManifoldArcs::Vector{MBD.CR3BPManifoldArc} = stopCrashes(CR3BPNegManifold)
CR3BPManifoldArcs::Vector{MBD.CR3BPManifoldArc} = append!(CR3BPPosManifoldArcs, CR3BPNegManifoldArcs)
pseudoManifold0::MBD.BCR4BPPseudoManifold = getPseudoManifoldByArclength(compOrbit1, dynamicsModel, pi*0/2, 100)
pseudoManifold0.TOF = pi*2
pseudoManifoldArcs0::Vector{MBD.BCR4BPPseudoManifoldArc} = stopCrashes(pseudoManifold0)
theta40::Vector{Float64} = collect(range(0, 2*pi, 21))
pseudoManifoldArcs::Vector{MBD.BCR4BPPseudoManifoldArc} = []
for t::Float64 in theta40[1:end-1]
    pseudoManifold::MBD.BCR4BPPseudoManifold = getPseudoManifoldByArclength(compOrbit1, dynamicsModel, t, 50)
    pseudoManifold.TOF = pi*2
    pseudoManifoldArcsSub::Vector{MBD.BCR4BPPseudoManifoldArc} = stopCrashes(pseudoManifold)
    append!(pseudoManifoldArcs, pseudoManifoldArcsSub)
end
pseudoManifold0::MBD.BCR4BPPseudoManifold = getPseudoManifoldByArclength(compOrbit1, dynamicsModel, 0.0, 50)

mf = MATLAB.MatFile("Output/ManifoldComparison.mat", "w")
exportBCR4BP12Manifold(BCR4BPPosManifold, BCR4BPManifoldArcs, mf, :BCR4BPManifold)
exportCR3BPManifold(CR3BPPosManifold, CR3BPManifoldArcs, mf, :CR3BPManifold)
exportPseudoManifold(pseudoManifold0, pseudoManifoldArcs, mf, :pseudoManifold)
MATLAB.close(mf)

println()
end
