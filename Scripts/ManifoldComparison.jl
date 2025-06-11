"""
Script for comparing CR3BP and BCR4BP invariant and pseudo-manifolds

Author: Jonathan Richmond
C: 6/10/25
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
p1::Int64, q1::Int64 = 3, 2
guessOrbit1::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", CR3BPP*q1/p1)
orbit1::MBD.BCR4BP12PeriodicOrbit = getResonantOrbit(targeter, guessOrbit1, p1, q1, Deltaeps = 0.001, tol = 1E-9, refTol = 1E-9, JTol = 1E-2)
compOrbit1::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", targetP*q1/p1)
propTime::Float64 = pi*3
BCR4BPManifold::MBD.BCR4BP12Manifold = getManifoldByArclength(orbit1, "Unstable", "Positive", 25/get12CharLength(systemData), 100*p1)
BCR4BPManifold.TOF = propTime
CR3BPManifold::MBD.CR3BPManifold = getManifoldByArclength(compOrbit1, "Unstable", "Positive", 25/getCharLength(CR3BPSystemData), 100)
CR3BPManifold.TOF = propTime
pseudoManifold::MBD.BCR4BPPseudoManifold = getPseudoManifoldByArclength(compOrbit1, dynamicsModel, pi*1/2, 100)
pseudoManifold.TOF = pi*2

mf = MATLAB.MatFile("Output/ManifoldComparison.mat", "w")
exportBCR4BP12Manifold(BCR4BPManifold, mf, :BCR4BPManifold)
exportCR3BPManifold(CR3BPManifold, mf, :CR3BPManifold)
exportPseudoManifold(pseudoManifold, mf, :pseudoManifold)
MATLAB.close(mf)

println()
end
