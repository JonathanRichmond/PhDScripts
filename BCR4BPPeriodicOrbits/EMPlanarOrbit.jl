"""
Script for BCR4BP Earth-Moon planar orbits

Author: Jonathan Richmond
C: 4/23/25
U: 6/10/25
"""
module EMPlanar
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

familyFile::String = "FamilyData/CR3BPEML2Lyapunovs.csv"
CR3BPP::Float64, targetP::Float64 = 6.2840049277323855, getSynodicPeriod(dynamicsModel)

p1::Int64, q1::Int64 = 1, 1
guessOrbit1::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", CR3BPP*q1/p1)
orbit1::MBD.BCR4BP12PeriodicOrbit = getResonantOrbit(targeter, guessOrbit1, p1, q1, Deltaeps = 0.01, tol = 1E-11, refTol = 1E-11, JTol = 2E-3)
compOrbit1::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", targetP*q1/p1)

p2::Int64, q2::Int64 = 3, 2
guessOrbit2::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", CR3BPP*q2/p2)
orbit2::MBD.BCR4BP12PeriodicOrbit = getResonantOrbit(targeter, guessOrbit2, p2, q2, Deltaeps = 0.001, tol = 1E-9, refTol = 1E-9, JTol = 1E-2)
compOrbit2::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", targetP*q2/p2)

mf = MATLAB.MatFile("Output/EMPlanarOrbit.mat", "w")
exportBCR4BP12Orbit(orbit1, mf, :BCR4BPOrbit1)
exportCR3BPOrbit(guessOrbit1, mf, :CR3BPGuessOrbit1)
exportCR3BPOrbit(compOrbit1, mf, :CR3BPCompOrbit1)
exportBCR4BP12Orbit(orbit2, mf, :BCR4BPOrbit2)
exportCR3BPOrbit(guessOrbit2, mf, :CR3BPGuessOrbit2)
exportCR3BPOrbit(compOrbit2, mf, :CR3BPCompOrbit2)
MATLAB.close(mf)

println()
end
