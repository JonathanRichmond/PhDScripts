"""
Script for BCR4BP Earth-Moon spatial orbits

Author: Jonathan Richmond
C: 6/11/25
"""
module EMSpatial
println()

using MBD, Logging, MATLAB

global_logger(ConsoleLogger(stderr, Logging.Info)) # Debug, Info, Warn, Error

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

familyFile::String = "FamilyData/CR3BPEML1Halos.csv"
CR3BPP::Float64, targetP::Float64 = 6.2840049277323855, getSynodicPeriod(dynamicsModel)

p1::Int64, q1::Int64 = 3, 1
guessOrbit1::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", CR3BPP*q1/p1; choiceIndex = 1)
orbit1::MBD.BCR4BP12PeriodicOrbit = getResonantOrbit(targeter, guessOrbit1, p1, q1, Deltaeps = 0.001, tol = 1E-10, refTol = 1E-10, JTol = 2E-3)
compOrbit1::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", targetP*q1/p1; choiceIndex = 1)

# homoSystemData::MBD.BCR4BPSystemData = MBD.shallowClone(systemData)
# homoDynamicsModel = MBD.BCR4BP12DynamicsModel(homoSystemData)
# homoTargeter = SpatialPerpP12Targeter(homoDynamicsModel)
# homoSystemData.primaryData[3].gravParam *= 0.089
# targetP= q1*getSynodicPeriod(homoDynamicsModel)
# solution::MBD.BCR4BP12MultipleShooterProblem = correct(homoTargeter, [0.9863068146479186, 0.0, 0.047384187137750715, 0.0, -0.6809168211880687, 0.0, 0.0], targetP, q1; tol = 1E-10, JTol = 1E-2)
# println(solution.nodes[1].state.data[1:7])

# p2::Int64, q2::Int64 = 3, 2
# guessOrbit2::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", CR3BPP*q2/p2)
# orbit2::MBD.BCR4BP12PeriodicOrbit = getResonantOrbit(targeter, guessOrbit2, p2, q2, Deltaeps = 0.001, tol = 1E-9, refTol = 1E-9, JTol = 1E-2)
# compOrbit2::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", targetP*q2/p2)

mf = MATLAB.MatFile("Output/EMSpatialOrbit.mat", "w")
exportBCR4BP12Orbit(orbit1, mf, :BCR4BPOrbit1)
exportCR3BPOrbit(guessOrbit1, mf, :CR3BPGuessOrbit1)
exportCR3BPOrbit(compOrbit1, mf, :CR3BPCompOrbit1)
# exportBCR4BP12Orbit(orbit2, mf, :BCR4BPOrbit2)
# exportCR3BPOrbit(guessOrbit2, mf, :CR3BPGuessOrbit2)
# exportCR3BPOrbit(compOrbit2, mf, :CR3BPCompOrbit2)
MATLAB.close(mf)

println()
end
