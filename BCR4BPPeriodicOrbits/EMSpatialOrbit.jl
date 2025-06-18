"""
Script for BCR4BP Earth-Moon spatial orbits

Author: Jonathan Richmond
C: 6/11/25
U: 6/18/25
"""
module EMSpatial
println()

using MBD, Logging, MATLAB

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../BCR4BPTargeters/SpatialPerpP.jl")
include("../CR3BPTargeters/SpatialPerpVy.jl")
include("../Utilities/Export.jl")

systemData = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
dynamicsModel = MBD.BCR4BP12DynamicsModel(systemData)
SB1DynamicsModel = MBD.BCR4BP41DynamicsModel(systemData)
CR3BPSystemData = MBD.CR3BPSystemData("Earth", "Moon")
CR3BPDynamicsModel = MBD.CR3BPDynamicsModel(CR3BPSystemData)
Earth::MBD.BodyData, Moon::MBD.BodyData = systemData.primaryData[1], systemData.primaryData[2]

propagator = MBD.Propagator()
targeter = SpatialPerpP12Targeter(dynamicsModel)
CR3BPTargeter = SpatialPerpVyTargeter(CR3BPDynamicsModel)

familyFile::String = "FamilyData/CR3BPEML2Halos.csv"
CR3BPP::Float64, targetP::Float64 = 6.2840049277323855, getSynodicPeriod(dynamicsModel)

p1::Int64, q1::Int64 = 3, 1
guessOrbit1::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", CR3BPP*q1/p1; choiceIndex = 1)
orbit1::MBD.BCR4BP12PeriodicOrbit = getResonantOrbit(targeter, guessOrbit1, 0.0, p1, q1, Deltaeps = 0.001, tol = 1E-10, refTol = 1E-10, JTol = 2E-3)
compOrbit1::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", targetP*q1/p1; choiceIndex = 1)
(orbit1ICSB1::Vector{Vector{Float64}}, ~) = rotating12ToRotating41(dynamicsModel, [orbit1.initialCondition], [0.0])

p2::Int64, q2::Int64 = 3, 1
guessOrbit2::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", CR3BPP*q2/p2)
orbit2::MBD.BCR4BP12PeriodicOrbit = getResonantOrbit(targeter, guessOrbit2, pi*1.0, p2, q2, Deltaeps = 0.001, tol = 1E-10, refTol = 1E-10, JTol = 2E-3)
compOrbit2::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", targetP*q2/p2)
(orbit2ICSB1::Vector{Vector{Float64}}, ~) = rotating12ToRotating41(dynamicsModel, [orbit2.initialCondition], [0.0])

mf = MATLAB.MatFile("Output/EMSpatialOrbit.mat", "w")
exportBCR4BP12Orbit(orbit1, mf, :BCR4BPOrbit1)
exportCR3BPOrbit(guessOrbit1, mf, :CR3BPGuessOrbit1)
exportCR3BPOrbit(compOrbit1, mf, :CR3BPCompOrbit1)
exportBCR4BP41Trajectory(orbit1ICSB1[1]..., orbit1.period*get12CharTime(systemData)/get41CharTime(systemData), SB1DynamicsModel, mf, :BCR4BPSB1Orbit1)
exportBCR4BP12Orbit(orbit2, mf, :BCR4BPOrbit2)
exportCR3BPOrbit(guessOrbit2, mf, :CR3BPGuessOrbit2)
exportCR3BPOrbit(compOrbit2, mf, :CR3BPCompOrbit2)
exportBCR4BP41Trajectory(orbit2ICSB1[1]..., orbit2.period*get12CharTime(systemData)/get41CharTime(systemData), SB1DynamicsModel, mf, :BCR4BPSB1Orbit2)
MATLAB.close(mf)

println()
end
