"""
Script for BCR4BP Earth-Moon planar orbits

Author: Jonathan Richmond
C: 4/23/25
U: 6/19/25
"""
module EMPlanar
println()

using MBD, Logging, MATLAB

global_logger(ConsoleLogger(stderr, Logging.Info)) # Debug, Info, Warn, Error

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

systemData0 = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
systemData0.P4Mass = 0.0
dynamicsModel0 = MBD.BCR4BP12DynamicsModel(systemData0)
targeter0 = PlanarPerpP12Targeter(dynamicsModel0)
systemData1 = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
systemData1.P4Mass = 0.001*systemData.P4Mass
dynamicsModel1 = MBD.BCR4BP12DynamicsModel(systemData1)
targeter1 = PlanarPerpP12Targeter(dynamicsModel1)

familyFile::String = "FamilyData/CR3BPEML2Lyapunovs.csv"
targetP::Float64 = getSynodicPeriod(dynamicsModel)

p1::Int64, q1::Int64 = 1, 1
compOrbit1::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", targetP*q1/p1)
println("Converged CR3BP Orbit:\n\tIC:\t$(compOrbit1.initialCondition)\n\tP:\t$(compOrbit1.period)\n")
solution0::MBD.BCR4BP12MultipleShooterProblem = correct(targeter0, push!(copy(compOrbit1.initialCondition), 0.0), targetP, tol = 1E-11, JTol = 2E-3)
println("Converged Homotopy Orbit 0:\n\tIC:\t$(solution0.nodes[1].state.data[1:7])\n\tP:\t$(getPeriod(targeter0, solution0))\n")
solution1::MBD.BCR4BP12MultipleShooterProblem = correct(targeter1, copy(solution0.nodes[1].state.data[1:7]), targetP, tol = 1E-11, JTol = 2E-3)
println("Converged Homotopy Orbit 1:\n\tIC:\t$(solution1.nodes[1].state.data[1:7])\n\tP:\t$(getPeriod(targeter1, solution1))\n")
continuationEngine = MBD.P4MassContinuationEngine(solution0, solution1, 0.001, 0.1)
q0JumpCheck = MBD.BoundingBoxJumpCheck("IntialState", [0.9 1.5; 1.0 2.0])
addJumpCheck!(continuationEngine, q0JumpCheck)
numStepsEndCheck = MBD.NumberStepsContinuationEndCheck(5)
addEndCheck!(continuationEngine, numStepsEndCheck)
solutions::MBD.BCR4BP12ContinuationFamily = doContinuation!(continuationEngine, solution0, solution1)

# p2::Int64, q2::Int64 = 2, 1
# guessOrbit2::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", CR3BPP*q2/p2)
# orbit2::MBD.BCR4BP12PeriodicOrbit = getResonantOrbit(targeter, guessOrbit2, 0.0, p2, q2, Deltaeps = 0.01, tol = 1E-11, refTol = 1E-11, JTol = 2E-3)
# compOrbit2::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", targetP*q2/p2)

mf = MATLAB.MatFile("Output/EMPlanarOrbit.mat", "w")
# exportBCR4BP12Orbit(orbit1, mf, :BCR4BPOrbit1)
exportCR3BPOrbit(compOrbit1, mf, :CR3BPCompOrbit1)
# exportBCR4BP12Orbit(orbit2, mf, :BCR4BPOrbit2)
# exportCR3BPOrbit(guessOrbit2, mf, :CR3BPGuessOrbit2)
# exportCR3BPOrbit(compOrbit2, mf, :CR3BPCompOrbit2)
MATLAB.close(mf)

println()
end
