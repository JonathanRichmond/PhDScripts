"""
Script for BCR4BP Earth-Moon spatial orbits

Author: Jonathan Richmond
C: 6/11/25
U: 7/8/25
"""
module EMSpatial
println()

using MBD, Logging, MATLAB

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../BCR4BPTargeters/SpatialAxP.jl")
include("../BCR4BPTargeters/SpatialContP.jl")
include("../BCR4BPTargeters/SpatialPerpP.jl")
include("../CR3BPTargeters/SpatialPerpVy.jl")
include("../CR3BPTargeters/SpatialAxialJC.jl")
include("../Utilities/Export.jl")

systemData = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
dynamicsModel = MBD.BCR4BP12DynamicsModel(systemData)
SB1DynamicsModel = MBD.BCR4BP41DynamicsModel(systemData)
CR3BPSystemData = MBD.CR3BPSystemData("Earth", "Moon")
CR3BPDynamicsModel = MBD.CR3BPDynamicsModel(CR3BPSystemData)
Earth::MBD.BodyData, Moon::MBD.BodyData = systemData.primaryData[1], systemData.primaryData[2]

propagator = MBD.Propagator()
axTargeter = SpatialAxialP12Targeter(dynamicsModel)
contTargeter = SpatialContP12Targeter(dynamicsModel)
perpTargeter = SpatialPerpP12Targeter(dynamicsModel)
CR3BPVyTargeter = SpatialPerpVyTargeter(CR3BPDynamicsModel)
CR3BPAxTargeter = SpatialAxialJCTargeter(CR3BPDynamicsModel)

familyFile::String = "FamilyData/CR3BPEML2Axials.csv"
p::Int64, q::Int64 = 11, 7
numSegs::Int64 = 4*q
compOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPAxTargeter, familyFile, "Period", getSynodicPeriod(dynamicsModel)*q/p; choiceIndex = 1)
println("Converged $p:$q CR3BP Orbit:\n\tIC:\t$(compOrbit.initialCondition)\n\tP:\t$(compOrbit.period)\n")
# orbitArc::MBD.CR3BPArc = propagate(propagator, compOrbit.initialCondition, [0, compOrbit.period/2], CR3BPDynamicsModel)
# halfState::Vector{Float64} = getStateByIndex(orbitArc, -1)
# compOrbit = MBD.CR3BPPeriodicOrbit(CR3BPDynamicsModel, halfState, compOrbit.period, Matrix{Float64}(compOrbit.monodromy))
# println(compOrbit.initialCondition)
# q0JumpCheck = MBD.BoundingBoxJumpCheck("Node 1 State", [0.7 1.3; -0.2 0; 0 0.5])
q0JumpCheck = MBD.BoundingBoxJumpCheck("Node 1 State", [1.0 1.2; NaN NaN; 0.2 0.6])
orbit::MBD.BCR4BP12PeriodicOrbit = getResonantOrbit(axTargeter, compOrbit, numSegs, 0.0, p, q, q0JumpCheck, tol = 1E-11, refTol = 1E-11)
(orbitSB1::Vector{Vector{Float64}}, ~) = rotating12ToRotating41(dynamicsModel, [orbit.initialCondition], [0.0])

mf = MATLAB.MatFile("Output/EMSpatialOrbit.mat", "w")
exportCR3BPOrbit(compOrbit, mf, :CR3BPCompOrbit)
exportBCR4BP12Orbit(orbit, mf, :BCR4BPOrbit)
exportBCR4BP41Trajectory(orbitSB1[1]..., orbit.period*get12CharTime(systemData)/get41CharTime(systemData), SB1DynamicsModel, mf, :BCR4BPSB1Orbit)
MATLAB.close(mf)

println()
end
