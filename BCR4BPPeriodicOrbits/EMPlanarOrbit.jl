"""
Script for BCR4BP Earth-Moon planar orbits

Author: Jonathan Richmond
C: 4/23/25
U: 6/23/25
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
SB1DynamicsModel = MBD.BCR4BP41DynamicsModel(systemData)
CR3BPSystemData = MBD.CR3BPSystemData("Earth", "Moon")
CR3BPDynamicsModel = MBD.CR3BPDynamicsModel(CR3BPSystemData)
Earth::MBD.BodyData, Moon::MBD.BodyData = systemData.primaryData[1], systemData.primaryData[2]

propagator = MBD.Propagator()
targeter = PlanarPerpP12Targeter(dynamicsModel)
CR3BPTargeter = PlanarPerpJCTargeter(CR3BPDynamicsModel)

familyFile::String = "FamilyData/CR3BPEML2Lyapunovs.csv"
p::Int64, q::Int64 = 2, 1
compOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", getSynodicPeriod(dynamicsModel)*q/p; choiceIndex = 1)
println("Converged $p:$q CR3BP Orbit:\n\tIC:\t$(compOrbit.initialCondition)\n\tP:\t$(compOrbit.period)\n")
q0JumpCheck = MBD.BoundingBoxJumpCheck("IntialState", [0.9 1.3; 0.0 0.5])
orbit::MBD.BCR4BP12PeriodicOrbit = getResonantOrbit(targeter, compOrbit, 0.0, p, q, q0JumpCheck)
(orbitSB1::Vector{Vector{Float64}}, ~) = rotating12ToRotating41(dynamicsModel, [orbit.initialCondition], [0.0])

mf = MATLAB.MatFile("Output/EMPlanarOrbit.mat", "w")
exportCR3BPOrbit(compOrbit, mf, :CR3BPCompOrbit)
exportBCR4BP12Orbit(orbit, mf, :BCR4BPOrbit)
exportBCR4BP41Trajectory(orbitSB1[1]..., orbit.period*get12CharTime(systemData)/get41CharTime(systemData), SB1DynamicsModel, mf, :BCR4BPSB1Orbit)
MATLAB.close(mf)

println()
end
