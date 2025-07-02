"""
Script for computing orbit metrics

Author: Jonathan Richmond
C: 6/30/25
U: 7/2/25
"""
module OrbMetrics
println()

using MBD, MATLAB, Logging

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../../CR3BPTargeters/PlanarPerpJC.jl")
include("../../CR3BPTargeters/PlanarPerpP.jl")
include("../../CR3BPTargeters/SpatialPerpJCMS.jl")
include("../../CR3BPTargeters/SpatialPerpVy.jl")
include("../../Utilities/Export.jl")

systemData = MBD.CR3BPSystemData("Earth", "Moon")
dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
BCR4BPSystemData = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
BCR4BPDynamicsModel = MBD.BCR4BP12DynamicsModel(BCR4BPSystemData)

propagator = MBD.Propagator()
targeter = SpatialPerpJCMSTargeter(dynamicsModel)
familyFile::String = "FamilyData/CR3BPEML1P7HO1s.csv"
numNodes::Int64 = 15

targetP::Float64 = getSynodicPeriod(BCR4BPDynamicsModel)
values::Vector{Float64} = [targetP*7/2]
orbits::Vector{MBD.CR3BPPeriodicOrbit} = []

mf = MATLAB.MatFile("RA/Katalyst/OrbitMetrics.mat", "w")
for v::Int64 in 1:length(values)
    # orbit::MBD.CR3BPPeriodicOrbit = interpOrbit(targeter, familyFile, "Period", values[v])
    orbit::MBD.CR3BPMSPeriodicOrbit = interpOrbit(targeter, familyFile, "Period", values[v], numNodes)
    println("Converged Orbit:\n\tIC:\t$(orbit.initialCondition)\n\tP:\t$(orbit.period) ($(orbit.period*getCharTime(systemData)/3600/24) d.)\n\tJC:\t$(getJacobiConstant(orbit))\n\tStab.:\t$(getStabilityIndex(orbit))\n")
    orbitArc::MBD.CR3BPArc = propagate(propagator, orbit.initialCondition, [0, orbit.period], dynamicsModel)
    inertialStates::Vector{Vector{Float64}} = rotatingToPrimaryInertial(dynamicsModel, 1, orbitArc.states, orbitArc.times)
    exportCR3BPOrbit(orbit, mf, Symbol("orbit"*string(v)))
    exportInertialTrajectory(inertialStates, orbitArc.times, mf, Symbol("inertOrbit"*string(v)))
end
MATLAB.close(mf)

println()
end
