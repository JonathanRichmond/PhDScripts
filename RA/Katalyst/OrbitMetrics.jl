"""
Script for computing orbit metrics

Author: Jonathan Richmond
C: 6/30/25
"""
module OrbMetrics
println()

using MBD, MATLAB

include("../../CR3BPTargeters/PlanarPerpJC.jl")
include("../../CR3BPTargeters/PlanarPerpP.jl")
include("../../CR3BPTargeters/SpatialPerpVy.jl")
include("../../Utilities/Export.jl")

systemData = MBD.CR3BPSystemData("Earth", "Moon")
dynamicsModel = MBD.CR3BPDynamicsModel(systemData)

propagator = MBD.Propagator()
targeter = PlanarPerpPTargeter(dynamicsModel)
familyFile::String = "FamilyData/CR3BPEMCyclers.csv"

values::Vector{Float64} = [5.0, 7.25, 9.5]
orbits::Vector{MBD.CR3BPPeriodicOrbit} = []

mf = MATLAB.MatFile("RA/Katalyst/OrbitMetrics.mat", "w")
for v::Int64 in 1:length(values)
    orbit::MBD.CR3BPPeriodicOrbit = interpOrbit(targeter, familyFile, "Period", values[v])
    println("Converged Orbit:\n\tIC:\t$(orbit.initialCondition)\n\tP:\t$(orbit.period) ($(orbit.period*getCharTime(systemData)/3600/24) d.)\n\tJC:\t$(getJacobiConstant(orbit))\n\tStab.:\t$(getStabilityIndex(orbit))\n")
    orbitArc::MBD.CR3BPArc = propagate(propagator, orbit.initialCondition, [0, orbit.period], dynamicsModel)
    inertialStates::Vector{Vector{Float64}} = rotatingToPrimaryInertial(dynamicsModel, 1, orbitArc.states, orbitArc.times)
    exportCR3BPOrbit(orbit, mf, Symbol("orbit"*string(v)))
    exportInertialTrajectory(inertialStates, orbitArc.times, mf, Symbol("inertOrbit"*string(v)))
end
MATLAB.close(mf)

println()
end
