"""
Script for computing visibility for constellations

Author: Jonathan Richmond
C: 7/21/25
U: 7/22/25
"""
module Const
println()

using MBD, Logging, MATLAB

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../../CR3BPTargeters/PlanarPerpJC.jl")
include("../../Utilities/Export.jl")

systemData = MBD.CR3BPSystemData("Earth", "Moon")
dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
BCR4BPSystemData = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
BCR4BPDynamicsModel = MBD.BCR4BP12DynamicsModel(BCR4BPSystemData)

targeter = PlanarPerpJCTargeter(dynamicsModel)
propagator = MBD.Propagator()
tSyn::Float64 = getSynodicPeriod(BCR4BPDynamicsModel)
tTarg::Float64 = 28.07*24*3600/getCharTime(systemData)

numSC::Int64 = 4
times::Vector{Float64} = collect(range(0, tSyn, 501))
constICs::Vector{Vector{Float64}} = []
constArcs::Vector{MBD.CR3BPArc} = []

p1::Int64, q1::Int64 = 1, 1
constOrbit1::MBD.CR3BPPeriodicOrbit = interpOrbit(targeter, "FamilyData/CR3BPEMResonant2_1Pros.csv", "Period", tTarg*q1/p1, choiceIndex = 2)
println("Converged Constellation Orbit 1:\n\tIC:\t$(constOrbit1.initialCondition)\n\tP:\t$(constOrbit1.period) ($(constOrbit1.period*getCharTime(systemData)/3600/24) d.)\n\tJC:\t$(getJacobiConstant(constOrbit1))\n\tStab.:\t$(getStabilityIndex(constOrbit1))\n")
xi01::Float64 = 0.0
for o::Int64 = 1:q1
    xi::Float64 = xi01+2*pi*(o-1)/q1
    orbitArc::MBD.CR3BPArc = propagate(propagator, constOrbit1.initialCondition, [0, xi*constOrbit1.period/(2*pi)], dynamicsModel)
    orbitIC::Vector{Float64} = getStateByIndex(orbitArc, -1)[1:6]
    push!(constICs, orbitIC)
    push!(constArcs, propagate(propagator, orbitIC, times, dynamicsModel))
end

p2::Int64, q2::Int64 = 1, 1
constOrbit2::MBD.CR3BPPeriodicOrbit = interpOrbit(targeter, "FamilyData/CR3BPEMResonant3_1Pros.csv", "Period", tTarg*q2/p2)
println("Converged Constellation Orbit 2:\n\tIC:\t$(constOrbit2.initialCondition)\n\tP:\t$(constOrbit2.period) ($(constOrbit2.period*getCharTime(systemData)/3600/24) d.)\n\tJC:\t$(getJacobiConstant(constOrbit2))\n\tStab.:\t$(getStabilityIndex(constOrbit2))\n")
xi02::Float64 = pi
for o::Int64 = 1:q2
    xi::Float64 = xi02+2*pi*(o-1)/q2
    orbitArc::MBD.CR3BPArc = propagate(propagator, constOrbit2.initialCondition, [0, xi*constOrbit2.period/(2*pi)], dynamicsModel)
    orbitIC::Vector{Float64} = getStateByIndex(orbitArc, -1)[1:6]
    push!(constICs, orbitIC)
    push!(constArcs, propagate(propagator, orbitIC, times, dynamicsModel))
end

p3::Int64, q3::Int64 = 1, 1
constOrbit3::MBD.CR3BPPeriodicOrbit = interpOrbit(targeter, "FamilyData/CR3BPEML1Lyapunovs.csv", "Period", tTarg*q3/p3)
println("Converged Constellation Orbit 3:\n\tIC:\t$(constOrbit3.initialCondition)\n\tP:\t$(constOrbit3.period) ($(constOrbit3.period*getCharTime(systemData)/3600/24) d.)\n\tJC:\t$(getJacobiConstant(constOrbit3))\n\tStab.:\t$(getStabilityIndex(constOrbit3))\n")
xi03::Float64 = 0.0
for o::Int64 = 1:q3
    xi::Float64 = xi03+2*pi*(o-1)/q3
    orbitArc::MBD.CR3BPArc = propagate(propagator, constOrbit3.initialCondition, [0, xi*constOrbit3.period/(2*pi)], dynamicsModel)
    orbitIC::Vector{Float64} = getStateByIndex(orbitArc, -1)[1:6]
    push!(constICs, orbitIC)
    push!(constArcs, propagate(propagator, orbitIC, times, dynamicsModel))
end

p4::Int64, q4::Int64 = 1, 1
constOrbit4::MBD.CR3BPPeriodicOrbit = interpOrbit(targeter, "FamilyData/CR3BPEML2Lyapunovs.csv", "Period", tTarg*q4/p4)
println("Converged Constellation Orbit 4:\n\tIC:\t$(constOrbit4.initialCondition)\n\tP:\t$(constOrbit4.period) ($(constOrbit4.period*getCharTime(systemData)/3600/24) d.)\n\tJC:\t$(getJacobiConstant(constOrbit4))\n\tStab.:\t$(getStabilityIndex(constOrbit4))\n")
xi04::Float64 = 0.0
for o::Int64 = 1:q4
    xi::Float64 = xi04+2*pi*(o-1)/q4
    orbitArc::MBD.CR3BPArc = propagate(propagator, constOrbit4.initialCondition, [0, xi*constOrbit4.period/(2*pi)], dynamicsModel)
    orbitIC::Vector{Float64} = getStateByIndex(orbitArc, -1)[1:6]
    push!(constICs, orbitIC)
    push!(constArcs, propagate(propagator, orbitIC, times, dynamicsModel))
end

mf = MATLAB.MatFile("RA/Katalyst/Constellation.mat", "w")
exportCR3BPOrbit(constOrbit1, mf, :ConstellationOrbit1)
exportCR3BPOrbit(constOrbit2, mf, :ConstellationOrbit2)
exportCR3BPOrbit(constOrbit3, mf, :ConstellationOrbit3)
exportCR3BPOrbit(constOrbit4, mf, :ConstellationOrbit4)
exportArrays(constICs, mf, :ObserverICs)
MATLAB.close(mf)

println()
end
