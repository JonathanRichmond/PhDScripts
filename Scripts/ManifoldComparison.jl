"""
Script for comparing CR3BP and BCR4BP invariant and pseudo-manifolds

Author: Jonathan Richmond
C: 6/10/25
U: 7/10/25
"""
module ManComp
println()

using MBD, DifferentialEquations, LinearAlgebra, Logging, MATLAB
using Base.Threads

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
p::Int64, q::Int64 = 2, 1
compOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", getSynodicPeriod(dynamicsModel)*q/p; choiceIndex = 1)
println("Converged $p:$q CR3BP Orbit:\n\tIC:\t$(compOrbit.initialCondition)\n\tP:\t$(compOrbit.period)\n\tStab.:\t$(getStabilityIndex(compOrbit))\n\tJC:\t$(getJacobiConstant(compOrbit))\n")
# orbitArc::MBD.CR3BPArc = propagate(propagator, compOrbit.initialCondition, [0, compOrbit.period/2], CR3BPDynamicsModel)
# halfState::Vector{Float64} = getStateByIndex(orbitArc, -1)
# compOrbit = MBD.CR3BPPeriodicOrbit(CR3BPDynamicsModel, halfState, compOrbit.period, Matrix{Float64}(compOrbit.monodromy))
q0JumpCheck = MBD.BoundingBoxJumpCheck("Initial State", [0.9 1.0; -1.5 -1.0])
orbit::MBD.BCR4BP12PeriodicOrbit = getResonantOrbit(targeter, compOrbit, 4*q, 0.0, p, q, q0JumpCheck)

propTime::Float64 = pi*10
R_H::Float64 = get41CharLength(systemData)*(systemData.primaryData[1].mass/(3*(systemData.primaryData[3].mass+systemData.primaryData[1].mass)))^(1/3)/get12CharLength(systemData)
HillsSphereEventBCR4BP = DifferentialEquations.ContinuousCallback(MBD.p1BCR4BP12DistanceCondition, terminateAffect!)
HillsSphereEventCR3BP = DifferentialEquations.ContinuousCallback(MBD.p1CR3BPDistanceCondition, terminateAffect!)

BCR4BPPosManifold::MBD.BCR4BP12Manifold = getManifoldByArclength(orbit, "Unstable", "Positive", 25/get12CharLength(systemData), 200*p)
BCR4BPPosManifold.TOF = propTime
BCR4BPPosManifoldArcs::Vector{MBD.BCR4BP12ManifoldArc} = stopCrashes(BCR4BPPosManifold)
BCR4BPNegManifold::MBD.BCR4BP12Manifold = getManifoldByArclength(orbit, "Unstable", "Negative", 25/get12CharLength(systemData), 200*p)
BCR4BPNegManifold.TOF = propTime
BCR4BPNegManifoldArcs::Vector{MBD.BCR4BP12ManifoldArc} = stopCrashes(BCR4BPNegManifold)
BCR4BPManifoldArcs::Vector{MBD.BCR4BP12ManifoldArc} = append!(BCR4BPPosManifoldArcs, BCR4BPNegManifoldArcs)
HillsTOFBCR4BP::Vector{Float64} = Vector{Float64}(undef, length(BCR4BPManifoldArcs))
vEscBCR4BP::Vector{Float64} = Vector{Float64}(undef, length(BCR4BPManifoldArcs))
Threads.@threads for a::Int64 in 1:length(BCR4BPManifoldArcs)
    arc::MBD.BCR4BP12Arc = propagateWithEvent(propagator, HillsSphereEventBCR4BP, real(BCR4BPManifoldArcs[a].initialCondition), [0, BCR4BPManifoldArcs[a].TOF], dynamicsModel, [R_H])
    BCR4BPManifoldArcs[a].TOF = getTimeByIndex(arc, -1)
    HillsTOFBCR4BP[a] = BCR4BPManifoldArcs[a].TOF
    (HillsStates::Vector{Vector{Float64}}, ~) = rotating12ToPrimaryInertial(dynamicsModel, 1, [getStateByIndex(arc, -1)], [BCR4BPManifoldArcs[a].TOF])
    vEscBCR4BP[a] = LinearAlgebra.norm(HillsStates[1][4:6])
end

CR3BPPosManifold::MBD.CR3BPManifold = getManifoldByArclength(compOrbit, "Unstable", "Positive", 25/getCharLength(CR3BPSystemData), 200)
CR3BPPosManifold.TOF = propTime
CR3BPPosManifoldArcs::Vector{MBD.CR3BPManifoldArc} = stopCrashes(CR3BPPosManifold)
CR3BPNegManifold::MBD.CR3BPManifold = getManifoldByArclength(compOrbit, "Unstable", "Negative", 25/getCharLength(CR3BPSystemData), 200)
CR3BPNegManifold.TOF = propTime
CR3BPNegManifoldArcs::Vector{MBD.CR3BPManifoldArc} = stopCrashes(CR3BPNegManifold)
CR3BPManifoldArcs::Vector{MBD.CR3BPManifoldArc} = append!(CR3BPPosManifoldArcs, CR3BPNegManifoldArcs)
HillsTOFCR3BP::Vector{Float64} = Vector{Float64}(undef, length(CR3BPManifoldArcs))
vEscCR3BP::Vector{Float64} = Vector{Float64}(undef, length(CR3BPManifoldArcs))
Threads.@threads for a::Int64 in 1:length(CR3BPManifoldArcs)
    arc::MBD.CR3BPArc = propagateWithEvent(propagator, HillsSphereEventCR3BP, real(CR3BPManifoldArcs[a].initialCondition), [0, CR3BPManifoldArcs[a].TOF], CR3BPDynamicsModel, [R_H])
    CR3BPManifoldArcs[a].TOF = getTimeByIndex(arc, -1)
    HillsTOFCR3BP[a] = CR3BPManifoldArcs[a].TOF
    HillsStates::Vector{Vector{Float64}} = rotatingToPrimaryInertial(CR3BPDynamicsModel, 1, [getStateByIndex(arc, -1)], [CR3BPManifoldArcs[a].TOF])
    vEscCR3BP[a] = LinearAlgebra.norm(HillsStates[1][4:6])
end

pseudoManifold0::MBD.BCR4BPPseudoManifold = getPseudoManifoldByArclength(compOrbit, dynamicsModel, 0.0, 50)
pseudoManifold0.TOF = propTime
theta40::Vector{Float64} = collect(range(0, 2*pi, 21))
pseudoManifoldArcs::Vector{MBD.BCR4BPPseudoManifoldArc} = []
for t::Float64 in theta40[1:end-1]
    pseudoManifold::MBD.BCR4BPPseudoManifold = getPseudoManifoldByArclength(compOrbit, dynamicsModel, t, 50)
    pseudoManifold.TOF = pseudoManifold0.TOF
    pseudoManifoldArcsSub::Vector{MBD.BCR4BPPseudoManifoldArc} = stopCrashes(pseudoManifold)
    append!(pseudoManifoldArcs, pseudoManifoldArcsSub)
end
HillsTOFPseudo::Vector{Float64} = Vector{Float64}(undef, length(pseudoManifoldArcs))
vEscPseudo::Vector{Float64} = Vector{Float64}(undef, length(pseudoManifoldArcs))
Threads.@threads for a::Int64 in 1:length(pseudoManifoldArcs)
    arc::MBD.BCR4BP12Arc = propagateWithEvent(propagator, HillsSphereEventBCR4BP, real(pseudoManifoldArcs[a].initialCondition), [0, pseudoManifoldArcs[a].TOF], dynamicsModel, [R_H])
    pseudoManifoldArcs[a].TOF = getTimeByIndex(arc, -1)
    HillsTOFPseudo[a] = pseudoManifoldArcs[a].TOF
    (HillsStates::Vector{Vector{Float64}}, ~) = rotating12ToPrimaryInertial(dynamicsModel, 1, [getStateByIndex(arc, -1)], [pseudoManifoldArcs[a].TOF])
    vEscPseudo[a] = LinearAlgebra.norm(HillsStates[1][4:6])
end

mf = MATLAB.MatFile("Output/ManifoldComparison.mat", "w")
exportCR3BPOrbit(compOrbit, mf, :CR3BPCompOrbit)
exportBCR4BP12Orbit(orbit, mf, :BCR4BPOrbit)
exportBCR4BP12Manifold(BCR4BPPosManifold, BCR4BPManifoldArcs, mf, :BCR4BPManifold)
exportCR3BPManifold(CR3BPPosManifold, CR3BPManifoldArcs, mf, :CR3BPManifold)
exportPseudoManifold(pseudoManifold0, pseudoManifoldArcs, mf, :pseudoManifold)
MATLAB.put_variable(mf, :HillsTOFBCR4BP, HillsTOFBCR4BP)
MATLAB.put_variable(mf, :vEscBCR4BP, vEscBCR4BP)
MATLAB.put_variable(mf, :HillsTOFCR3BP, HillsTOFCR3BP)
MATLAB.put_variable(mf, :vEscCR3BP, vEscCR3BP)
MATLAB.put_variable(mf, :HillsTOFPseudo, HillsTOFPseudo)
MATLAB.put_variable(mf, :vEscPseudo, vEscPseudo)
MATLAB.close(mf)

println()
end
