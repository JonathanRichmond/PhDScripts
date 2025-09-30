"""
Script for analyzing escape parameters and characteristics

Author: Jonathan Richmond
C: 9/22/25
U: 9/30/25
"""
module EscAn
println()

using MBD, Logging, MATLAB, SPICE

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../BCR4BPTargeters/PlanarPerpP.jl")
include("../CR3BPTargeters/PlanarPerpJC.jl")
include("../Utilities/Export.jl")

SPICE.furnsh("SPICEKernels/naif0012.tls", "SPICEKernels/de430.bsp", "SPICEKernels/de440.bsp")

ESystemData = MBD.KSystemData("Earth")
SSystemData = MBD.KSystemData("Sun")
EMSystemData = MBD.CR3BPSystemData("Earth", "Moon")
EMSSystemData = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
EDynamicsModel = MBD.KDynamicsModel(ESystemData)
SDynamicsModel = MBD.KDynamicsModel(SSystemData)
EMDynamicsModel = MBD.CR3BPDynamicsModel(EMSystemData)
EMSDynamicsModel = MBD.BCR4BP12DynamicsModel(EMSSystemData)
SB1DynamicsModel = MBD.BCR4BP41DynamicsModel(EMSSystemData)
Earth::MBD.BodyData, Moon::MBD.BodyData, Sun::MBD.BodyData = EMSSystemData.primaryData[1], EMSSystemData.primaryData[2], EMSSystemData.primaryData[3]

propagator = MBD.Propagator()
CR3BPTargeter = PlanarPerpJCTargeter(EMDynamicsModel)
BCR4BPTargeter = PlanarPerpP12Targeter(EMSDynamicsModel)

familyFile::String = "FamilyData/CR3BPEML1Lyapunovs.csv"
p::Int64, q::Int64 = 2, 1
compOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, familyFile, "Period", getSynodicPeriod(EMSDynamicsModel)*q/p; choiceIndex = 1)
println("Converged $p:$q CR3BP Orbit:\n\tIC:\t$(compOrbit.initialCondition)\n\tP:\t$(compOrbit.period)\n\tStab.:\t$(getStabilityIndex(compOrbit))\n\tJC:\t$(getJacobiConstant(compOrbit))\n")
q0JumpCheck = MBD.BoundingBoxJumpCheck("Initial State", [0.9 1.0; -1.5 -1.0])
orbit::MBD.BCR4BP12PeriodicOrbit = getResonantOrbit(BCR4BPTargeter, compOrbit, 4*q, 0.0, p, q, q0JumpCheck)

E1_SB1::Vector{Float64} = getInstantaneousEquilibriumPoint(SB1DynamicsModel, 2, 0.0)
H_SB1_E1::Float64 = getHamiltonian(SB1DynamicsModel, append!(E1_SB1, [0.0, 0.0, 0.0, 0.0]))

propTime::Float64 = pi*15
R_H::Float64 = get41CharLength(EMSSystemData)*(EMSSystemData.primaryData[1].mass/(3*(EMSSystemData.primaryData[3].mass+EMSSystemData.primaryData[1].mass)))^(1/3)/get12CharLength(EMSSystemData)
HillsSphereEventCR3BP = DifferentialEquations.ContinuousCallback(MBD.p1CR3BPDistanceCondition, terminateAffect!)
HillsSphereEventBCR4BP = DifferentialEquations.ContinuousCallback(MBD.p1BCR4BP12DistanceCondition, terminateAffect!)

# CR3BPPosManifold::MBD.CR3BPManifold = getManifoldByArclength(compOrbit, "Unstable", "Positive", 25/getCharLength(CR3BPSystemData), 200)
# CR3BPPosManifold.TOF = propTime
# CR3BPPosManifoldArcs::Vector{MBD.CR3BPManifoldArc} = stopCrashes(CR3BPPosManifold)
# CR3BPNegManifold::MBD.CR3BPManifold = getManifoldByArclength(compOrbit, "Unstable", "Negative", 25/getCharLength(CR3BPSystemData), 200)
# CR3BPNegManifold.TOF = propTime
# CR3BPNegManifoldArcs::Vector{MBD.CR3BPManifoldArc} = stopCrashes(CR3BPNegManifold)
# CR3BPManifoldArcs::Vector{MBD.CR3BPManifoldArc} = append!(CR3BPPosManifoldArcs, CR3BPNegManifoldArcs)
# HillsTOFCR3BP::Vector{Float64} = Vector{Float64}(undef, length(CR3BPManifoldArcs))
# vEscCR3BP::Vector{Float64} = Vector{Float64}(undef, length(CR3BPManifoldArcs))
# Threads.@threads for a::Int64 in 1:length(CR3BPManifoldArcs)
#     arc::MBD.CR3BPArc = propagateWithEvent(propagator, HillsSphereEventCR3BP, real(CR3BPManifoldArcs[a].initialCondition), [0, CR3BPManifoldArcs[a].TOF], CR3BPDynamicsModel, [R_H])
#     CR3BPManifoldArcs[a].TOF = getTimeByIndex(arc, -1)
#     HillsTOFCR3BP[a] = CR3BPManifoldArcs[a].TOF
#     HillsStates::Vector{Vector{Float64}} = rotatingToPrimaryInertial(CR3BPDynamicsModel, 1, [getStateByIndex(arc, -1)], [CR3BPManifoldArcs[a].TOF])
#     vEscCR3BP[a] = LinearAlgebra.norm(HillsStates[1][4:6])
# end

# BCR4BPPosManifold::MBD.BCR4BP12Manifold = getManifoldByArclength(orbit, "Unstable", "Positive", 25/get12CharLength(systemData), 200*p)
# BCR4BPPosManifold.TOF = propTime
# BCR4BPPosManifoldArcs::Vector{MBD.BCR4BP12ManifoldArc} = stopCrashes(BCR4BPPosManifold)
# BCR4BPNegManifold::MBD.BCR4BP12Manifold = getManifoldByArclength(orbit, "Unstable", "Negative", 25/get12CharLength(systemData), 200*p)
# BCR4BPNegManifold.TOF = propTime
# BCR4BPNegManifoldArcs::Vector{MBD.BCR4BP12ManifoldArc} = stopCrashes(BCR4BPNegManifold)
# BCR4BPManifoldArcs::Vector{MBD.BCR4BP12ManifoldArc} = append!(BCR4BPPosManifoldArcs, BCR4BPNegManifoldArcs)
# HillsTOFBCR4BP::Vector{Float64} = Vector{Float64}(undef, length(BCR4BPManifoldArcs))
# vEscBCR4BP::Vector{Float64} = Vector{Float64}(undef, length(BCR4BPManifoldArcs))
# Threads.@threads for a::Int64 in 1:length(BCR4BPManifoldArcs)
#     arc::MBD.BCR4BP12Arc = propagateWithEvent(propagator, HillsSphereEventBCR4BP, real(BCR4BPManifoldArcs[a].initialCondition), [0, BCR4BPManifoldArcs[a].TOF], dynamicsModel, [R_H])
#     BCR4BPManifoldArcs[a].TOF = getTimeByIndex(arc, -1)
#     HillsTOFBCR4BP[a] = BCR4BPManifoldArcs[a].TOF
#     (HillsStates::Vector{Vector{Float64}}, ~) = rotating12ToPrimaryInertial(dynamicsModel, 1, [getStateByIndex(arc, -1)], [BCR4BPManifoldArcs[a].TOF])
#     vEscBCR4BP[a] = LinearAlgebra.norm(HillsStates[1][4:6])
# end

pseudoManifold::MBD.BCR4BPPseudoManifold = getPseudoManifoldByArclength(compOrbit, EMSDynamicsModel, 0.0, 50)
pseudoManifold.TOF = propTime
pseudoManifoldArcs::Vector{MBD.BCR4BPPseudoManifoldArc} = stopCrashes(pseudoManifold)
manifoldArcs::Vector{MBD.BCR4BPPseudoManifoldArc} = pseudoManifoldArcs[4:10]
q0_SB1::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(manifoldArcs))
r::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(manifoldArcs))
t::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(manifoldArcs))
a_E::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(manifoldArcs))
a_S::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(manifoldArcs))
e_E::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(manifoldArcs))
e_S::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(manifoldArcs))
H_SB1::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(manifoldArcs))
quad_SB1::Vector{Vector{Int64}} = Vector{Vector{Int64}}(undef, length(manifoldArcs))
vEsc::Vector{Float64} = zeros(Float64, length(manifoldArcs))
for m::Int64 in 1:length(manifoldArcs)
    arcHills::MBD.BCR4BP12Arc = propagateWithEvent(propagator, HillsSphereEventBCR4BP, real(manifoldArcs[m].initialCondition), [0, manifoldArcs[m].TOF], EMSDynamicsModel, [R_H])
    arcTime::MBD.BCR4BP12Arc = propagate(propagator, real(manifoldArcs[m].initialCondition), [0, manifoldArcs[m].TOF], EMSDynamicsModel)
    manifoldArcs[m].TOF = getTimeByIndex(arcTime, -1)
    arcHills_EarthEJ2000::Vector{Vector{Float64}} = rotating12ToPrimaryEcliptic(EMSDynamicsModel, "ECLIPJ2000", 1, "JAN 1 2030", arcHills.states, arcHills.times)[1]
    arcTime_EarthEJ2000::Vector{Vector{Float64}} = rotating12ToPrimaryEcliptic(EMSDynamicsModel, "ECLIPJ2000", 1, "JAN 1 2030", arcTime.states, arcTime.times)[1]
    arcTime_SunEJ2000::Vector{Vector{Float64}} = rotating12ToPrimaryEcliptic(EMSDynamicsModel, "ECLIPJ2000", 4, "JAN 1 2030", arcTime.states, arcTime.times)[1]
    arcTime_SB1::Vector{Vector{Float64}} = rotating12ToRotating41(EMSDynamicsModel, arcTime.states, arcTime.times)[1]
    q0_SB1[m] = arcTime_SB1[1]
    vEsc[m] = LinearAlgebra.norm(arcHills_EarthEJ2000[end][4:6])
    r[m] = zeros(Float64, getStateCount(arcTime))
    t[m] = zeros(Float64, getStateCount(arcTime))
    a_E[m] = zeros(Float64, getStateCount(arcTime))
    a_S[m] = zeros(Float64, getStateCount(arcTime))
    e_E[m] = zeros(Float64, getStateCount(arcTime))
    e_S[m] = zeros(Float64, getStateCount(arcTime))
    H_SB1[m] = zeros(Float64, getStateCount(arcTime))
    quad_SB1[m] = zeros(Int64, getStateCount(arcTime))
    for s::Int64 in 1:getStateCount(arcTime)
        r[m][s] = getExcursion(EMSDynamicsModel, 1, getStateByIndex(arcTime, s))
        t[m][s] = getTimeByIndex(arcTime, s)
        elementStates_E = getOrbitalElements(EDynamicsModel, append!(arcTime_EarthEJ2000[s][1:3]*get12CharLength(EMSDynamicsModel), arcTime_EarthEJ2000[s][4:6]*get12CharLength(EMSDynamicsModel)./get12CharTime(EMSDynamicsModel)))
        elementStates_S = getOrbitalElements(SDynamicsModel, append!(arcTime_SunEJ2000[s][1:3]*get12CharLength(EMSDynamicsModel), arcTime_SunEJ2000[s][4:6]*get12CharLength(EMSDynamicsModel)./get12CharTime(EMSDynamicsModel)))
        a_E[m][s] = elementStates_E[1]
        a_S[m][s] = elementStates_S[1]
        e_E[m][s] = elementStates_E[2]
        e_S[m][s] = elementStates_S[2]
        H_SB1[m][s] = getHamiltonian(SB1DynamicsModel, arcTime_SB1[s])
        quad_SB1[m][s] = getQuadrant(SB1DynamicsModel, arcTime_SB1[s])
    end
end

mf = MATLAB.MatFile("Output/EscapeAnalysis.mat", "w")
exportCR3BPOrbit(compOrbit, mf, :CR3BPCompOrbit)
exportBCR4BP12Orbit(orbit, mf, :BCR4BPOrbit)
for t::Int64 in 1:length(manifoldArcs)
    exportBCR4BP12Trajectory(manifoldArcs[t].initialCondition..., manifoldArcs[t].TOF, EMSDynamicsModel, mf, Symbol("Trajectory", t))
    exportBCR4BP41Trajectory(q0_SB1[t]..., manifoldArcs[t].TOF*get12CharTime(EMSSystemData)/get41CharTime(EMSSystemData), SB1DynamicsModel, mf, Symbol("SB1Trajectory", t))
end
MATLAB.put_variable(mf, :d, r)
MATLAB.put_variable(mf, :t, t)
MATLAB.put_variable(mf, :a_osc_E, a_E)
MATLAB.put_variable(mf, :a_osc_S, a_S)
MATLAB.put_variable(mf, :e_osc_E, e_E)
MATLAB.put_variable(mf, :e_osc_S, e_S)
MATLAB.put_variable(mf, :H_SB1, H_SB1)
MATLAB.put_variable(mf, :quad_SB1, quad_SB1)
MATLAB.put_variable(mf, :v_esc, vEsc)
MATLAB.put_variable(mf, :H_E1, H_SB1_E1)
MATLAB.close(mf)

SPICE.kclear()

println()
end
