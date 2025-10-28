"""
Script for analyzing escape parameters and characteristics

Author: Jonathan Richmond
C: 9/22/25
U: 10/24/25
"""
# module EscAn
println("Running EscapeAnalysis.jl...\n")

using MBD, Logging, MATLAB, SPICE

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../BCR4BPTargeters/PlanarPerpP.jl")
include("../CR3BPTargeters/PlanarPerpJC.jl")
include("../Utilities/Export.jl")

mf = MATLAB.MatFile("Output/EscapeAnalysis.mat", "w")

SPICE.furnsh("SPICEKernels/naif0012.tls", "SPICEKernels/de430.bsp", "SPICEKernels/de440.bsp")

ESystemData = MBD.KSystemData("Earth")
SSystemData = MBD.KSystemData("Sun")
EMSystemData = MBD.CR3BPSystemData("Earth", "Moon")
SB1SystemData = MBD.CR3BPSystemData("Sun", "Earth_Barycenter")
EMSSystemData = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
EDynamicsModel = MBD.KDynamicsModel(ESystemData)
SDynamicsModel = MBD.KDynamicsModel(SSystemData)
EMDynamicsModel = MBD.CR3BPDynamicsModel(EMSystemData)
SB1CR3BPDynamicsModel = MBD.CR3BPDynamicsModel(SB1SystemData)
EMSDynamicsModel = MBD.BCR4BP12DynamicsModel(EMSSystemData)
SB1DynamicsModel = MBD.BCR4BP41DynamicsModel(EMSSystemData)
Earth::MBD.BodyData, Moon::MBD.BodyData, Sun::MBD.BodyData = EMSSystemData.primaryData[1], EMSSystemData.primaryData[2], EMSSystemData.primaryData[3]

propagator = MBD.Propagator()
EMCR3BPTargeter = PlanarPerpJCTargeter(EMDynamicsModel)
SB1CR3BPTargeter = PlanarPerpJCTargeter(SB1CR3BPDynamicsModel)
BCR4BPTargeter = PlanarPerpP12Targeter(EMSDynamicsModel)

familyFile::String = "FamilyData/CR3BPEML1Lyapunovs.csv"
p::Int64, q::Int64 = 2, 1
compOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(EMCR3BPTargeter, familyFile, "Period", getSynodicPeriod(EMSDynamicsModel)*q/p; choiceIndex = 1)
println("Converged $p:$q CR3BP Orbit:\n\tIC:\t$(compOrbit.initialCondition)\n\tP:\t$(compOrbit.period)\n\tStab.:\t$(getStabilityIndex(compOrbit))\n\tJC:\t$(getJacobiConstant(compOrbit))\n")
# q0JumpCheck = MBD.BoundingBoxJumpCheck("Node 1 State", [0.9 1.0; 1.0 2.0])
q0JumpCheck = MBD.BoundingBoxJumpCheck("Node 1 State", [0.9 1.0; -1.0 0])
orbit::MBD.BCR4BP12PeriodicOrbit = getResonantOrbit(BCR4BPTargeter, compOrbit, 4*q, 0.0, p, q, q0JumpCheck)

L2_EM::Vector{Float64} = getEquilibriumPoint(EMDynamicsModel, 2)
L2_SB1::Vector{Float64} = getEquilibriumPoint(SB1CR3BPDynamicsModel, 2)
E2_SB1::Vector{Float64} = getInstantaneousEquilibriumPoint(SB1DynamicsModel, 2, 0.0)
JC_EM_L2::Float64 = getJacobiConstant(EMDynamicsModel, append!(L2_EM, 0, 0, 0))
JC_SB1_L2::Float64 = getJacobiConstant(SB1CR3BPDynamicsModel, append!(L2_SB1, 0, 0, 0))
H_SB1_E2::Float64 = getHamiltonian(SB1DynamicsModel, append!(E2_SB1, [0.0, 0.0, 0.0, 0.0]))

propTime::Float64 = pi*15
R_H::Float64 = get41CharLength(EMSSystemData)*(EMSSystemData.primaryData[1].mass/(3*(EMSSystemData.primaryData[3].mass+EMSSystemData.primaryData[1].mass)))^(1/3)/get12CharLength(EMSSystemData)
sphereEventCR3BP = DifferentialEquations.ContinuousCallback(MBD.p1CR3BPDistanceCondition, terminateAffect!)
P2SphereEventCR3BP = DifferentialEquations.ContinuousCallback(MBD.p2DistanceCondition, terminateAffect!)
sphereEventBCR4BP = DifferentialEquations.ContinuousCallback(MBD.b1BCR4BP12DistanceCondition, terminateAffect!)

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
numArcs::Int64 = length(manifoldArcs)
r::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, numArcs)
t::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, numArcs)
a_E::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, numArcs)
a_S::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, numArcs)
e_E::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, numArcs)
e_S::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, numArcs)
quad_SB1::Vector{Vector{Int64}} = Vector{Vector{Int64}}(undef, numArcs)
JC_EM::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, numArcs)
JC_SB1::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, numArcs)
H_EM::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, numArcs)
vEsc::Vector{Float64} = zeros(Float64, numArcs)
hz::Vector{Float64} = zeros(Float64, numArcs)
for m::Int64 in 1:numArcs
    arcEsc::MBD.BCR4BP12Arc = propagateWithEvent(propagator, sphereEventBCR4BP, real(manifoldArcs[m].initialCondition), [0, manifoldArcs[m].TOF], EMSDynamicsModel, [750000/get12CharLength(EMSSystemData)])
    arcHills::MBD.BCR4BP12Arc = propagateWithEvent(propagator, sphereEventBCR4BP, real(manifoldArcs[m].initialCondition), [0, manifoldArcs[m].TOF], EMSDynamicsModel, [R_H])
    arcTime::MBD.BCR4BP12Arc = propagate(propagator, real(manifoldArcs[m].initialCondition), [0, manifoldArcs[m].TOF], EMSDynamicsModel)
    manifoldArcs[m].TOF = getTimeByIndex(arcTime, -1)
    arcHills_EarthEJ2000::Vector{Vector{Float64}} = rotating12ToPrimaryEcliptic(EMSDynamicsModel, "ECLIPJ2000", 1, "JAN 1 2030", arcHills.states, arcHills.times)[1]
    arcTime_EarthEJ2000::Vector{Vector{Float64}} = rotating12ToPrimaryEcliptic(EMSDynamicsModel, "ECLIPJ2000", 1, "JAN 1 2030", arcTime.states, arcTime.times)[1]
    arcTime_SunEJ2000::Vector{Vector{Float64}} = rotating12ToPrimaryEcliptic(EMSDynamicsModel, "ECLIPJ2000", 4, "JAN 1 2030", arcTime.states, arcTime.times)[1]
    (arcTime_SB1::Vector{Vector{Float64}}, t_SB1::Vector{Float64}) = rotating12ToRotating41(EMSDynamicsModel, arcTime.states, arcTime.times)
    vEsc[m] = LinearAlgebra.norm(arcHills_EarthEJ2000[end][4:6])
    escState::Vector{Float64} = getStateByIndex(arcEsc, -1)
    hz[m] = (abs(getTimeByIndex(arcEsc, -1)) < abs(manifoldArcs[m].TOF)) ? escState[1]*escState[5]-escState[2]*escState[4] : NaN
    numStates::Int64 = getStateCount(arcTime)
    r[m] = zeros(Float64, numStates)
    t[m] = zeros(Float64, numStates)
    a_E[m] = zeros(Float64, numStates)
    a_S[m] = zeros(Float64, numStates)
    e_E[m] = zeros(Float64, numStates)
    e_S[m] = zeros(Float64, numStates)
    x_SB1 = zeros(Float64, numStates)
    y_SB1 = zeros(Float64, numStates)
    z_SB1 = zeros(Float64, numStates)
    xdot_SB1 = zeros(Float64, numStates)
    ydot_SB1 = zeros(Float64, numStates)
    zdot_SB1 = zeros(Float64, numStates)
    theta2_SB1 = zeros(Float64, numStates)
    H_SB1 = zeros(Float64, numStates)
    quad_SB1[m] = zeros(Int64, numStates)
    JC_EM[m] = zeros(Float64, numStates)
    JC_SB1[m] = zeros(Float64, numStates)
    H_EM[m] = zeros(Float64, numStates)
    for s::Int64 in 1:numStates
        r[m][s] = getExcursion(EMSDynamicsModel, 1, getStateByIndex(arcTime, s))
        t[m][s] = getTimeByIndex(arcTime, s)
        elementStates_E = getOrbitalElements(EDynamicsModel, append!(arcTime_EarthEJ2000[s][1:3]*get12CharLength(EMSDynamicsModel), arcTime_EarthEJ2000[s][4:6]*get12CharLength(EMSDynamicsModel)./get12CharTime(EMSDynamicsModel)))
        elementStates_S = getOrbitalElements(SDynamicsModel, append!(arcTime_SunEJ2000[s][1:3]*get12CharLength(EMSDynamicsModel), arcTime_SunEJ2000[s][4:6]*get12CharLength(EMSDynamicsModel)./get12CharTime(EMSDynamicsModel)))
        a_E[m][s] = elementStates_E[1]
        a_S[m][s] = elementStates_S[1]
        e_E[m][s] = elementStates_E[2]
        e_S[m][s] = elementStates_S[2]
        x_SB1[s] = arcTime_SB1[s][1]
        y_SB1[s] = arcTime_SB1[s][2]
        z_SB1[s] = arcTime_SB1[s][3]
        xdot_SB1[s] = arcTime_SB1[s][4]
        ydot_SB1[s] = arcTime_SB1[s][5]
        zdot_SB1[s] = arcTime_SB1[s][6]
        theta2_SB1[s] = arcTime_SB1[s][7]
        H_SB1[s] = getHamiltonian(SB1DynamicsModel, arcTime_SB1[s])
        quad_SB1[m][s] = getQuadrant(SB1DynamicsModel, arcTime_SB1[s])
        JC_EM[m][s] = getJacobiConstant(EMDynamicsModel, getStateByIndex(arcTime, s)[1:6])
        JC_SB1[m][s] = getJacobiConstant(SB1CR3BPDynamicsModel, arcTime_SB1[s])
        H_EM[m][s] = getHamiltonian(EMSDynamicsModel, getStateByIndex(arcTime, s)[1:7])
    end
    exportBCR4BP41Trajectory(x_SB1, y_SB1, z_SB1, xdot_SB1, ydot_SB1, zdot_SB1, theta2_SB1, t_SB1, H_SB1, t_SB1[end], mf, Symbol("SB1Trajectory", m))
end

L1EscapeOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(SB1CR3BPTargeter, "FamilyData/CR3BPSB1L1Lyapunovs.csv", "JC", 3.0008; choiceIndex = 1)
L1StableManifold::MBD.CR3BPManifold = getManifoldByArclength(L1EscapeOrbit, "Stable", "Positive", 1000/getCharLength(SB1SystemData), 200)
L1StableManifold.TOF = pi
L1StableManifoldArcs::Vector{MBD.CR3BPManifoldArc} = stopCrashes(L1StableManifold)
Threads.@threads for a::Int64 in 1:length(L1StableManifoldArcs)
    arc::MBD.CR3BPArc = propagate(propagator, real(L1StableManifoldArcs[a].initialCondition), [0, -L1StableManifoldArcs[a].TOF], SB1CR3BPDynamicsModel)
    L1StableManifoldArcs[a].TOF = getTimeByIndex(arc, -1)
end
exportCR3BPManifold(L1StableManifold, L1StableManifoldArcs, mf, :L1StableMan)

L2EscapeOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(SB1CR3BPTargeter, "FamilyData/CR3BPSB1L2Lyapunovs.csv", "JC", 3.0008; choiceIndex = 1)
L2StableManifold::MBD.CR3BPManifold = getManifoldByArclength(L2EscapeOrbit, "Stable", "Negative", 1000/getCharLength(SB1SystemData), 200)
L2StableManifold.TOF = pi
L2StableManifoldArcs::Vector{MBD.CR3BPManifoldArc} = stopCrashes(L2StableManifold)
Threads.@threads for a::Int64 in 1:length(L2StableManifoldArcs)
    arc::MBD.CR3BPArc = propagate(propagator, real(L2StableManifoldArcs[a].initialCondition), [0, -L2StableManifoldArcs[a].TOF], SB1CR3BPDynamicsModel)
    L2StableManifoldArcs[a].TOF = getTimeByIndex(arc, -1)
end
exportCR3BPManifold(L2StableManifold, L2StableManifoldArcs, mf, :L2StableMan)

exportCR3BPOrbit(compOrbit, mf, :CR3BPCompOrbit)
exportBCR4BP12Orbit(orbit, mf, :BCR4BPOrbit)
for t::Int64 in 1:length(manifoldArcs)
    exportBCR4BP12Trajectory(manifoldArcs[t].initialCondition..., manifoldArcs[t].TOF, EMSDynamicsModel, mf, Symbol("Trajectory", t))
end
MATLAB.put_variable(mf, :d, r)
MATLAB.put_variable(mf, :t, t)
MATLAB.put_variable(mf, :a_osc_E, a_E)
MATLAB.put_variable(mf, :a_osc_S, a_S)
MATLAB.put_variable(mf, :e_osc_E, e_E)
MATLAB.put_variable(mf, :e_osc_S, e_S)
MATLAB.put_variable(mf, :quad_SB1, quad_SB1)
MATLAB.put_variable(mf, :JC_EM, JC_EM)
MATLAB.put_variable(mf, :JC_SB1, JC_SB1)
MATLAB.put_variable(mf, :v_esc, vEsc)
MATLAB.put_variable(mf, :h_z, hz)
MATLAB.put_variable(mf, :JC_EM_L2, JC_EM_L2)
MATLAB.put_variable(mf, :JC_SB1_L2, JC_SB1_L2)
MATLAB.put_variable(mf, :H_EM, H_EM)
MATLAB.put_variable(mf, :H_E2, H_SB1_E2)
MATLAB.close(mf)

SPICE.kclear()

println()
# end
