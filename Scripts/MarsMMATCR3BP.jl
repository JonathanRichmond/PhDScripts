"""
Script for computing CR3BP MMATs between Earth-Moon and Sun-Mars systems

Author: Jonathan LeFevre Richmond
C: 1/24/26
"""
# module MMMATCR3BP
println("Running MarsMMATCR3BP.jl...\n")

using MBD, DifferentialEquations, Logging, MATLAB, SPICE

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../CR3BPTargeters/MMATArrivalPhase.jl")
include("../CR3BPTargeters/SpatialPerpJC.jl")
include("../Utilities/Export.jl")

function findIntersection(targeter::MMATArrivalPhaseTargeter, KeplerianDynamicsModel::MBD.KDynamicsModel, arrivalDynamicsModel::MBD.CR3BPDynamicsModel, oe_arrivalBody::Vector{Float64}, initialEpoch::String, oe_dep_SoI::Vector{Float64}, oe_bridge_peri::Vector{Float64}, oe_arr_SoI_guess::Vector{Float64}, q_arr_SoI::Vector{Float64}, r_int::Float64, theta_bridge_int_guess::Float64, u_arr_guess::Float64, n::Int64, o::Int64)
    theta_arr_int_guess::Float64 = 2*o*pi+(-1)^o*acos((oe_arr_SoI_guess[1]*(1-oe_arr_SoI_guess[2]^2)/r_int-1)/oe_arr_SoI_guess[2])
    omega_arr_guess::Float64 = (u_arr_guess-(theta_arr_int_guess+n*pi) > 0) ? u_arr_guess-(theta_arr_int_guess+n*pi) : 2*pi+u_arr_guess-(theta_arr_int_guess+n*pi)
    sigma_guess::Float64 = atan(cos(oe_arrivalBody[3])*tan(oe_arrivalBody[4]))+omega_arr_guess+oe_arr_SoI_guess[6]
    delta::Float64 = atan(q_arr_SoI[2], q_arr_SoI[1])
    X_guess::Vector{Float64} = [theta_bridge_int_guess, sigma_guess-(oe_arrivalBody[4]+oe_arrivalBody[5]+oe_arrivalBody[6])-delta, theta_arr_int_guess]
    X::Vector{Float64} = correct(targeter, X_guess, oe_bridge_peri, arrivalDynamicsModel, initialEpoch, q_arr_SoI)
    Q_bridge_int_SI::Vector{Float64} = getCartesianState(KeplerianDynamicsModel, append!(oe_bridge_peri[1:5], X[1]))
    t_bridgeConic::Float64 = solveKeplersEquation(KeplerianDynamicsModel, append!(oe_bridge_peri[1:5], X[1]))
    q_arr_SoI_SI::Vector{Float64} = rotatingToSunEclipJ2000(arrivalDynamicsModel, initialEpoch, [q_arr_SoI], [X[2]])[1]
    Q_arr_SoI_SI::Vector{Float64} = append!(q_arr_SoI_SI[1:3].*getCharLength(arrivalDynamicsModel), q_arr_SoI_SI[4:6].*getCharLength(arrivalDynamicsModel)./getCharTime(arrivalDynamicsModel))
    oe_arr_SoI::Vector{Float64} = getOrbitalElements(KeplerianDynamicsModel, Q_arr_SoI_SI)
    oe_arr_int::Vector{Float64} = append!(oe_arr_SoI[1:5], X[3])
    Q_arr_int_SI::Vector{Float64} = getCartesianState(KeplerianDynamicsModel, oe_arr_int)
    Deltav_2::Float64 = LinearAlgebra.norm(Q_arr_int_SI[4:6]-Q_bridge_int_SI[4:6])
    t_arrConic::Float64 = solveKeplersEquation(KeplerianDynamicsModel, oe_arr_SoI)-solveKeplersEquation(KeplerianDynamicsModel, oe_arr_int)
    (t_arrConic < 0) && (t_arrConic += MBD.getPeriod(KeplerianDynamicsModel, oe_arr_SoI))

    return [t_bridgeConic, Deltav_2, oe_arr_int..., t_arrConic, X[2]]
end

function run_MarsMMATCR3BP()
    mf = MATLAB.MatFile("Output/MarsMMATCR3BP.mat", "w")

    SPICE.furnsh("SPICEKernels/naif0012.tls", "SPICEKernels/de430.bsp", "SPICEKernels/de440.bsp", "SPICEKernels/mar097.bsp")

    EMSystemData = MBD.CR3BPSystemData("Earth", "Moon")
    SESystemData = MBD.CR3BPSystemData("Sun", "Earth")
    SSystemData = MBD.KSystemData("Sun")
    SMSystemData = MBD.CR3BPSystemData("Sun", "Mars")
    EMDynamicsModel = MBD.CR3BPDynamicsModel(EMSystemData)
    SEDynamicsModel = MBD.CR3BPDynamicsModel(SESystemData)
    SDynamicsModel = MBD.KDynamicsModel(SSystemData)
    SMDynamicsModel = MBD.CR3BPDynamicsModel(SMSystemData)
    Earth::MBD.BodyData, Moon::MBD.BodyData = EMSystemData.primaryData[1], EMSystemData.primaryData[2]
    Sun::MBD.BodyData, Mars::MBD.BodyData = SMSystemData.primaryData[1], SMSystemData.primaryData[2]

    lstarEM::Float64 = getCharLength(EMSystemData)
    lstarSE::Float64 = getCharLength(SESystemData)
    lstarSM::Float64 = getCharLength(SMSystemData)
    tstarEM::Float64 = getCharTime(EMSystemData)
    tstarSE::Float64 = getCharTime(SESystemData)
    tstarSM::Float64 = getCharTime(SMSystemData)

    propagator = MBD.Propagator()
    momentumPropagator = MBD.Propagator(equationType = MBD.MOMENTUM)
    EMTargeter = SpatialPerpJCTargeter(EMDynamicsModel)
    SMTargeter = SpatialPerpJCTargeter(SMDynamicsModel)
    MMATTargeter = MMATArrivalPhaseTargeter(SDynamicsModel, sqrt(eps(Float64)))

    orbitDepartureEvent = DifferentialEquations.ContinuousCallback(momentumDifferenceConditionCR3BP, terminateAffect!)
    P1DistanceEvent = DifferentialEquations.ContinuousCallback(p1CR3BPDistanceCondition, terminateAffect!)
    P2DistanceEvent = DifferentialEquations.ContinuousCallback(p2CR3BPDistanceCondition, terminateAffect!)

    MoonSoI::Float64 = lstarSE*(Moon.mass/Sun.mass)^(2/5)/lstarEM
    EarthSoI::Float64 = 0.09877
    MarsSoI::Float64 = 0.05375
    EMMomentumDiff::Float64 = 1E-3
    SMMomentumDiff::Float64 = 1E-5

    EMJC::Float64 = 3.03
    SMJC::Float64 = 3.0001857
    initialEpoch::String = "Jan 1 2030"
    initialEpochTime::Float64 = SPICE.str2et(initialEpoch)
    days::Vector{Float64} = [0.0, 15.0]

    q_E_0::Vector{Float64} = getEphemerides(initialEpoch, [0.0], "Earth", "Sun", "ECLIPJ2000")[1][1]
    q_M_0::Vector{Float64} = getEphemerides(initialEpoch, [0.0], "Mars", "Sun", "ECLIPJ2000")[1][1]
    oe_E_0::Vector{Float64} = SPICE.oscltx(q_E_0, initialEpochTime, Sun.gravParam)
    oe_M_0::Vector{Float64} = SPICE.oscltx(q_M_0, initialEpochTime, Sun.gravParam)
    theta_E_0::Float64 = oe_E_0[6]
    i_M::Float64 = oe_M_0[3]
    Omega_M::Float64 = oe_M_0[4]
    omega_M::Float64 = oe_M_0[5]
    theta_M_0::Float64 = oe_M_0[6]

    coarseEMOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(EMTargeter, "FamilyData/CR3BPEML1Halos.csv", "JC", EMJC)
    coarseEMOrbit.initialCondition[3] *= -1 # Flip to northern halo
    EMSolution::MBD.CR3BPMultipleShooterProblem = correct(EMTargeter, coarseEMOrbit.initialCondition, [0, coarseEMOrbit.period], EMJC)
    EMOrbit = MBD.CR3BPPeriodicOrbit(EMDynamicsModel, EMSolution.nodes[1].state.data[1:6], getPeriod(EMTargeter, EMSolution), getMonodromy(EMTargeter, EMSolution))
    println("Converged Earth-Moon orbit:\n\tIC:\t$(EMOrbit.initialCondition)\n\tP:\t$(EMOrbit.period)\n\tJC:\t$(getJacobiConstant(EMOrbit))\n")
    EMPosUnstableManifold::MBD.CR3BPManifold = getManifoldByArclength(EMOrbit, "Unstable", "Positive", 25/lstarEM, 50)
    EMNegUnstableManifold::MBD.CR3BPManifold = getManifoldByArclength(EMOrbit, "Unstable", "Negative", 25/lstarEM, 50)
    EMPosUnstableManifold.TOF, EMNegUnstableManifold.TOF = 40.0, 40.0
    EMPosUnstableArcs::Vector{MBD.CR3BPManifoldArc} = stopCrashes(EMPosUnstableManifold)
    EMNegUnstableArcs::Vector{MBD.CR3BPManifoldArc} = stopCrashes(EMNegUnstableManifold)
    departureArcs::Vector{MBD.CR3BPManifoldArc} = []
    for manifoldArc::MBD.CR3BPManifoldArc in append!(copy(EMPosUnstableArcs), copy(EMNegUnstableArcs))
        arc::MBD.CR3BPArc = propagateWithEvent(propagator, P2DistanceEvent, real(manifoldArc.initialCondition), [0, manifoldArc.TOF], EMDynamicsModel, [MoonSoI])
        (abs(getTimeByIndex(arc, -1)) == abs(manifoldArc.TOF)) && continue
        manifoldArc.TOF = getTimeByIndex(arc, -1)
        push!(departureArcs, manifoldArc)
    end
    println("Available departure arcs: $(length(departureArcs))\n")

    coarseSMOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(SMTargeter, "FamilyData/CR3BPSML1Halos.csv", "JC", SMJC)
    coarseSMOrbit.initialCondition[3] *= -1 # Flip to northern halo
    SMSolution::MBD.CR3BPMultipleShooterProblem = correct(SMTargeter, coarseSMOrbit.initialCondition, [0, coarseSMOrbit.period], SMJC)
    SMOrbit = MBD.CR3BPPeriodicOrbit(SMDynamicsModel, SMSolution.nodes[1].state.data[1:6], getPeriod(SMTargeter, SMSolution), getMonodromy(SMTargeter, SMSolution))
    println("Converged Sun-Mars orbit:\n\tIC:\t$(SMOrbit.initialCondition)\n\tP:\t$(SMOrbit.period)\n\tJC:\t$(getJacobiConstant(SMOrbit))\n")
    SMPosStableManifold::MBD.CR3BPManifold = getManifoldByArclength(SMOrbit, "Stable", "Positive", 1000/lstarSM, 50)
    SMNegStableManifold::MBD.CR3BPManifold = getManifoldByArclength(SMOrbit, "Stable", "Negative", 1000/lstarSM, 50)
    SMPosStableManifold.TOF, SMNegStableManifold.TOF = -20.0, -20.0
    SMPosStableArcs::Vector{MBD.CR3BPManifoldArc} = stopCrashes(SMPosStableManifold)
    SMNegStableArcs::Vector{MBD.CR3BPManifoldArc} = stopCrashes(SMNegStableManifold)
    arrivalArcs::Vector{MBD.CR3BPManifoldArc} = []
    arrivalPeriapses::Vector{Float64} = []
    for manifoldArc::MBD.CR3BPManifoldArc in append!(copy(SMPosStableArcs), copy(SMNegStableArcs))
        arc::MBD.CR3BPArc = propagateWithEvent(propagator, P2DistanceEvent, real(manifoldArc.initialCondition), [0, manifoldArc.TOF], SMDynamicsModel, [MarsSoI])
        (abs(getTimeByIndex(arc, -1)) == abs(manifoldArc.TOF)) && continue
        manifoldArc.TOF = getTimeByIndex(arc, -1)
        q_SoI_SI::Vector{Float64} = rotatingToPrimaryInertial(SMDynamicsModel, 1, [getStateByIndex(arc, -1)], [0.0])[1]
        Q_SoI_SI::Vector{Float64} = append!(copy(q_SoI_SI[1:3]).*lstarSM, copy(q_SoI_SI[4:6]).*lstarSM/tstarSM)
        oe::Vector{Float64} = getOrbitalElements(SDynamicsModel, Q_SoI_SI)
        if oe[2] < 1.0
            push!(arrivalArcs, manifoldArc)
            push!(arrivalPeriapses, oe[1]*(1-oe[2]))
        end
    end
    println("Available arrival arcs: $(length(arrivalArcs))\n")
    minPeriapsisIndex::Int64 = argmin(arrivalPeriapses)
    arrivalArc::MBD.CR3BPManifoldArc = arrivalArcs[minPeriapsisIndex]
    println("Choosing minimum periapsis arrival arc ($minPeriapsisIndex): r_{p} = $(arrivalPeriapses[minPeriapsisIndex]) km\n")

    SMOrbitArc::MBD.CR3BPArc = propagate(propagator, SMOrbit.initialCondition, [0, arrivalArc.orbitTime*SMOrbit.period], SMDynamicsModel)
    q_arr::Vector{Float64} = getStateByIndex(SMOrbitArc, -1)
    orbitArrivalArc::MBD.CR3BPArc = propagateWithEvent(momentumPropagator, orbitDepartureEvent, appendExtraInitialConditions(SMDynamicsModel, real(arrivalArc.initialCondition), MBD.MOMENTUM), [0, arrivalArc.TOF], SMDynamicsModel, [momentumPropagator, SMDynamicsModel, q_arr, SMMomentumDiff])
    t_orbitArrival::Float64 = -1*getTimeByIndex(orbitArrivalArc, -1)*tstarSM
    arrivalArcProp::MBD.CR3BPArc = propagate(propagator, real(arrivalArc.initialCondition), [0, arrivalArc.TOF], SMDynamicsModel)
    t_arrival::Float64 = -1*arrivalArc.TOF*tstarSM

    q_arr_SoI_SI_guess::Vector{Float64} = rotatingToSunEclipJ2000(SMDynamicsModel, initialEpoch, [getStateByIndex(arrivalArcProp, -1)], [0.0])[1]
    Q_arr_SoI_SI_guess::Vector{Float64} = append!(copy(q_arr_SoI_SI_guess[1:3]).*lstarSM, copy(q_arr_SoI_SI_guess[4:6]).*lstarSM/tstarSM)
    oe_arr_SoI_guess::Vector{Float64} = getOrbitalElements(SDynamicsModel, Q_arr_SoI_SI_guess)
    bridgeApoapsis::Float64 = oe_arr_SoI_guess[1]*(1+oe_arr_SoI_guess[2])
    r_arr_lower::Float64 = oe_arr_SoI_guess[1]*(1-oe_arr_SoI_guess[2])
    r_arr_upper::Float64 = oe_arr_SoI_guess[1]*(1+oe_arr_SoI_guess[2])
    println("Arrival radius bounds: [$r_arr_lower, $r_arr_upper] km\n")

    println("\nComputing MMAT transfers...")
    for day::Int64 in 1:length(days)
        count::Int64 = 0
        println("\tComputing transfers for $(SPICE.et2utc(initialEpochTime+3600*24*days[day], "C", 0))...")
        for d::Int64 in 1:length(departureArcs)
            theta_E_dep::Float64 = theta_E_0+days[day]*3600*24/lstarSE
            departureArc::MBD.CR3BPManifoldArc = departureArcs[d]
            EMOrbitArc::MBD.CR3BPArc = propagate(propagator, EMOrbit.initialCondition, [0, departureArc.orbitTime*EMOrbit.period], EMDynamicsModel)
            q_dep::Vector{Float64} = getStateByIndex(EMOrbitArc, -1)
            orbitDepartureArc::MBD.CR3BPArc = propagateWithEvent(momentumPropagator, orbitDepartureEvent, appendExtraInitialConditions(EMDynamicsModel, real(departureArc.initialCondition), MBD.MOMENTUM), [0, departureArc.TOF], EMDynamicsModel, [momentumPropagator, EMDynamicsModel, q_dep, EMMomentumDiff])
            t_orbitDeparture::Float64 = getTimeByIndex(orbitDepartureArc, -1)*tstarEM
            departureArcEMProp::MBD.CR3BPArc = propagate(propagator, real(departureArc.initialCondition), [0, departureArc.TOF], EMDynamicsModel)
            q_dep_blend_EI_EM::Vector{Float64} = rotatingToPrimaryEclipJ2000(EMDynamicsModel, initialEpoch, [getStateByIndex(departureArcEMProp, -1)], [getTimeByIndex(departureArcEMProp, -1)+days[day]*3600*24/tstarEM])[1]
            Q_dep_blend_EI::Vector{Float64} = append!(q_dep_blend_EI_EM[1:3].*lstarEM, q_dep_blend_EI_EM[4:6].*lstarEM./tstarEM)
            q_dep_blend_EI_SE::Vector{Float64} = append!(Q_dep_blend_EI[1:3]./lstarSE, Q_dep_blend_EI[4:6].*tstarSE./lstarSE)
            e_blend_SE::Float64 = (getTimeByIndex(departureArcEMProp, -1)*tstarEM+days[day]*3600*24)/tstarSE
            q_dep_blend_SE::Vector{Float64} = secondaryEclipJ2000ToRotating(SEDynamicsModel, initialEpoch, [q_dep_blend_EI_SE], [e_blend_SE])[1]
            testProp::MBD.CR3BPArc = propagateWithEvent(propagator, P2DistanceEvent, q_dep_blend_SE, [0, 4.0*pi], SEDynamicsModel, [Earth.bodyRadius/lstarSE])
            testTime::Float64 = getTimeByIndex(testProp, -1)
            departureArcSEProp::MBD.CR3BPArc = propagateWithEvent(propagator, P2DistanceEvent, q_dep_blend_SE, [0, testTime], SEDynamicsModel, [EarthSoI])
            (abs(getTimeByIndex(departureArcSEProp, -1)) == abs(testTime)) && continue
            t_departure::Float64 = departureArc.TOF*tstarEM+getTimeByIndex(departureArcSEProp, -1)*tstarSE

            q_dep_SoI_SI::Vector{Float64} = rotatingToSunEclipJ2000(SEDynamicsModel, initialEpoch, [getStateByIndex(departureArcSEProp, -1)], [(days[day]*3600*24+t_departure)/tstarSE])[1]
            Q_dep_SoI_SI::Vector{Float64} = append!(q_dep_SoI_SI[1:3].*lstarSE, q_dep_SoI_SI[4:6].*lstarSE./tstarSE)
            oe_dep_SoI::Vector{Float64} = getOrbitalElements(SDynamicsModel, Q_dep_SoI_SI)
            t_depConic::Float64 = MBD.getPeriod(SDynamicsModel, oe_dep_SoI)-solveKeplersEquation(SDynamicsModel, oe_dep_SoI)
            bridgePeriapsis::Float64 = oe_dep_SoI[1]*(1-oe_dep_SoI[2])
            bridgeRatio::Float64 = bridgePeriapsis/bridgeApoapsis
            e_bridge::Float64 = (1-bridgeRatio)/(1+bridgeRatio)
            a_bridge::Float64 = bridgePeriapsis/(1-e_bridge)
            Q_dep_peri_SI::Vector{Float64} = getCartesianState(SDynamicsModel, append!(oe_dep_SoI[1:5], 0.0))

            DeltaOmega_guess::Float64 = Omega_M-oe_dep_SoI[4]
            psi_guess::Float64 = acos(cos(i_M)*cos(oe_dep_SoI[3])+sin(i_M)*sin(oe_dep_SoI[3])*cos(DeltaOmega_guess))
            u_bridge_sin_guess::Float64 = sin(pi-i_M)*sin(DeltaOmega_guess)/sin(psi_guess)
            u_bridge_cos_guess::Float64 = (sin(i_M)*cos(DeltaOmega_guess)/cos(oe_dep_SoI[3])-cos(psi_guess)*tan(oe_dep_SoI[3]))/sin(psi_guess)
            u_bridge_guess::Float64 = atan(u_bridge_sin_guess, u_bridge_cos_guess)
            u_arr_sin_guess::Float64 = sin(oe_dep_SoI[3])*sin(DeltaOmega_guess)/sin(psi_guess)
            u_arr_cos_guess::Float64 = cos(DeltaOmega_guess)*cos(u_bridge_guess)+sin(DeltaOmega_guess)*sin(u_bridge_guess)*cos(oe_dep_SoI[3])
            u_arr_guess::Float64 = atan(u_arr_sin_guess, u_arr_cos_guess)
            oe_bridge_peri::Vector{Float64} = [a_bridge, e_bridge, oe_dep_SoI[3], oe_dep_SoI[4], oe_dep_SoI[5], 0.0]
            Q_bridge_peri_SI::Vector{Float64} = getCartesianState(SDynamicsModel, oe_bridge_peri)
            Deltav_1::Float64 = LinearAlgebra.norm(Q_bridge_peri_SI[4:6]-Q_dep_peri_SI[4:6])

            theta_bridge_int_guess::Float64 = u_bridge_guess-oe_dep_SoI[5]
            r_int_0::Float64 = a_bridge*(1-e_bridge^2)/(1+e_bridge*cos(theta_bridge_int_guess))
            r_int_1::Float64 = a_bridge*(1-e_bridge^2)/(1+e_bridge*cos(theta_bridge_int_guess+pi))

            if (r_arr_lower < r_int_0 < r_arr_upper)
                try
                    intersect_a::Vector{Float64} = findIntersection(MMATTargeter, SDynamicsModel, SMDynamicsModel, oe_M_0, initialEpoch, oe_dep_SoI, oe_bridge_peri, oe_arr_SoI_guess, getStateByIndex(arrivalArcProp, -1), r_int_0, theta_bridge_int_guess, u_arr_guess, 0, 0)
                    theta_M_arr_a::Float64 = intersect_a[10]+t_arrival/tstarSM
                    theta_M_dep_a::Float64 = intersect_a[10]-(intersect_a[9]+intersect_a[1]+t_depConic+t_departure)/tstarSM
                    while theta_M_dep_a < 0.0
                        theta_M_dep_a += 2*pi
                    end
                    TOF_a::Float64 = t_departure-t_orbitDeparture+t_depConic+intersect_a[1]+intersect_a[9]+t_arrival-t_orbitArrival
                    count += 1
                    exportCR3BPMMAT(initialEpoch, days[day]*3600*24, theta_E_dep, theta_M_arr_a, theta_M_dep_a, departureArc, departureArcSEProp, SDynamicsModel, oe_dep_SoI, t_depConic, oe_bridge_peri, intersect_a[1], intersect_a[3:8], intersect_a[9], arrivalArc, Deltav_1, intersect_a[2], TOF_a, mf, Symbol("Transfer_", Int(days[day]), "_", count, "a"))
                catch
                    continue
                end
                intersect_b::Vector{Float64} = findIntersection(MMATTargeter, SDynamicsModel, SMDynamicsModel, oe_M_0, initialEpoch, oe_dep_SoI, oe_bridge_peri, oe_arr_SoI_guess, getStateByIndex(arrivalArcProp, -1), r_int_0, theta_bridge_int_guess, u_arr_guess, 0, 1)
                theta_M_arr_b::Float64 = intersect_b[10]+t_arrival/tstarSM
                theta_M_dep_b::Float64 = intersect_b[10]-(intersect_b[9]+intersect_b[1]+t_depConic+t_departure)/tstarSM
                while theta_M_dep_b < 0.0
                    theta_M_dep_b += 2*pi
                end
                TOF_b::Float64 = t_departure-t_orbitDeparture+t_depConic+intersect_b[1]+intersect_b[9]+t_arrival-t_orbitArrival
                exportCR3BPMMAT(initialEpoch, days[day]*3600*24, theta_E_dep, theta_M_arr_b, theta_M_dep_b, departureArc, departureArcSEProp, SDynamicsModel, oe_dep_SoI, t_depConic, oe_bridge_peri, intersect_b[1], intersect_b[3:8], intersect_b[9], arrivalArc, Deltav_1, intersect_b[2], TOF_b, mf, Symbol("Transfer_", Int(days[day]), "_", count, "b"))
            elseif (r_arr_lower < r_int_1 < r_arr_upper)
                try
                    intersect_a = findIntersection(MMATTargeter, SDynamicsModel, SMDynamicsModel, oe_M_0, initialEpoch, oe_dep_SoI, oe_bridge_peri, oe_arr_SoI_guess, getStateByIndex(arrivalArcProp, -1), r_int_1, theta_bridge_int_guess+pi, u_arr_guess, 1, 0)
                    theta_M_arr_a = intersect_a[10]+t_arrival/tstarSM
                    theta_M_dep_a = intersect_a[10]-(intersect_a[9]+intersect_a[1]+t_depConic+t_departure)/tstarSM
                    while theta_M_dep_a < 0.0
                        theta_M_dep_a += 2*pi
                    end
                    TOF_a = t_departure-t_orbitDeparture+t_depConic+intersect_a[1]+intersect_a[9]+t_arrival-t_orbitArrival
                    count += 1
                    exportCR3BPMMAT(initialEpoch, days[day]*3600*24, theta_E_dep, theta_M_arr_a, theta_M_dep_a, departureArc, departureArcSEProp, SDynamicsModel, oe_dep_SoI, t_depConic, oe_bridge_peri, intersect_a[1], intersect_a[3:8], intersect_a[9], arrivalArc, Deltav_1, intersect_a[2], TOF_a, mf, Symbol("Transfer_", Int(days[day]), "_", count, "a"))
                catch
                    continue
                end
                intersect_b = findIntersection(MMATTargeter, SDynamicsModel, SMDynamicsModel, oe_M_0, initialEpoch, oe_dep_SoI, oe_bridge_peri, oe_arr_SoI_guess, getStateByIndex(arrivalArcProp, -1), r_int_1, theta_bridge_int_guess+pi, u_arr_guess, 1, 1)
                theta_M_arr_b = intersect_b[10]+t_arrival/tstarSM
                theta_M_dep_b = intersect_b[10]-(intersect_b[9]+intersect_b[1]+t_depConic+t_departure)/tstarSM
                while theta_M_dep_b < 0.0
                    theta_M_dep_b += 2*pi
                end
                TOF_b = t_departure-t_orbitDeparture+t_depConic+intersect_b[1]+intersect_b[9]+t_arrival-t_orbitArrival
                exportCR3BPMMAT(initialEpoch, days[day]*3600*24, theta_E_dep, theta_M_arr_b, theta_M_dep_b, departureArc, departureArcSEProp, SDynamicsModel, oe_dep_SoI, t_depConic, oe_bridge_peri, intersect_b[1], intersect_b[3:8], intersect_b[9], arrivalArc, Deltav_1, intersect_b[2], TOF_b, mf, Symbol("Transfer_", Int(days[day]), "_", count, "b"))
            end
        end
        println("\t\tTransfers found: $(count*2)")
    end

    MATLAB.close(mf)
    SPICE.kclear()
end

# end
