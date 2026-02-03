"""
Script for computing CR3BP MMATs between Earth-Moon and Sun-Mars systems

Author: Jonathan LeFevre Richmond
C: 1/24/26
U: 2/3/26
"""
# module MMMATCR3BP

using MBD, DifferentialEquations, Ephemerides, Logging, MATLAB, SPICE

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../CR3BPTargeters/MMATArrivalPhase.jl")
include("../CR3BPTargeters/SpatialPerpJC.jl")
include("../Utilities/Export.jl")

struct ArrCache
    arc::MBD.CR3BPArc
    oe_arr_SoI_guess::Vector{Float64}
    r_arr_l::Float64
    r_arr_u::Float64
    t_arrival::Float64
    t_orbitArrival::Float64
end

struct CharValues
    lstar::Float64
    tstar::Float64
end

struct DepCache
    arc::MBD.CR3BPArc
    t_orbitDeparture::Float64
end

function setupEnvironment()::MMATEnv
    EMSystemData = MBD.CR3BPSystemData("Earth", "Moon")
    SESystemData = MBD.CR3BPSystemData("Sun", "Earth")
    SSystemData = MBD.KSystemData("Sun")
    SMSystemData = MBD.CR3BPSystemData("Sun", "Mars")
    EMDynamicsModel = MBD.CR3BPDynamicsModel(EMSystemData)
    SEDynamicsModel = MBD.CR3BPDynamicsModel(SESystemData)
    SDynamicsModel = MBD.KDynamicsModel(SSystemData)
    SMDynamicsModel = MBD.CR3BPDynamicsModel(SMSystemData)

    Earth::MBD.BodyData, Moon::MBD.BodyData = EMSystemData.primaryData
    Sun::MBD.BodyData, Mars::MBD.BodyData = SMSystemData.primaryData

    charValues::NamedTuple = (EM = CharValues(getCharLength(EMSystemData), getCharTime(EMSystemData)),
                              SE = CharValues(getCharLength(SESystemData), getCharTime(SESystemData)),
                              SM = CharValues(getCharLength(SMSystemData), getCharTime(SMSystemData)))
    
    propagator = MBD.Propagator()
    momentumPropagator = MBD.Propagator(equationType = MBD.MOMENTUM)
    orbitDepartureEvent = DifferentialEquations.ContinuousCallback(momentumDifferenceConditionCR3BP, terminateAffect!)
    P1DistanceEvent = DifferentialEquations.ContinuousCallback(p1CR3BPDistanceCondition, terminateAffect!)
    P2DistanceEvent = DifferentialEquations.ContinuousCallback(p2CR3BPDistanceCondition, terminateAffect!)

    MoonSoI::Float64 = charValues.SE.lstar*(Moon.mass/Sun.mass)^(2/5)/charValues.EM.lstar
    EarthSoI::Float64 = 0.09877
    MarsSoI::Float64 = 0.05375
    EMMomentumDiff::Float64 = 1E-3
    SMMomentumDiff::Float64 = 1E-5

    initialEpoch::String = "Jan 1 2030"
    initialEpochTime::Float64 = SPICE.str2et(initialEpoch)
    # days::Vector{Float64} = [0.0]
    days::Vector{Float64} = collect(0.0:36.0:360.0)
    epochTimes::Vector{Float64} = initialEpochTime .+ days .* 3600 .* 24
    epochs::Vector{String} = [SPICE.et2utc(et, "C", 0) for et in epochTimes]

    q_E_0::Vector{Float64} = getEphemerides(initialEpoch, [0.0], "Earth", "Sun", "ECLIPJ2000")[1][1]
    q_M_0::Vector{Float64} = getEphemerides(initialEpoch, [0.0], "Mars", "Sun", "ECLIPJ2000")[1][1]
    oe_E::Vector{Float64} = SPICE.oscltx(q_E_0, initialEpochTime, Sun.gravParam)
    oe_M::Vector{Float64} = SPICE.oscltx(q_M_0, initialEpochTime, Sun.gravParam)

    return MMATEnv(EMDynamicsModel, SDynamicsModel, SEDynamicsModel, SMDynamicsModel, Earth, Mars, Moon, Sun, charValues, momentumPropagator, orbitDepartureEvent, propagator, P1DistanceEvent, P2DistanceEvent, EarthSoI, EMMomentumDiff, MarsSoI, MoonSoI, SMMomentumDiff, days, epochs, initialEpoch, initialEpochTime, oe_E, oe_M)
end

function computeDepartureArcs(env::MMATEnv, targeter, JC::Float64)
    dynamicsModel::MBD.CR3BPDynamicsModel = env.EMDynamicsModel
    lstar::Float64 = env.charValues.EM.lstar
    tstar::Float64 = env.charValues.EM.tstar

    coarseOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(targeter, "FamilyData/CR3BPEML1Halos.csv", "JC", JC)
    coarseOrbit.initialCondition[3] *= -1 # Flip to northern halo
    solution::MBD.CR3BPMultipleShooterProblem = correct(targeter, coarseOrbit.initialCondition, [0, coarseOrbit.period], JC)
    orbit = MBD.CR3BPPeriodicOrbit(dynamicsModel, solution.nodes[1].state.data[1:6], getPeriod(targeter, solution), getMonodromy(targeter, solution))
    println("Converged Earth-Moon orbit:\n\tIC:\t$(orbit.initialCondition)\n\tP:\t$(orbit.period)\n\tJC:\t$(getJacobiConstant(orbit))\n")

    posUnstableManifold::MBD.CR3BPManifold = getManifoldByArclength(orbit, "Unstable", "Positive", 25/lstar, 100)
    negUnstableManifold::MBD.CR3BPManifold = getManifoldByArclength(orbit, "Unstable", "Negative", 25/lstar, 100)
    posUnstableManifold.TOF, negUnstableManifold.TOF = 40.0, 40.0
    unstableManifoldArcs::Vector{MBD.CR3BPManifoldArc} = vcat(stopCrashes(posUnstableManifold), stopCrashes(negUnstableManifold))
    arcs::Vector{MBD.CR3BPManifoldArc} = []
    cache::Vector{DepCache} = []
    for manifoldArc::MBD.CR3BPManifoldArc in unstableManifoldArcs
        arc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.P2DistanceEvent, real(manifoldArc.initialCondition), [0, manifoldArc.TOF], dynamicsModel, [env.MoonSoI])
        (abs(getTimeByIndex(arc, -1)) == abs(manifoldArc.TOF)) && continue
        manifoldArc.TOF = getTimeByIndex(arc, -1)
        push!(arcs, manifoldArc)
        orbitArc::MBD.CR3BPArc = propagate(env.propagator, orbit.initialCondition, [0, manifoldArc.orbitTime*orbit.period], dynamicsModel)
        q_dep::Vector{Float64} = getStateByIndex(orbitArc, -1)
        orbitDepartureArc::MBD.CR3BPArc = propagateWithEvent(env.momentumPropagator, env.orbitDepartureEvent, appendExtraInitialConditions(dynamicsModel, real(manifoldArc.initialCondition), MBD.MOMENTUM), [0, manifoldArc.TOF], dynamicsModel, [env.momentumPropagator, dynamicsModel, q_dep, env.EMMomentumDiff])
        t_orbitDeparture::Float64 = getTimeByIndex(orbitDepartureArc, -1)*tstar
        push!(cache, DepCache(arc, t_orbitDeparture))
    end
    println("Available departure arcs: $(length(arcs))\n")

    return arcs, cache
end

function computeArrivalArc(env::MMATEnv, targeter, JC::Float64)
    dynamicsModel::MBD.CR3BPDynamicsModel = env.SMDynamicsModel
    lstar::Float64 = env.charValues.SM.lstar
    tstar::Float64 = env.charValues.SM.tstar

    coarseOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(targeter, "FamilyData/CR3BPSML1Halos.csv", "JC", JC)
    coarseOrbit.initialCondition[3] *= -1 # Flip to northern halo
    solution::MBD.CR3BPMultipleShooterProblem = correct(targeter, coarseOrbit.initialCondition, [0, coarseOrbit.period], JC)
    orbit = MBD.CR3BPPeriodicOrbit(dynamicsModel, solution.nodes[1].state.data[1:6], getPeriod(targeter, solution), getMonodromy(targeter, solution))
    println("Converged Sun-Mars orbit:\n\tIC:\t$(orbit.initialCondition)\n\tP:\t$(orbit.period)\n\tJC:\t$(getJacobiConstant(orbit))\n")

    posStableManifold::MBD.CR3BPManifold = getManifoldByArclength(orbit, "Stable", "Positive", 1000/lstar, 50)
    negStableManifold::MBD.CR3BPManifold = getManifoldByArclength(orbit, "Stable", "Negative", 1000/lstar, 50)
    posStableManifold.TOF, negStableManifold.TOF = -40.0, -40.0
    stableManifoldArcs::Vector{MBD.CR3BPManifoldArc} = vcat(stopCrashes(posStableManifold), stopCrashes(negStableManifold))
    arcs::Vector{MBD.CR3BPManifoldArc} = []
    periapses::Vector{Float64} = []
    for manifoldArc::MBD.CR3BPManifoldArc in stableManifoldArcs
        arc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.P2DistanceEvent, real(manifoldArc.initialCondition), [0, manifoldArc.TOF], dynamicsModel, [env.MarsSoI])
        (abs(getTimeByIndex(arc, -1)) == abs(manifoldArc.TOF)) && continue
        manifoldArc.TOF = getTimeByIndex(arc, -1)
        q_SoI_SI::Vector{Float64} = rotatingToPrimaryInertial(dynamicsModel, 1, [getStateByIndex(arc, -1)], [0.0])[1]
        @views Q_SoI_SI = similar(q_SoI_SI)
        Q_SoI_SI[1:3] .= q_SoI_SI[1:3].*lstar
        Q_SoI_SI[4:6] .= q_SoI_SI[4:6].*lstar./tstar
        oe::Vector{Float64} = getOrbitalElements(env.SDynamicsModel, Q_SoI_SI)
        if oe[2] < 1.0
            push!(arcs, manifoldArc)
            push!(periapses, oe[1]*(1-oe[2]))
        end
    end
    println("Available arrival arcs: $(length(arcs))")

    minPeriapsisIndex::Int64 = argmin(periapses)
    arrivalArc::MBD.CR3BPManifoldArc = arcs[minPeriapsisIndex]
    println("Choosing minimum periapsis arrival arc ($minPeriapsisIndex): r_{p} = $(periapses[minPeriapsisIndex]) km")

    orbitArc::MBD.CR3BPArc = propagate(env.propagator, orbit.initialCondition, [0, arrivalArc.orbitTime*orbit.period], dynamicsModel)
    q_arr::Vector{Float64} = getStateByIndex(orbitArc, -1)
    orbitArrivalArc::MBD.CR3BPArc = propagateWithEvent(env.momentumPropagator, env.orbitDepartureEvent, appendExtraInitialConditions(dynamicsModel, real(arrivalArc.initialCondition), MBD.MOMENTUM), [0, arrivalArc.TOF], dynamicsModel, [env.momentumPropagator, dynamicsModel, q_arr, env.SMMomentumDiff])
    t_orbitArrival::Float64 = -1*getTimeByIndex(orbitArrivalArc, -1)*tstar
    arrivalArcProp::MBD.CR3BPArc = propagate(env.propagator, real(arrivalArc.initialCondition), [0, arrivalArc.TOF], dynamicsModel)
    t_arrival::Float64 = -1*arrivalArc.TOF*tstar
    q_arr_SoI_SI_guess::Vector{Float64} = rotatingToSunEclipJ2000(dynamicsModel, env.initialEpoch, [getStateByIndex(arrivalArcProp, -1)], [0.0])[1]
    @views Q_arr_SoI_SI_guess = similar(q_arr_SoI_SI_guess)
    Q_arr_SoI_SI_guess[1:3] .= q_arr_SoI_SI_guess[1:3].*lstar
    Q_arr_SoI_SI_guess[4:6] .= q_arr_SoI_SI_guess[4:6].*lstar./tstar
    oe_arr_SoI_guess::Vector{Float64} = getOrbitalElements(env.SDynamicsModel, Q_arr_SoI_SI_guess)
    r_arr_l::Float64 = oe_arr_SoI_guess[1]*(1-oe_arr_SoI_guess[2])
    r_arr_u::Float64 = oe_arr_SoI_guess[1]*(1+oe_arr_SoI_guess[2])
    println("Arrival radius bounds: [$r_arr_l, $r_arr_u] km\n")

    return arrivalArc, ArrCache(arrivalArcProp, oe_arr_SoI_guess, r_arr_l, r_arr_u, t_arrival, t_orbitArrival)
end

function findIntersection(targeter::MMATArrivalPhaseTargeter, env::MMATEnv, arrCache::ArrCache, oe_dep_SoI::Vector{Float64}, oe_bridge_peri::Vector{Float64}, q_arr_SoI::Vector{Float64}, r_int::Float64, theta_bridge_int_guess::Float64, u_arr_guess::Float64, n::Int64, o::Int64)
    theta_arr_int_guess::Float64 = 2*o*pi+(-1)^o*acos((arrCache.oe_arr_SoI_guess[1]*(1-arrCache.oe_arr_SoI_guess[2]^2)/r_int-1)/arrCache.oe_arr_SoI_guess[2])
    omega_arr_guess::Float64 = (u_arr_guess-(theta_arr_int_guess+n*pi) > 0) ? u_arr_guess-(theta_arr_int_guess+n*pi) : 2*pi+u_arr_guess-(theta_arr_int_guess+n*pi)
    sigma_guess::Float64 = atan(cos(env.oe_M[3])*tan(env.oe_M[4]))+omega_arr_guess+arrCache.oe_arr_SoI_guess[6]
    delta::Float64 = atan(q_arr_SoI[2], q_arr_SoI[1])
    X_guess::Vector{Float64} = [theta_bridge_int_guess, sigma_guess-(env.oe_M[4]+env.oe_M[5]+env.oe_M[6])-delta, theta_arr_int_guess]
    X::Vector{Float64} = correct(targeter, env, X_guess, oe_bridge_peri, q_arr_SoI)
    Q_bridge_int_SI::Vector{Float64} = getCartesianState(env.SDynamicsModel, append!(oe_bridge_peri[1:5], X[1]))
    t_bridgeConic::Float64 = solveKeplersEquation(env.SDynamicsModel, append!(oe_bridge_peri[1:5], X[1]))
    q_arr_SoI_SI::Vector{Float64} = rotatingToSunEclipJ2000(env.SMDynamicsModel, env.initialEpoch, [q_arr_SoI], [X[2]])[1]
    @views Q_arr_SoI_SI = similar(q_arr_SoI_SI)
    Q_arr_SoI_SI[1:3] .= q_arr_SoI_SI[1:3].*env.charValues.SM.lstar
    Q_arr_SoI_SI[4:6] .= q_arr_SoI_SI[4:6].*env.charValues.SM.lstar./env.charValues.SM.tstar
    oe_arr_SoI::Vector{Float64} = getOrbitalElements(env.SDynamicsModel, Q_arr_SoI_SI)
    oe_arr_int::Vector{Float64} = append!(oe_arr_SoI[1:5], X[3])
    Q_arr_int_SI::Vector{Float64} = getCartesianState(env.SDynamicsModel, oe_arr_int)
    Deltav_2::Float64 = LinearAlgebra.norm(Q_arr_int_SI[4:6]-Q_bridge_int_SI[4:6])
    t_arrConic::Float64 = solveKeplersEquation(env.SDynamicsModel, oe_arr_SoI)-solveKeplersEquation(env.SDynamicsModel, oe_arr_int)
    (t_arrConic < 0) && (t_arrConic += MBD.getPeriod(env.SDynamicsModel, oe_arr_SoI))

    return [t_bridgeConic, Deltav_2, oe_arr_int..., t_arrConic, X[2]]
end

function computeMMAT(targeter::MMATArrivalPhaseTargeter, env::MMATEnv, depCache::DepCache, arrCache::ArrCache, oe_dep_SoI::Vector{Float64}, oe_bridge_peri::Vector{Float64}, r_int::Float64, theta_bridge_int_guess::Float64, u_arr_guess::Float64, n::Int64, o::Int64, t_depConic::Float64, t_departure::Float64)
    (n == 1) && (theta_bridge_int_guess += pi)
    intersect::Vector{Float64} = findIntersection(targeter, env, arrCache, oe_dep_SoI, oe_bridge_peri, getStateByIndex(arrCache.arc, -1), r_int, theta_bridge_int_guess, u_arr_guess, n, o)
    theta_M_arr::Float64 = intersect[10]+arrCache.t_arrival/env.charValues.SM.tstar
    theta_M_dep::Float64 = intersect[10]-(intersect[9]+intersect[1]+t_depConic+t_departure)/env.charValues.SM.tstar
    while theta_M_dep < 0.0
        theta_M_dep += 2*pi
    end
    TOF::Float64 = t_departure-depCache.t_orbitDeparture+t_depConic+intersect[1]+intersect[9]+arrCache.t_arrival-arrCache.t_orbitArrival
    
    return [intersect..., theta_M_arr, theta_M_dep, TOF]
end

function computeMMATs(env::MMATEnv, targeter::MMATArrivalPhaseTargeter, depArcs::Vector{MBD.CR3BPManifoldArc}, depCache::Vector{DepCache}, arrArc::MBD.CR3BPManifoldArc, arrCache::ArrCache, mf::MATLAB.MatFile)
    println("Computing MMAT transfers...")
    for (idx::Int64, day::Float64) in pairs(env.days)
        e_dep::Float64 = 3600*24*day
        println("\tComputing transfers for $(env.epochs[idx])...")
        theta_E_dep::Float64 = env.oe_E[6]+e_dep/env.charValues.SE.tstar
        count0::Int64, count1::Int64 = 0, 0
        fails::Int64 = 0
        for d::Int64 in 1:length(depArcs)
            depArc::MBD.CR3BPManifoldArc = depArcs[d]
            cache::DepCache = depCache[d]
            q_dep_blend_EI_EM::Vector{Float64} = rotatingToPrimaryEclipJ2000(env.EMDynamicsModel, env.initialEpoch, [getStateByIndex(cache.arc, -1)], [getTimeByIndex(cache.arc, -1)+e_dep/env.charValues.EM.tstar])[1]
            @views Q_dep_blend_EI = similar(q_dep_blend_EI_EM)
            Q_dep_blend_EI[1:3] .= q_dep_blend_EI_EM[1:3].*env.charValues.EM.lstar
            Q_dep_blend_EI[4:6] .= q_dep_blend_EI_EM[4:6].*env.charValues.EM.lstar./env.charValues.EM.tstar
            @views q_dep_blend_EI_SE = similar(Q_dep_blend_EI)
            q_dep_blend_EI_SE[1:3] .= Q_dep_blend_EI[1:3]./env.charValues.SE.lstar
            q_dep_blend_EI_SE[4:6] .= Q_dep_blend_EI[4:6].*env.charValues.SE.tstar./env.charValues.SE.lstar
            e_blend_SE::Float64 = (getTimeByIndex(cache.arc, -1)*env.charValues.EM.tstar+e_dep)/env.charValues.SE.tstar
            q_dep_blend_SE::Vector{Float64} = secondaryEclipJ2000ToRotating(env.SEDynamicsModel, env.initialEpoch, [q_dep_blend_EI_SE], [e_blend_SE])[1]
            testProp::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.P2DistanceEvent, q_dep_blend_SE, [0, 4.0*pi], env.SEDynamicsModel, [env.Earth.bodyRadius/env.charValues.SE.lstar])
            testTime::Float64 = getTimeByIndex(testProp, -1)
            depArcSEProp::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.P2DistanceEvent, q_dep_blend_SE, [0, testTime], env.SEDynamicsModel, [env.EarthSoI])
            (abs(getTimeByIndex(depArcSEProp, -1)) == abs(testTime)) && continue
            t_departure::Float64 = depArc.TOF*env.charValues.EM.tstar+getTimeByIndex(depArcSEProp, -1)*env.charValues.SE.tstar

            q_dep_SoI_SI::Vector{Float64} = rotatingToSunEclipJ2000(env.SEDynamicsModel, env.initialEpoch, [getStateByIndex(depArcSEProp, -1)], [(e_dep+t_departure)/env.charValues.SE.tstar])[1]
            @views Q_dep_SoI_SI = similar(q_dep_SoI_SI)
            Q_dep_SoI_SI[1:3] .= q_dep_SoI_SI[1:3].*env.charValues.SE.lstar
            Q_dep_SoI_SI[4:6] .= q_dep_SoI_SI[4:6].*env.charValues.SE.lstar./env.charValues.SE.tstar
            oe_dep_SoI::Vector{Float64} = getOrbitalElements(env.SDynamicsModel, Q_dep_SoI_SI)
            t_depConic::Float64 = MBD.getPeriod(env.SDynamicsModel, oe_dep_SoI)-solveKeplersEquation(env.SDynamicsModel, oe_dep_SoI)
            bridgePeriapsis::Float64 = oe_dep_SoI[1]*(1-oe_dep_SoI[2])
            bridgeRatio::Float64 = bridgePeriapsis/arrCache.r_arr_u
            e_bridge::Float64 = (1-bridgeRatio)/(1+bridgeRatio)
            a_bridge::Float64 = bridgePeriapsis/(1-e_bridge)
            Q_dep_peri_SI::Vector{Float64} = getCartesianState(env.SDynamicsModel, append!(oe_dep_SoI[1:5], 0.0))

            DeltaOmega_guess::Float64 = env.oe_M[4]-oe_dep_SoI[4]
            s_i_M::Float64, c_i_M::Float64 = sincos(env.oe_M[3])
            s_i_dep::Float64, c_i_dep::Float64 = sincos(oe_dep_SoI[3])
            s_DeltaOmega_guess::Float64, c_DeltaOmega_guess::Float64 = sincos(DeltaOmega_guess)
            psi_guess::Float64 = acos(c_i_M*c_i_dep+s_i_M*s_i_dep*c_DeltaOmega_guess)
            s_psi_guess::Float64, c_psi_guess::Float64 = sincos(psi_guess)
            s_u_bridge_guess::Float64 = sin(pi-env.oe_M[3])*s_DeltaOmega_guess/s_psi_guess
            c_u_bridge_guess::Float64 = (s_i_M*c_DeltaOmega_guess/c_i_dep-c_psi_guess*(s_i_dep/c_i_dep))/s_psi_guess
            u_bridge_guess::Float64 = atan(s_u_bridge_guess, c_u_bridge_guess)
            s_u_arr_guess::Float64 = s_i_dep*s_DeltaOmega_guess/s_psi_guess
            c_u_arr_guess::Float64 = c_DeltaOmega_guess*c_u_bridge_guess+s_DeltaOmega_guess*s_u_bridge_guess*c_i_dep
            u_arr_guess::Float64 = atan(s_u_arr_guess, c_u_arr_guess)
            oe_bridge_peri::Vector{Float64} = [a_bridge, e_bridge, oe_dep_SoI[3], oe_dep_SoI[4], oe_dep_SoI[5], 0.0]
            Q_bridge_peri_SI::Vector{Float64} = getCartesianState(env.SDynamicsModel, oe_bridge_peri)
            Deltav_1::Float64 = LinearAlgebra.norm(Q_bridge_peri_SI[4:6]-Q_dep_peri_SI[4:6])

            theta_bridge_int_guess::Float64 = u_bridge_guess-oe_dep_SoI[5]
            r_int_0::Float64 = a_bridge*(1-e_bridge^2)/(1+e_bridge*cos(theta_bridge_int_guess))
            r_int_1::Float64 = a_bridge*(1-e_bridge^2)/(1+e_bridge*cos(theta_bridge_int_guess+pi))

            if (arrCache.r_arr_l < r_int_0 < arrCache.r_arr_u)
                try
                    transfer_a::Vector{Float64} = computeMMAT(targeter, env, cache, arrCache, oe_dep_SoI, oe_bridge_peri, r_int_0, theta_bridge_int_guess, u_arr_guess, 0, 0, t_depConic, t_departure)
                    count0 += 1
                    exportCR3BPMMAT(env, e_dep, transfer_a[1:10], theta_E_dep, transfer_a[11], transfer_a[12], depArc, depArcSEProp, oe_dep_SoI, t_depConic, oe_bridge_peri, arrArc, Deltav_1, transfer_a[13], mf, Symbol("Transfer_Day", Int(day), "_Type0_", count0, "a"))
                catch
                    fails += 2
                    continue
                end
                transfer_b::Vector{Float64} = computeMMAT(targeter, env, cache, arrCache, oe_dep_SoI, oe_bridge_peri, r_int_0, theta_bridge_int_guess, u_arr_guess, 0, 1, t_depConic, t_departure)
                exportCR3BPMMAT(env, e_dep, transfer_b[1:10], theta_E_dep, transfer_b[11], transfer_b[12], depArc, depArcSEProp, oe_dep_SoI, t_depConic, oe_bridge_peri, arrArc, Deltav_1, transfer_b[13], mf, Symbol("Transfer_Day", Int(day), "_Type0_", count0, "b"))
            elseif (arrCache.r_arr_l < r_int_1 < arrCache.r_arr_u)
                try
                    transfer_a = computeMMAT(targeter, env, cache, arrCache, oe_dep_SoI, oe_bridge_peri, r_int_1, theta_bridge_int_guess, u_arr_guess, 1, 0, t_depConic, t_departure)
                    count1 += 1
                    exportCR3BPMMAT(env, e_dep, transfer_a[1:10], theta_E_dep, transfer_a[11], transfer_a[12], depArc, depArcSEProp, oe_dep_SoI, t_depConic, oe_bridge_peri, arrArc, Deltav_1, transfer_a[13], mf, Symbol("Transfer_Day", Int(day), "_Type1_", count1, "a"))
                catch
                    fails += 2
                    continue
                end
                transfer_b = computeMMAT(targeter, env, cache, arrCache, oe_dep_SoI, oe_bridge_peri, r_int_1, theta_bridge_int_guess, u_arr_guess, 1, 1, t_depConic, t_departure)
                exportCR3BPMMAT(env, e_dep, transfer_b[1:10], theta_E_dep, transfer_b[11], transfer_b[12], depArc, depArcSEProp, oe_dep_SoI, t_depConic, oe_bridge_peri, arrArc, Deltav_1, transfer_b[13], mf, Symbol("Transfer_Day", Int(day), "_Type1_", count1, "b"))
            end
        end
        println("\t\tTransfers found: $(count0*2+count1*2) ($fails failed)")
    end
end

function run_MarsMMATCR3BP()
    mf = MATLAB.MatFile("Output/MarsMMATCR3BP.mat", "w")
    SPICE.furnsh("SPICEKernels/naif0012.tls", "SPICEKernels/de430.bsp", "SPICEKernels/de440.bsp", "SPICEKernels/mar097.bsp")

    env = setupEnvironment()

    EMTargeter = SpatialPerpJCTargeter(env.EMDynamicsModel)
    EMJC::Float64 = 3.03
    depArcs::Vector{MBD.CR3BPManifoldArc}, depCache::Vector{DepCache} = computeDepartureArcs(env, EMTargeter, EMJC)

    SMTargeter = SpatialPerpJCTargeter(env.SMDynamicsModel)
    SMJC::Float64 = 3.0001857
    arrArc::MBD.CR3BPManifoldArc, arrCache::ArrCache = computeArrivalArc(env, SMTargeter, SMJC)

    # eph = Ephemerides.EphemerisProvider(["SPICEKernels/de430.bsp", "SPICEKernels/de440.bsp", "SPICEKernels/mar097.bsp"])

    MMATTargeter = MMATArrivalPhaseTargeter(env.SDynamicsModel, sqrt(eps(Float64)))
    computeMMATs(env, MMATTargeter, depArcs, depCache, arrArc, arrCache, mf)

    MATLAB.close(mf)
    SPICE.kclear()
end

# end
