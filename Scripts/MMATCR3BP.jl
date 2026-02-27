"""
Script for computing CR3BP MMATs between Earth-Moon and Sun-planet systems

Author: Jonathan LeFevre Richmond
C: 1/24/26
U: 2/27/26
"""

module MMATCR3BP

using MBD, DifferentialEquations, Logging, MATLAB, SPICE
using ..EphemeridesLoader: Ephemerides, FrameTransformations

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../CR3BPTargeters/MMATArrivalPhase.jl")
include("../CR3BPTargeters/PlanarPerpJC.jl")
include("../CR3BPTargeters/SpatialPerpJC.jl")
include("../Utilities/Export.jl")
include("../Utilities/Frames.jl")

# (SoI, momentum diff)
const sysParams = Dict{String, Tuple{Float64, Float64}}("Venus" => (0.1107, 1E-5),
                                                        "Earth" => (0.09877, 1E-5),
                                                        "Mars" => (0.05375, 1E-5))
const targeterMap = Dict{}("L1Lyapunov" => PlanarPerpJCTargeter,
                           "L2Lyapunov" => PlanarPerpJCTargeter,
                           "L1Halo" => SpatialPerpJCTargeter,
                           "L2Halo" => SpatialPerpJCTargeter)

struct ArrCache
    arc::MBD.CR3BPArc
    oe_arr_SoI_guess::Vector{Float64}
    r_arr_l::Float64
    r_arr_u::Float64
    t_arrival::Float64
    t_orbitArrival::Float64
end

struct ArrivalCase
    JCs::Vector{Vector{Float64}}
    mode::Symbol
    planet::String
end

struct CharValues
    lstar::Float64
    mstar::Float64
    tstar::Float64
end

struct DepCache
    arc::MBD.CR3BPArc
    t_orbitDeparture::Float64
end

struct MMAT_ext
    branch::Int64
    day::Float64
    type::Int64

    Deltav_1::Float64
    depArc::MBD.CR3BPManifoldArc
    depArcSEProp::MBD.CR3BPArc
    oe_bridge_peri::Vector{Float64}
    oe_dep_SoI::Vector{Float64}
    transfer::Vector{Float64}
    t_depConic::Float64
end

struct MMAT_int
    branch::Int64
    day::Float64
    type::Int64

    depArc::MBD.CR3BPManifoldArc
    depArcSEProp::MBD.CR3BPArc
    oe_dep_SoI::Vector{Float64}
    transfer::Vector{Float64}
end

function buildFrame(eph::Ephemerides.EphemerisProvider, Sun::MBD.BodyData, Earth::MBD.BodyData, Moon::MBD.BodyData, arrData::MBD.BodyData)
    frame = FrameTransformations.FrameSystem{2, Float64}()
    FrameTransformations.add_axes_icrf!(frame)
    FrameTransformations.add_axes_ecl2000!(frame)
    FrameTransformations.add_point!(frame, :SSB, 0, 1)
    FrameTransformations.add_point_ephemeris!(frame, eph, :Sun, Int64(Sun.spiceID))
    FrameTransformations.add_point_ephemeris!(frame, eph, :B1, firstdigit(Earth.spiceID))
    FrameTransformations.add_point_ephemeris!(frame, eph, :Earth, Int64(Earth.spiceID))
    FrameTransformations.add_point_ephemeris!(frame, eph, :Moon, Int64(Moon.spiceID))
    FrameTransformations.add_point_ephemeris!(frame, eph, Symbol(arrData.name, "_B"), firstdigit(arrData.spiceID))
    FrameTransformations.add_point_ephemeris!(frame, eph, Symbol(arrData.name), Int64(arrData.spiceID))

    return frame
end

function setupEnvironment(eph::Ephemerides.EphemerisProvider, arrBody::String)::MMATEnv
    EMSystemData = MBD.CR3BPSystemData("Earth", "Moon")
    SESystemData = MBD.CR3BPSystemData("Sun", "Earth")
    SSystemData = MBD.KSystemData("Sun")
    SArrSystemData = MBD.CR3BPSystemData("Sun", arrBody)
    EMDynamicsModel = MBD.CR3BPDynamicsModel(EMSystemData)
    SEDynamicsModel = MBD.CR3BPDynamicsModel(SESystemData)
    SDynamicsModel = MBD.KDynamicsModel(SSystemData)
    SArrDynamicsModel = MBD.CR3BPDynamicsModel(SArrSystemData)

    Earth::MBD.BodyData, Moon::MBD.BodyData = EMSystemData.primaryData
    Sun::MBD.BodyData, arrData::MBD.BodyData = SArrSystemData.primaryData

    charValues::NamedTuple = (EM = CharValues(getCharLength(EMSystemData), getCharMass(EMSystemData), getCharTime(EMSystemData)),
                              SE = CharValues(getCharLength(SESystemData), getCharMass(SESystemData), getCharTime(SESystemData)),
                              SArr = CharValues(getCharLength(SArrSystemData), getCharMass(SArrSystemData), getCharTime(SArrSystemData)))
    
    propagator = MBD.Propagator()
    momentumPropagator = MBD.Propagator(equationType = MBD.MOMENTUM)
    orbitDepartureEvent = DifferentialEquations.ContinuousCallback(momentumDifferenceConditionCR3BP, terminateAffect!)
    P1DistanceEvent = DifferentialEquations.ContinuousCallback(p1CR3BPDistanceCondition, terminateAffect!)
    P2DistanceEvent = DifferentialEquations.ContinuousCallback(p2CR3BPDistanceCondition, terminateAffect!)

    MoonSoI::Float64 = charValues.SE.lstar*(Moon.mass/Sun.mass)^(2/5)/charValues.EM.lstar
    EarthSoI::Float64 = sysParams["Earth"][1]
    arrSoI::Float64 = sysParams[arrBody][1]
    EMMomentumDiff::Float64 = 1E-3
    SArrMomentumDiff::Float64 = sysParams[arrBody][2]

    initialEpoch::String = "Jan 1 2030"
    initialEpochTime::Float64 = SPICE.str2et(initialEpoch)
    # days::Vector{Float64} = collect(0.0:1.0:5.0)
    days::Vector{Float64} = collect(0.0:1.0:364.0)
    # days::Vector{Float64} = collect((10*365.0):1.0:(11*365.0-1.0)) # 2040
    epochTimes::Vector{Float64} = initialEpochTime .+ days .* 3600 .* 24
    epochs::Vector{String} = [SPICE.et2utc(et, "C", 0) for et in epochTimes]

    frame::FrameTransformations.FrameSystem = buildFrame(eph, Sun, Earth, Moon, arrData)

    Q_E_0::Vector{Float64} = FrameTransformations.vector6(frame, Int64(Sun.spiceID), Int64(Earth.spiceID), 17, initialEpochTime)
    Q_arrBody_0::Vector{Float64} = FrameTransformations.vector6(frame, Int64(Sun.spiceID), Int64(arrData.spiceID), 17, initialEpochTime)
    oe_E::Vector{Float64} = getOrbitalElements(SDynamicsModel, Q_E_0)
    oe_arrBody::Vector{Float64} = getOrbitalElements(SDynamicsModel, Q_arrBody_0)

    return MMATEnv(EMDynamicsModel, SDynamicsModel, SEDynamicsModel, SArrDynamicsModel, Earth, arrData, Moon, Sun, charValues, momentumPropagator, orbitDepartureEvent, propagator, P1DistanceEvent, P2DistanceEvent, EarthSoI, EMMomentumDiff, arrSoI, MoonSoI, SArrMomentumDiff, days, epochs, frame, initialEpoch, initialEpochTime, oe_E, oe_arrBody)
end

function computeDepartureArcs(env::MMATEnv, targeter, family::String, JC::Float64, flip::Bool)
    dynamicsModel::MBD.CR3BPDynamicsModel = env.EMDynamicsModel
    lstar::Float64 = env.charValues.EM.lstar
    tstar::Float64 = env.charValues.EM.tstar

    coarseOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(targeter, "FamilyData/CR3BPEM$(family)s.csv", "JC", JC)
    if flip
        if occursin("Halo", family)
            (coarseOrbit.initialCondition[3] *= -1) # Flip to northern halo
        end
    end
    solution::MBD.CR3BPMultipleShooterProblem = correct(targeter, coarseOrbit.initialCondition, [0, coarseOrbit.period], JC)
    orbit = MBD.CR3BPPeriodicOrbit(dynamicsModel, solution.nodes[1].state.data[1:6], getPeriod(targeter, solution), getMonodromy(targeter, solution))
    println("Converged Earth-Moon orbit:\n\tIC:\t$(orbit.initialCondition)\n\tP:\t$(orbit.period)\n\tJC:\t$(getJacobiConstant(orbit))\n")

    posUnstableManifold::MBD.CR3BPManifold = getManifoldByArclength(orbit, "Unstable", "Positive", 25/lstar, 50)
    negUnstableManifold::MBD.CR3BPManifold = getManifoldByArclength(orbit, "Unstable", "Negative", 25/lstar, 50)
    posUnstableManifold.TOF, negUnstableManifold.TOF = 40.0, 40.0
    unstableManifoldArcs::Vector{MBD.CR3BPManifoldArc} = vcat(stopCrashes(posUnstableManifold), stopCrashes(negUnstableManifold))
    nThreads::Int64 = Threads.maxthreadid()
    arcs_tls::Vector{Vector{MBD.CR3BPManifoldArc}} = [Vector{MBD.CR3BPManifoldArc}(undef, 0) for _ in 1:nThreads]
    cache_tls::Vector{Vector{DepCache}} = [Vector{DepCache}(undef, 0) for _ in 1:nThreads]
    Threads.@threads for a::Int64 in eachindex(unstableManifoldArcs)
        id::Int64 = Threads.threadid()
        manifoldArc::MBD.CR3BPManifoldArc = unstableManifoldArcs[a]
        arc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.P2DistanceEvent, real(manifoldArc.initialCondition), [0, manifoldArc.TOF], dynamicsModel, [env.MoonSoI])
        (abs(getTimeByIndex(arc, -1)) == abs(manifoldArc.TOF)) && continue
        localArc = MBD.CR3BPManifoldArc(manifoldArc.periodicOrbit, manifoldArc.orbitTime, manifoldArc.direction, manifoldArc.d, manifoldArc.initialCondition, getTimeByIndex(arc, -1))
        orbitArc::MBD.CR3BPArc = propagate(env.propagator, orbit.initialCondition, [0, localArc.orbitTime*orbit.period], dynamicsModel)
        q_dep::Vector{Float64} = getStateByIndex(orbitArc, -1)
        orbitDepartureArc::MBD.CR3BPArc = propagateWithEvent(env.momentumPropagator, env.orbitDepartureEvent, appendExtraInitialConditions(dynamicsModel, real(localArc.initialCondition), MBD.MOMENTUM), [0, localArc.TOF], dynamicsModel, [env.momentumPropagator, dynamicsModel, q_dep, env.EMMomentumDiff])
        t_orbitDeparture::Float64 = getTimeByIndex(orbitDepartureArc, -1)*tstar
        push!(arcs_tls[id], localArc)
        push!(cache_tls[id], DepCache(arc, t_orbitDeparture))
    end
    arcs::Vector{MBD.CR3BPManifoldArc} = reduce(vcat, arcs_tls)
    cache::Vector{DepCache} = reduce(vcat, cache_tls)
    println("Available departure arcs: $(length(arcs))\n")

    return arcs, cache
end

function computeArrivalArc_ext(env::MMATEnv, targeter, arrBody::String, family::String, JC::Float64, flip::Bool)
    dynamicsModel::MBD.CR3BPDynamicsModel = env.SArrDynamicsModel
    lstar::Float64 = env.charValues.SArr.lstar
    tstar::Float64 = env.charValues.SArr.tstar

    coarseOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(targeter, "FamilyData/CR3BPS$(arrBody[firstindex(arrBody)])$(family)s.csv", "JC", JC; JTol = 0.01)
    if flip
        if occursin("Halo", family)
            (coarseOrbit.initialCondition[3] *= -1) # Flip to northern halo
        end
    end
    solution::MBD.CR3BPMultipleShooterProblem = correct(targeter, coarseOrbit.initialCondition, [0, coarseOrbit.period], JC; JTol = 0.01)
    orbit = MBD.CR3BPPeriodicOrbit(dynamicsModel, solution.nodes[1].state.data[1:6], getPeriod(targeter, solution), getMonodromy(targeter, solution))
    println("Converged Sun-$arrBody orbit:\n\tIC:\t$(orbit.initialCondition)\n\tP:\t$(orbit.period)\n\tJC:\t$(getJacobiConstant(orbit))\n")

    posStableManifold::MBD.CR3BPManifold = getManifoldByArclength(orbit, "Stable", "Positive", 1000/lstar, 100)
    negStableManifold::MBD.CR3BPManifold = getManifoldByArclength(orbit, "Stable", "Negative", 1000/lstar, 100)
    posStableManifold.TOF, negStableManifold.TOF = -40.0, -40.0
    stableManifoldArcs::Vector{MBD.CR3BPManifoldArc} = vcat(stopCrashes(posStableManifold), stopCrashes(negStableManifold))
    nThreads::Int64 = Threads.maxthreadid()
    arcs_tls::Vector{Vector{MBD.CR3BPManifoldArc}} = [Vector{MBD.CR3BPManifoldArc}(undef, 0) for _ in 1:nThreads]
    periapses_tls::Vector{Vector{Float64}} = [Vector{Float64}(undef, 0) for _ in 1:nThreads]
    Threads.@threads for a::Int64 in eachindex(stableManifoldArcs)
        id::Int64 = Threads.threadid()
        manifoldArc::MBD.CR3BPManifoldArc = stableManifoldArcs[a]
        arc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.P2DistanceEvent, real(manifoldArc.initialCondition), [0, manifoldArc.TOF], dynamicsModel, [env.arrSoI])
        (abs(getTimeByIndex(arc, -1)) == abs(manifoldArc.TOF)) && continue
        localArc = MBD.CR3BPManifoldArc(manifoldArc.periodicOrbit, manifoldArc.orbitTime, manifoldArc.direction, manifoldArc.d, manifoldArc.initialCondition, getTimeByIndex(arc, -1))
        q_SoI_SI::Vector{Float64} = rotatingToPrimaryInertial(dynamicsModel, 1, [getStateByIndex(arc, -1)], [0.0])[1]
        @views Q_SoI_SI = similar(q_SoI_SI)
        Q_SoI_SI[1:3] .= q_SoI_SI[1:3].*lstar
        Q_SoI_SI[4:6] .= q_SoI_SI[4:6].*lstar./tstar
        oe::Vector{Float64} = getOrbitalElements(env.SDynamicsModel, Q_SoI_SI)
        if oe[2] < 1.0
            push!(arcs_tls[id], localArc)
            push!(periapses_tls[id], oe[1]*(1-oe[2]))
        end
    end
    arcs::Vector{MBD.CR3BPManifoldArc} = reduce(vcat, arcs_tls)
    periapses::Vector{Float64} = reduce(vcat, periapses_tls)
    println("Available arrival arcs: $(length(arcs))")

    minPeriapsisIndex::Int64 = argmin(periapses)
    arrivalArc::MBD.CR3BPManifoldArc = arcs[minPeriapsisIndex]
    println("Choosing minimum periapsis arrival arc ($minPeriapsisIndex): r_{p} = $(periapses[minPeriapsisIndex]) km")

    orbitArc::MBD.CR3BPArc = propagate(env.propagator, orbit.initialCondition, [0, arrivalArc.orbitTime*orbit.period], dynamicsModel)
    q_arr::Vector{Float64} = getStateByIndex(orbitArc, -1)
    orbitArrivalArc::MBD.CR3BPArc = propagateWithEvent(env.momentumPropagator, env.orbitDepartureEvent, appendExtraInitialConditions(dynamicsModel, real(arrivalArc.initialCondition), MBD.MOMENTUM), [0, arrivalArc.TOF], dynamicsModel, [env.momentumPropagator, dynamicsModel, q_arr, env.SArrMomentumDiff])
    t_orbitArrival::Float64 = -1*getTimeByIndex(orbitArrivalArc, -1)*tstar
    arrivalArcProp::MBD.CR3BPArc = propagate(env.propagator, real(arrivalArc.initialCondition), [0, arrivalArc.TOF], dynamicsModel)
    t_arrival::Float64 = -1*arrivalArc.TOF*tstar
    q_arr_SoI_SI_guess::Vector{Float64} = rotatingToSunEclipJ2000(dynamicsModel, env.frame, env.initialEpochTime, [getStateByIndex(arrivalArcProp, -1)], [0.0])[1]
    @views Q_arr_SoI_SI_guess = similar(q_arr_SoI_SI_guess)
    Q_arr_SoI_SI_guess[1:3] .= q_arr_SoI_SI_guess[1:3].*lstar
    Q_arr_SoI_SI_guess[4:6] .= q_arr_SoI_SI_guess[4:6].*lstar./tstar
    oe_arr_SoI_guess::Vector{Float64} = getOrbitalElements(env.SDynamicsModel, Q_arr_SoI_SI_guess)
    r_arr_l::Float64 = oe_arr_SoI_guess[1]*(1-oe_arr_SoI_guess[2])
    r_arr_u::Float64 = oe_arr_SoI_guess[1]*(1+oe_arr_SoI_guess[2])
    println("Arrival radius bounds: [$r_arr_l, $r_arr_u] km\n")

    return arrivalArc, ArrCache(arrivalArcProp, oe_arr_SoI_guess, r_arr_l, r_arr_u, t_arrival, t_orbitArrival)
end

function computeArrivalArc_int(env::MMATEnv, targeter, arrBody::String, family::String, JC::Float64, flip::Bool)
    dynamicsModel::MBD.CR3BPDynamicsModel = env.SArrDynamicsModel
    lstar::Float64 = env.charValues.SArr.lstar
    tstar::Float64 = env.charValues.SArr.tstar

    coarseOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(targeter, "FamilyData/CR3BPS$(arrBody[firstindex(arrBody)])$(family)s.csv", "JC", JC)
    if flip
        if occursin("Halo", family)
            (coarseOrbit.initialCondition[3] *= -1) # Flip to northern halo
        end
    end
    solution::MBD.CR3BPMultipleShooterProblem = correct(targeter, coarseOrbit.initialCondition, [0, coarseOrbit.period], JC)
    orbit = MBD.CR3BPPeriodicOrbit(dynamicsModel, solution.nodes[1].state.data[1:6], getPeriod(targeter, solution), getMonodromy(targeter, solution))
    println("Converged Sun-$arrBody orbit:\n\tIC:\t$(orbit.initialCondition)\n\tP:\t$(orbit.period)\n\tJC:\t$(getJacobiConstant(orbit))\n")

    posStableManifold::MBD.CR3BPManifold = getManifoldByArclength(orbit, "Stable", "Positive", 1000/lstar, 100)
    negStableManifold::MBD.CR3BPManifold = getManifoldByArclength(orbit, "Stable", "Negative", 1000/lstar, 100)
    posStableManifold.TOF, negStableManifold.TOF = -40.0, -40.0
    stableManifoldArcs::Vector{MBD.CR3BPManifoldArc} = vcat(stopCrashes(posStableManifold), stopCrashes(negStableManifold))
    nThreads::Int64 = Threads.maxthreadid()
    arcs_tls::Vector{Vector{MBD.CR3BPManifoldArc}} = [Vector{MBD.CR3BPManifoldArc}(undef, 0) for _ in 1:nThreads]
    apoapses_tls::Vector{Vector{Float64}} = [Vector{Float64}(undef, 0) for _ in 1:nThreads]
    Threads.@threads for a::Int64 in eachindex(stableManifoldArcs)
        id::Int64 = Threads.threadid()
        manifoldArc::MBD.CR3BPManifoldArc = stableManifoldArcs[a]
        arc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.P2DistanceEvent, real(manifoldArc.initialCondition), [0, manifoldArc.TOF], dynamicsModel, [env.arrSoI])
        (abs(getTimeByIndex(arc, -1)) == abs(manifoldArc.TOF)) && continue
        localArc = MBD.CR3BPManifoldArc(manifoldArc.periodicOrbit, manifoldArc.orbitTime, manifoldArc.direction, manifoldArc.d, manifoldArc.initialCondition, getTimeByIndex(arc, -1))
        q_SoI_SI::Vector{Float64} = rotatingToPrimaryInertial(dynamicsModel, 1, [getStateByIndex(arc, -1)], [0.0])[1]
        @views Q_SoI_SI = similar(q_SoI_SI)
        Q_SoI_SI[1:3] .= q_SoI_SI[1:3].*lstar
        Q_SoI_SI[4:6] .= q_SoI_SI[4:6].*lstar./tstar
        oe::Vector{Float64} = getOrbitalElements(env.SDynamicsModel, Q_SoI_SI)
        if oe[2] < 1.0
            push!(arcs_tls[id], localArc)
            push!(apoapses_tls[id], oe[1]*(1+oe[2]))
        end
    end
    arcs::Vector{MBD.CR3BPManifoldArc} = reduce(vcat, arcs_tls)
    apoapses::Vector{Float64} = reduce(vcat, apoapses_tls)
    println("Available arrival arcs: $(length(arcs))")

    maxApoapsisIndex::Int64 = argmax(apoapses)
    arrivalArc::MBD.CR3BPManifoldArc = arcs[maxApoapsisIndex]
    println("Choosing maximum apoapsis arrival arc ($maxApoapsisIndex): r_{a} = $(apoapses[maxApoapsisIndex]) km")

    orbitArc::MBD.CR3BPArc = propagate(env.propagator, orbit.initialCondition, [0, arrivalArc.orbitTime*orbit.period], dynamicsModel)
    q_arr::Vector{Float64} = getStateByIndex(orbitArc, -1)
    orbitArrivalArc::MBD.CR3BPArc = propagateWithEvent(env.momentumPropagator, env.orbitDepartureEvent, appendExtraInitialConditions(dynamicsModel, real(arrivalArc.initialCondition), MBD.MOMENTUM), [0, arrivalArc.TOF], dynamicsModel, [env.momentumPropagator, dynamicsModel, q_arr, env.SArrMomentumDiff])
    t_orbitArrival::Float64 = -1*getTimeByIndex(orbitArrivalArc, -1)*tstar
    arrivalArcProp::MBD.CR3BPArc = propagate(env.propagator, real(arrivalArc.initialCondition), [0, arrivalArc.TOF], dynamicsModel)
    t_arrival::Float64 = -1*arrivalArc.TOF*tstar
    q_arr_SoI_SI_guess::Vector{Float64} = rotatingToSunEclipJ2000(dynamicsModel, env.frame, env.initialEpochTime, [getStateByIndex(arrivalArcProp, -1)], [0.0])[1]
    @views Q_arr_SoI_SI_guess = similar(q_arr_SoI_SI_guess)
    Q_arr_SoI_SI_guess[1:3] .= q_arr_SoI_SI_guess[1:3].*lstar
    Q_arr_SoI_SI_guess[4:6] .= q_arr_SoI_SI_guess[4:6].*lstar./tstar
    oe_arr_SoI_guess::Vector{Float64} = getOrbitalElements(env.SDynamicsModel, Q_arr_SoI_SI_guess)
    r_arr_l::Float64 = oe_arr_SoI_guess[1]*(1-oe_arr_SoI_guess[2])
    r_arr_u::Float64 = oe_arr_SoI_guess[1]*(1+oe_arr_SoI_guess[2])
    println("Arrival radius bounds: [$r_arr_l, $r_arr_u] km\n")

    return arrivalArc, ArrCache(arrivalArcProp, oe_arr_SoI_guess, r_arr_l, r_arr_u, t_arrival, t_orbitArrival)
end

function getPhaseDiff(env::MMATEnv, oe_bridge_peri::Vector{Float64}, arrArcGuess::MBD.CR3BPManifoldArc, theta_arrBody_bridge::Float64, X::Vector{Float64})
    oe_bridge_int::Vector{Float64} = append!(oe_bridge_peri[1:5], X[1])
    t_bridgeConic::Float64 = solveKeplersEquation(env.SDynamicsModel, oe_bridge_int)
    arrivalArc::MBD.CR3BPManifoldArc = getManifoldArcByTime(arrArcGuess.periodicOrbit, "Stable", arrArcGuess.direction, arrArcGuess.d, mod(X[3], 1.0))
    arrivalArcProp::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.P2DistanceEvent, real(arrivalArc.initialCondition), [0, -40.0], env.SArrDynamicsModel, [env.arrSoI])
    q_arr_SoI::Vector{Float64} = getStateByIndex(arrivalArcProp, -1)
    q_arr_SoI_SI::Vector{Float64} = rotatingToSunEclipJ2000(env.SArrDynamicsModel, env.frame, env.initialEpochTime, [q_arr_SoI], [X[4]])[1]
    @views Q_arr_SoI_SI = similar(q_arr_SoI_SI)
    Q_arr_SoI_SI[1:3] .= q_arr_SoI_SI[1:3].*env.charValues.SArr.lstar
    Q_arr_SoI_SI[4:6] .= q_arr_SoI_SI[4:6].*env.charValues.SArr.lstar./env.charValues.SArr.tstar
    oe_arr_SoI::Vector{Float64} = getOrbitalElements(env.SDynamicsModel, Q_arr_SoI_SI)
    oe_arr_int::Vector{Float64} = append!(oe_arr_SoI[1:5], X[2])
    t_arrConic::Float64 = solveKeplersEquation(env.SDynamicsModel, oe_arr_SoI)-solveKeplersEquation(env.SDynamicsModel, oe_arr_int)
    theta_arrBody_SoI::Float64 = theta_arrBody_bridge+(t_bridgeConic+t_arrConic)/env.charValues.SArr.tstar

    return mod(X[4]-theta_arrBody_SoI+pi, 2.0*pi)-pi
end

function findIntersection_ext(targeter::MMATArrivalPhaseTargeter, env::MMATEnv, arrCache::ArrCache, oe_bridge_peri::Vector{Float64}, q_arr_SoI::Vector{Float64}, r_int::Float64, theta_bridge_int_guess::Float64, u_arr_guess::Float64, n::Int64, o::Int64)
    theta_arr_int_guess::Float64 = 2*o*pi+(-1)^o*acos((arrCache.oe_arr_SoI_guess[1]*(1-arrCache.oe_arr_SoI_guess[2]^2)/r_int-1)/arrCache.oe_arr_SoI_guess[2])
    omega_arr_guess::Float64 = (u_arr_guess-(theta_arr_int_guess+n*pi) > 0) ? u_arr_guess-(theta_arr_int_guess+n*pi) : 2*pi+u_arr_guess-(theta_arr_int_guess+n*pi)
    sigma_guess::Float64 = env.oe_arrBody[4]+omega_arr_guess+arrCache.oe_arr_SoI_guess[6]
    delta::Float64 = atan(q_arr_SoI[2], q_arr_SoI[1])
    X_guess::Vector{Float64} = [theta_bridge_int_guess, sigma_guess-(env.oe_arrBody[4]+env.oe_arrBody[5]+env.oe_arrBody[6])-delta, theta_arr_int_guess]
    X::Vector{Float64} = correct_ext(targeter, env, X_guess, oe_bridge_peri, q_arr_SoI)
    oe_bridge_int::Vector{Float64} = append!(oe_bridge_peri[1:5], X[1])
    Q_bridge_int_SI::Vector{Float64} = getCartesianState(env.SDynamicsModel, oe_bridge_int)
    t_bridgeConic::Float64 = solveKeplersEquation(env.SDynamicsModel, oe_bridge_int)
    q_arr_SoI_SI::Vector{Float64} = rotatingToSunEclipJ2000(env.SArrDynamicsModel, env.frame, env.initialEpochTime, [q_arr_SoI], [X[2]])[1]
    @views Q_arr_SoI_SI = similar(q_arr_SoI_SI)
    Q_arr_SoI_SI[1:3] .= q_arr_SoI_SI[1:3].*env.charValues.SArr.lstar
    Q_arr_SoI_SI[4:6] .= q_arr_SoI_SI[4:6].*env.charValues.SArr.lstar./env.charValues.SArr.tstar
    oe_arr_SoI::Vector{Float64} = getOrbitalElements(env.SDynamicsModel, Q_arr_SoI_SI)
    oe_arr_int::Vector{Float64} = append!(oe_arr_SoI[1:5], X[3])
    Q_arr_int_SI::Vector{Float64} = getCartesianState(env.SDynamicsModel, oe_arr_int)
    Deltav_2::Float64 = LinearAlgebra.norm(Q_arr_int_SI[4:6]-Q_bridge_int_SI[4:6])
    t_arrConic::Float64 = solveKeplersEquation(env.SDynamicsModel, oe_arr_SoI)-solveKeplersEquation(env.SDynamicsModel, oe_arr_int)
    (t_arrConic < 0) && (t_arrConic += MBD.getPeriod(env.SDynamicsModel, oe_arr_SoI))

    return [t_bridgeConic, Deltav_2, oe_arr_int..., t_arrConic, X[2]]
end

function findIntersection_int(targeter::MMATArrivalPhaseTargeter, env::MMATEnv, a_bridge::Float64, e_bridge::Float64, arrCache::ArrCache, oe_dep_SoI::Vector{Float64}, q_arr_SoI::Vector{Float64}, r_int::Float64, theta_dep_int_guess::Float64, u_bridge_guess::Float64, n::Int64, o::Int64)
    theta_bridge_int_guess::Float64 = 2*o*pi+(-1)^o*acos((a_bridge*(1-e_bridge^2)/r_int-1)/e_bridge)
    omega_bridge_guess::Float64 = (u_bridge_guess-(theta_bridge_int_guess+n*pi) > 0) ? u_bridge_guess-(theta_bridge_int_guess+n*pi) : 2*pi+u_bridge_guess-(theta_bridge_int_guess+n*pi)
    sigma_guess::Float64 = env.oe_arrBody[4]+omega_bridge_guess+arrCache.oe_arr_SoI_guess[6]
    delta::Float64 = atan(q_arr_SoI[2], q_arr_SoI[1])
    X_guess::Vector{Float64} = [theta_dep_int_guess, sigma_guess-(env.oe_arrBody[4]+env.oe_arrBody[5]+env.oe_arrBody[6])-delta, theta_bridge_int_guess]
    X::Vector{Float64} = correct_int(targeter, env, X_guess, oe_dep_SoI, q_arr_SoI, a_bridge, e_bridge)
    oe_dep_int::Vector{Float64} = append!(oe_dep_SoI[1:5], X[1])
    Q_dep_int_SI::Vector{Float64} = getCartesianState(env.SDynamicsModel, oe_dep_int)
    t_depConic::Float64 = solveKeplersEquation(env.SDynamicsModel, oe_dep_int)-solveKeplersEquation(env.SDynamicsModel, oe_dep_SoI)
    (t_depConic < 0) && (t_depConic += MBD.getPeriod(env.SDynamicsModel, oe_dep_SoI))
    q_arr_SoI_SI::Vector{Float64} = rotatingToSunEclipJ2000(env.SArrDynamicsModel, env.frame, env.initialEpochTime, [q_arr_SoI], [X[2]])[1]
    @views Q_arr_SoI_SI = similar(q_arr_SoI_SI)
    Q_arr_SoI_SI[1:3] .= q_arr_SoI_SI[1:3].*env.charValues.SArr.lstar
    Q_arr_SoI_SI[4:6] .= q_arr_SoI_SI[4:6].*env.charValues.SArr.lstar./env.charValues.SArr.tstar
    oe_arr_SoI::Vector{Float64} = getOrbitalElements(env.SDynamicsModel, Q_arr_SoI_SI)
    oe_bridge_peri::Vector{Float64} = [a_bridge, e_bridge, oe_arr_SoI[3], oe_arr_SoI[4], oe_arr_SoI[5], 0.0]
    oe_bridge_int::Vector{Float64} = append!(oe_bridge_peri[1:5], X[3])
    Q_bridge_int_SI::Vector{Float64} = getCartesianState(env.SDynamicsModel, oe_bridge_int)
    Deltav_1::Float64 = LinearAlgebra.norm(Q_bridge_int_SI[4:6]-Q_dep_int_SI[4:6])
    t_bridgeConic::Float64 = MBD.getPeriod(env.SDynamicsModel, oe_bridge_int)-solveKeplersEquation(env.SDynamicsModel, oe_bridge_int)
    Q_bridge_peri_SI::Vector{Float64} = getCartesianState(env.SDynamicsModel, oe_bridge_peri)
    oe_arr_peri::Vector{Float64} = append!(oe_arr_SoI[1:5], 0.0)
    Q_arr_peri_SI::Vector{Float64} = getCartesianState(env.SDynamicsModel, oe_arr_peri)
    Deltav_2::Float64 = LinearAlgebra.norm(Q_arr_peri_SI[4:6]-Q_bridge_peri_SI[4:6])
    t_arrConic::Float64 = solveKeplersEquation(env.SDynamicsModel, oe_arr_SoI)

    return [t_depConic, Deltav_1, oe_bridge_int..., t_bridgeConic, Deltav_2, oe_arr_peri..., t_arrConic, X[2]]
end

function computeMMAT_ext(targeter::MMATArrivalPhaseTargeter, env::MMATEnv, depCache::DepCache, arrCache::ArrCache, oe_bridge_peri::Vector{Float64}, r_int::Float64, theta_bridge_int_guess::Float64, u_arr_guess::Float64, n::Int64, o::Int64, t_depConic::Float64, t_departure::Float64)
    (n == 1) && (theta_bridge_int_guess += pi)
    intersect::Vector{Float64} = findIntersection_ext(targeter, env, arrCache, oe_bridge_peri, getStateByIndex(arrCache.arc, -1), r_int, theta_bridge_int_guess, u_arr_guess, n, o)
    theta_arrBody_arr::Float64 = intersect[10]+arrCache.t_arrival/env.charValues.SArr.tstar
    theta_arrBody_dep::Float64 = intersect[10]-(intersect[9]+intersect[1]+t_depConic+t_departure)/env.charValues.SArr.tstar
    while theta_arrBody_arr < 0.0
        theta_arrBody_arr += 2*pi
    end
    while theta_arrBody_dep < 0.0
        theta_arrBody_dep += 2*pi
    end
    TOF::Float64 = t_departure-depCache.t_orbitDeparture+t_depConic+intersect[1]+intersect[9]+arrCache.t_arrival-arrCache.t_orbitArrival
    
    return [intersect..., theta_arrBody_arr, theta_arrBody_dep, TOF]
end

function computeMMAT_int(targeter::MMATArrivalPhaseTargeter, env::MMATEnv, depCache::DepCache, a_bridge::Float64, e_bridge::Float64, arrCache::ArrCache, oe_dep_SoI::Vector{Float64}, r_int::Float64, theta_dep_int_guess::Float64, u_bridge_guess::Float64, n::Int64, o::Int64, t_departure::Float64)
    (n == 1) && (theta_dep_int_guess += pi)
    intersect::Vector{Float64} = findIntersection_int(targeter, env, a_bridge, e_bridge, arrCache, oe_dep_SoI, getStateByIndex(arrCache.arc, -1), r_int, theta_dep_int_guess, u_bridge_guess, n, o)
    theta_arrBody_arr::Float64 = intersect[18]+arrCache.t_arrival/env.charValues.SArr.tstar
    theta_arrBody_dep::Float64 = intersect[18]-(intersect[17]+intersect[9]+intersect[1]+t_departure)/env.charValues.SArr.tstar
    while theta_arrBody_arr < 0.0
        theta_arrBody_arr += 2*pi
    end
    while theta_arrBody_dep < 0.0
        theta_arrBody_dep += 2*pi
    end
    TOF::Float64 = t_departure-depCache.t_orbitDeparture+intersect[1]+intersect[9]+intersect[17]+arrCache.t_arrival-arrCache.t_orbitArrival
    
    return [intersect..., theta_arrBody_arr, theta_arrBody_dep, TOF]
end

function computeMMATs_ext(env::MMATEnv, targeter::MMATArrivalPhaseTargeter, depArcs::Vector{MBD.CR3BPManifoldArc}, depCache::Vector{DepCache}, arrArc::MBD.CR3BPManifoldArc, arrCache::ArrCache, mf::MATLAB.MatFile)
    println("Computing MMAT transfers...")
    for (idx::Int64, day::Float64) in pairs(env.days)
        e_dep::Float64 = 3600*24*day
        println("\tComputing transfers for $(env.epochs[idx])...")
        theta_E_dep::Float64 = env.oe_E[6]+e_dep/env.charValues.SE.tstar
        nThreads::Int64 = Threads.maxthreadid()
        results_tls::Vector{Vector{MMAT_ext}} = [Vector{MMAT_ext}() for _ in 1:nThreads]
        fails::Int64 = 0
        Threads.@threads for d::Int64 in eachindex(depArcs)
            id::Int64 = Threads.threadid()
            depArc::MBD.CR3BPManifoldArc = depArcs[d]
            cache::DepCache = depCache[d]
            q_dep_blend_EI_EM::Vector{Float64} = rotatingToPrimaryEclipJ2000(env.EMDynamicsModel, env.frame, env.initialEpochTime, [getStateByIndex(cache.arc, -1)], [getTimeByIndex(cache.arc, -1)+e_dep/env.charValues.EM.tstar])[1]
            @views Q_dep_blend_EI = similar(q_dep_blend_EI_EM)
            Q_dep_blend_EI[1:3] .= q_dep_blend_EI_EM[1:3].*env.charValues.EM.lstar
            Q_dep_blend_EI[4:6] .= q_dep_blend_EI_EM[4:6].*env.charValues.EM.lstar./env.charValues.EM.tstar
            @views q_dep_blend_EI_SE = similar(Q_dep_blend_EI)
            q_dep_blend_EI_SE[1:3] .= Q_dep_blend_EI[1:3]./env.charValues.SE.lstar
            q_dep_blend_EI_SE[4:6] .= Q_dep_blend_EI[4:6].*env.charValues.SE.tstar./env.charValues.SE.lstar
            e_blend_SE::Float64 = (getTimeByIndex(cache.arc, -1)*env.charValues.EM.tstar+e_dep)/env.charValues.SE.tstar
            q_dep_blend_SE::Vector{Float64} = secondaryEclipJ2000ToRotating(env.SEDynamicsModel, env.frame, env.initialEpochTime, [q_dep_blend_EI_SE], [e_blend_SE])[1]
            testProp::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.P2DistanceEvent, q_dep_blend_SE, [0, 4.0*pi], env.SEDynamicsModel, [env.Earth.bodyRadius/env.charValues.SE.lstar])
            testTime::Float64 = getTimeByIndex(testProp, -1)
            depArcSEProp::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.P2DistanceEvent, q_dep_blend_SE, [0, testTime], env.SEDynamicsModel, [env.EarthSoI])
            (abs(getTimeByIndex(depArcSEProp, -1)) == abs(testTime)) && continue
            t_departure::Float64 = depArc.TOF*env.charValues.EM.tstar+getTimeByIndex(depArcSEProp, -1)*env.charValues.SE.tstar

            q_dep_SoI_SI::Vector{Float64} = rotatingToSunEclipJ2000(env.SEDynamicsModel, env.frame, env.initialEpochTime, [getStateByIndex(depArcSEProp, -1)], [(e_dep+t_departure)/env.charValues.SE.tstar])[1]
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

            DeltaOmega_guess::Float64 = env.oe_arrBody[4]-oe_dep_SoI[4]
            s_i_arrBody::Float64, c_i_arrBody::Float64 = sincos(env.oe_arrBody[3])
            s_i_dep::Float64, c_i_dep::Float64 = sincos(oe_dep_SoI[3])
            s_DeltaOmega_guess::Float64, c_DeltaOmega_guess::Float64 = sincos(DeltaOmega_guess)
            psi_guess::Float64 = acos(c_i_arrBody*c_i_dep+s_i_arrBody*s_i_dep*c_DeltaOmega_guess)
            s_psi_guess::Float64, c_psi_guess::Float64 = sincos(psi_guess)
            s_u_bridge_guess::Float64 = sin(pi-env.oe_arrBody[3])*s_DeltaOmega_guess/s_psi_guess
            c_u_bridge_guess::Float64 = (s_i_arrBody*c_DeltaOmega_guess/c_i_dep-c_psi_guess*(s_i_dep/c_i_dep))/s_psi_guess
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
                    transfer_a::Vector{Float64} = computeMMAT_ext(targeter, env, cache, arrCache, oe_bridge_peri, r_int_0, theta_bridge_int_guess, u_arr_guess, 0, 0, t_depConic, t_departure)
                    push!(results_tls[id], MMAT_ext(0, day, 0, Deltav_1, depArc, depArcSEProp, oe_bridge_peri, oe_dep_SoI, transfer_a, t_depConic))
                catch
                    fails += 1
                end
                try
                    transfer_b::Vector{Float64} = computeMMAT_ext(targeter, env, cache, arrCache, oe_bridge_peri, r_int_0, theta_bridge_int_guess, u_arr_guess, 0, 1, t_depConic, t_departure)
                    push!(results_tls[id], MMAT_ext(1, day, 0, Deltav_1, depArc, depArcSEProp, oe_bridge_peri, oe_dep_SoI, transfer_b, t_depConic))
                catch
                    fails += 1
                    continue
                end
            end
            if (arrCache.r_arr_l < r_int_1 < arrCache.r_arr_u)
                try
                    transfer_a = computeMMAT_ext(targeter, env, cache, arrCache, oe_bridge_peri, r_int_1, theta_bridge_int_guess, u_arr_guess, 1, 0, t_depConic, t_departure)
                    push!(results_tls[id], MMAT_ext(0, day, 1, Deltav_1, depArc, depArcSEProp, oe_bridge_peri, oe_dep_SoI, transfer_a, t_depConic))
                catch
                    fails += 1
                end
                try
                    transfer_b = computeMMAT_ext(targeter, env, cache, arrCache, oe_bridge_peri, r_int_1, theta_bridge_int_guess, u_arr_guess, 1, 1, t_depConic, t_departure)
                    push!(results_tls[id], MMAT_ext(1, day, 1, Deltav_1, depArc, depArcSEProp, oe_bridge_peri, oe_dep_SoI, transfer_b, t_depConic))
                catch
                    fails += 1
                    continue
                end
            end
        end
        results::Vector{MMAT_ext} = reduce(vcat, results_tls)
        empty!.(results_tls)
        count0::Int64 = count1::Int64 = 0
        for r::MMAT_ext in results
            if r.type == 0
                count0 += 1
                label = Symbol("transfer_day", Int(day), "_type0", r.branch == 0 ? "a_" : "b_", count0)
            else
                count1 += 1
                label = Symbol("transfer_day", Int(day), "_type1", r.branch == 0 ? "a_" : "b_", count1)
            end
            exportCR3BPMMAT_ext(env, e_dep, r.transfer[1:10], theta_E_dep, r.transfer[11], r.transfer[12], r.depArc, r.depArcSEProp, r.oe_dep_SoI, r.t_depConic, r.oe_bridge_peri, arrArc, r.Deltav_1, r.transfer[13], mf, label)
        end
        println("\t\tTransfers found: $(count0+count1) ($fails failed)")
        empty!(results)
    end
end

function computeMMATs_int(env::MMATEnv, targeter::MMATArrivalPhaseTargeter, depArcs::Vector{MBD.CR3BPManifoldArc}, depCache::Vector{DepCache}, arrArc::MBD.CR3BPManifoldArc, arrCache::ArrCache, mf::MATLAB.MatFile)
    println("Computing MMAT transfers...")
    for (idx::Int64, day::Float64) in pairs(env.days)
        e_dep::Float64 = 3600*24*day
        println("\tComputing transfers for $(env.epochs[idx])...")
        theta_E_dep::Float64 = env.oe_E[6]+e_dep/env.charValues.SE.tstar
        nThreads::Int64 = Threads.maxthreadid()
        results_tls::Vector{Vector{MMAT_int}} = [Vector{MMAT_int}() for _ in 1:nThreads]
        fails::Int64 = 0
        Threads.@threads for d::Int64 in collect(1:5:length(depArcs))
            id::Int64 = Threads.threadid()
            depArc::MBD.CR3BPManifoldArc = depArcs[d]
            cache::DepCache = depCache[d]
            q_dep_blend_EI_EM::Vector{Float64} = rotatingToPrimaryEclipJ2000(env.EMDynamicsModel, env.frame, env.initialEpochTime, [getStateByIndex(cache.arc, -1)], [getTimeByIndex(cache.arc, -1)+e_dep/env.charValues.EM.tstar])[1]
            @views Q_dep_blend_EI = similar(q_dep_blend_EI_EM)
            Q_dep_blend_EI[1:3] .= q_dep_blend_EI_EM[1:3].*env.charValues.EM.lstar
            Q_dep_blend_EI[4:6] .= q_dep_blend_EI_EM[4:6].*env.charValues.EM.lstar./env.charValues.EM.tstar
            @views q_dep_blend_EI_SE = similar(Q_dep_blend_EI)
            q_dep_blend_EI_SE[1:3] .= Q_dep_blend_EI[1:3]./env.charValues.SE.lstar
            q_dep_blend_EI_SE[4:6] .= Q_dep_blend_EI[4:6].*env.charValues.SE.tstar./env.charValues.SE.lstar
            e_blend_SE::Float64 = (getTimeByIndex(cache.arc, -1)*env.charValues.EM.tstar+e_dep)/env.charValues.SE.tstar
            q_dep_blend_SE::Vector{Float64} = secondaryEclipJ2000ToRotating(env.SEDynamicsModel, env.frame, env.initialEpochTime, [q_dep_blend_EI_SE], [e_blend_SE])[1]
            testProp::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.P2DistanceEvent, q_dep_blend_SE, [0, 4.0*pi], env.SEDynamicsModel, [env.Earth.bodyRadius/env.charValues.SE.lstar])
            testTime::Float64 = getTimeByIndex(testProp, -1)
            depArcSEProp::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.P2DistanceEvent, q_dep_blend_SE, [0, testTime], env.SEDynamicsModel, [env.EarthSoI])
            (abs(getTimeByIndex(depArcSEProp, -1)) == abs(testTime)) && continue
            t_departure::Float64 = depArc.TOF*env.charValues.EM.tstar+getTimeByIndex(depArcSEProp, -1)*env.charValues.SE.tstar

            q_dep_SoI_SI::Vector{Float64} = rotatingToSunEclipJ2000(env.SEDynamicsModel, env.frame, env.initialEpochTime, [getStateByIndex(depArcSEProp, -1)], [(e_dep+t_departure)/env.charValues.SE.tstar])[1]
            @views Q_dep_SoI_SI = similar(q_dep_SoI_SI)
            Q_dep_SoI_SI[1:3] .= q_dep_SoI_SI[1:3].*env.charValues.SE.lstar
            Q_dep_SoI_SI[4:6] .= q_dep_SoI_SI[4:6].*env.charValues.SE.lstar./env.charValues.SE.tstar
            oe_dep_SoI::Vector{Float64} = getOrbitalElements(env.SDynamicsModel, Q_dep_SoI_SI)
            bridgeApoapsis::Float64 = oe_dep_SoI[1]*(1+oe_dep_SoI[2])
            bridgeRatio::Float64 = arrCache.r_arr_l/bridgeApoapsis
            e_bridge::Float64 = (1-bridgeRatio)/(1+bridgeRatio)
            a_bridge::Float64 = arrCache.r_arr_l/(1-e_bridge)

            DeltaOmega_guess::Float64 = env.oe_arrBody[4]-oe_dep_SoI[4]
            s_i_arrBody::Float64, c_i_arrBody::Float64 = sincos(env.oe_arrBody[3])
            s_i_dep::Float64, c_i_dep::Float64 = sincos(oe_dep_SoI[3])
            s_DeltaOmega_guess::Float64, c_DeltaOmega_guess::Float64 = sincos(DeltaOmega_guess)
            psi_guess::Float64 = acos(c_i_arrBody*c_i_dep+s_i_arrBody*s_i_dep*c_DeltaOmega_guess)
            s_psi_guess::Float64, c_psi_guess::Float64 = sincos(psi_guess)
            s_u_dep_guess::Float64 = sin(pi-env.oe_arrBody[3])*s_DeltaOmega_guess/s_psi_guess
            c_u_dep_guess::Float64 = (s_i_arrBody*c_DeltaOmega_guess/c_i_dep-c_psi_guess*(s_i_dep/c_i_dep))/s_psi_guess
            u_dep_guess::Float64 = atan(s_u_dep_guess, c_u_dep_guess)
            s_u_bridge_guess::Float64 = s_i_dep*s_DeltaOmega_guess/s_psi_guess
            c_u_bridge_guess::Float64 = c_DeltaOmega_guess*c_u_dep_guess+s_DeltaOmega_guess*s_u_dep_guess*c_i_dep
            u_bridge_guess::Float64 = atan(s_u_bridge_guess, c_u_bridge_guess)

            theta_dep_int_guess::Float64 = u_dep_guess-oe_dep_SoI[5]
            r_int_0::Float64 = oe_dep_SoI[1]*(1-oe_dep_SoI[2]^2)/(1+oe_dep_SoI[2]*cos(theta_dep_int_guess))
            r_int_1::Float64 = oe_dep_SoI[1]*(1-oe_dep_SoI[2]^2)/(1+oe_dep_SoI[2]*cos(theta_dep_int_guess+pi))

            if (a_bridge*(1-e_bridge) < r_int_0 < a_bridge*(1+e_bridge))
                try
                    transfer_a::Vector{Float64} = computeMMAT_int(targeter, env, cache, a_bridge, e_bridge, arrCache, oe_dep_SoI, r_int_0, theta_dep_int_guess, u_bridge_guess, 0, 0, t_departure)
                    push!(results_tls[id], MMAT_int(0, day, 0, depArc, depArcSEProp, oe_dep_SoI, transfer_a))
                catch
                    fails += 1
                end
                try
                    transfer_b::Vector{Float64} = computeMMAT_int(targeter, env, cache, a_bridge, e_bridge, arrCache, oe_dep_SoI, r_int_0, theta_dep_int_guess, u_bridge_guess, 0, 1, t_departure)
                    push!(results_tls[id], MMAT_int(1, day, 0, depArc, depArcSEProp, oe_dep_SoI, transfer_b))
                catch
                    fails += 1
                    continue
                end
            end
            if (a_bridge*(1-e_bridge) < r_int_1 < a_bridge*(1+e_bridge))
                try
                    transfer_a = computeMMAT_int(targeter, env, cache, a_bridge, e_bridge, arrCache, oe_dep_SoI, r_int_1, theta_dep_int_guess, u_bridge_guess, 1, 0, t_departure)
                    push!(results_tls[id], MMAT_int(0, day, 1, depArc, depArcSEProp, oe_dep_SoI, transfer_a))
                catch
                    fails += 1
                end
                try
                    transfer_b = computeMMAT_int(targeter, env, cache, a_bridge, e_bridge, arrCache, oe_dep_SoI, r_int_1, theta_dep_int_guess, u_bridge_guess, 1, 1, t_departure)
                    push!(results_tls[id], MMAT_int(1, day, 1, depArc, depArcSEProp, oe_dep_SoI, transfer_b))
                catch
                    fails += 1
                    continue
                end
            end
        end
        results::Vector{MMAT_int} = reduce(vcat, results_tls)
        empty!.(results_tls)
        count0::Int64 = count1::Int64 = 0
        for r::MMAT_int in results
            if r.type == 0
                count0 += 1
                label = Symbol("transfer_day", Int(day), "_type0", r.branch == 0 ? "a_" : "b_", count0)
            else
                count1 += 1
                label = Symbol("transfer_day", Int(day), "_type1", r.branch == 0 ? "a_" : "b_", count1)
            end
            exportCR3BPMMAT_int(env, e_dep, r.transfer[1:18], theta_E_dep, r.transfer[19], r.transfer[20], r.depArc, r.depArcSEProp, r.oe_dep_SoI, arrArc, r.transfer[21], mf, label)
        end
        println("\t\tTransfers found: $(count0+count1) ($fails failed)")
        empty!(results)
    end
end

function run_MMATCR3BP_ext(arrBody::String, EMFamily::String, EMJC::Float64, SArrFamily::String, SArrJC::Float64; EMFlip::Bool = false, SArrFlip::Bool = false)
    mf = MATLAB.MatFile("Output/MMAT/$(arrBody)MMATCR3BP_$(replace(string(EMJC), '.' => '_'))_$(EMFamily)_flip$(EMFlip)_$(replace(string(SArrJC), '.' => '_'))_$(SArrFamily)_flip$(SArrFlip).mat", "w")
    eph = Ephemerides.EphemerisProvider(["SPICEKernels/de430.bsp", "SPICEKernels/de440.bsp", "SPICEKernels/mar099s.bsp"])
    
    SPICE.furnsh("SPICEKernels/naif0012.tls")
    env::MMATEnv = setupEnvironment(eph, arrBody)
    SPICE.kclear()

    EMTargeter = targeterMap[EMFamily](env.EMDynamicsModel)
    depArcs::Vector{MBD.CR3BPManifoldArc}, depCache::Vector{DepCache} = computeDepartureArcs(env, EMTargeter, EMFamily, EMJC, EMFlip)

    SArrTargeter = targeterMap[SArrFamily](env.SArrDynamicsModel)
    arrArc::MBD.CR3BPManifoldArc, arrCache::ArrCache = computeArrivalArc_ext(env, SArrTargeter, arrBody, SArrFamily, SArrJC, SArrFlip)

    MMATTargeter = MMATArrivalPhaseTargeter(env.SDynamicsModel, sqrt(eps(Float64)))
    computeMMATs_ext(env, MMATTargeter, depArcs, depCache, arrArc, arrCache, mf)

    MATLAB.close(mf)
end

function run_MMATCR3BP_int(arrBody::String, EMFamily::String, EMJC::Float64, SArrFamily::String, SArrJC::Float64; EMFlip::Bool = false, SArrFlip::Bool = false)
    mf = MATLAB.MatFile("Output/MMAT/$(arrBody)MMATCR3BP_$(replace(string(EMJC), '.' => '_'))_$(EMFamily)_flip$(EMFlip)_$(replace(string(SArrJC), '.' => '_'))_$(SArrFamily)_flip$(SArrFlip).mat", "w")
    eph = Ephemerides.EphemerisProvider(["SPICEKernels/de430.bsp", "SPICEKernels/de440.bsp", "SPICEKernels/mar099s.bsp"])
    
    SPICE.furnsh("SPICEKernels/naif0012.tls")
    env::MMATEnv = setupEnvironment(eph, arrBody)
    SPICE.kclear()

    EMTargeter = targeterMap[EMFamily](env.EMDynamicsModel)
    depArcs::Vector{MBD.CR3BPManifoldArc}, depCache::Vector{DepCache} = computeDepartureArcs(env, EMTargeter, EMFamily, EMJC, EMFlip)

    SArrTargeter = targeterMap[SArrFamily](env.SArrDynamicsModel)
    arrArc::MBD.CR3BPManifoldArc, arrCache::ArrCache = computeArrivalArc_int(env, SArrTargeter, arrBody, SArrFamily, SArrJC, SArrFlip)

    MMATTargeter = MMATArrivalPhaseTargeter(env.SDynamicsModel, sqrt(eps(Float64)))
    computeMMATs_int(env, MMATTargeter, depArcs, depCache, arrArc, arrCache, mf)

    MATLAB.close(mf)
end

function populateJobs(depJCs::Vector{Float64})
    arrOrbits::Vector{String} = ["L1Lyapunov", "L1Halo", "L2Lyapunov", "L2Halo"]
    arrJCsMars::Vector{Vector{Float64}} = [[3.0002, 3.0001, 3.0], [3.000187, 3.00012, 3.00006], [3.0002, 3.0001, 3.0], [3.000186, 3.00012, 3.00007]]
    arrJCsVenus::Vector{Vector{Float64}} = [[3.0007, 3.00035, 3.0], [3.000714, 3.0005, 3.0003], [3.0007, 3.00035, 3.0], [3.000713, 3.0005, 3.0003]]
    cases::Vector{ArrivalCase} = [ArrivalCase(arrJCsMars, :ext, "Mars"), ArrivalCase(arrJCsVenus, :int, "Venus")]

    jobs = Vector{Tuple{String, Symbol, Float64, String, Float64, Bool}}()
    for depJC::Float64 in depJCs
        for (o::Int64, arrOrbit::String) in enumerate(arrOrbits)
            for case::ArrivalCase in cases
                for arrJC::Float64 in case.JCs[o]
                    push!(jobs, (case.planet, case.mode, depJC, arrOrbit, arrJC, false))
                    if occursin("Halo", arrOrbit)
                        push!(jobs, (case.planet, case.mode, depJC, arrOrbit, arrJC, true))
                    end
                end
            end
        end
    end

    return jobs
end

function MMATAnalysis(depOrbit::String, depFlip::Bool, jobs)
    failedJobs = Vector{Tuple{String, Symbol, Float64, String, Float64, Bool, Any}}()
    for (j::Int64, (planet::String, mode::Symbol, depJC::Float64, arrOrbit::String, arrJC::Float64, arrFlip::Bool)) in enumerate(jobs)
        if mode == :ext
            println("\nJob $j/$(length(jobs)): Running MMAT for exterior transfer from Earth-Moon $depJC $depOrbit (flip = $depFlip) to Sun-$planet $arrJC $arrOrbit (flip = $arrFlip):\n")
            try
                MMATCR3BP.run_MMATCR3BP_ext(planet, depOrbit, depJC, arrOrbit, arrJC; EMFlip = depFlip, SArrFlip = arrFlip)
            catch err
                println("\nJob failed")
                push!(failedJobs, (planet, mode, depJC, arrOrbit, arrJC, arrFlip, err))
            end
        else
            println("\nRunning MMAT for interior transfer from Earth-Moon $depJC $depOrbit (flip = $depFlip) to Sun-$planet $arrJC $arrOrbit (flip = $arrFlip):\n")
            try
                MMATCR3BP.run_MMATCR3BP_int(planet, depOrbit, depJC, arrOrbit, arrJC; EMFlip = depFlip, SArrFlip = arrFlip)
            catch err
                println("\nJob failed")
                push!(failedJobs, (planet, mode, depJC, arrOrbit, arrJC, arrFlip, err))
            end
        end
    end
    println("Jobs failed: $(length(failedJobs))")

    return failedJobs
end

end # MMATCR3BP
