"""
Script for computing CR3BP escape trajectories in the Earth-Moon system

Author: Jonathan LeFevre Richmond
C: 6/16/26
U: 9/24/26
"""

module EscCR3BP

using MBD, Clustering, DifferentialEquations, Distances, HDF5, JMcDM, LinearAlgebra, Logging, MATLAB, StaticArrays, Statistics

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../CR3BPTargeters/PlanarPerpJC.jl")
include("../Utilities/Export.jl")

struct CharValues
    lstar::Float64
    mstar::Float64
    tstar::Float64
end

struct EscEnv
    EDynamicsModel::MBD.KDynamicsModel
    EMDynamicsModel::MBD.CR3BPDynamicsModel
    EMEoMs::MBD.CR3BPEquationsOfMotion
    SEDynamicsModel::MBD.CR3BPDynamicsModel
    SMDynamicsModel::MBD.CR3BPDynamicsModel

    primaries::Vector{MBD.BodyData}
    Sun::MBD.BodyData

    charValues::NamedTuple

    apoapsisEvent
    arclengthEvent
    endEvents
    escapeEvent
    flybyEvent
    MoonEvent
    orbitTargeter
    periapsisEvent
    propagator::MBD.Propagator
    propagator_AL::MBD.Propagator
    propagator_STM::MBD.Propagator

    circHillv::Float64
    EarthHill_EM::Float64
    EarthRadius_EM::Float64
    MoonHill_EM::Float64
    MoonRadius_EM::Float64
end

function apoapsisCondition(state::Vector{Float64}, time::Float64, integrator)
    time-integrator.sol.prob.tspan[1] < 1E-6 ? -0.01 : LinearAlgebra.dot(state[1:3]-integrator.p[2], state[4:6])
end

function apseEndAffect!(integrator) # Need to find better way to do this
    if ((integrator.p[3] != 2) ? (integrator.u[1] < 0.74) : (integrator.u[1] > 0.92))
        DifferentialEquations.terminate!(integrator)
    end
end

function endAffect!(integrator, index)
    idx = findfirst(!=(0), index)
    if idx == 1
        if (index[idx] == 1)
            integrator.p[2][1].count += 1
            if integrator.p[10]
                push!(integrator.p[2][1].states, copy(integrator.u))
                push!(integrator.p[2][1].times, copy(integrator.t))
            end
        elseif index[idx] == -1
            integrator.p[2][2].count += 1
            if integrator.p[10]
                push!(integrator.p[2][2].states, copy(integrator.u))
                push!(integrator.p[2][2].times, copy(integrator.t))
            end
        end
    else
        integrator.p[2][idx+1].count += 1
        push!(integrator.p[2][idx+1].states, copy(integrator.u))
        push!(integrator.p[2][idx+1].times, copy(integrator.t))
        DifferentialEquations.terminate!(integrator)
    end
end

function endConditions(out::SubArray{Float64}, state::Vector{Float64}, time::Float64, integrator)
    if time-integrator.sol.prob.tspan[1] < 1E-6
        out[1] = (integrator.p[3] == :peri ? 0.01 : -0.01)
    else
        out[1] = LinearAlgebra.dot(state[1:3]-integrator.p[4], state[4:6])
    end
    out[2] = LinearAlgebra.norm(state[1:3])-integrator.p[7]
    out[3] = LinearAlgebra.norm(state[1:3]-integrator.p[5])-integrator.p[8]
    out[4] = LinearAlgebra.norm(state[1:3]-integrator.p[6])-integrator.p[9]
end

function escapeCondition(state::Vector{Float64}, time::Float64, integrator)
    LinearAlgebra.norm(state[1:3])-integrator.p[2]
end

function periapsisCondition(state::Vector{Float64}, time::Float64, integrator)
    time-integrator.sol.prob.tspan[1] < 1E-6 ? 0.01 : LinearAlgebra.dot(state[1:3]-integrator.p[2], state[4:6])
end

function xValueCondition(state::Vector{Float64}, time::Float64, integrator)
    state[1]-integrator.p[2]
end

function setupEnvironment()::EscEnv
    ESystemData = MBD.KSystemData("Earth")
    EMSystemData = MBD.CR3BPSystemData("Earth", "Moon")
    SESystemData = MBD.CR3BPSystemData("Sun", "Earth")
    SMSystemData = MBD.CR3BPSystemData("Sun", "Mars")
    EDynamicsModel = MBD.KDynamicsModel(ESystemData)
    EMDynamicsModel = MBD.CR3BPDynamicsModel(EMSystemData)
    SEDynamicsModel = MBD.CR3BPDynamicsModel(SESystemData)
    SMDynamicsModel = MBD.CR3BPDynamicsModel(SMSystemData)
    EMEoMs = MBD.CR3BPEquationsOfMotion(MBD.SIMPLE, EMDynamicsModel)

    primaries::Vector{MBD.BodyData} = EMSystemData.primaryData
    Sun::MBD.BodyData = SESystemData.primaryData[1]

    charValues::NamedTuple = (EM = CharValues(getCharLength(EMSystemData), getCharMass(EMSystemData), getCharTime(EMSystemData)),
                              SE = CharValues(getCharLength(SESystemData), getCharMass(SESystemData), getCharTime(SESystemData)),
                              SM = CharValues(getCharLength(SMSystemData), getCharMass(SMSystemData), getCharTime(SMSystemData)))
    
    propagator = MBD.Propagator()
    propagator_STM = MBD.Propagator(equationType = MBD.STM)
    propagator_AL = MBD.Propagator(equationType = MBD.ARCLENGTH)
    apoapsisEvent = DifferentialEquations.ContinuousCallback(apoapsisCondition, nothing, apseEndAffect!)
    arclengthEvent = DifferentialEquations.ContinuousCallback(arclengthCondition, terminateAffect!)
    endEvents = DifferentialEquations.VectorContinuousCallback(endConditions, endAffect!, 4)
    escapeEvent = DifferentialEquations.ContinuousCallback(escapeCondition, terminateAffect!)
    flybyEvent = DifferentialEquations.ContinuousCallback(p1CR3BPDistanceCondition, terminateAffect!)
    MoonEvent = DifferentialEquations.ContinuousCallback(xValueCondition, terminateAffect!)
    periapsisEvent = DifferentialEquations.ContinuousCallback(periapsisCondition, apseEndAffect!, nothing)
    orbitTargeter = PlanarPerpJCTargeter(EMDynamicsModel)

    EarthRadius_EM::Float64 = primaries[1].bodyRadius/charValues.EM.lstar
    MoonRadius_EM::Float64 = primaries[2].bodyRadius/charValues.EM.lstar
    EarthHill_EM::Float64 = charValues.SE.lstar*cbrt(getMassRatio(SESystemData)/3)/charValues.EM.lstar
    MoonHill_EM::Float64 = cbrt(getMassRatio(EMSystemData)/3)
    circHillv::Float64 = sqrt(2*primaries[1].gravParam/(EarthHill_EM*charValues.EM.lstar))

    return EscEnv(EDynamicsModel, EMDynamicsModel, EMEoMs, SEDynamicsModel, SMDynamicsModel, primaries, Sun, charValues, apoapsisEvent, arclengthEvent, endEvents, escapeEvent, flybyEvent, MoonEvent, orbitTargeter, periapsisEvent, propagator, propagator_AL, propagator_STM, circHillv, EarthHill_EM, EarthRadius_EM, MoonHill_EM, MoonRadius_EM)
end

function isInterior(q::AbstractVector{Float64}, mu::Float64, rE::AbstractVector{Float64}, rM::AbstractVector{Float64})
    r1::Float64 = LinearAlgebra.norm(q[1:2]-rE[1:2])
    r2::Float64 = LinearAlgebra.norm(q[1:2]-rM[1:2])
    dOmegadx::Float64 = q[1]-(1-mu)*(q[1]+mu)/r1^3-mu*(q[1]-1+mu)/r2^3
    dOmegady::Float64 = q[2]-(1-mu)*q[2]/r1^3-mu*q[2]/r2^3

    return q[1]*dOmegadx+q[2]*dOmegady <= 0
end

function getGrid(env::EscEnv, n::Int64, primary::Int64)
    if primary == 2
        radius::Float64 = 0.2
        center::Vector{Float64} = getPrimaryState(env.EMDynamicsModel, primary)[1:2]
    else
        radius = 1.0
        center = [0.0, 0.0]
    end

    xGrid::Vector{Float64} = collect(range(center[1]-radius, center[1]+radius, n))
    yGrid::Vector{Float64} = collect(range(center[2]-radius, center[2]+radius, n))
    rRect::Matrix{StaticArrays.SVector{2, Float64}} = [StaticArrays.SA[x, y] for x in xGrid, y in yGrid]

    mu::Float64 = getMassRatio(env.EMDynamicsModel)
    rE::StaticArrays.SVector{2, Float64} = StaticArrays.SVector{2, Float64}(getPrimaryState(env.EMDynamicsModel, 1)[1:2])
    rM::StaticArrays.SVector{2, Float64} = StaticArrays.SVector{2, Float64}(getPrimaryState(env.EMDynamicsModel, 2)[1:2])
    
    lunarMask::BitMatrix = LinearAlgebra.norm.(rRect .- Ref(rM)) .> env.MoonHill_EM
    mask::Matrix{Bool} = (primary == 2) ? .~lunarMask : (lunarMask .& isInterior.(rRect, mu, Ref(rE), Ref(rM)))

    rGrid::Vector{StaticArrays.SVector{2, Float64}} = rRect[mask]

    return rGrid
end

function isPeriapsis(env::EscEnv, primary::Int64, r::AbstractVector{Float64}, v::Vector{Float64})
    center::Vector{Float64} = (primary == 0 ? [0.0, 0.0] : getPrimaryState(env.EMDynamicsModel, primary)[1:2])
    d::Vector{Float64} = r-center
    mu::Float64 = getMassRatio(env.EMDynamicsModel)
    r1::Float64 = LinearAlgebra.norm(r-getPrimaryState(env.EMDynamicsModel, 1)[1:2])
    r2::Float64 = LinearAlgebra.norm(r-getPrimaryState(env.EMDynamicsModel, 2)[1:2])
    dOmegadx::Float64 = r[1]-(1-mu)*(r[1]+mu)/r1^3-mu*(r[1]-1+mu)/r2^3
    dOmegady::Float64 = r[2]-(1-mu)*r[2]/r1^3-mu*r[2]/r2^3
    xddot::Float64 = 2*v[2]+dOmegadx
    yddot::Float64 = -2*v[1]+dOmegady

    return v[1]^2+v[2]^2+LinearAlgebra.dot(d, [xddot, yddot]) > 0
end

function computeApseVelocities(env::EscEnv, JC::Float64, rGrid::Vector{StaticArrays.SVector{2, Float64}})
    OmegaGrid::Vector{Float64} = map(q -> getPseudopotential(env.EMDynamicsModel, Vector(push(q, 0.0))), rGrid)
    v2Grid::Vector{Float64} = 2 .* OmegaGrid .- JC

    return map(v2 -> (v2 < 0 ? NaN : sqrt(v2)), v2Grid)
end

function computeApseStates(env::EscEnv, primary::Int64, JC::Float64, apse::Symbol, grade::Symbol, rGrid::Vector{StaticArrays.SVector{2, Float64}})
    center::StaticArrays.SVector{2, Float64} = (primary == 0 ? StaticArrays.SA[0.0, 0.0] : StaticArrays.SVector{2, Float64}(getPrimaryState(env.EMDynamicsModel, primary)[1:2]))
    dGrid::Vector{StaticArrays.SVector{2, Float64}} = rGrid .- Ref(center)
    thatGrid::Vector{StaticArrays.SVector{2, Float64}} = map(r -> StaticArrays.SA[-r[2], r[1]] ./ norm(r), dGrid)
    vMagGrid::Vector{Float64} = computeApseVelocities(env, JC, rGrid)
    qGrid::Vector{StaticArrays.MVector{6, Float64}} = map(r -> StaticArrays.MVector{6, Float64}([r[1], r[2], 0.0, NaN, NaN, NaN]), rGrid)
    gradeSign::Float64 = (grade == :retro) ? -1.0 : 1.0
    Threads.@threads for j::Int64 in eachindex(qGrid)
        isnan(vMagGrid[j]) && continue
        v::Vector{Float64} = gradeSign*vMagGrid[j] .* thatGrid[j]
        peri::Bool = isPeriapsis(env, primary, rGrid[j], v)
        if ((apse == :peri) && !peri) || ((apse == :apo) && peri)
            qGrid[j][4] = Inf; qGrid[j][5] = Inf; qGrid[j][6] = Inf
            continue
        else
            qGrid[j][4] = v[1]; qGrid[j][5] = v[2]; qGrid[j][6] = 0.0
        end
    end

    return qGrid
end

function countApses(env::EscEnv, primary::Int64, apse::Symbol, IC::AbstractVector{Float64}, apses::Bool)
    peri = MBD.EventTracker(0, :peri, [], [])
    apo = MBD.EventTracker(0, :apo, [], [])
    escape = MBD.EventTracker(0, :escape, [], [])
    crashEarth = MBD.EventTracker(0, :earth, [], [])
    crashMoon = MBD.EventTracker(0, :moon, [], [])
    center::Vector{Float64} = (primary == 0 ? zeros(Float64, 3) : getPrimaryState(env.EMDynamicsModel, primary)[1:3])
    r_Earth::Vector{Float64} = getPrimaryState(env.EMDynamicsModel, 1)[1:3]
    r_Moon::Vector{Float64} = getPrimaryState(env.EMDynamicsModel, 2)[1:3]
    params::Vector{Any} = [apse, center, r_Earth, r_Moon, env.EarthHill_EM, env.EarthRadius_EM, env.MoonRadius_EM, apses]
    (_, eventTrackers::Vector{EventTracker}) = propagateWithEvents(env.propagator, env.endEvents, Vector(IC), [0, 12.0*pi], env.EMDynamicsModel, [peri, apo, escape, crashEarth, crashMoon], params)
    
    return eventTrackers
end

function apseMapCR3BP(env::EscEnv, JC::Float64, primary::Int64, rGrid::Vector{StaticArrays.SVector{2, Float64}}, mf::MATLAB.MatFile; apse::Symbol = :peri, grade::Symbol = :pro, apses::Bool = true)
    qGrid::Vector{StaticArrays.MVector{6, Float64}} = computeApseStates(env, primary, JC, apse, grade, rGrid)

    flags::Vector{Int64} = fill(9, size(qGrid))
    flags[map(q -> isinf(q[4]), qGrid)] .= 8
    valid::Vector{Bool} = map(q -> (!isnan(q[4]) && !isinf(q[4])), qGrid) 
    qProp::Vector{StaticArrays.MVector{6, Float64}} = qGrid[valid]
    qMap::Vector{Int64} = findall(valid)
    counts::Vector{Int64} = zeros(Int64, size(qGrid))
    periapses::Vector{Vector{StaticArrays.SVector{6, Float64}}} = [StaticArrays.SVector{6, Float64}[] for _ in qGrid]
    apoapses::Vector{Vector{StaticArrays.SVector{6, Float64}}} = [StaticArrays.SVector{6, Float64}[] for _ in qGrid]
    println("Propagating $(length(qProp)) CR3BP trajectories with $(Threads.nthreads()) threads...")
    Threads.@threads for j in eachindex(qProp)
        eventTrackers::Vector{MBD.EventTracker} = countApses(env, primary, apse, qProp[j], apses)
        apsesCount::Int64 = (apse == :peri ? eventTrackers[1].count : eventTrackers[2].count)
        if (eventTrackers[4].count != 0) || (eventTrackers[5].count != 0)
            flags[qMap[j]] = 7
        elseif eventTrackers[3].count != 0
            flags[qMap[j]] = min(apsesCount, 5)
        else
            flags[qMap[j]] = 6
        end
        counts[qMap[j]] = apsesCount
        if apses
            periapses[qMap[j]] = map(q -> StaticArrays.SVector{6, Float64}(q), eventTrackers[1].states)
            apoapses[qMap[j]] = map(q -> StaticArrays.SVector{6, Float64}(q), eventTrackers[2].states)
        end
    end

    flagCounts::Vector{Int64} = [count(==(f), flags) for f in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]]
    println("\tDirect escapes:\t$(flagCounts[1])")
    println("\tLonger escapes:\t$(sum(flagCounts[2:6]))")
    println("\tCaptures:\t$(flagCounts[7])")
    println("\tCrashes:\t$(flagCounts[8])")
    println("\tInvalid apses:\t$(flagCounts[9])")
    println("\tZVCs:\t\t$(flagCounts[10])")

    if apses
        println("Exporting map with apses...")
        validPeri::Vector{Int64} = findall(!isempty, periapses)
        validApo::Vector{Int64} = findall(!isempty, apoapses)
        periStates::Vector{StaticArrays.SVector{6, Float64}} = [StaticArrays.SVector{6, Float64}(state) for idx in validPeri for state in periapses[idx]]
        periIndices::Vector{Int64} = [idx for idx in validPeri for _ in 1:length(periapses[idx])]
        apoStates::Vector{StaticArrays.SVector{6, Float64}} = [StaticArrays.SVector{6, Float64}(state) for idx in validApo for state in apoapses[idx]]
        apoIndices::Vector{Int64} = [idx for idx in validApo for _ in 1:length(apoapses[idx])]
        exportCR3BPApseMap(env.EMDynamicsModel, primary, apse, grade, JC, qGrid, flags, counts, periStates, periIndices, apoStates, apoIndices, mf, Symbol("map_", replace(string(JC), "." => "_")))
    else
        println("Exporting map without apses...")
        exportCR3BPApseMap(env.EMDynamicsModel, primary, apse, grade, JC, qGrid, flags, counts, Vector{StaticArrays.SVector{6, Float64}}(), Vector{Int64}(), Vector{StaticArrays.SVector{6, Float64}}(), Vector{Int64}(), mf, Symbol("map_", replace(string(JC), "." => "_")))
    end
end

function pruneVolumeData(JCRange::Vector{Float64}, indices::Vector{Int64}, volFileName::String)
    HDF5.h5open(volFileName, "r") do file
        volJC = vec(read(file["JCVolume"]))
        JCStartIdx::Int64 = findfirst(j -> j <= JCRange[1], volJC)
        JCEndIdx::Int64 = findlast(j -> j >= JCRange[2], volJC)
        newVolJC::Vector{Float64} = volJC[JCStartIdx:JCEndIdx,1]

        nJC::Int64 = length(newVolJC)
        nIdx::Int64 = length(indices)

        volFlags = file["flagsVolume"]
        newVolFlags = Matrix{Int64}(undef, (nJC,nIdx))
        volqs = file["qVolume"]
        newVolqs = Array{Float64, 3}(undef, (6,nIdx,nJC))
        for (i::Int64, idx::Int64) in enumerate(indices)
            newVolFlags[:,i] = volFlags[JCStartIdx:JCEndIdx,idx]
            newVolqs[:,i,:] = volqs[:,idx,JCStartIdx:JCEndIdx]
        end

        return (newVolJC, newVolFlags, newVolqs)
    end
end

function getHohmannCost(env::EscEnv, q0::Vector{Float64})
    R2::Float64 = getEquilibriumPoint(env.SMDynamicsModel, 1)[1]*env.charValues.SM.lstar
    escArc::MBD.CR3BPArc = propagateWithEvent(env.propagator_STM, env.escapeEvent, appendExtraInitialConditions(env.EMDynamicsModel, q0, MBD.STM), [0, 12.0*pi], env.EMDynamicsModel, [env.EarthHill_EM])
    qf::Vector{Float64} = getStateByIndex(escArc, -1)
    (Eesc::Float64, _) = getMaxHeliocentricEnergy(env, qf)

    return abs(sqrt((8*Eesc^2*R2)/(env.Sun.gravParam-2*Eesc*R2))-sqrt(-2*Eesc))+abs(sqrt(env.Sun.gravParam/R2)-sqrt((2*env.Sun.gravParam^2)/(R2*(env.Sun.gravParam-2*Eesc*R2))))
end

function getEscapeEnergy(env::EscEnv, q0::Vector{Float64})
    arc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.escapeEvent, q0, [0, 12.0*pi], env.EMDynamicsModel, [env.EarthHill_EM])
    qf::Vector{Float64} = getStateByIndex(arc, -1)
    qf_I::Vector{Float64} = rotatingToPrimaryInertial(env.EMDynamicsModel, 1, [qf], [0.0])[1]
    Qf_I::Vector{Float64} = append!(qf_I[1:3].*env.charValues.EM.lstar, qf_I[4:6].*env.charValues.EM.lstar./env.charValues.EM.tstar)
    oef::Vector{Float64} = getOrbitalElements(env.primaries[1].gravParam, Qf_I)
    E::Float64 = getEnergy(env.EDynamicsModel, oef)
    v::Float64 = sqrt(2*(E+env.primaries[1].gravParam/(env.EarthHill_EM*env.charValues.EM.lstar)))

    return (E, v)
end

function getEscapeManeuverEfficiencyMetric(v1::Float64, v2::Float64, Deltavs::Vector{Float64})
    DeltavMags::Vector{Float64} = abs.(Deltavs)
    totalDeltav::Float64 = sum(DeltavMags)

    return (v2-v1)/totalDeltav
end

function getCaptureManeuverEfficiencyMetric(v1::Float64, v2::Float64, Deltavs::Vector{Float64})
    DeltavMags::Vector{Float64} = abs.(Deltavs)
    totalDeltav::Float64 = sum(DeltavMags)

    return v1/((v1-v2)*totalDeltav*1000)
end

function getEnergyGradient(env::EscEnv, q0::Vector{Float64})
    mu::Float64 = env.primaries[1].gravParam
    arc::MBD.CR3BPArc = propagateWithEvent(env.propagator_STM, env.escapeEvent, appendExtraInitialConditions(env.EMDynamicsModel, q0, MBD.STM), [0, 12.0*pi], env.EMDynamicsModel, [env.EarthHill_EM])
    qf::Vector{Float64} = getStateByIndex(arc, -1)
    qf_I::Vector{Float64} = rotatingToPrimaryInertial(env.EMDynamicsModel, 1, [qf[1:6]], [0.0])[1]
    Qf_I::Vector{Float64} = append!(qf_I[1:3].*env.charValues.EM.lstar, qf_I[4:6].*env.charValues.EM.lstar./env.charValues.EM.tstar)
    dEdQI::Matrix{Float64} = [mu*Qf_I[1:3]'./LinearAlgebra.norm(Qf_I[1:3])^3 Qf_I[4:6]';]
    lstar::Float64 = env.charValues.EM.lstar
    tstar::Float64 = env.charValues.EM.tstar
    dQIdqI::Matrix{Float64} = lstar.*LinearAlgebra.diagm([1.0, 1.0, 1.0, 1.0/tstar, 1.0/tstar, 1.0/tstar])
    dqIdqR::Matrix{Float64} = [LinearAlgebra.I zeros(Float64, (3,3)); [0 -1.0 0; 1.0 0 0; 0 0 0] LinearAlgebra.I]
    dqRdq0::Matrix{Float64} = getStateTransitionMatrix(env.EMDynamicsModel, qf)
    dq0dv0::Matrix{Float64} = [zeros(Float64, (3,3)); LinearAlgebra.I]
    qfdot::Vector{Float64} = zeros(Float64, 6)
    computeDerivatives!(qfdot, qf[1:6], (env.EMEoMs,), 0.0)
    dqRdtau::Matrix{Float64} = reshape(qfdot, (6,1))
    dgdqR::Matrix{Float64} = [Qf_I[1:3]'./LinearAlgebra.norm(Qf_I[1:3]) zeros(Float64, (1,3))]
    dtaudv0::Matrix{Float64} = -dgdqR*dqRdq0*dq0dv0./(dgdqR*dqRdtau)
    dv0dalpha::Matrix{Float64} = reshape(q0[4:6], (3,1))./LinearAlgebra.norm(q0[4:6])

    return only(dEdQI*dQIdqI*dqIdqR*(dqRdq0*dq0dv0+dqRdtau*dtaudv0)*dv0dalpha)
end

# function getEnergyGradientFull(env::EscEnv, q0::Vector{Float64})
#     mu::Float64 = env.primaries[1].gravParam
#     arc::MBD.CR3BPArc = propagateWithEvent(env.propagator_STM, env.escapeEvent, appendExtraInitialConditions(env.EMDynamicsModel, q0, MBD.STM), [0, 12.0*pi], env.EMDynamicsModel, [env.EarthHill_EM])
#     qf::Vector{Float64} = getStateByIndex(arc, -1)
#     qf_I::Vector{Float64} = rotatingToPrimaryInertial(env.EMDynamicsModel, 1, [qf[1:6]], [0.0])[1]
#     Qf_I::Vector{Float64} = append!(qf_I[1:3].*env.charValues.EM.lstar, qf_I[4:6].*env.charValues.EM.lstar./env.charValues.EM.tstar)
#     dEdQI::Matrix{Float64} = [mu*Qf_I[1:3]'./LinearAlgebra.norm(Qf_I[1:3])^3 Qf_I[4:6]';]
#     lstar::Float64 = env.charValues.EM.lstar
#     tstar::Float64 = env.charValues.EM.tstar
#     dQIdqI::Matrix{Float64} = lstar.*LinearAlgebra.diagm([1.0, 1.0, 1.0, 1.0/tstar, 1.0/tstar, 1.0/tstar])
#     dqIdqR::Matrix{Float64} = [LinearAlgebra.I zeros(Float64, (3,3)); [0 -1.0 0; 1.0 0 0; 0 0 0] LinearAlgebra.I]
#     dqRdq0::Matrix{Float64} = getStateTransitionMatrix(env.EMDynamicsModel, qf)
#     qfdot::Vector{Float64} = zeros(Float64, 6)
#     computeDerivatives!(qfdot, qf[1:6], (env.EMEoMs,), 0.0)
#     dqRdtau::Matrix{Float64} = reshape(qfdot, (6,1))
#     dgdqR::Matrix{Float64} = [Qf_I[1:3]'./LinearAlgebra.norm(Qf_I[1:3]) zeros(Float64, (1,3))]
#     dtaudq0::Matrix{Float64} = -dgdqR*dqRdq0./(dgdqR*dqRdtau)

#     return dEdQI*dQIdqI*dqIdqR*(dqRdq0+dqRdtau*dtaudq0)
# end

# function getEnergyGradientApo(env::EscEnv, primary::Int64, q0::Vector{Float64})
#     center::Vector{Float64} = (primary == 0 ? zeros(Float64, 3) : getPrimaryState(env.EMDynamicsModel, primary)[1:3])
#     periArc::MBD.CR3BPArc = propagateWithEvent(env.propagator_STM, env.periapsisEvent, appendExtraInitialConditions(env.EMDynamicsModel, q0, MBD.STM), [0, 4*pi], env.EMDynamicsModel, [center, primary])
#     qPeri::Vector{Float64} = getStateByIndex(periArc, -1)
#     grad2::Matrix{Float64} = getEnergyGradientFull(env, qPeri[1:6])
#     dqPeridq0::Matrix{Float64} = getStateTransitionMatrix(env.EMDynamicsModel, qPeri)
#     dq0dv0::Matrix{Float64} = [zeros(Float64, (3,3)); LinearAlgebra.I]
#     qPeridot::Vector{Float64} = zeros(Float64, 6)
#     computeDerivatives!(qPeridot, qPeri[1:6], (env.EMEoMs,), 0.0)
#     dqPeridtau::Matrix{Float64} = reshape(qPeridot, (6,1))
#     dgdqPeri::Matrix{Float64} = [qPeri[4:6]' qPeri[1:3]';]
#     dtaudv0::Matrix{Float64} = -dgdqPeri*dqPeridq0*dq0dv0./(dgdqPeri*dqPeridtau)
#     dv0dalpha::Matrix{Float64} = reshape(q0[4:6], (3,1))./LinearAlgebra.norm(q0[4:6])
#     dqPeridalpha::Matrix{Float64} = (dqPeridq0*dq0dv0+dqPeridtau*dtaudv0)*dv0dalpha

#     return only(grad2*dqPeridalpha)
# end

function getPeriluneDistance(env::EscEnv, q0::Vector{Float64})
    arc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.flybyEvent, q0, [0, 3.0*pi], env.EMDynamicsModel, [1.2])
    dMin::Float64 = 10.0*env.charValues.EM.lstar
    for s::Int64 in 1:getStateCount(arc)
        d::Float64 = getExcursion(env.EMDynamicsModel, 2, getStateByIndex(arc, s))
        (d < dMin) && (dMin = d)
    end

    return dMin
end

function discreteFrechet(traj1::Vector{Vector{Float64}}, traj2::Vector{Vector{Float64}})
    n1::Int64 = length(traj1)
    n2::Int64 = length(traj2)
    cost::Matrix{Float64} = zeros(Float64, (n1,n2))
    cost[1,1] = LinearAlgebra.norm(traj1[1]-traj2[1])
    for i::Int64 in 2:n1
        cost[i,1] = max(cost[i-1,1], LinearAlgebra.norm(traj1[i]-traj2[1]))
    end
    for j::Int64 in 2:n2
        cost[1,j] = max(cost[1,j-1], LinearAlgebra.norm(traj1[1]-traj2[j]))
    end
    for i::Int64 in 2:n1
        for j::Int64 in 2:n2
            minVal::Float64 = min(cost[i-1,j], cost[i,j-1], cost[i-1,j-1])
            cost[i,j] = max(minVal, LinearAlgebra.norm(traj1[i]-traj2[j]))
        end
    end

    return cost[n1,n2]
end

function optimizeForTransit(env::EscEnv, JC::Float64, primary::Int64, q0::Vector{Float64}, volJCs::Vector{Float64}, flags::Vector{Int64}, qs::Matrix{Float64})
    v::Float64 = sqrt(q0[4]^2+q0[5]^2)
    grad::Float64 = getEnergyGradient(env, q0)
    qPrev::Vector{Float64} = copy(q0)
    qNew::Vector{Float64} = copy(q0)
    Deltav::Float64 = 0.0
    JCNew::Float64 = JC
    JCIdx::Int64 = 0
    if primary == 2
        Deltav = 100.0
        JCNew = 0.0

        return (Deltav, JCNew)
    end
    if grad >= 0
        JCIdx = findfirst(x -> x < JC, volJCs)
        while (grad > 0) && (JCIdx <= length(volJCs)) && (flags[JCIdx] == 0)
            qPrev = copy(qNew)
            qNew = qs[:,JCIdx]
            grad = getEnergyGradient(env, qNew)
            JCIdx += 1
        end
        if JCIdx > length(volJCs)
            Deltav = 100.0
            JCNew = 0.0

            return (Deltav, JCNew)
        end
    else
        JCIdx = findlast(x -> x > JC, volJCs)
        while (grad < 0) && (JCIdx >= 1) && (flags[JCIdx] == 0)
            qPrev = copy(qNew)
            qNew = qs[:,JCIdx]
            grad = getEnergyGradient(env, qNew)
            JCIdx -= 1
        end
    end
    (EPrev::Float64, _) = getEscapeEnergy(env, qPrev)
    (ENew::Float64, _) = getEscapeEnergy(env, qNew)
    qOpt::Vector{Float64} = (EPrev > ENew) ? copy(qPrev) : copy(qNew)
    Deltav = sqrt(qOpt[4]^2+qOpt[5]^2)-v
    JCNew = getJacobiConstant(env.EMDynamicsModel, qOpt)

    return (Deltav, JCNew)
end

function feasibleEscape(JC::Float64, q0::Vector{Float64}, volJCs::Vector{Float64}, flags::Vector{Int64}, qs::Matrix{Float64})
    v::Float64 = sqrt(q0[4]^2+q0[5]^2)
    JCIdx::Int64 = findfirst(x -> x < JC, volJCs)
    qNew::Vector{Float64} = zeros(Float64, 6)
    Deltav::Float64 = 0.0
    JCNew::Float64 = 0.0
    prevIdx = findprev(f -> f < 6, vec(flags), JCIdx-1)
    if !isnothing(prevIdx)
        qNew = qs[:,prevIdx]
        Deltav = sqrt(qNew[4]^2+qNew[5]^2)-v
        JCNew = volJCs[prevIdx]
    else
        nextIdx = findnext(f -> f < 6, vec(flags), JCIdx)
        if !isnothing(nextIdx)
            qNew = qs[:,nextIdx]
            Deltav = sqrt(qNew[4]^2+qNew[5]^2)-v
            JCNew = volJCs[nextIdx]
        else
            Deltav = 100.0
            JCNew = 0.0
        end
    end

    return (Deltav, JCNew)
end

function findClosest(q::Vector{Float64}, flag::Int64, qs::Matrix{Float64}, flags::Vector{Int64})
    xSort::Vector{Float64} = sort(qs[1,:])
    ySort::Vector{Float64} = sort(qs[2,:])
    xIdx::Int64 = searchsortedlast(xSort, q[1])
    xIdx = min(xIdx, length(xSort)-1)
    yIdx::Float64 = searchsortedlast(ySort, q[2])
    yIdx = min(yIdx, length(ySort)-1)
    corners::Vector{Tuple{Int64, Int64}} = [(xIdx, yIdx), (xIdx+1, yIdx), (xIdx, yIdx+1), (xIdx+1, yIdx+1)]
    bestDist::Float64 = Inf
    bestIdx::Int64 = 0
    for (i::Int64, j::Int64) in corners
        idx = findfirst(k -> ((abs(qs[1,k]-xSort[i]) < 1E-5) && (abs(qs[2,k]-ySort[j]) < 1E-5)), eachindex(qs[1,:]))
        if flags[idx] == flag
            dx::Float64 = xSort[i]-q[1]
            dy::Float64 = ySort[j]-q[2]
            dist::Float64 = sqrt(dx^2+dy^2)
            if dist < bestDist
                bestDist = dist
                bestIdx = idx
            end
        end
    end
    if bestIdx == 0
        throw(ErrorException("No matching points"))
    end

    return (qs[:,bestIdx], bestIdx)
end

function findClosestNoFlag(q::Vector{Float64}, qs::Matrix{Float64})
    xSort::Vector{Float64} = sort(qs[1,:])
    ySort::Vector{Float64} = sort(qs[2,:])
    xIdx::Int64 = searchsortedlast(xSort, q[1])
    xIdx = min(xIdx, length(xSort)-1)
    yIdx::Float64 = searchsortedlast(ySort, q[2])
    yIdx = min(yIdx, length(ySort)-1)
    corners::Vector{Tuple{Int64, Int64}} = [(xIdx, yIdx), (xIdx+1, yIdx), (xIdx, yIdx+1), (xIdx+1, yIdx+1)]
    bestDist::Float64 = Inf
    bestIdx::Int64 = 0
    for (i::Int64, j::Int64) in corners
        idx = findfirst(k -> ((abs(qs[1,k]-xSort[i]) < 1E-5) && (abs(qs[2,k]-ySort[j]) < 1E-5)), eachindex(qs[1,:]))
        dx::Float64 = xSort[i]-q[1]
        dy::Float64 = ySort[j]-q[2]
        dist::Float64 = sqrt(dx^2+dy^2)
        if dist < bestDist
            bestDist = dist
            bestIdx = idx
        end
    end

    return (qs[:,bestIdx], bestIdx)
end

function trajOptimizeForTransit(env::EscEnv, JC::Float64, primary::Int64, q::Vector{Float64}, flag::Int64, qs::Matrix{Float64}, flags::Vector{Int64}, volFileName::String)
    if primary == 2
        Deltav = 100.0
        JCNew = 0.0

        return (Deltav, JCNew)
    end
    (qClosest::Vector{Float64}, idx::Int64) = findClosest(q, flag, qs, flags)
    JCRange::Vector{Float64} = [3.18, 2.8]
    (volJCs::Vector{Float64}, volFlags::Matrix{Int64}, volqs::Array{Float64, 3}) = pruneVolumeData(JCRange, [idx], volFileName)
    (_, JCNew::Float64) = optimizeForTransit(env, JC, primary, qClosest, volJCs, vec(volFlags), dropdims(volqs, dims = 2))
    r::Vector{Float64} = q[1:2]
    vMag::Float64 = computeApseVelocities(env, JCNew, [StaticArrays.SVector{2, Float64}(r)])[1]
    Deltav::Float64 = vMag-sqrt(q[4]^2+q[5]^2)

    return (Deltav, JCNew)
end

function updateTrajectory(env::EscEnv, JC::Float64, primary::Int64, q::Vector{Float64})
    vMag::Float64 = computeApseVelocities(env, JC, [StaticArrays.SVector{2, Float64}(q[1:2])])[1]
    Deltav::Float64 = vMag-sqrt(q[4]^2+q[5]^2)
    vhat::Vector{Float64} = q[4:6]./LinearAlgebra.norm(q[4:5])
    qNew::Vector{Float64} = append!(q[1:3], vMag .* vhat)
    eventTrackers::Vector{MBD.EventTracker} = countApses(env, primary, :peri, qNew, true)
    apsesCount::Int64 = eventTrackers[1].count
    flag::Int64 = 9
    if (eventTrackers[4].count != 0) || (eventTrackers[5].count != 0)
        flag = 7
    elseif eventTrackers[3].count != 0
        flag = min(apsesCount, 5)
    else
        flag = 6
    end

    return (Deltav, flag)
end

function trajFeasibleEscape(env::EscEnv, JC::Float64, primary::Int64, q::Vector{Float64}, flag::Int64, qs::Matrix{Float64}, flags::Vector{Int64}, volFileName::String)
    (qClosest::Vector{Float64}, idx::Int64) = findClosestNoFlag(q, qs)
    JCRange::Vector{Float64} = [3.18, 2.8]
    (volJCs::Vector{Float64}, volFlags::Matrix{Int64}, volqs::Array{Float64, 3}) = pruneVolumeData(JCRange, [idx], volFileName)
    (_, JCNew::Float64) = feasibleEscape(JC, qClosest, volJCs, vec(volFlags), dropdims(volqs, dims = 2))
    (Deltav::Float64, flag::Int64) = updateTrajectory(env, JCNew, primary, q)
    JCNewIdx::Int64 = findfirst(x -> x == JCNew, volJCs)
    iter::Int64 = 0
    step::Int64 = JCNew < JC ? 1 : -1
    while (flag > 5) && (JCNewIdx > 1)
        JCNewIdx += step
        !(1 <= JCNewIdx <= length(volJCs)) && break
        JCNew = volJCs[JCNewIdx]
        (Deltav, flag) = updateTrajectory(env, JCNew, primary, q)
        iter += 1
    end
    if flag > 5
        JCNewIdx = findfirst(x -> x < JC, volJCs)
        while (JCNewIdx < length(volJCs)) && (flag > 5)
            JCNewIdx += 1
            JCNew = volJCs[JCNewIdx]
            (Deltav, flag) = updateTrajectory(env, JCNew, primary, q)
        end
    end
    if flag > 5
        return (100.0, 0.0)
    else
        return (Deltav, JCNew)
    end
end

function escapeAnalysisCR3BP(env::EscEnv, JC::Float64, primary::Int64, apse::Symbol, flags::Vector{Int64}, qs::Matrix{Float64}, periapses::Matrix{Float64}, periapsesIndices::Vector{Int64}, volFileName::String, idx::Int64, indirectIdx::Int64, indirectFlag::Int64, mf::MATLAB.MatFile)
    esc0Indices::Vector{Int64} =  findall(flags .== 0)
    n_esc0::Int64 = length(esc0Indices)
    esc0q_0s::Matrix{Float64} = qs[:,esc0Indices]
    esc0t_fs::Vector{Float64} = Vector{Float64}(undef, n_esc0)
    esc0Es::Vector{Float64} = Vector{Float64}(undef, n_esc0)
    esc0dMins::Vector{Float64} = Vector{Float64}(undef, n_esc0)
    Threads.@threads for idx::Int64 in eachindex(esc0Indices)
        arc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.escapeEvent, qs[:,esc0Indices[idx]], [0, 12.0*pi], env.EMDynamicsModel, [env.EarthHill_EM])
        esc0t_fs[idx] = getTimeByIndex(arc, -1)
        (esc0Es[idx], _) = getEscapeEnergy(env, qs[:,esc0Indices[idx]])
        esc0dMins[idx] = getPeriluneDistance(env, qs[:,esc0Indices[idx]])
    end

    MATLAB.put_variable(mf, :esc0q0, esc0q_0s)
    MATLAB.put_variable(mf, :esc0tf, esc0t_fs)
    MATLAB.put_variable(mf, :esc0E, esc0Es)
    MATLAB.put_variable(mf, :esc0d, esc0dMins)

    esc1Indices::Vector{Int64} = findall(flags .== 1)
    n_esc1::Int64 = length(esc1Indices)
    esc1q_0s::Matrix{Float64} = qs[:,esc1Indices]
    esc1t_fs::Vector{Float64} = Vector{Float64}(undef, n_esc1)
    esc1Es::Vector{Float64} = Vector{Float64}(undef, n_esc1)
    Threads.@threads for idx::Int64 in eachindex(esc1Indices)
        arc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.escapeEvent, qs[:,esc1Indices[idx]], [0, 12.0*pi], env.EMDynamicsModel, [env.EarthHill_EM])
        esc1t_fs[idx] = getTimeByIndex(arc, -1)
        (esc1Es[idx], _) = getEscapeEnergy(env, qs[:,esc1Indices[idx]])
    end

    MATLAB.put_variable(mf, :esc1q0, esc1q_0s)
    MATLAB.put_variable(mf, :esc1tf, esc1t_fs)
    MATLAB.put_variable(mf, :esc1E, esc1Es)

    JCRange::Vector{Float64} = [3.18, 2.8]
    (volJCs::Vector{Float64}, volFlags::Matrix{Int64}, volqs::Array{Float64, 3}) = pruneVolumeData(JCRange, [idx], volFileName)

    @views vs = sqrt.(qs[4,:].^2 .+ qs[5,:].^2)
    Deltav2s::Vector{Float64} = fill(NaN, (length(vec(volFlags))))
    escvs::Vector{Float64} = copy(Deltav2s)
    # grads::Vector{Float64} = copy(Deltav2s)
    gammas::Vector{Float64} = copy(Deltav2s)
    (_, volIndirectFlags::Matrix{Int64}, volIndirectqs::Array{Float64, 3}) = pruneVolumeData(JCRange, [indirectIdx], volFileName)
    Deltav2bs::Vector{Float64} = fill(NaN, (length(vec(volIndirectFlags))))
    escvbs::Vector{Float64} = copy(Deltav2bs)
    gammabs::Vector{Float64} = copy(Deltav2bs)

    if flags[idx] < 6
        (_, escvOld::Float64) = getEscapeEnergy(env, qs[:,idx])
        
        if flags[idx] == 0
            escIndices::Vector{Int64} = findall(f -> f == 0, vec(volFlags))
            Threads.@threads for escIdx::Int64 in escIndices
                q::Vector{Float64} = volqs[:,1,escIdx]
                Deltav2s[escIdx] = sqrt(q[4]^2+q[5]^2)-vs[idx]
                (_, escvs[escIdx]) = getEscapeEnergy(env, q)
                # grads[escIdx] = getEnergyGradient(env, q)
                gammas[escIdx] = getEscapeManeuverEfficiencyMetric(escvOld, escvs[escIdx], [Deltav2s[escIdx]] .* env.charValues.EM.lstar ./ env.charValues.EM.tstar)
            end
        end

        escbIndices::Vector{Int64} = findall(f -> f == indirectFlag, vec(volIndirectFlags))
        Threads.@threads for escbIdx::Int64 in escbIndices
            q::Vector{Float64} = volIndirectqs[:,1,escbIdx]
            Deltav2bs[escbIdx] = sqrt(q[4]^2+q[5]^2)-vs[indirectIdx]
            (_, escvbs[escbIdx]) = getEscapeEnergy(env, q)
            gammabs[escbIdx] = getEscapeManeuverEfficiencyMetric(escvOld, escvbs[escbIdx], [Deltav2bs[escbIdx]] .* env.charValues.EM.lstar ./ env.charValues.EM.tstar)
        end
    end

    MATLAB.put_variable(mf, :Deltav2s, Deltav2s)
    MATLAB.put_variable(mf, :Escapevs, escvs)
    # MATLAB.put_variable(mf, :DeltaEs, grads)
    MATLAB.put_variable(mf, :metrics, gammas)
    MATLAB.put_variable(mf, :Deltav2bs, Deltav2bs)
    MATLAB.put_variable(mf, :Escapevbs, escvbs)
    MATLAB.put_variable(mf, :metricbs, gammabs)

    """Maneuver Sequencing"""
    qTraj::Vector{Float64} = qs[:,idx]
    vMag::Float64 = 0
    vhat::Vector{Float64} = zeros(Float64, 3)
    escE::Float64 = NaN
    escv::Float64 = NaN
    Deltav1::Float64 = 0
    Deltav2::Float64 = 0
    qNew::Vector{Float64} = zeros(Float64, 6)
    JCNew::Float64 = JC
    escENew::Float64 = NaN
    escvNew::Float64 = NaN
    if apse == :peri
        if flags[idx] == 0
            println("Direct: $(flags[idx])")
            (escE, escv) = getEscapeEnergy(env, qTraj)
            (Deltav2, JCNew) = optimizeForTransit(env, JC, primary, qTraj, volJCs, volFlags[:,1], volqs[:,1,:])
            vMag = LinearAlgebra.norm(qTraj[4:6])
            vhat = qTraj[4:6]./vMag
            qNew = append!(qTraj[1:3], (vMag+Deltav2) .* vhat)
            (escENew, escvNew) = getEscapeEnergy(env, qNew)
        elseif (0 < flags[idx] < 6)
            println("Indirect: $(flags[idx])")
            (escE, escv) = getEscapeEnergy(env, qTraj)
            periIndices::Vector{Int64} = findall(i -> i == idx, periapsesIndices)
            qPeri::Vector{Float64} = periapses[:,periIndices[end]]
            (Deltav2, JCNew) = trajOptimizeForTransit(env, primary, JC, qPeri, 0, qs, flags, volFileName)
            vMag = LinearAlgebra.norm(qPeri[4:6])
            vhat = qPeri[4:6]./vMag
            qNew = append!(qPeri[1:3], (vMag+Deltav2) .* vhat)
            (escENew, escvNew) = getEscapeEnergy(env, qNew)
        elseif (flags[idx] == 6) || (flags[idx] == 7)
            println("Failure: $(flags[idx])")
            (Deltav1, JCNew) = feasibleEscape(JC, qTraj, volJCs, volFlags[:,1], volqs[:,1,:])
            (Deltav1 == 100.0) && println("High energy maneuver")
            vMag = LinearAlgebra.norm(qTraj[4:6])
            vhat = qTraj[4:6]./vMag
            qNew = append!(qTraj[1:3], (vMag+Deltav1) .* vhat)
            (escENew, escvNew) = getEscapeEnergy(env, qNew)
        end
    elseif apse == :apo
        if flags[idx] == 0
            println("Direct: $(flags[idx])")
        elseif (0 < flags[idx] < 6)
            println("Indirect: $(flags[idx])")
        elseif (flags[idx] == 6) || (flags[idx] == 7)
            println("Failure: $(flags[idx])")
            (Deltav1, JCNew) = feasibleEscape(JC, qTraj, volJCs, volFlags[:,1], volqs[:,1,:])
            (Deltav1 == 100.0) && println("High energy maneuver")
            vMag = LinearAlgebra.norm(qTraj[4:6])
            vhat = qTraj[4:6]./vMag
            qNew = append!(qTraj[1:3], (vMag+Deltav1) .* vhat)
            (escENew, escvNew) = getEscapeEnergy(env, qNew)
        end
    end

    println("Old energy: $escE")
    println("Old velocity: $escv")
    println("Delta-v 1: $Deltav1")
    println("Delta-v 2: $Deltav2")
    println("New JC: $JCNew")
    println("New energy: $escENew")
    println("New velocity: $escvNew")

    MATLAB.put_variable(mf, :escE, escE)
    MATLAB.put_variable(mf, :escv, escv)
    MATLAB.put_variable(mf, :Deltav1, Deltav1)
    MATLAB.put_variable(mf, :Deltav2, Deltav2)
    MATLAB.put_variable(mf, :qMan, qNew)
    MATLAB.put_variable(mf, :newJC, JCNew)
    MATLAB.put_variable(mf, :newEscE, escENew)
    MATLAB.put_variable(mf, :newEscv, escvNew)
end

function assistedEscapeAnalysisCR3BP(env::EscEnv, JC::Float64, primary::Int64, apse::Symbol, maneuver::Int64, flags::Vector{Int64}, qs::Matrix{Float64}, volFileName::String, mf::MATLAB.MatFile)
    indices::Vector{Int64} = collect(1:length(flags))
    JCRange::Vector{Float64} = [3.18, 2.8]
    (volJCs::Vector{Float64}, volFlags::Matrix{Int64}, volqs::Array{Float64, 3}) = pruneVolumeData(JCRange, indices, volFileName)
    Deltav1s::Vector{Float64} = zeros(Float64, size(qs, 2))
    Deltav2s::Vector{Float64} = copy(Deltav1s)
    JCNews::Vector{Float64} = JC .* ones(Float64, size(qs, 2))
    Threads.@threads for j::Int64 in eachindex(indices)
        idx::Int64 = indices[j]
        qTraj::Vector{Float64} = qs[:,idx]
        if apse == :peri
            if (flags[idx] == 0)
                if maneuver == 2
                    (Deltav2s[idx], JCNews[idx]) = optimizeForTransit(env, JC, primary, qTraj, volJCs, volFlags[:,j], volqs[:,j,:])
                end
            elseif (0 < flags[idx] < 6)
            elseif (flags[idx] == 6) || (flags[idx] == 7)
                if maneuver == 1
                    (Deltav1s[idx], JCNews[idx]) = feasibleEscape(JC, qTraj, volJCs, volFlags[:,j], volqs[:,j,:])
                end
            else
                Deltav1s[idx] = NaN
                Deltav2s[idx] = NaN
                JCNews[idx] = NaN
            end
        elseif apse == :apo
            if flags[idx] == 0
            elseif (0 < flags[idx] < 6)
            elseif (flags[idx] == 6) || (flags[idx] == 7)
                if maneuver == 1
                    (Deltav1s[idx], JCNews[idx]) = feasibleEscape(JC, qTraj, volJCs, volFlags[:,j], volqs[:,j,:])
                end
            else
                Deltav1s[idx] = NaN
                Deltav2s[idx] = NaN
                JCNews[idx] = NaN
            end
        end
    end

    MATLAB.put_variable(mf, :Deltav1s, Deltav1s)
    MATLAB.put_variable(mf, :Deltav2s, Deltav2s)
    MATLAB.put_variable(mf, :newJCs, JCNews)
end

function getPeriapsisStates(env::EscEnv, primary::Int64, orbit::MBD.CR3BPPeriodicOrbit, mf::MATLAB.MatFile)
    posUnstableManifold::MBD.CR3BPManifold = getManifoldByArclength(orbit, "Unstable", "Positive", 25/env.charValues.EM.lstar, 100)
    negUnstableManifold::MBD.CR3BPManifold = getManifoldByArclength(orbit, "Unstable", "Negative", 25/env.charValues.EM.lstar, 100)
    posUnstableManifold.TOF, negUnstableManifold.TOF = 4.0*pi, 4.0*pi
    unstableManifoldArcs::Vector{MBD.CR3BPManifoldArc} = vcat(stopCrashes(posUnstableManifold), stopCrashes(negUnstableManifold))
    center::Vector{Float64} = (primary == 0 ? zeros(Float64, 3) : getPrimaryState(env.EMDynamicsModel, primary)[1:3])
    q0s::Matrix{Float64} = Matrix{Float64}(undef, (6,length(unstableManifoldArcs)))
    qPeris::Matrix{Float64} = Matrix{Float64}(undef, (6, length(unstableManifoldArcs)))
    tPeris::Vector{Float64} = Vector{Float64}(undef, length(unstableManifoldArcs))
    proPeris::Vector{Bool} = Vector{Bool}(undef, length(unstableManifoldArcs))
    Threads.@threads for m::Int64 in eachindex(unstableManifoldArcs)
        q0::Vector{Float64} = real(unstableManifoldArcs[m].initialCondition)
        q0s[:,m] = q0
        periArc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.periapsisEvent, q0, [0, unstableManifoldArcs[m].TOF], env.EMDynamicsModel, [center, primary])
        qPeri::Vector{Float64} = getStateByIndex(periArc, -1)
        tPeri::Float64 = getTimeByIndex(periArc, -1)
        h::Vector{Float64} = cross(qPeri[1:3], qPeri[4:6])
        if (tPeri < unstableManifoldArcs[m].TOF)
            qPeris[:,m] = qPeri
            tPeris[m] = tPeri
            proPeris[m] = (h[3] > 0)
        else
            qPeris[:,m] = fill(NaN, 6)
            tPeris[m] = NaN
        end
    end

    MATLAB.put_variable(mf, :initialStates, q0s)
    MATLAB.put_variable(mf, :periStates, qPeris)
    MATLAB.put_variable(mf, :periTimes, tPeris)
    MATLAB.put_variable(mf, :periProgrades, proPeris)

    return (q0s, qPeris, tPeris)
end

function trajAssistedEscapeAnalysisCR3BP(env::EscEnv, JC::Float64, primary::Int64, qPeri::Vector{Float64}, qs::Matrix{Float64}, flags::Vector{Int64}, volFileName::String, mf::MATLAB.MatFile)
    eventTrackers::Vector{MBD.EventTracker} = countApses(env, primary, :peri, qPeri, true)
    apsesCount::Int64 = eventTrackers[1].count
    flag::Int64 = 9
    if (eventTrackers[4].count != 0) || (eventTrackers[5].count != 0)
        flag = 7
    elseif eventTrackers[3].count != 0
        flag = min(apsesCount, 5)
    else
        flag = 6
    end
    periapses::Vector{Vector{Float64}} = append!([qPeri], eventTrackers[1].states)
    periapsesTimes::Vector{Float64} = append!([0.0], eventTrackers[1].times)
    apoapses::Vector{Vector{Float64}} = eventTrackers[2].states

    q0::Vector{Float64} = copy(qPeri)
    newFlag::Int64 = 10
    vMag::Float64 = 0.0
    vhat::Vector{Float64} = zeros(Float64, 3)
    qEsc::Vector{Float64} = copy(q0)
    escE::Float64 = NaN
    escv::Float64 = NaN
    Deltav1::Float64 = 0.0
    Deltav2::Float64 = 0.0
    q1::Vector{Float64} = copy(q0)
    t1::Float64 = 0.0
    q2::Vector{Float64} = Vector{Float64}(undef, 6)
    t2::Float64 = 0.0
    JCNew::Float64 = NaN
    escENew::Float64 = NaN
    escvNew::Float64 = NaN
    gammaNew::Float64 = NaN
    if flag == 0
        println("Direct: $flag")
        (escE, escv) = getEscapeEnergy(env, q0)
        (Deltav2, JCNew) = trajOptimizeForTransit(env, JC, primary, q1, 0, qs, flags, volFileName)
        vMag = LinearAlgebra.norm(qPeri[4:6])
        vhat = qPeri[4:6]./vMag
        q2 = append!(q1[1:3], (vMag+Deltav2) .* vhat)
        (escENew, escvNew) = getEscapeEnergy(env, q2)
        gammaNew = getEscapeManeuverEfficiencyMetric(escv, escvNew, [Deltav2] .* env.charValues.EM.lstar ./ env.charValues.EM.tstar)
    elseif (0 < flag < 6)
        println("Indirect: $flag")
        (escE, escv) = getEscapeEnergy(env, q0)
        qEsc = periapses[end]
        t2 = periapsesTimes[end]
        (Deltav2, JCNew) = trajOptimizeForTransit(env, JC, primary, qEsc, 0, qs, flags, volFileName)
        vMag = LinearAlgebra.norm(qEsc[4:6])
        vhat = qEsc[4:6]./vMag
        q2 = append!(qEsc[1:3], (vMag+Deltav2) .* vhat)
        (escENew, escvNew) = getEscapeEnergy(env, q2)
        gammaNew = getEscapeManeuverEfficiencyMetric(escv, escvNew, [Deltav2] .* env.charValues.EM.lstar ./ env.charValues.EM.tstar)
    elseif (flag == 6) || (flag == 7)
        println("Failure: $flag")
        t1s::Vector{Float64} = Vector{Float64}(undef, length(periapses))
        Deltav1s::Vector{Float64} = Vector{Float64}(undef, length(periapses))
        JC1s::Vector{Float64} = Vector{Float64}(undef, length(periapses))
        q1s::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(periapses))
        escE1s::Vector{Float64} = Vector{Float64}(undef, length(periapses))
        escv1s::Vector{Float64} = Vector{Float64}(undef, length(periapses))
        gammas::Vector{Float64} = Vector{Float64}(undef, length(periapses))
        midFlags::Vector{Int64} = Vector{Int64}(undef, length(periapses))
        for p::Int64 in eachindex(periapses)
            qPeri::Vector{Float64} = periapses[p]
            t1s[p] = periapsesTimes[p]
            if (primary == 2) && (LinearAlgebra.norm(qPeri[1:2]-getPrimaryState(env.EMDynamicsModel, 2)[1:2]) > env.MoonHill_EM)
                Deltav1s[p] = Inf
                JC1s[p] = 0.0
                q1s[p] = append!(qPeri[1:3], (vMag+Deltav1s[p]) .* vhat)
                escE1s[p] = NaN
                escv1s[p] = NaN
                gammas[p] = NaN
                midFlags[p] = flag
                continue
            end
            (Deltav1s[p], JC1s[p]) = trajFeasibleEscape(env, JC, primary, qPeri, flag, qs, flags, volFileName)
            (Deltav1s[p] == 100.0) && println("Periapsis $p is a high-energy maneuver")
            vMag = LinearAlgebra.norm(qPeri[4:6])
            vhat = qPeri[4:6]./vMag
            q1s[p] = append!(qPeri[1:3], (vMag+Deltav1s[p]) .* vhat)
            (escE1s[p], escv1s[p]) = getEscapeEnergy(env, q1s[p])
            gammas[p] = getCaptureManeuverEfficiencyMetric(env.circHillv, escv1s[p], [Deltav1s[p]] .* env.charValues.EM.lstar ./ env.charValues.EM.tstar)
            (_, midFlags[p]) = updateTrajectory(env, JC1s[p], primary, q1s[p])
        end
        filt::Vector{Int64} = filter(p -> Deltav1s[p] < 100.0, eachindex(periapses))
        for p::Int64 in eachindex(periapses)[filt]
            println("Peri $p: Deltav = $(Deltav1s[p]) | gamma = $(gammas[p]) | t = $(t1s[p]) | TOF = $(t1s[p]+pi*midFlags[p])")
        end
        isempty(filt) && throw(ErrorException("No low-energy maneuvers available"))
        decision::Matrix{Float64} = [gammas[filt] t1s[filt]+pi .* midFlags[filt]]
        weights::Vector{Float64} = [0.5, 0.5]
        functions::Vector{Function} = [maximum, minimum]
        opt::JMcDM.SawResult = JMcDM.saw(decision, weights, functions)
        bestIdx::Int64 = opt.bestIndex
        println("Chose periapsis $(filt[bestIdx])")
        Deltav1 = Deltav1s[filt[bestIdx]]
        q1 = q1s[filt[bestIdx]]
        t1 = t1s[filt[bestIdx]]
        JC1::Float64 = JC1s[filt[bestIdx]]
        escE1::Float64 = escE1s[filt[bestIdx]]
        escv1::Float64 = escv1s[filt[bestIdx]]
        gamma::Float64 = gammas[filt[bestIdx]]
        midFlag::Int64 = midFlags[filt[bestIdx]]
        println("Intermediate JC: $JC1")
        println("Intermediate energy: $escE1")
        println("Intermediate velocity: $escv1")
        println("Intermediate gamma: $gamma")
        println("Intermediate flag: $midFlag")
        (_, volFlags::Matrix{Int64}, volqs::Array{Float64, 3}) = pruneVolumeData([JC1, JC1], collect(eachindex(flags)), volFileName)
        newqs::Matrix{Float64} = dropdims(volqs, dims = 3)
        newFlags::Vector{Int64} = vec(volFlags)
        newEventTrackers::Vector{MBD.EventTracker} = countApses(env, primary, :peri, q1, true)
        newPeriapses::Vector{Vector{Float64}} = newEventTrackers[1].states
        newPeriapsesTimes::Vector{Float64} = newEventTrackers[1].times
        if !isempty(newPeriapses)
            if primary == 2
                lastIdx = findlast(q -> LinearAlgebra.norm(q[1:2]-getPrimaryState(env.EMDynamicsModel, 2)[1:2]) < env.MoonHill_EM, newPeriapses)
                if !isnothing(lastIdx)
                    qEsc = newPeriapses[lastIdx]
                    t2 = newPeriapsesTimes[lastIdx]
                    newFlag = length(newPeriapses)-lastIdx
                else
                    qEsc = q1
                end
            else
                if isInterior(newPeriapses[end], getMassRatio(env.EMDynamicsModel), getPrimaryState(env.EMDynamicsModel, 1), getPrimaryState(env.EMDynamicsModel, 2)) && (LinearAlgebra.norm(newPeriapses[end][1:2]-getPrimaryState(env.EMDynamicsModel, 2)[1:2]) > env.MoonHill_EM)
                    qEsc = newPeriapses[end]
                    t2 = newPeriapsesTimes[end]
                    newFlag = 0
                else
                    qEsc = newPeriapses[end-1]
                    t2 = newPeriapsesTimes[end-1]
                    newFlag = 1
                end
            end
        else
            qEsc = q1
        end
        (Deltav2, JCNew) = try
            trajOptimizeForTransit(env, JC1, primary, qEsc, newFlag, newqs, newFlags, volFileName)
        catch err
            println(err)
            (0.0, JC1)
        end
        vMag = LinearAlgebra.norm(qEsc[4:6])
        vhat = qEsc[4:6]./vMag
        q2 = append!(qEsc[1:3], (vMag+Deltav2) .* vhat)
        (escENew, escvNew) = getEscapeEnergy(env, q2)
        gammaNew = getEscapeManeuverEfficiencyMetric(escv1, escvNew, [Deltav2] .* env.charValues.EM.lstar ./ env.charValues.EM.tstar)
    end

    println("Old energy: $escE")
    println("Old velocity: $escv")
    println("Delta-v 1: $Deltav1")
    println("Delta-v 2: $Deltav2")
    println("New JC: $JCNew")
    println("New energy: $escENew")
    println("New velocity: $escvNew")
    println("New gamma: $gammaNew")

    MATLAB.put_variable(mf, :escE, escE)
    MATLAB.put_variable(mf, :escv, escv)
    MATLAB.put_variable(mf, :q0, q0)
    MATLAB.put_variable(mf, :t1, t1)
    MATLAB.put_variable(mf, :Deltav1, Deltav1)
    MATLAB.put_variable(mf, :q1, q1)
    MATLAB.put_variable(mf, :t2, t2)
    MATLAB.put_variable(mf, :Deltav2, Deltav2)
    MATLAB.put_variable(mf, :q2, q2)
    MATLAB.put_variable(mf, :newJC, JCNew)
    MATLAB.put_variable(mf, :newEscE, escENew)
    MATLAB.put_variable(mf, :newEscv, escvNew)
    MATLAB.put_variable(mf, :peris, reduce(hcat, periapses))
    MATLAB.put_variable(mf, :apos, reduce(hcat, apoapses))
end

function clusterTrajectoriesCR3BP(env::EscEnv, flags::Vector{Int64}, qs::Matrix{Float64}, mf::MATLAB.MatFile)
    esc0Indices::Vector{Int64} = findall(flags .== 0)
    n_esc0::Int64 = length(esc0Indices)
    features::Matrix{Float64} = Matrix{Float64}(undef, (n_esc0,6))
    rE::Vector{Float64} = getPrimaryState(env.EMDynamicsModel, 1)[1:2]
    rM::Vector{Float64} = getPrimaryState(env.EMDynamicsModel, 2)[1:2]
    Threads.@threads for idx::Int64 in eachindex(esc0Indices)
        q::Vector{Float64} = qs[:,esc0Indices[idx]]
        r::Float64 = LinearAlgebra.norm(q[1:2]-rE)
        features[idx,1] = r
        features[idx,2] = (q[1]-rE[1])/r
        features[idx,3] = q[2]/r
        arc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.MoonEvent, q, [0, 4.0*pi], env.EMDynamicsModel, [rM[1]])
        qf::Vector{Float64} = getStateByIndex(arc, -1)
        features[idx,4] = qf[2]
        features[idx,5] = qf[4]
        features[idx,6] = qf[5]
    end
    means::Matrix{Float64} = Statistics.mean(features, dims = 1)
    stds::Matrix{Float64} = Statistics.std(features, dims = 1)
    features_stand::Matrix{Float64} = (features .- means) ./ stds
    obs::Matrix{Float64} = copy(transpose(features_stand))
    dist::Matrix{Float64} = Distances.pairwise(Distances.Euclidean(), obs, dims = 2)
    tree::Clustering.Hclust = Clustering.hclust(dist, linkage = :ward)
    kBest::Int64 = 0
    scoreBest::Float64 = -1.0
    clustersBest::Vector{Int64} = []
    println("Evaluating cluster configurations...")
    for kCurrent::Int64 in 5:12
        clustersCurrent::Vector{Int64} = Clustering.cutree(tree, k = kCurrent)
        silScores::Vector{Float64} = Clustering.silhouettes(clustersCurrent, dist)
        avgSilScore::Float64 = Statistics.mean(silScores)
        println("\tClusters (k) = $kCurrent: Mean Silhouette Score = $(round(avgSilScore, digits = 4))")
        if avgSilScore > scoreBest
            scoreBest = avgSilScore
            kBest = kCurrent
            clustersBest = clustersCurrent
        end
    end
    println("Optimal number of clusters: $kBest")

    MATLAB.put_variable(mf, :clusters, clustersBest)
end

function run_apseMapCR3BP(JC::Float64, n::Int64, primary::Int64; apse::Symbol = :peri, grade::Symbol = :pro, apses::Bool = true)
    mf = MATLAB.MatFile("Output/ApseMaps/CR3BP_$(string(primary))_$(string(apse))_$(string(grade))_$(string(n))_$(string(JC)).mat", "w")
        
    env::EscEnv = setupEnvironment()

    rGrid::Vector{StaticArrays.SVector{2, Float64}} = getGrid(env, n, primary)
    apseMapCR3BP(env, JC, primary, rGrid, mf; apse = apse, grade = grade, apses = apses)
    
    MATLAB.close(mf)
end

function run_apseMapsCR3BP(JCs::Vector{Float64}, n::Int64, primary::Int64; apse::Symbol = :peri, grade::Symbol = :pro, apses::Bool = false)
    mf = MATLAB.MatFile("Output/ApseMaps/CR3BPJCVolume_$(string(primary))_$(string(apse))_$(string(grade))_$(string(n))_$(string(minimum(JCs)))_$(string(maximum(JCs))).mat", "w")

    env::EscEnv = setupEnvironment()

    rGrid::Vector{StaticArrays.SVector{2, Float64}} = getGrid(env, n, primary)
    o::Int64 = length(JCs)
    for j::Int64 in eachindex(JCs)
        println("\nProducing map $j / $o: JC = $(JCs[j])...")
        apseMapCR3BP(env, JCs[j], primary, rGrid, mf; apse = apse, grade = grade, apses = apses)
    end

    MATLAB.close(mf)
end

function run_escapeAnalysisCR3BP(fileName::String, mapName::String, primary::Int64, apse::Symbol, volFileName::String, directIdx::Int64, indirectIdx::Int64, indirectFlag::Int64)
    mf_in = MATLAB.MatFile(fileName, "r")

    map::Dict{String, Any} = get_variable(mf_in, mapName)

    MATLAB.close(mf_in)

    JC::Float64 = map["JC"]
    qs::Matrix{Float64} = map["q"]
    flags::Vector{Int64} = map["flags"]
    periapses::Matrix{Float64} = map["periapses"]
    periapsesIndices::Vector{Int64} = map["periapsesIndices"]

    mf_out = MATLAB.MatFile("Output/EscapeAnalysisCR3BP.mat", "w")    

    env::EscEnv = setupEnvironment()

    escapeAnalysisCR3BP(env, JC, primary, apse, flags, qs, periapses, periapsesIndices, volFileName, directIdx, indirectIdx, indirectFlag, mf_out)

    MATLAB.close(mf_out)
end

function run_assistedEscapeAnalysisCR3BP(fileName::String, mapName::String, primary::Int64, apse::Symbol, maneuver::Int64, volFileName::String)
    mf_in = MATLAB.MatFile(fileName, "r")

    map::Dict{String, Any} = get_variable(mf_in, mapName)

    MATLAB.close(mf_in)

    JC::Float64 = map["JC"]
    qs::Matrix{Float64} = map["q"]
    flags::Vector{Int64} = map["flags"]

    mf_out = MATLAB.MatFile("Output/AssistedEscapeAnalysisCR3BP.mat", "w")

    env::EscEnv = setupEnvironment()

    assistedEscapeAnalysisCR3BP(env, JC, primary, apse, maneuver, flags, qs, volFileName, mf_out)

    MATLAB.close(mf_out)
end

function run_assistApseStatesCR3BP(fileName::String, mapName::String, primary::Int64, family::String, idx::Int64, volFileName::String)
    mf_in = MATLAB.MatFile(fileName, "r")

    map::Dict{String, Any} = get_variable(mf_in, mapName)

    MATLAB.close(mf_in)

    JC::Float64 = map["JC"]
    qs::Matrix{Float64} = map["q"]
    flags::Vector{Int64} = map["flags"]

    mf_out = MATLAB.MatFile("Output/ApseStatesCR3BP.mat", "w")

    env::EscEnv = setupEnvironment()

    orbit::MBD.CR3BPPeriodicOrbit = interpOrbit(env.orbitTargeter, "FamilyData/CR3BPEM$(family)s.csv", "JC", JC)

    MATLAB.put_variable(mf_out, :qOrbit, orbit.initialCondition)
    MATLAB.put_variable(mf_out, :orbitP, orbit.period)

    (_, qPeris::Matrix{Float64}, _) = getPeriapsisStates(env, primary, orbit, mf_out)

    trajAssistedEscapeAnalysisCR3BP(env, JC, primary, qPeris[:,idx], qs, flags, volFileName, mf_out)
    
    MATLAB.close(mf_out)
end

function run_clusterTrajectoriesCR3BP(fileName::String, mapName::String)
    mf_in = MATLAB.MatFile(fileName, "r")

    map::Dict{String, Any} = get_variable(mf_in, mapName)

    MATLAB.close(mf_in)

    qs::Matrix{Float64} = map["q"]
    flags::Vector{Int64} = map["flags"]

    mf_out = MATLAB.MatFile("Output/ClusterTrajectoriesCR3BP.mat", "w")

    env::EscEnv = setupEnvironment()

    clusterTrajectoriesCR3BP(env, flags, qs, mf_out)

    MATLAB.close(mf_out)
end

end # EscCR3BP
