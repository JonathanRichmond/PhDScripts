"""
Script for computing CR3BP escape trajectories in the Earth-Moon system

Author: Jonathan LeFevre Richmond
C: 6/16/26
U: 8/5/26
"""

module EscCR3BP

using MBD, Clustering, DifferentialEquations, GLMakie, HDF5, LinearAlgebra, Logging, MATLAB, StaticArrays, Statistics

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../Utilities/Export.jl")
include("../Utilities/Plot.jl")

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
    periapsisEvent
    propagator::MBD.Propagator
    propagator_AL::MBD.Propagator
    propagator_STM::MBD.Propagator

    EarthHill_EM::Float64
    EarthRadius_EM::Float64
    MoonRadius_EM::Float64
end

function apoapsisCondition(state::Vector{Float64}, time::Float64, integrator)
    time-integrator.sol.prob.tspan[1] < 1E-6 ? 0.01 : LinearAlgebra.dot(state[1:3]-integrator.p[2], state[4:6])
end

function endAffect!(integrator, index)
    idx = findfirst(!=(0), index)
    if idx == 1
        if index[idx] == 1
            integrator.p[2][1].count += 1
            push!(integrator.p[2][1].states, integrator.u)
        else
            integrator.p[2][2].count += 1
            push!(integrator.p[2][2].states, integrator.u)
        end
    else
        integrator.p[2][idx+1].count += 1
        push!(integrator.p[2][idx+1].states, integrator.u)
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
    time-integrator.sol.prob.tspan[1] < 1E-6 ? -0.01 : LinearAlgebra.dot(state[1:3]-integrator.p[2], state[4:6])
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
    apoapsisEvent = DifferentialEquations.ContinuousCallback(apoapsisCondition, terminateAffect!)
    arclengthEvent = DifferentialEquations.ContinuousCallback(arclengthCondition, terminateAffect!)
    endEvents = DifferentialEquations.VectorContinuousCallback(endConditions, endAffect!, 4)
    escapeEvent = DifferentialEquations.ContinuousCallback(escapeCondition, terminateAffect!)
    flybyEvent = DifferentialEquations.ContinuousCallback(p1CR3BPDistanceCondition, terminateAffect!)
    MoonEvent = DifferentialEquations.ContinuousCallback(xValueCondition, terminateAffect!)
    periapsisEvent = DifferentialEquations.ContinuousCallback(periapsisCondition, terminateAffect!)

    EarthRadius_EM::Float64 = primaries[1].bodyRadius/charValues.EM.lstar
    MoonRadius_EM::Float64 = primaries[2].bodyRadius/charValues.EM.lstar
    EarthHill_EM::Float64 = charValues.SE.lstar*cbrt(getMassRatio(SESystemData)/3)/charValues.EM.lstar

    return EscEnv(EDynamicsModel, EMDynamicsModel, EMEoMs, SEDynamicsModel, SMDynamicsModel, primaries, Sun, charValues, apoapsisEvent, arclengthEvent, endEvents, escapeEvent, flybyEvent, MoonEvent, periapsisEvent, propagator, propagator_AL, propagator_STM, EarthHill_EM, EarthRadius_EM, MoonRadius_EM)
end

function getGrid(env::EscEnv, n::Int64, primary::Int64)
    if primary == 2
        radius::Float64 = 0.3
        center::Vector{Float64} = getPrimaryState(env.EMDynamicsModel, primary)[1:2]
    else
        radius = 1.25
        center = [0.0, 0.0]
    end

    xGrid::Vector{Float64} = collect(range(center[1]-radius, center[1]+radius, n))
    yGrid::Vector{Float64} = collect(range(center[2]-radius, center[2]+radius, n))
    rGrid::Matrix{StaticArrays.SVector{2, Float64}} = [StaticArrays.SA[x, y] for x in xGrid, y in yGrid]

    return rGrid
end

function isPeriapsis(dynamicsModel::MBD.CR3BPDynamicsModel, primary::Int64, r::StaticArrays.SVector{2, Float64}, dMag::Float64, v::Vector{Float64})
    mu::Float64 = (primary == 0 ? (dynamicsModel.systemData.primaryData[1].gravParam+dynamicsModel.systemData.primaryData[2].gravParam) : dynamicsModel.systemData.primaryData[primary].gravParam)
    vInert::StaticArrays.SVector{2, Float64} = StaticArrays.SA[v[1]-r[2], v[2]+r[1]]
    vInertMag::Float64 = LinearAlgebra.norm(vInert)
    vCirc::Float64 = sqrt(mu/(dMag*getCharLength(dynamicsModel)))*getCharTime(dynamicsModel)/getCharLength(dynamicsModel)

    return vInertMag > vCirc
end

function computeApseVelocities(env::EscEnv, JC::Float64, rGrid::Matrix{StaticArrays.SVector{2, Float64}})
    OmegaGrid::Matrix{Float64} = map(q -> getPseudopotential(env.EMDynamicsModel, Vector(push(q, 0.0))), rGrid)
    v2Grid::Matrix{Float64} = 2 .* OmegaGrid .- JC

    return map(v2 -> (v2 < 0 ? NaN : sqrt(v2)), v2Grid)
end

function computeApseStates(env::EscEnv, primary::Int64, JC::Float64, apse::Symbol, grade::Symbol, rGrid::Matrix{StaticArrays.SVector{2, Float64}})
    center::StaticArrays.SVector{2, Float64} = (primary == 0 ? StaticArrays.SA[0.0, 0.0] : StaticArrays.SVector{2, Float64}(getPrimaryState(env.EMDynamicsModel, primary)[1:2]))
    dGrid::Matrix{StaticArrays.SVector{2, Float64}} = rGrid .- Ref(center)
    dMagGrid::Matrix{Float64} = LinearAlgebra.norm.(dGrid)
    thatGrid::Matrix{StaticArrays.SVector{2, Float64}} = map(r -> StaticArrays.SA[-r[2], r[1]] ./ norm(r), dGrid)
    vMagGrid::Matrix{Float64} = computeApseVelocities(env, JC, rGrid)
    qGrid::Matrix{StaticArrays.MVector{6, Float64}} = map(r -> StaticArrays.MVector{6, Float64}([r[1], r[2], 0.0, NaN, NaN, NaN]), rGrid)
    gradeSign::Float64 = (grade == :retro) ? -1.0 : 1.0
    Threads.@threads for j::Int64 in eachindex(qGrid)
        isnan(vMagGrid[j]) && continue
        v::Vector{Float64} = gradeSign*vMagGrid[j] .* thatGrid[j]
        peri::Bool = isPeriapsis(env.EMDynamicsModel, primary, rGrid[j], dMagGrid[j], v)
        if ((apse == :peri) && !peri) || ((apse == :apo) && peri)
            qGrid[j][4] = Inf; qGrid[j][5] = Inf; qGrid[j][6] = Inf
            continue
        else
            qGrid[j][4] = v[1]; qGrid[j][5] = v[2]; qGrid[j][6] = 0.0
        end
    end

    return qGrid
end

function countApses(env::EscEnv, primary::Int64, apse::Symbol, IC::AbstractVector{Float64})
    peri = MBD.EventTracker(0, :peri, [])
    apo = MBD.EventTracker(0, :apo, [])
    escape = MBD.EventTracker(0, :escape, [])
    crashEarth = MBD.EventTracker(0, :earth, [])
    crashMoon = MBD.EventTracker(0, :moon, [])
    center::Vector{Float64} = (primary == 0 ? zeros(Float64, 3) : getPrimaryState(env.EMDynamicsModel, primary)[1:3])
    r_Earth::Vector{Float64} = getPrimaryState(env.EMDynamicsModel, 1)[1:3]
    r_Moon::Vector{Float64} = getPrimaryState(env.EMDynamicsModel, 2)[1:3]
    params::Vector{Any} = [apse, center, r_Earth, r_Moon, env.EarthHill_EM, env.EarthRadius_EM, env.MoonRadius_EM]
    (_, eventTrackers::Vector{EventTracker}) = propagateWithEvents(env.propagator, env.endEvents, Vector(IC), [0, 12.0*pi], env.EMDynamicsModel, [peri, apo, escape, crashEarth, crashMoon], params)
    
    return eventTrackers
end

function apseMapCR3BP(env::EscEnv, JC::Float64, n::Int64, primary::Int64, rGrid::Matrix{StaticArrays.SVector{2, Float64}}, mf::MATLAB.MatFile; apse::Symbol = :peri, grade::Symbol = :pro)
    qGrid::Matrix{StaticArrays.MVector{6, Float64}} = computeApseStates(env, primary, JC, apse, grade, rGrid)

    flags::Matrix{Int64} = fill(9, size(qGrid))
    flags[map(q -> isinf(q[4]), qGrid)] .= 8
    valid::Matrix{Bool} = map(q -> (!isnan(q[4]) && !isinf(q[4])), qGrid) 
    qProp::Vector{StaticArrays.SVector{6, Float64}} = qGrid[valid]
    qMap::Vector{CartesianIndex{2}} = findall(valid)
    counts::Matrix{Int64} = zeros(Int64, size(qGrid))
    periapses::Matrix{Vector{StaticArrays.SVector{6, Float64}}} = [StaticArrays.SVector{6, Float64}[] for _ in qGrid]
    apoapses::Matrix{Vector{StaticArrays.SVector{6, Float64}}} = [StaticArrays.SVector{6, Float64}[] for _ in qGrid]
    println("Propagating $(length(qProp)) / $(n^2) CR3BP trajectories with $(Threads.nthreads()) threads...")
    Threads.@threads for j in eachindex(qProp)
        eventTrackers::Vector{MBD.EventTracker} = countApses(env, primary, apse, qProp[j])
        apsesCount::Int64 = (apse == :peri ? eventTrackers[1].count : eventTrackers[2].count)
        if (eventTrackers[4].count != 0) || (eventTrackers[5].count != 0)
            flags[qMap[j]] = 7
        elseif eventTrackers[3].count != 0
            flags[qMap[j]] = min(apsesCount, 5)
        else
            flags[qMap[j]] = 6
        end
        counts[qMap[j]] = apsesCount
        periapses[qMap[j]] = map(q -> StaticArrays.SVector{6, Float64}(q), eventTrackers[1].states)
        apoapses[qMap[j]] = map(q -> StaticArrays.SVector{6, Float64}(q), eventTrackers[2].states)
    end

    flagCounts::Vector{Int64} = [count(==(f), flags) for f in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]]
    println("\tDirect escapes:\t$(flagCounts[1])")
    println("\tLonger escapes:\t$(sum(flagCounts[2:6]))")
    println("\tCaptures:\t$(flagCounts[7])")
    println("\tCrashes:\t$(flagCounts[8])")
    println("\tInvalid apses:\t$(flagCounts[9])")
    println("\tZVCs:\t\t$(flagCounts[10])")

    println("Exporting map...")
    validPeri::Vector{CartesianIndex{2}} = findall(!isempty, periapses)
    validApo::Vector{CartesianIndex{2}} = findall(!isempty, apoapses)
    periStates::Vector{StaticArrays.SVector{6, Float64}} = reduce(vcat, periapses[validPeri])
    periIndices::Vector{CartesianIndex{2}} = reduce(vcat, [fill(idx, length(periapses[idx])) for idx in validPeri])
    periLinear::Vector{Int64} = [LinearIndices(periapses)[idx] for idx in periIndices]
    apoStates::Vector{StaticArrays.SVector{6, Float64}} = reduce(vcat, apoapses[validApo])
    apoIndices::Vector{CartesianIndex{2}} = reduce(vcat, [fill(idx, length(apoapses[idx])) for idx in validApo])
    apoLinear::Vector{Int64} = [LinearIndices(apoapses)[idx] for idx in apoIndices]
    exportCR3BPApseMap(env.EMDynamicsModel, primary, apse, grade, JC, qGrid, flags, counts, periStates, periLinear, apoStates, apoLinear, mf, Symbol("map_", replace(string(JC), "." => "_")))
end

function getMaxHeliocentricEnergy(env::EscEnv, q::Vector{Float64})
    n_theta::Int64 = 1001
    theta_S::Vector{Float64} = collect(range(0.0, 2.0*pi, n_theta))
    r_B1::Float64 = env.charValues.SE.lstar
    v_B1::Float64 = sqrt(env.Sun.gravParam/r_B1)
    R_B1::Matrix{Float64} = r_B1.*hcat(cos.(theta_S), sin.(theta_S), zeros(Float64, n_theta))
    V_B1::Matrix{Float64} = v_B1.*hcat(-sin.(theta_S), cos.(theta_S), zeros(Float64, n_theta))
    Q_B1::Matrix{Float64} = hcat(R_B1, V_B1)
    Q_EI::Vector{Float64} = [q[1], q[2], q[3], (q[4]-q[2])/env.charValues.EM.tstar, (q[5]+q[1])/env.charValues.EM.tstar, q[6]/env.charValues.EM.tstar].*env.charValues.EM.lstar
    Q_SI::Matrix{Float64} = Q_EI' .+ Q_B1
    E::Vector{Float64} = Vector{Float64}(undef, n_theta)
    Threads.@threads for j::Int64 in 1:n_theta
        E[j] = 0.5*dot(Q_SI[j,4:6], Q_SI[j,4:6])-env.Sun.gravParam/LinearAlgebra.norm(Q_SI[j,1:3])
    end
    idx::Int64 = argmax(E)

    return (E[idx], Q_SI[idx,:])
end

function pruneVolumeData(JCRange::Vector{Float64}, indices::Vector{Int64}, volFileName::String)
    HDF5.h5open(volFileName, "r") do file
        volJC = vec(read(file["JCVolume"]))
        JCStartIdx::Int64 = findfirst(j -> j < JCRange[1], volJC)
        JCEndIdx::Int64 = findlast(j -> j > JCRange[2], volJC)
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

function getPartials_HohmannWRTEnergy(env::EscEnv, R2::Float64, E::Float64)
    muS::Float64 = env.Sun.gravParam
    a::Float64 = muS-2*E*R2
    A1::Float64 = ((2^(3/2)*R2^2*abs(E))/(a^2*sqrt(R2/a)))+((2^(3/2)*E*sqrt(R2/a))/abs(E))-(1/sqrt(-2*E))
    A2::Float64 = 2^(3/2)*sqrt(R2/a)*abs(E)-sqrt(-2*E)
    A::Float64 = A1*A2/abs(A2)
    B1::Float64 = muS*sqrt(2)*sqrt(1/(R2*a))-sqrt(muS/R2)
    B2::Float64 = a^2*sqrt(1/(R2*a))
    B::Float64 = muS*sqrt(2)*B1/(B2*abs(B1))
    
    return A+B
end

function getSmoothManeuverGradient(env::EscEnv, q0::Vector{Float64})
    R2::Float64 = getEquilibriumPoint(env.SMDynamicsModel, 1)[1]*env.charValues.SM.lstar
    escArc::MBD.CR3BPArc = propagateWithEvent(env.propagator_STM, env.escapeEvent, appendExtraInitialConditions(env.EMDynamicsModel, q0, MBD.STM), [0, 12.0*pi], env.EMDynamicsModel, [env.EarthHill_EM])
    qf::Vector{Float64} = getStateByIndex(escArc, -1)
    (Eesc::Float64, Qf_SI::Vector{Float64}) = getMaxHeliocentricEnergy(env, qf)
    dHdE::Float64 = getPartials_HohmannWRTEnergy(env, R2, Eesc)
    dEdqfI::Vector{Float64} = vcat(-env.Sun.gravParam .* Qf_SI[1:3] ./ LinearAlgebra.norm(Qf_SI[1:3])^3, Qf_SI[4:6]./env.charValues.EM.tstar).*env.charValues.EM.lstar
    Omega::Matrix{Float64} = [0 -1 0; 1 0 0; 0 0 0]
    dqfIdqf::Matrix{Float64} = [LinearAlgebra.I zeros(Float64, (3,3)); Omega LinearAlgebra.I]
    Phi::Matrix{Float64} = getStateTransitionMatrix(env.EMDynamicsModel, qf)
    qfdot::Vector{Float64} = zeros(Float64, 6)
    computeDerivatives!(qfdot, qf[1:6], (env.EMEoMs,), 0.0)
    gradg::Vector{Float64} = [qf[1], qf[2], qf[3], 0, 0, 0]./LinearAlgebra.norm(qf[1:3])
    dtfdalpha::Float64 = -gradg'*Phi[:,4:6]*(q0[4:6]./LinearAlgebra.norm(q0[4:6]))./(gradg'*qfdot)

    return dHdE*dEdqfI'*dqfIdqf*(Phi[:,4:6]*(q0[4:6]./LinearAlgebra.norm(q0[4:6]))+qfdot*dtfdalpha)
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

    return getEnergy(env.EDynamicsModel, oef)
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

function optimizeForFlyby(env::EscEnv, JC::Float64, q0::Vector{Float64}, idx::Int64, volJCs::Vector{Float64}, volFlags::Matrix{Int64}, volqs::Array{Float64, 3})
    v::Float64 = sqrt(q0[4]^2+q0[5]^2)
    grad::Float64 = getEnergyGradient(env, q0)
    qPrev::Vector{Float64} = copy(q0)
    qNew::Vector{Float64} = copy(q0)
    JCIdx::Int64 = 0
    if grad >= 0
        JCIdx = findfirst(x -> x < JC, volJCs)
        while (grad > 0) && (JCIdx <= length(volJCs)) && (volFlags[JCIdx,idx] == 0)
            qPrev = copy(qNew)
            qNew = volqs[:,idx,JCIdx]
            grad = getEnergyGradient(env, qNew)
            JCIdx += 1
        end
        (JCIdx > length(volJCs)) && throw(ErrorException("Gradient continues to increase: $idx"))
    else
        JCIdx = findlast(x -> x > JC, volJCs)
        while (grad < 0) && (JCIdx >= 1) && (volFlags[JCIdx,idx] == 0)
            qPrev = copy(qNew)
            qNew = volqs[:,idx,JCIdx]
            grad = getEnergyGradient(env, qNew)
            JCIdx -= 1
        end
        (JCIdx < 1) && throw(ErrorException("Gradient continues to decrease: $idx"))
    end
    EPrev::Float64 = getEscapeEnergy(env, qPrev)
    ENew::Float64 = getEscapeEnergy(env, qNew)
    qOpt::Vector{Float64} = (EPrev > ENew) ? copy(qPrev) : copy(qNew)
    Deltav::Float64 = sqrt(qOpt[4]^2+qOpt[5]^2)-v

    return Deltav
end

function escapeAnalysisCR3BP(env::EscEnv, JC::Float64, primary::Int64, flags::Matrix{Int64}, qs::Matrix{Float64}, radius::Float64, type::Symbol, volFileName::String, idx::Int64, mf::MATLAB.MatFile)
    @views rs = sqrt.(qs[1,:].^2 .+ qs[2,:].^2)
    mask::Vector{Bool} = rs .<= radius
    combinedMask::Vector{Bool} = (vec(flags) .== 0) .& mask
    esc0Indices::Vector{Int64} =  findall(combinedMask)
    n_esc0::Int64 = length(esc0Indices)
    esc0q_0s::Matrix{Float64} = qs[:,esc0Indices]
    esc0t_fs::Vector{Float64} = Vector{Float64}(undef, n_esc0)
    esc0Es::Vector{Float64} = Vector{Float64}(undef, n_esc0)
    Threads.@threads for idx::Int64 in eachindex(esc0Indices)
        arc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.escapeEvent, qs[:,esc0Indices[idx]], [0, 12.0*pi], env.EMDynamicsModel, [env.EarthHill_EM])
        esc0t_fs[idx] = getTimeByIndex(arc, -1)
        esc0Es[idx] = getEscapeEnergy(env, qs[:,esc0Indices[idx]])
    end

    MATLAB.put_variable(mf, :esc0q0, esc0q_0s)
    MATLAB.put_variable(mf, :esc0tf, esc0t_fs)
    MATLAB.put_variable(mf, :esc0E, esc0Es)

    JCRange::Vector{Float64} = (type == :flyby) ? [3.16, 2.5] : [2.5, 0.89]
    (volJCs::Vector{Float64}, volFlags::Matrix{Int64}, volqs::Array{Float64, 3}) = pruneVolumeData(JCRange, [idx], volFileName)
    @views vs = sqrt.(qs[4,:].^2 .+ qs[5,:].^2)
    Deltav2s::Vector{Float64} = fill(NaN, (length(vec(volFlags))))
    escEs::Vector{Float64} = copy(Deltav2s)
    grads::Vector{Float64} = copy(Deltav2s)
    escIndices::Vector{Int64} = findall(f -> f == 0, vec(volFlags))
    Threads.@threads for escIdx::Int64 in escIndices
        q::Vector{Float64} = volqs[:,1,escIdx]
        Deltav2s[escIdx] = sqrt(q[4]^2+q[5]^2)-vs[idx]
        arc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.escapeEvent, q, [0, 12.0*pi], env.EMDynamicsModel, [env.EarthHill_EM])
        qf::Vector{Float64} = getStateByIndex(arc, -1)
        escEs[escIdx] = getEscapeEnergy(env, qf)
        grads[escIdx] = getEnergyGradient(env, q)
    end

    MATLAB.put_variable(mf, :Deltav2s, Deltav2s)
    MATLAB.put_variable(mf, :EscapeEs, escEs)
    MATLAB.put_variable(mf, :DeltaEs, grads)

    """Maneuver Sequencing"""
    qTraj::Vector{Float64} = qs[:,idx]
    # Deltav1::Float64 = 0
    Deltav2::Float64 = 0
    if flags[idx] == 0
        Deltav2 = optimizeForFlyby(env, JC, qTraj, 1, volJCs, volFlags, volqs)
    elseif (0 < flags[idx] < 6)
    elseif (flags[idx] == 6) || (flags[idx] == 7)
    end

    MATLAB.put_variable(mf, :Deltav2, Deltav2)

    # smoothGrad::Float64 = getSmoothManeuverGradient(env, qs[:,idx])
    # Deltav1::Float64 = NaN
    # dMin::Float64 = Inf
    # if abs(smoothGrad) <= 1.0
    #     Deltav1 = 0
    #     flybyArc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.flybyEvent, qs[:,idx], [0, 3.0*pi], env.EMDynamicsModel, [1.1])
    #     for s::Vector{Float64} in flybyArc.states
    #         d::Float64 = getExcursion(env.EMDynamicsModel, 2, s)
    #         (d < dMin) && (dMin = d)
    #     end
    # else
    #     dir::Float64 = -sign(smoothGrad)
    #     grad::Float64 = dir+smoothGrad
    #     if grad > 0
    #         JCIdx::Int64 = findfirst(j -> j <= JC, volJC)
    #         if JCIdx === nothing
    #             Deltav1 = 0
    #             flybyArc = propagateWithEvent(env.propagator, env.flybyEvent, qs[:,idx], [0, 3.0*pi], env.EMDynamicsModel, [1.1])
    #             for s::Vector{Float64} in flybyArc.states
    #                 d::Float64 = getExcursion(env.EMDynamicsModel, 2, s)
    #                 (d < dMin) && (dMin = d)
    #             end
    #         end
    #         totalCosts::Vector{Float64} = []
    #         while (grad >= 0) && (JCIdx > 1)
    #             JCIdx -= 1
    #             qGrad::Vector{Float64} = volqs[:,1,JCIdx]
    #             grad = -1+getSmoothManeuverGradient(env, volqs[:,1,JCIdx])
    #             totalCost::Float64 = getHohmannCost(env, qGrad)+abs(sqrt(qGrad[4]^2+qGrad[5]^2)-vs[idx])*env.charValues.EM.lstar/env.charValues.EM.tstar
    #             push!(totalCosts, totalCost)
    #         end
    #         adjust::Int64 = argmin(totalCosts)
    #         optIdx::Int64 = findfirst(j -> j <= JC, volJC)-adjust
    #         qOpt::Vector{Float64} = volqs[:,1,optIdx]
    #         Deltav1 = sqrt(qOpt[4]^2+qOpt[5]^2)-vs[idx]
    #         flybyArc = propagateWithEvent(env.propagator, env.flybyEvent, qOpt, [0, 3.0*pi], env.EMDynamicsModel, [1.1])
    #         for s::Vector{Float64} in flybyArc.states
    #             d::Float64 = getExcursion(env.EMDynamicsModel, 2, s)
    #             (d < dMin) && (dMin = d)
    #         end
    #     elseif grad < 0
    #         JCIdx = findlast(j -> j >= JC, volJC)
    #         if JCIdx === nothing
    #             Deltav1 = 0
    #             flybyArc = propagateWithEvent(env.propagator, env.flybyEvent, qs[:,idx], [0, 3.0*pi], env.EMDynamicsModel, [1.1])
    #             for s::Vector{Float64} in flybyArc.states
    #                 d::Float64 = getExcursion(env.EMDynamicsModel, 2, s)
    #                 (d < dMin) && (dMin = d)
    #             end
    #         end
    #         totalCosts = []
    #         while (grad <= 0) && (JCIdx < length(volJC))
    #             JCIdx += 1
    #             qGrad::Vector{Float64} = volqs[:,1,JCIdx]
    #             grad = 1+getSmoothManeuverGradient(env, volqs[:,1,JCIdx])
    #             totalCost::Float64 = getHohmannCost(env, qGrad)+abs(sqrt(qGrad[4]^2+qGrad[5]^2)-vs[idx])*env.charValues.EM.lstar/env.charValues.EM.tstar
    #             push!(totalCosts, totalCost)
    #         end
    #         adjust = argmin(totalCosts)
    #         optIdx = findlast(j -> j >= JC, volJC)+adjust
    #         qOpt = volqs[:,1,optIdx]
    #         Deltav1 = sqrt(qOpt[4]^2+qOpt[5]^2)-vs[idx]
    #         flybyArc = propagateWithEvent(env.propagator, env.flybyEvent, qOpt, [0, 3.0*pi], env.EMDynamicsModel, [1.1])
    #         for s::Vector{Float64} in flybyArc.states
    #             d::Float64 = getExcursion(env.EMDynamicsModel, 2, s)
    #             (d < dMin) && (dMin = d)
    #         end
    #     end
    # end

    # MATLAB.put_variable(mf, :Deltav1, Deltav1)
    # MATLAB.put_variable(mf, :flybyDistance, dMin-env.MoonRadius_EM*env.charValues.EM.lstar)

    # esc1Indices::Vector{Int64} = findall(f -> f == 1, vec(flags))
    # n_esc1::Int64 = length(esc1Indices)
    # esc1q_0s::Matrix{Float64} = qs[:,esc1Indices]
    # esc1t_fs::Vector{Float64} = Vector{Float64}(undef, n_esc1)
    # esc1Es::Vector{Float64} = Vector{Float64}(undef, n_esc1)
    # Threads.@threads for idx::Int64 in eachindex(esc1Indices)
    #     arc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.escapeEvent, qs[:,esc1Indices[idx]], [0, 12.0*pi], env.EMDynamicsModel, [env.EarthHill_EM])
    #     q_f::Vector{Float64} = getStateByIndex(arc, -1)
    #     esc1t_fs[idx] = getTimeByIndex(arc, -1)
    #     (esc1Es[idx], _) = getMaxHeliocentricEnergy(env, q_f)
    # end

    # MATLAB.put_variable(mf, :esc1q0, esc1q_0s)
    # MATLAB.put_variable(mf, :esc1tf, esc1t_fs)
    # MATLAB.put_variable(mf, :esc1E, esc1Es)

    # (volJC::Vector{Float64}, _, _) = pruneVolumeData(3.16, [idx], volFileName)
    # @views vs = sqrt.(qs[4,:].^2 .+ qs[5,:].^2)
    # q::Vector{Float64} = qs[:,idx]
    # center::Vector{Float64} = (primary == 0 ? zeros(Float64, 3) : getPrimaryState(env.EMDynamicsModel, primary)[1:3])
    # apoapsisArc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.apoapsisEvent, q, [0, 4.0*pi], env.EMDynamicsModel, [center])
    # q_apo::Vector{Float64} = getStateByIndex(apoapsisArc, -1)
    # println(q_apo[1:2])
    # vMag::Float64 = LinearAlgebra.norm(q_apo[4:6])
    # vhat::Vector{Float64} = q_apo[4:6]./vMag
    # Deltav1s::Vector{Float64} = fill(NaN, (2*length(volJC)))
    # flagsNew::Vector{Int64} = fill(9, (2*length(volJC)))
    # Threads.@threads for JCIdx::Int64 in eachindex(volJC)
    #     vMagNew::Float64 = computeApseVelocities(env, volJC[JCIdx], fill(StaticArrays.SVector{2, Float64}(q_apo[1:2]), 1, 1))[1,1]
    #     isnan(vMagNew) && continue
    #     vNew::Vector{Float64} = vhat.*vMagNew
    #     peri::Bool = isPeriapsis(env.EMDynamicsModel, primary, StaticArrays.SVector{2, Float64}(q_apo[1:2]), LinearAlgebra.norm(q_apo[1:2]), vNew[1:2])
    #     q_apoNew::Vector{Float64} = [q_apo[1], q_apo[2], 0.0, NaN, NaN, NaN]
    #     if peri
    #         flagsNew[JCIdx] = 8
    #         continue
    #     else
    #         q_apoNew[4:6] = vNew
    #     end
    #     Deltav1s[JCIdx] = vMagNew-vMag
    #     eventTrackers::Vector{MBD.EventTracker} = countApses(env, primary, :peri, q_apoNew)
    #     apsesCount::Int64 = eventTrackers[1].count-1
    #     if (eventTrackers[4].count != 0) || (eventTrackers[5].count != 0)
    #         flagsNew[JCIdx] = 7
    #     elseif eventTrackers[3].count != 0
    #         flagsNew[JCIdx] = min(apsesCount, 5)
    #     else
    #         flagsNew[JCIdx] = 6
    #     end

    #     vNewFlip::Vector{Float64} = vhat.*(-vMagNew)
    #     periFlip::Bool = isPeriapsis(env.EMDynamicsModel, primary, StaticArrays.SVector{2, Float64}(q_apo[1:2]), LinearAlgebra.norm(q_apo[1:2]), vNewFlip[1:2])
    #     q_apoNewFlip::Vector{Float64} = [q_apo[1], q_apo[2], 0.0, NaN, NaN, NaN]
    #     if periFlip
    #         flagsNew[JCIdx+length(volJC)] = 8
    #         continue
    #     else
    #         q_apoNewFlip[4:6] = vNewFlip
    #     end
    #     Deltav1s[JCIdx+length(volJC)] = -(vMagNew+vMag)
    #     eventTrackersFlip::Vector{MBD.EventTracker} = countApses(env, primary, :peri, q_apoNewFlip)
    #     apsesCountFlip::Int64 = eventTrackersFlip[1].count-1
    #     if (eventTrackersFlip[4].count != 0) || (eventTrackersFlip[5].count != 0)
    #         flagsNew[JCIdx+length(volJC)] = 7
    #     elseif eventTrackers[3].count != 0
    #         flagsNew[JCIdx+length(volJC)] = min(apsesCountFlip, 5)
    #     else
    #         flagsNew[JCIdx+length(volJC)] = 6
    #     end
    # end

    # MATLAB.put_variable(mf, :ApoState, q_apo)
    # MATLAB.put_variable(mf, :Deltav1s, Deltav1s)
    # MATLAB.put_variable(mf, :NewFlags, flagsNew)

    # esc2Indices::Vector{Int64} = findall(f -> f == 2, vec(flags))
    # n_esc2::Int64 = length(esc2Indices)
    # esc2q_0s::Matrix{Float64} = qs[:,esc2Indices]
    # esc2t_fs::Vector{Float64} = Vector{Float64}(undef, n_esc2)
    # esc2Es::Vector{Float64} = Vector{Float64}(undef, n_esc2)
    # Threads.@threads for idx::Int64 in eachindex(esc2Indices)
    #     arc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.escapeEvent, qs[:,esc2Indices[idx]], [0, 12.0*pi], env.EMDynamicsModel, [env.EarthHill_EM])
    #     q_f::Vector{Float64} = getStateByIndex(arc, -1)
    #     esc2t_fs[idx] = getTimeByIndex(arc, -1)
    #     (esc2Es[idx], _) = getMaxHeliocentricEnergy(env, q_f)
    # end

    # MATLAB.put_variable(mf, :esc2q0, esc2q_0s)
    # MATLAB.put_variable(mf, :esc2tf, esc2t_fs)
    # MATLAB.put_variable(mf, :esc2E, esc2Es)
end

function assisted0EscapeAnalysisCR3BP(flags::Matrix{Int64}, qs::Matrix{Float64}, volFileName::String, mf::MATLAB.MatFile)
    mask::Vector{Bool} = (vec(flags) .== 6) .| (vec(flags) .== 7)
    nonEscIndices::Vector{Int64} = findall(mask)
    (_, volFlags::Matrix{Int64}, volqs::Array{Float64, 3}) = pruneVolumeData(3.16, nonEscIndices, volFileName)
    @views vs = sqrt.(qs[4,:].^2 .+ qs[5,:].^2)
    Deltavs::Vector{Float64} = zeros(Float64, length(vs))
    Threads.@threads for j::Int64 in eachindex(nonEscIndices)
        idx::Int64 = nonEscIndices[j]
        escIndices::Vector{Int64} = findall(f -> f < 6, @view volFlags[:,j])
        if isempty(escIndices)
            Deltavs[idx] = NaN
        else
            minDeltav::Float64 = Inf
            for escIdx::Int64 in escIndices
                Deltav::Float64 = sqrt(volqs[4,j,escIdx]^2+volqs[5,j,escIdx]^2)-vs[idx]
                (abs(Deltav) < abs(minDeltav)) && (minDeltav = Deltav)
            end
            Deltavs[idx] = minDeltav
        end
    end
    Threads.@threads for j::Int64 in eachindex(flags)
        if (flags[j] == 8) || (flags[j] == 9)
            Deltavs[j] = NaN
        end
    end

    MATLAB.put_variable(mf, :Deltav0s, Deltavs)
end

function assisted1EscapeAnalysisCR3BP(env::EscEnv, JC::Float64, flags::Matrix{Int64}, qs::Matrix{Float64}, radius::Float64, volFileName::String, mf::MATLAB.MatFile)
    @views rs = sqrt.(qs[1,:].^2 .+ qs[2,:].^2)
    mask::Vector{Bool} = rs .<= radius
    combinedMask::Vector{Bool} = (vec(flags) .== 0) .& mask
    directIndices::Vector{Int64} = findall(combinedMask)
    (volJC::Vector{Float64}, volFlags::Matrix{Int64}, volqs::Array{Float64, 3}) = pruneVolumeData(3.16, directIndices, volFileName)
    @views vs = sqrt.(qs[4,:].^2 .+ qs[5,:].^2)
    Deltavs::Vector{Float64} = zeros(Float64, length(vs))
    maneuverMetrics::Vector{Float64} = zeros(Float64, length(vs))
    Threads.@threads for j::Int64 in eachindex(directIndices)
        idx::Int64 = directIndices[j]
        JCIndices::Vector{Int64} = findall(f -> f == 0, @view volFlags[:,j])
        isempty(JCIndices) && continue
        smoothGrad::Float64 = getSmoothManeuverGradient(env, qs[:,idx])
        (abs(smoothGrad) <= 1.0) && continue
        totalCosts::Vector{Float64} = []
        cost0::Float64 = getHohmannCost(env, qs[:,idx])
        push!(totalCosts, cost0)
        dir::Float64 = -sign(smoothGrad)
        grad::Float64 = dir+smoothGrad
        if grad > 0
            JCIdx = findlast(j -> j > JC, volJC[JCIndices])
            (JCIdx === nothing) && continue
            currentJC::Int64 = JCIdx
            while (grad >= 0) && (currentJC >= 1)
                qGrad::Vector{Float64} = volqs[:,j,JCIndices[currentJC]]
                grad = -1+getSmoothManeuverGradient(env, qGrad)
                v_q::Float64 = sqrt(qGrad[4]^2+qGrad[5]^2)
                totalCost::Float64 = getHohmannCost(env, qGrad)+abs(v_q-vs[idx])*env.charValues.EM.lstar/env.charValues.EM.tstar
                push!(totalCosts, totalCost)
                currentJC -= 1
            end
            adjust::Int64 = argmin(totalCosts)-1
            (adjust == 0) && continue
            optIdx::Int64 = JCIdx-adjust+1
            qOpt::Vector{Float64} = volqs[:,j,JCIndices[optIdx]]
            Deltavs[idx] = sqrt(qOpt[4]^2+qOpt[5]^2)-vs[idx]
            costOpt::Float64 = getHohmannCost(env, qOpt)+abs(Deltavs[idx])*env.charValues.EM.lstar/env.charValues.EM.tstar
            maneuverMetrics[idx] = cost0-costOpt
        elseif grad < 0
            JCIdx = findfirst(j -> j < JC, volJC[JCIndices])
            (JCIdx === nothing) && continue
            currentJC = JCIdx
            while (grad <= 0) && (currentJC <= length(JCIndices))
                qGrad::Vector{Float64} = volqs[:,j,JCIndices[currentJC]]
                grad = 1+getSmoothManeuverGradient(env, qGrad)
                v_q::Float64 = sqrt(qGrad[4]^2+qGrad[5]^2)
                totalCost::Float64 = getHohmannCost(env, qGrad)+abs(v_q-vs[idx])*env.charValues.EM.lstar/env.charValues.EM.tstar
                push!(totalCosts, totalCost)
                currentJC += 1
            end
            adjust = argmin(totalCosts)-1
            (adjust == 0) && continue
            optIdx = JCIdx+adjust-1
            qOpt = volqs[:,j,JCIndices[optIdx]]
            Deltavs[idx] = sqrt(qOpt[4]^2+qOpt[5]^2)-vs[idx]
            costOpt = getHohmannCost(env, qOpt)+abs(Deltavs[idx])*env.charValues.EM.lstar/env.charValues.EM.tstar
            maneuverMetrics[idx] = cost0-costOpt
        end
    end
    Threads.@threads for j::Int64 in eachindex(flags)
        if (flags[j] == 8) || (flags[j] == 9)
            Deltavs[j] = NaN
            maneuverMetrics[j] = NaN
        elseif (flags[j] == 0) && (rs[j] > radius)
            Deltavs[j] = Inf
            maneuverMetrics[j] = Inf
        end
    end

    MATLAB.put_variable(mf, :Deltav1s, Deltavs)
    MATLAB.put_variable(mf, :maneuverMetrics, maneuverMetrics)
end

function assistedEscapeAnalysisCR3BP(env::EscEnv, JC::Float64, flags::Matrix{Int64}, qs::Matrix{Float64}, maneuvers::Int64, radius::Float64, volFileName::String, mf::MATLAB.MatFile)
    # if maneuvers == 0
    #     assisted0EscapeAnalysisCR3BP(flags, qs, volFileName, mf)
    # elseif maneuvers == 1
    #     assisted1EscapeAnalysisCR3BP(env, JC, flags, qs, radius, volFileName, mf)
    # end
    @views rs = sqrt.(qs[1,:].^2 .+ qs[2,:].^2)
    mask::Vector{Bool} = rs .<= radius
    indices::Vector{Int64} = findall(mask)
    JCRange::Vector{Float64} = [3.16, 2.5]
    (volJCs::Vector{Float64}, volFlags::Matrix{Int64}, volqs::Array{Float64, 3}) = pruneVolumeData(JCRange, indices, volFileName)
    Deltav1s::Vector{Float64} = zeros(Float64, size(qs, 2))
    Deltav2s::Vector{Float64} = copy(Deltav1s)
    Threads.@threads for j::Int64 in eachindex(indices)
        idx::Int64 = indices[j]
        qTraj::Vector{Float64} = qs[:,idx]
        if flags[idx] == 0
            Deltav2s[idx] = optimizeForFlyby(env, JC, qTraj, j, volJCs, volFlags, volqs)
        elseif (0 < flags[idx] < 6)
            # throw(ErrorException("Not defined yet"))
        elseif (flags[idx] == 6) || (flags[idx] == 7)
            # throw(ErrorException("Not defined yet"))
        else
            Deltav1s[idx] = NaN
            Deltav2s[idx] = NaN
        end
    end

    MATLAB.put_variable(mf, :Deltav1s, Deltav1s)
    MATLAB.put_variable(mf, :Deltav2s, Deltav2s)
end

function clusterTrajectoriesCR3BP(env::EscEnv, flags::Matrix{Int64}, qs::Matrix{Float64}, radius::Float64, mf::MATLAB.MatFile)
    @views rs = sqrt.(qs[1,:].^2 .+ qs[2,:].^2)
    mask::Vector{Bool} = rs .<= radius
    combinedMask::Vector{Bool} = (vec(flags) .== 0) .& mask
    esc0Indices::Vector{Int64} =  findall(combinedMask)
    n_esc0::Int64 = length(esc0Indices)
    steps::Int64 = 51
    trajectories::Vector{Vector{Vector{Float64}}} = Vector{Vector{Vector{Float64}}}(undef, n_esc0)
    yMoons::Vector{Float64} = Vector{Float64}(undef, n_esc0)
    Threads.@threads for idx::Int64 in eachindex(esc0Indices)
        arc::MBD.CR3BPArc = propagateWithEvent(env.propagator_AL, env.flybyEvent, appendExtraInitialConditions(env.EMDynamicsModel, qs[:,esc0Indices[idx]], MBD.ARCLENGTH), [0, 12.0*pi], env.EMDynamicsModel, [1.2])
        arclength::Float64 = getStateByIndex(arc, -1)[end]   
        arclengths::Vector{Float64} = collect(LinRange(0.0, arclength, steps))
        path::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, steps)
        path[1] = qs[1:2,esc0Indices[idx]]
        for a::Int64 in eachindex(arclengths[2:end])
            shortArc::MBD.CR3BPArc = propagateWithEvent(env.propagator_AL, env.arclengthEvent, appendExtraInitialConditions(env.EMDynamicsModel, qs[:,esc0Indices[idx]], MBD.ARCLENGTH), [0, 12.0*pi], env.EMDynamicsModel, [env.EMDynamicsModel, arclengths[a+1]])
            qf::Vector{Float64} = getStateByIndex(shortArc, -1)
            path[a+1] = [qf[1], qf[2]]
        end
        trajectories[idx] = path
        MoonArc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.MoonEvent, qs[:,esc0Indices[idx]], [0, 12.0*pi], env.EMDynamicsModel, [1-getMassRatio(env.EMDynamicsModel)])
        yMoons[idx] = getStateByIndex(MoonArc, -1)[2]
    end
    dist::Matrix{Float64} = zeros(Float64, (n_esc0,n_esc0))
    Threads.@threads for i::Int64 in 1:n_esc0
        for j::Int64 in (i+1):n_esc0
            d::Float64 = discreteFrechet(trajectories[i], trajectories[j])
            dist[i,j] = d
            dist[j,i] = d
        end
    end
    weight::Float64 = 5.0
    augDist::Matrix{Float64} = copy(dist)
    Threads.@threads for i::Int64 in 1:n_esc0
        for j::Int64 in (i+1):n_esc0
            yPenalty::Float64 = abs(yMoons[i]-yMoons[j])
            augDist[i,j] += weight*yPenalty
            augDist[j,i] = augDist[i,j]
        end
    end
    tree::Clustering.Hclust = Clustering.hclust(dist, linkage = :average)
    kBest::Int64 = 3
    scoreBest::Float64 = -1.0
    clustersBest::Vector{Int64} = []
    println("Evaluating cluster configurations...")
    for kCurrent::Int64 in 3:10
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
    cluster1LocalIndices::Vector{Int64} = findall(x -> x == 1, clustersBest)
    subDist::Matrix{Float64} = dist[cluster1LocalIndices,cluster1LocalIndices]
    subTree::Clustering.Hclust = Clustering.hclust(subDist, linkage = :average)
    subkBest::Int64 = 2
    subScoreBest::Float64 = -1.0
    subClustersBest::Vector{Int64} = []
    println("Evaluating sub-cluster configurations...")
    for kCurrent::Int64 in 3:10
        clustersCurrent::Vector{Int64} = Clustering.cutree(subTree, k = kCurrent)
        silScores::Vector{Float64} = Clustering.silhouettes(clustersCurrent, subDist)
        avgSilScore::Float64 = Statistics.mean(silScores)
        println("\tSub-clusters (k) = $kCurrent: Mean Silhouette Score = $(round(avgSilScore, digits = 4))")
        if avgSilScore > subScoreBest
            subScoreBest = avgSilScore
            subkBest = kCurrent
            subClustersBest = clustersCurrent
        end
    end
    println("Optimal number of sub-clusters: $subkBest")
    clusters::Vector{Int64} = zeros(Int64, n_esc0)
    for (i::Int64, idx::Int64) in enumerate(cluster1LocalIndices)
        clusters[idx] = subClustersBest[i]
    end
    for j::Int64 in 1:n_esc0
        originalId::Int64 = clustersBest[j]
        (originalId != 1) && (clusters[j] = originalId+subkBest-1)
    end
    println("Total number of clusters after refinement: $(maximum(clusters))")

    MATLAB.put_variable(mf, :clusters, clusters)
end

function run_apseMapCR3BP(JC::Float64, n::Int64, primary::Int64; apse::Symbol = :peri, grade::Symbol = :pro)
    mf = MATLAB.MatFile("Output/ApseMaps/CR3BP_$(string(primary))_$(string(apse))_$(string(grade))_$(string(n))_$(string(JC)).mat", "w")
        
    env::EscEnv = setupEnvironment()

    rGrid::Matrix{StaticArrays.SVector{2, Float64}} = getGrid(env, n, primary)
    apseMapCR3BP(env, JC, n, primary, rGrid, mf; apse = apse, grade = grade)
    
    MATLAB.close(mf)
end

function run_apseMapsCR3BP(JCs::Vector{Float64}, n::Int64, primary::Int64; apse::Symbol = :peri, grade::Symbol = :pro)
    mf = MATLAB.MatFile("Output/ApseMaps/CR3BPJCVolume_$(string(primary))_$(string(apse))_$(string(grade))_$(string(n))_$(string(minimum(JCs)))_$(string(maximum(JCs))).mat", "w")

    env::EscEnv = setupEnvironment()

    rGrid::Matrix{StaticArrays.SVector{2, Float64}} = getGrid(env, n, primary)
    o::Int64 = length(JCs)
    for j::Int64 in eachindex(JCs)
        println("\nProducing map $j / $o: JC = $(JCs[j])...")
        apseMapCR3BP(env, JCs[j], n, primary, rGrid, mf; apse = apse, grade = grade)
    end

    MATLAB.close(mf)
end

function run_escapeAnalysisCR3BP(fileName::String, mapName::String, primary::Int64, radius::Float64, type::Symbol, volFileName::String, idx::Int64)
    mf_in = MATLAB.MatFile(fileName, "r")

    map::Dict{String, Any} = get_variable(mf_in, mapName)

    MATLAB.close(mf_in)

    JC::Float64 = map["JC"]
    qs::Matrix{Float64} = map["q"]
    flags::Matrix{Int64} = map["flags"]

    mf_out = MATLAB.MatFile("Output/EscapeAnalysisCR3BP.mat", "w")

    env::EscEnv = setupEnvironment()

    escapeAnalysisCR3BP(env, JC, primary, flags, qs, radius, type, volFileName, idx, mf_out)

    MATLAB.close(mf_out)
end

function run_assistedEscapeAnalysisCR3BP(fileName::String, mapName::String, maneuvers::Int64, radius::Float64, volFileName::String)
    mf_in = MATLAB.MatFile(fileName, "r")

    map::Dict{String, Any} = get_variable(mf_in, mapName)

    MATLAB.close(mf_in)

    JC::Float64 = map["JC"]
    qs::Matrix{Float64} = map["q"]
    flags::Matrix{Int64} = map["flags"]

    mf_out = MATLAB.MatFile("Output/AssistedEscapeAnalysisCR3BP.mat", "w")

    env::EscEnv = setupEnvironment()

    assistedEscapeAnalysisCR3BP(env, JC, flags, qs, maneuvers, radius, volFileName, mf_out)

    MATLAB.close(mf_out)
end

function run_clusterTrajectoriesCR3BP(fileName::String, mapName::String, radius::Float64)
    mf_in = MATLAB.MatFile(fileName, "r")

    map::Dict{String, Any} = get_variable(mf_in, mapName)

    MATLAB.close(mf_in)

    qs::Matrix{Float64} = map["q"]
    flags::Matrix{Int64} = map["flags"]

    mf_out = MATLAB.MatFile("Output/ClusterTrajectoriesCR3BP.mat", "w")

    env::EscEnv = setupEnvironment()

    clusterTrajectoriesCR3BP(env, flags, qs, radius, mf_out)

    MATLAB.close(mf_out)
end

end # EscCR3BP
