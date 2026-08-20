"""
Script for computing BCR4BP escape trajectories in the Earth-Moon system

Author: Jonathan LeFevre Richmond
C: 8/12/26
U: 8/20/26
"""

module EscBCR4BP

using MBD, Clustering, DifferentialEquations, Distances, HDF5, LinearAlgebra, Logging, MATLAB, StaticArrays, Statistics

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../Utilities/Export.jl")

struct CharValues
    lstar::Float64
    mstar::Float64
    tstar::Float64
end

struct EscEnv
    EDynamicsModel::MBD.KDynamicsModel
    EMDynamicsModel::MBD.CR3BPDynamicsModel
    EMEoMs::MBD.BCR4BP12EquationsOfMotion
    EMSDynamicsModel::MBD.BCR4BP12DynamicsModel
    SEDynamicsModel::MBD.CR3BPDynamicsModel

    primaries::Vector{MBD.BodyData}

    charValues::NamedTuple

    endEvents
    propagator::MBD.Propagator

    EarthHill_EM::Float64
    EarthRadius_EM::Float64
    MoonHill_EM::Float64
    MoonRadius_EM::Float64
end

function endAffect!(integrator, index)
    idx = findfirst(!=(0), index)
    if idx == 1
        if (index[idx] == 1) && (isInterior(integrator.u, get12MassRatio(integrator.p[1]), integrator.p[5], integrator.p[6])) && ((integrator.p[10] != 2) ? (integrator.u[1] < 0.84) : true)
            integrator.p[2][1].count += 1
            push!(integrator.p[2][1].states, integrator.u)
        elseif index[idx] == -1
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

function setupEnvironment()::EscEnv
    ESystemData = MBD.KSystemData("Earth")
    EMSystemData = MBD.CR3BPSystemData("Earth", "Moon")
    SESystemData = MBD.CR3BPSystemData("Sun", "Earth")
    EMSSystemData = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
    EDynamicsModel = MBD.KDynamicsModel(ESystemData)
    EMDynamicsModel = MBD.CR3BPDynamicsModel(EMSystemData)
    SEDynamicsModel = MBD.CR3BPDynamicsModel(SESystemData)
    EMSDynamicsModel = MBD.BCR4BP12DynamicsModel(EMSSystemData)
    EMEoMs = MBD.BCR4BP12EquationsOfMotion(MBD.SIMPLE, EMSDynamicsModel)

    primaries::Vector{MBD.BodyData} = EMSSystemData.primaryData

    charValues::NamedTuple = (EM = CharValues(getCharLength(EMSystemData), getCharMass(EMSystemData), getCharTime(EMSystemData)),
                              SE = CharValues(getCharLength(SESystemData), getCharMass(SESystemData), getCharTime(SESystemData)))
    
    propagator = MBD.Propagator()
    endEvents = DifferentialEquations.VectorContinuousCallback(endConditions, endAffect!, 4)

    EarthRadius_EM::Float64 = primaries[1].bodyRadius/charValues.EM.lstar
    MoonRadius_EM::Float64 = primaries[2].bodyRadius/charValues.EM.lstar
    EarthHill_EM::Float64 = charValues.SE.lstar*cbrt(getMassRatio(SESystemData)/3)/charValues.EM.lstar
    MoonHill_EM::Float64 = cbrt(getMassRatio(EMSystemData)/3)

    return EscEnv(EDynamicsModel, EMDynamicsModel, EMEoMs, EMSDynamicsModel, SEDynamicsModel, primaries, charValues, endEvents, propagator, EarthHill_EM, EarthRadius_EM, MoonHill_EM, MoonRadius_EM)
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
        radius::Float64 = 0.3
        center::Vector{Float64} = getPrimaryState(env.EMDynamicsModel, primary)[1:2]
    else
        radius = 1.25
        center = [0.0, 0.0]
    end

    xGrid::Vector{Float64} = collect(range(center[1]-radius, center[1]+radius, n))
    yGrid::Vector{Float64} = collect(range(center[2]-radius, center[2]+radius, n))
    rRect::Matrix{StaticArrays.SVector{2, Float64}} = [StaticArrays.SA[x, y] for x in xGrid, y in yGrid]

    mu::Float64 = get12MassRatio(env.EMSDynamicsModel)
    rE::StaticArrays.SVector{2, Float64} = StaticArrays.SVector{2, Float64}(getPrimaryState(env.EMDynamicsModel, 1)[1:2])
    rM::StaticArrays.SVector{2, Float64} = StaticArrays.SVector{2, Float64}(getPrimaryState(env.EMDynamicsModel, 2)[1:2])
    
    lunarMask::BitMatrix = LinearAlgebra.norm.(rRect .- Ref(rM)) .> 1.25*env.MoonHill_EM
    mask::Matrix{Bool} = (primary == 2) ? .~lunarMask : (lunarMask .& isInterior.(rRect, mu, Ref(rE), Ref(rM)))

    rGrid::Vector{StaticArrays.SVector{2, Float64}} = rRect[mask]

    return rGrid
end

function isPeriapsis(dynamicsModel::MBD.CR3BPDynamicsModel, primary::Int64, r::StaticArrays.SVector{2, Float64}, dMag::Float64, v::Vector{Float64})
    mu::Float64 = (primary == 0 ? (dynamicsModel.systemData.primaryData[1].gravParam+dynamicsModel.systemData.primaryData[2].gravParam) : dynamicsModel.systemData.primaryData[primary].gravParam)
    vInert::StaticArrays.SVector{2, Float64} = StaticArrays.SA[v[1]-r[2], v[2]+r[1]]
    vInertMag::Float64 = LinearAlgebra.norm(vInert)
    vCirc::Float64 = sqrt(mu/(dMag*getCharLength(dynamicsModel)))*getCharTime(dynamicsModel)/getCharLength(dynamicsModel)

    return vInertMag > vCirc
end

function computeApseVelocities(env::EscEnv, JC::Float64, rGrid::Vector{StaticArrays.SVector{2, Float64}})
    OmegaGrid::Vector{Float64} = map(q -> getPseudopotential(env.EMDynamicsModel, Vector(push(q, 0.0))), rGrid)
    v2Grid::Vector{Float64} = 2 .* OmegaGrid .- JC

    return map(v2 -> (v2 < 0 ? NaN : sqrt(v2)), v2Grid)
end

function computeApseStates(env::EscEnv, primary::Int64, JC::Float64, thetaS::Float64, apse::Symbol, grade::Symbol, rGrid::Vector{StaticArrays.SVector{2, Float64}})
    center::StaticArrays.SVector{2, Float64} = (primary == 0 ? StaticArrays.SA[0.0, 0.0] : StaticArrays.SVector{2, Float64}(getPrimaryState(env.EMDynamicsModel, primary)[1:2]))
    dGrid::Vector{StaticArrays.SVector{2, Float64}} = rGrid .- Ref(center)
    dMagGrid::Vector{Float64} = LinearAlgebra.norm.(dGrid)
    thatGrid::Vector{StaticArrays.SVector{2, Float64}} = map(r -> StaticArrays.SA[-r[2], r[1]] ./ norm(r), dGrid)
    vMagGrid::Vector{Float64} = computeApseVelocities(env, JC, rGrid)
    qGrid::Vector{StaticArrays.MVector{7, Float64}} = map(r -> StaticArrays.MVector{7, Float64}([r[1], r[2], 0.0, NaN, NaN, NaN, thetaS]), rGrid)
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
    params::Vector{Any} = [apse, center, r_Earth, r_Moon, env.EarthHill_EM, env.EarthRadius_EM, env.MoonRadius_EM, primary]
    (_, eventTrackers::Vector{EventTracker}) = propagateWithEvents(env.propagator, env.endEvents, Vector(IC), [0, 12.0*pi], env.EMSDynamicsModel, [peri, apo, escape, crashEarth, crashMoon], params)
    
    return eventTrackers
end

function apseMapBCR4BP(env::EscEnv, JC::Float64, thetaS::Float64, primary::Int64, rGrid::Vector{StaticArrays.SVector{2, Float64}}, mf::MATLAB.MatFile; apse::Symbol = :peri, grade::Symbol = :pro)
    qGrid::Vector{StaticArrays.MVector{7, Float64}} = computeApseStates(env, primary, JC, thetaS, apse, grade, rGrid)

    flags::Vector{Int64} = fill(9, size(qGrid))
    flags[map(q -> isinf(q[4]), qGrid)] .= 8
    valid::Vector{Bool} = map(q -> (!isnan(q[4]) && !isinf(q[4])), qGrid) 
    qProp::Vector{StaticArrays.MVector{7, Float64}} = qGrid[valid]
    qMap::Vector{Int64} = findall(valid)
    counts::Vector{Int64} = zeros(Int64, size(qGrid))
    periapses::Vector{Vector{StaticArrays.SVector{7, Float64}}} = [StaticArrays.SVector{7, Float64}[] for _ in qGrid]
    apoapses::Vector{Vector{StaticArrays.SVector{7, Float64}}} = [StaticArrays.SVector{7, Float64}[] for _ in qGrid]
    println("Propagating $(length(qProp)) CR3BP trajectories with $(Threads.nthreads()) threads...")
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
        periapses[qMap[j]] = map(q -> StaticArrays.SVector{7, Float64}(q), eventTrackers[1].states)
        apoapses[qMap[j]] = map(q -> StaticArrays.SVector{7, Float64}(q), eventTrackers[2].states)
    end

    flagCounts::Vector{Int64} = [count(==(f), flags) for f in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]]
    println("\tDirect escapes:\t$(flagCounts[1])")
    println("\tLonger escapes:\t$(sum(flagCounts[2:6]))")
    println("\tCaptures:\t$(flagCounts[7])")
    println("\tCrashes:\t$(flagCounts[8])")
    println("\tInvalid apses:\t$(flagCounts[9])")
    println("\tZVCs:\t\t$(flagCounts[10])")

    println("Exporting map...")
    validPeri::Vector{Int64} = findall(!isempty, periapses)
    validApo::Vector{Int64} = findall(!isempty, apoapses)
    periStates::Vector{StaticArrays.SVector{7, Float64}} = reduce(vcat, periapses[validPeri])
    periIndices::Vector{Int64} = reduce(vcat, [fill(idx, length(periapses[idx])) for idx in validPeri])
    apoStates::Vector{StaticArrays.SVector{7, Float64}} = reduce(vcat, apoapses[validApo])
    apoIndices::Vector{Int64} = reduce(vcat, [fill(idx, length(apoapses[idx])) for idx in validApo])
    exportBCR4BPApseMap(env.EMSDynamicsModel, primary, apse, grade, JC, thetaS, qGrid, flags, counts, periStates, periIndices, apoStates, apoIndices, mf, Symbol("map_", replace(string(JC), "." => "_")))
end

function run_apseMapBCR4BP(JC::Float64, thetaS::Float64, n::Int64, primary::Int64; apse::Symbol = :peri, grade::Symbol = :pro)
    mf = MATLAB.MatFile("Output/ApseMaps/BCR4BP_$(string(primary))_$(string(apse))_$(string(grade))_$(string(n))_$(string(JC))_$(string(round(thetaS, digits = 3))).mat", "w")
        
    env::EscEnv = setupEnvironment()

    rGrid::Vector{StaticArrays.SVector{2, Float64}} = getGrid(env, n, primary)
    apseMapBCR4BP(env, JC, thetaS, primary, rGrid, mf; apse = apse, grade = grade)
    
    MATLAB.close(mf)
end

function run_apseMapsthetaBCR4BP(JC::Float64, thetaS::Vector{Float64}, n::Int64, primary::Int64; apse::Symbol = :peri, grade::Symbol = :pro)
    mf = MATLAB.MatFile("Output/ApseMaps/BCR4BPthetaVolume_$(string(primary))_$(string(apse))_$(string(grade))_$(string(n))_$(string(JC)).mat", "w")

    env::EscEnv = setupEnvironment()

    rGrid::Vector{StaticArrays.SVector{2, Float64}} = getGrid(env, n, primary)
    o::Int64 = length(thetaS)
    for j::Int64 in eachindex(thetaS)
        println("\nProducing map $j / $o: theta = $(thetaS[j]*180.0/pi)...")
        apseMapBCR4BP(env, JC, thetaS[j], primary, rGrid, mf; apse = apse, grade = grade)
    end

    MATLAB.close(mf)
end

function run_apseMapsJCBCR4BP(JCs::Vector{Float64}, thetaS::Float64, n::Int64, primary::Int64; apse::Symbol = :peri, grade::Symbol = :pro)
    mf = MATLAB.MatFile("Output/ApseMaps/BCR4BPJCVolume_$(string(primary))_$(string(apse))_$(string(grade))_$(string(n))_$(string(round(thetaS, digits = 3)))_$(string(minimum(JCs)))_$(string(maximum(JCs))).mat", "w")

    env::EscEnv = setupEnvironment()

    rGrid::Vector{StaticArrays.SVector{2, Float64}} = getGrid(env, n, primary)
    o::Int64 = length(JCs)
    for j::Int64 in eachindex(JCs)
        println("\nProducing map $j / $o: JC = $(JCs[j])...")
        apseMapBCR4BP(env, JCs[j], thetaS, primary, rGrid, mf; apse = apse, grade = grade)
    end

    MATLAB.close(mf)
end

end # EscBCR4BP
