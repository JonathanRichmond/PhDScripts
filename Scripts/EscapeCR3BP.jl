"""
Script for computing CR3BP escape trajectories in the Earth-Moon system

Author: Jonathan LeFevre Richmond
C: 6/16/26
U: 6/18/26
"""

module EscCR3BP

using MBD, DifferentialEquations, LinearAlgebra, Logging, MATLAB, StaticArrays

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../Utilities/Export.jl")

struct CharValues
    lstar::Float64
    mstar::Float64
    tstar::Float64
end

struct EscEnv
    EMDynamicsModel::MBD.CR3BPDynamicsModel
    SEDynamicsModel::MBD.CR3BPDynamicsModel

    primaries::Vector{MBD.BodyData}
    Sun::MBD.BodyData

    charValues::NamedTuple

    endEvents
    propagator::MBD.Propagator

    EarthHill_EM::Float64
    EarthRadius_EM::Float64
    MoonRadius_EM::Float64
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

function setupEnvironment()::EscEnv
    EMSystemData = MBD.CR3BPSystemData("Earth", "Moon")
    SESystemData = MBD.CR3BPSystemData("Sun", "Earth")
    EMDynamicsModel = MBD.CR3BPDynamicsModel(EMSystemData)
    SEDynamicsModel = MBD.CR3BPDynamicsModel(SESystemData)

    primaries::Vector{MBD.BodyData} = EMSystemData.primaryData
    Sun::MBD.BodyData = SESystemData.primaryData[1]

    charValues::NamedTuple = (EM = CharValues(getCharLength(EMSystemData), getCharMass(EMSystemData), getCharTime(EMSystemData)),
                              SE = CharValues(getCharLength(SESystemData), getCharMass(SESystemData), getCharTime(SESystemData)))
    
    propagator = MBD.Propagator()
    endEvents = DifferentialEquations.VectorContinuousCallback(endConditions, endAffect!, 4)

    EarthRadius_EM::Float64 = primaries[1].bodyRadius/charValues.EM.lstar
    MoonRadius_EM::Float64 = primaries[2].bodyRadius/charValues.EM.lstar
    EarthHill_EM::Float64 = charValues.SE.lstar*cbrt(getMassRatio(SESystemData)/3)/charValues.EM.lstar

    return EscEnv(EMDynamicsModel, SEDynamicsModel, primaries, Sun, charValues, endEvents, propagator, EarthHill_EM, EarthRadius_EM, MoonRadius_EM)
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

function computeApseStates(env::EscEnv, primary::Int64, JC::Float64, apse::Symbol, grade::Symbol, rGrid::Matrix{StaticArrays.SVector{2, Float64}})
    center::StaticArrays.SVector{2, Float64} = (primary == 0 ? StaticArrays.SA[0.0, 0.0] : StaticArrays.SVector{2, Float64}(getPrimaryState(env.EMDynamicsModel, primary)[1:2]))

    dGrid::Matrix{StaticArrays.SVector{2, Float64}} = rGrid .- Ref(center)
    dMagGrid::Matrix{Float64} = LinearAlgebra.norm.(dGrid)
    thatGrid::Matrix{StaticArrays.SVector{2, Float64}} = map(r -> StaticArrays.SA[-r[2], r[1]] ./ norm(r), dGrid)
    OmegaGrid::Matrix{Float64} = map(q -> getPseudopotential(env.EMDynamicsModel, Vector(push(q, 0.0))), rGrid)
    v2Grid::Matrix{Float64} = 2 .* OmegaGrid .- JC
    vMagGrid::Matrix{Float64} = map(v2 -> (v2 < 0 ? NaN : sqrt(v2)), v2Grid)
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

function run_apseMapCR3BP(JC::Float64, n::Int64, primary::Int64; apse::Symbol = :peri, grade::Symbol = :pro)
    mf = MATLAB.MatFile("Output/ApseMaps/TestMap.mat", "w")
        
    env::EscEnv = setupEnvironment()

    rGrid::Matrix{StaticArrays.SVector{2, Float64}} = getGrid(env, n, primary)
    apseMapCR3BP(env, JC, n, primary, rGrid, mf; apse = apse, grade = grade)
    
    MATLAB.close(mf)
end

function run_apseMapsCR3BP(JCs::Vector{Float64}, n::Int64, primary::Int64; apse::Symbol = :peri, grade::Symbol = :pro)
    mf = MATLAB.MatFile("Output/ApseMaps/CR3BPJCVolume_$(string(primary))_$(string(apse))_$(string(grade))_$(string(n))_$(minimum(JCs))_$(maximum(JCs)).mat", "w")

    env::EscEnv = setupEnvironment()

    rGrid::Matrix{StaticArrays.SVector{2, Float64}} = getGrid(env, n, primary)
    o::Int64 = length(JCs)
    for j::Int64 in eachindex(JCs)
        println("\nProducing map $j / $o: JC = $(JCs[j])...")
        apseMapCR3BP(env, JCs[j], n, primary, rGrid, mf; apse = apse, grade = grade)
    end

    MATLAB.close(mf)
end

end # EscCR3BP
