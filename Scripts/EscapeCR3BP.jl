"""
Script for computing CR3BP escape trajectories in the Earth-Moon system

Author: Jonathan LeFevre Richmond
C: 6/16/26
U: 7/1/26
"""

module EscCR3BP

using MBD, DifferentialEquations, HDF5, LinearAlgebra, Logging, MATLAB, StaticArrays

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

    escapeEvent
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

function escapeCondition(state::Vector{Float64}, time::Float64, integrator)
    LinearAlgebra.norm(state[1:3])-integrator.p[2]
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
    escapeEvent = DifferentialEquations.ContinuousCallback(escapeCondition, terminateAffect!)

    EarthRadius_EM::Float64 = primaries[1].bodyRadius/charValues.EM.lstar
    MoonRadius_EM::Float64 = primaries[2].bodyRadius/charValues.EM.lstar
    EarthHill_EM::Float64 = charValues.SE.lstar*cbrt(getMassRatio(SESystemData)/3)/charValues.EM.lstar

    return EscEnv(EMDynamicsModel, SEDynamicsModel, primaries, Sun, charValues, escapeEvent, endEvents, propagator, EarthHill_EM, EarthRadius_EM, MoonRadius_EM)
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

    return maximum(E)
end

function pruneVolumeData(JC::Float64, indices::Vector{Int64}, volFileName::String)
    HDF5.h5open(volFileName, "r") do file
        volJC = vec(read(file["JCVolume"]))
        JCIdx::Int64 = findfirst(j -> j < JC, volJC)
        # newVolJC::Vector{Float64} = volJC[JCIdx:end,1]

        nJC::Int64 = length(volJC)-JCIdx+1
        nIdx::Int64 = length(indices)

        volFlags = file["flagsVolume"]
        newVolFlags = Matrix{Int64}(undef, (nJC,nIdx))
        volqs = file["qVolume"]
        newVolqs = Array{Float64, 3}(undef, (6,nIdx,nJC))
        for (i::Int64, idx::Int64) in enumerate(indices)
            newVolFlags[:,i] = volFlags[JCIdx:end,idx]
            newVolqs[:,i,:] = volqs[:,idx,JCIdx:end]
        end

        return (newVolFlags, newVolqs)
    end
end

function escapeAnalysisCR3BP(env::EscEnv, JC::Float64, flags::Matrix{Int64}, qs::Matrix{Float64}, volFileName::String, idx::Int64, mf::MATLAB.MatFile)
    esc0Indices::Vector{Int64} = findall(f -> f == 0, vec(flags))
    n_esc0::Int64 = length(esc0Indices)
    esc0q_0s::Matrix{Float64} = qs[:,esc0Indices]
    esc0t_fs::Vector{Float64} = Vector{Float64}(undef, n_esc0)
    esc0Es::Vector{Float64} = Vector{Float64}(undef, n_esc0)
    Threads.@threads for idx::Int64 in eachindex(esc0Indices)
        arc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.escapeEvent, qs[:,esc0Indices[idx]], [0, 12.0*pi], env.EMDynamicsModel, [env.EarthHill_EM])
        q_f::Vector{Float64} = getStateByIndex(arc, -1)
        esc0t_fs[idx] = getTimeByIndex(arc, -1)
        esc0Es[idx] = getMaxHeliocentricEnergy(env, q_f)
    end

    MATLAB.put_variable(mf, :esc0q0, esc0q_0s)
    MATLAB.put_variable(mf, :esc0tf, esc0t_fs)
    MATLAB.put_variable(mf, :esc0E, esc0Es)

    (volFlags::Matrix{Int64}, volqs::Array{Float64, 3}) = pruneVolumeData(JC, [idx], volFileName)
    @views vs = sqrt.(qs[4,:].^2 .+ qs[5,:].^2)
    Deltavs::Vector{Float64} = fill(NaN, (length(vec(volFlags))))
    escEs::Vector{Float64} = copy(Deltavs)
    escIndices::Vector{Int64} = findall(f -> f < 6, vec(volFlags))
    Threads.@threads for escIdx::Int64 in escIndices
        q::Vector{Float64} = volqs[:,1,escIdx]
        Deltavs[escIdx] = sqrt(q[4]^2+q[5]^2)-vs[idx]
        arc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.escapeEvent, q, [0, 12.0*pi], env.EMDynamicsModel, [env.EarthHill_EM])
        q_f::Vector{Float64} = getStateByIndex(arc, -1)
        escEs[escIdx] = getMaxHeliocentricEnergy(env, q_f)
    end

    MATLAB.put_variable(mf, :Deltav1s, Deltavs)
    MATLAB.put_variable(mf, :EscapeEs, escEs)

    esc1Indices::Vector{Int64} = findall(f -> f == 1, vec(flags))
    n_esc1::Int64 = length(esc1Indices)
    esc1q_0s::Matrix{Float64} = qs[:,esc1Indices]
    esc1t_fs::Vector{Float64} = Vector{Float64}(undef, n_esc1)
    esc1Es::Vector{Float64} = Vector{Float64}(undef, n_esc1)
    Threads.@threads for idx::Int64 in eachindex(esc1Indices)
        arc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.escapeEvent, qs[:,esc1Indices[idx]], [0, 12.0*pi], env.EMDynamicsModel, [env.EarthHill_EM])
        q_f::Vector{Float64} = getStateByIndex(arc, -1)
        esc1t_fs[idx] = getTimeByIndex(arc, -1)
        esc1Es[idx] = getMaxHeliocentricEnergy(env, q_f)
    end

    MATLAB.put_variable(mf, :esc1q0, esc1q_0s)
    MATLAB.put_variable(mf, :esc1tf, esc1t_fs)
    MATLAB.put_variable(mf, :esc1E, esc1Es)

    esc2Indices::Vector{Int64} = findall(f -> f == 2, vec(flags))
    n_esc2::Int64 = length(esc2Indices)
    esc2q_0s::Matrix{Float64} = qs[:,esc2Indices]
    esc2t_fs::Vector{Float64} = Vector{Float64}(undef, n_esc2)
    esc2Es::Vector{Float64} = Vector{Float64}(undef, n_esc2)
    Threads.@threads for idx::Int64 in eachindex(esc2Indices)
        arc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.escapeEvent, qs[:,esc2Indices[idx]], [0, 12.0*pi], env.EMDynamicsModel, [env.EarthHill_EM])
        q_f::Vector{Float64} = getStateByIndex(arc, -1)
        esc2t_fs[idx] = getTimeByIndex(arc, -1)
        esc2Es[idx] = getMaxHeliocentricEnergy(env, q_f)
    end

    MATLAB.put_variable(mf, :esc2q0, esc2q_0s)
    MATLAB.put_variable(mf, :esc2tf, esc2t_fs)
    MATLAB.put_variable(mf, :esc2E, esc2Es)
end

function assisted0EscapeAnalysisCR3BP(JC::Float64, flags::Matrix{Int64}, qs::Matrix{Float64}, volFileName::String, mf::MATLAB.MatFile)
    nonEscIndices::Vector{Int64} = findall(f -> (f == 6) || (f == 7), vec(flags))
    (volFlags::Matrix{Int64}, volqs::Array{Float64, 3}) = pruneVolumeData(JC, nonEscIndices, volFileName)
    @views vs = sqrt.(qs[4,:].^2 .+ qs[5,:].^2)
    Deltavs::Vector{Float64} = zeros(Float64, length(vs))
    Threads.@threads for j::Int64 in eachindex(nonEscIndices)
        idx::Int64 = nonEscIndices[j]
        escIdx = findfirst(f -> f < 6, @view volFlags[:,j])
        if escIdx !== nothing
            Deltavs[idx] = sqrt(volqs[4,j,escIdx]^2+volqs[5,j,escIdx]^2)-vs[idx]
        else
            Deltavs[idx] = NaN
        end
    end
    Threads.@threads for j::Int64 in eachindex(flags)
        if (flags[j] == 8) || (flags[j] == 9)
            Deltavs[j] = NaN
        end
    end

    MATLAB.put_variable(mf, :Deltav0s, Deltavs)
end

function assisted1EscapeAnalysisCR3BP(env::EscEnv, JC::Float64, flags::Matrix{Int64}, qs::Matrix{Float64}, volFileName::String, mf::MATLAB.MatFile)
    # directIndices::Vector{Int64} = findall(f -> f == 0, vec(flags))
    # (volFlags::Matrix{Int64}, volqs::Array{Float64, 3}) = pruneVolumeData(JC, directIndices, volFileName)
    # @views vs = sqrt.(qs[4,:].^2 .+ qs[5,:].^2)
    # Deltavs::Matrix{Float64} = fill(NaN, (size(volFlags)))
    # escEs::Matrix{Float64} = fill(NaN, (size(volFlags)))
    # Threads.@threads for j::Int64 in eachindex(directIndices)
    #     idx::Int64 = directIndices[j]
    #     escIndices::Vector{Int64} = findall(f -> f < 6, @view volFlags[:,j])
    #     for escIdx::Int64 in escIndices
    #         q::Vector{Float64} = volqs[:,j,escIdx]
    #         Deltavs[escIdx,j] = sqrt(q[4]^2+q[5]^2)-vs[idx]
    #         arc::MBD.CR3BPArc = propagateWithEvent(env.propagator, env.escapeEvent, q, [0, 12.0*pi], env.EMDynamicsModel, [env.EarthHill_EM])
    #         q_f::Vector{Float64} = getStateByIndex(arc, -1)
    #         escEs[escIdx,j] = getMaxHeliocentricEnergy(env, q_f)
    #     end
    # end

    # MATLAB.put_variable(mf, :Deltav1s, Deltavs)
    # MATLAB.put_variable(mf, :EscapeEs, escEs)
end

function assistedEscapeAnalysisCR3BP(env::EscEnv, JC::Float64, flags::Matrix{Int64}, qs::Matrix{Float64}, maneuvers::Int64, volFileName::String, mf::MATLAB.MatFile)
    if maneuvers == 0
        assisted0EscapeAnalysisCR3BP(JC, flags, qs, volFileName, mf)
    elseif maneuvers == 1
        assisted1EscapeAnalysisCR3BP(env, JC, flags, qs, volFileName, mf)
    end
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

function run_escapeAnalysisCR3BP(fileName::String, mapName::String, volFileName::String, idx::Int64)
    mf_in = MATLAB.MatFile(fileName, "r")

    map::Dict{String, Any} = get_variable(mf_in, mapName)

    MATLAB.close(mf_in)

    JC::Float64 = map["JC"]
    qs::Matrix{Float64} = map["q"]
    flags::Matrix{Int64} = map["flags"]

    mf_out = MATLAB.MatFile("Output/EscapeAnalysisCR3BP.mat", "w")

    env::EscEnv = setupEnvironment()

    escapeAnalysisCR3BP(env, JC, flags, qs, volFileName, idx, mf_out)

    MATLAB.close(mf_out)
end

function run_assistedEscapeAnalysisCR3BP(fileName::String, mapName::String, maneuvers::Int64, volFileName::String)
    mf_in = MATLAB.MatFile(fileName, "r")

    map::Dict{String, Any} = get_variable(mf_in, mapName)

    MATLAB.close(mf_in)

    JC::Float64 = map["JC"]
    qs::Matrix{Float64} = map["q"]
    flags::Matrix{Int64} = map["flags"]

    mf_out = MATLAB.MatFile("Output/AssistedEscapeAnalysisCR3BP.mat", "w")

    env::EscEnv = setupEnvironment()

    assistedEscapeAnalysisCR3BP(env, JC, flags, qs, maneuvers, volFileName, mf_out)

    MATLAB.close(mf_out)
end

end # EscCR3BP
