"""
Script for computing periapsis maps

Author: Jonathan Richmond
C: 10/9/25
U: 11/4/25
"""
# module PeriMap
println("Running PeriapsisMap.jl...\n")

using MBD, DifferentialEquations, LinearAlgebra, Logging, MATLAB, StaticArrays

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../CR3BPTargeters/PlanarPerpJC.jl")
include("../Utilities/Export.jl")

mutable struct Periapsis
    count::Int64
    states::Vector{Vector{Float64}}
end

function endConditionsBCR4BP(out::SubArray{Float64}, state::Vector{Float64}, time::Float64, integrator)
    if time < 1E-8
        out[1:4] = [1.0, 1.0, 1.0, 1.0]
    else
        q_B1::Vector{Float64} = [1-integrator.p[4], 0, 0, 0, 0, 0]
        q_E::Vector{Float64} = getPrimaryState(integrator.p[3], 1, state[7])
        q_M::Vector{Float64} = getPrimaryState(integrator.p[3], 2, state[7])
        q_2::Vector{Float64} = state[1:6]-q_B1
    
        out[1] = LinearAlgebra.dot(q_2[1:3], q_2[4:6])
        out[2] = LinearAlgebra.norm(q_2[1:3])-integrator.p[5]
        out[3] = integrator.p[6]-LinearAlgebra.norm(state[1:3]-q_E[1:3])
        out[4] = integrator.p[7]-LinearAlgebra.norm(state[1:3]-q_M[1:3])
    end
end

function endConditionsCR3BP(out::SubArray{Float64}, state::Vector{Float64}, time::Float64, integrator)
    if time < 1E-8
        out[1:4] = [1.0, 1.0, 1.0, 1.0]
    else
        q_E::Vector{Float64} = getPrimaryState(integrator.p[3], 1)
        q_M::Vector{Float64} = getPrimaryState(integrator.p[3], 2)
    
        out[1] = LinearAlgebra.dot(state[1:3], state[4:6])
        out[2] = LinearAlgebra.norm(state[1:3])-integrator.p[4]
        out[3] = integrator.p[5]-LinearAlgebra.norm(state[1:3]-q_E[1:3])
        out[4] = integrator.p[6]-LinearAlgebra.norm(state[1:3]-q_M[1:3])
    end
end

function periapsisCondition(state::Vector{Float64}, time::Float64, integrator)
    if abs(time) < 1E-8
        -1.0
    else
        LinearAlgebra.dot(state[1:3], state[4:6])
    end
end

function endAffectBCR4BP!(integrator, index)
    if index == 1
        integrator.p[8].count += 1
        if integrator.p[8].count >= 3
            integrator.p[2].flag = :capture
            DifferentialEquations.terminate!(integrator)
        end
    else
        if index == 2
            if integrator.p[8].count == 0
                integrator.p[2].flag = :escape
            elseif integrator.p[8].count == 1
                integrator.p[2].flag = :escape1
            elseif integrator.p[8].count == 2
                integrator.p[2].flag = :escape2
            end
        elseif index == 3
            integrator.p[2].flag = :earth
        elseif index == 4
            integrator.p[2].flag = :moon
        end
        DifferentialEquations.terminate!(integrator)
    end
end

function endAffectCR3BP!(integrator, index)
    if index == 1
        integrator.p[7].count += 1
        if integrator.p[7].count >= 3
            integrator.p[2].flag = :capture
            DifferentialEquations.terminate!(integrator)
        end
    else
        if index == 2
            if integrator.p[7].count == 0
                integrator.p[2].flag = :escape
            elseif integrator.p[7].count == 1
                integrator.p[2].flag = :escape1
            elseif integrator.p[7].count == 2
                integrator.p[2].flag = :escape2
            end
        elseif index == 3
            integrator.p[2].flag = :earth
        elseif index == 4
            integrator.p[2].flag = :moon
        end
        DifferentialEquations.terminate!(integrator)
    end
end

function endAffectManifold!(integrator)
    integrator.p[2].count += 1
    push!(integrator.p[2].states, copy(integrator.u))
    if integrator.p[2].count >= 6
        DifferentialEquations.terminate!(integrator)
    end
end

mf = MATLAB.MatFile("Output/PeriapsisMap.mat", "w")

EMCR3BPSystemData = MBD.CR3BPSystemData("Earth", "Moon")
SB1CR3BPSystemData = MBD.CR3BPSystemData("Sun", "Earth_Barycenter")
systemData = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
EMCR3BPDynamicsModel = MBD.CR3BPDynamicsModel(EMCR3BPSystemData)
SB1CR3BPDynamicsModel = MBD.CR3BPDynamicsModel(SB1CR3BPSystemData)
EMDynamicsModel = MBD.BCR4BP12DynamicsModel(systemData)
SB1DynamicsModel = MBD.BCR4BP41DynamicsModel(systemData)

propagator = MBD.Propagator()
endEventsBCR4BP = DifferentialEquations.VectorContinuousCallback(endConditionsBCR4BP, endAffectBCR4BP!, nothing, 4)
endEventsCR3BP = DifferentialEquations.VectorContinuousCallback(endConditionsCR3BP, endAffectCR3BP!, nothing, 4)
manifoldEvent = DifferentialEquations.ContinuousCallback(periapsisCondition, nothing, endAffectManifold!)

JCEM::Float64 = 3.0663
thetaM::Float64 = pi*0
R_H::Float64 = get41CharLength(systemData)*(systemData.primaryData[1].mass/(3*(systemData.primaryData[3].mass+systemData.primaryData[1].mass)))^(1/3)

radius::Float64 = 0.0035
# radius::Float64 = 0.00075
n::Int64 = 500

x_B1::Float64 = 1-get41MassRatio(systemData)
x_Moon::Float64 = 1-get41MassRatio(systemData)+(1-get12MassRatio(systemData))*get12CharLength(systemData)/get41CharLength(systemData)
xSB1::Vector{Float64} = range(x_B1-radius, x_B1+radius, n)
# xSB1::Vector{Float64} = range(x_Moon-radius, x_Moon+radius, n)
ySB1::Vector{Float64} = range(-radius, radius, n)
posSB1::Vector{StaticArrays.SVector{2, Float64}} = vec([StaticArrays.SVector(x, y) for x in xSB1, y in ySB1])
deltarSB1::Vector{StaticArrays.SVector{2, Float64}} = posSB1 .- Ref(StaticArrays.SVector{2, Float64}([1-get41MassRatio(systemData), 0]))
rSB1::Vector{Float64} = LinearAlgebra.norm.(deltarSB1)
rhatSB1::Vector{StaticArrays.SVector{2, Float64}} = deltarSB1 ./ rSB1
thatSB1::Vector{StaticArrays.SVector{2, Float64}} = [StaticArrays.SVector(-v[2], v[1]) for v in rhatSB1] # Prograde
# thatSB1::Vector{StaticArrays.SVector{2, Float64}} = [StaticArrays.SVector(v[2], -v[1]) for v in rhatSB1] # Retrograde
qGuessSB1::Vector{Vector{Float64}} = [[posSB1[j]..., 0, thatSB1[j]..., 0, thetaM] for j in eachindex(posSB1)]
qGuessEM::Vector{Vector{Float64}} = rotating41ToRotating12(SB1DynamicsModel, qGuessSB1, zeros(length(qGuessSB1)))[1]
OmegaEM::Vector{Float64} = map(q -> getPseudopotential(EMCR3BPDynamicsModel, q[1:3]), qGuessEM)
v2EM::Vector{Float64} = 2 .* OmegaEM .- JCEM
v2EM[v2EM .< 0] .= NaN
vEM::Vector{Float64} = sqrt.(v2EM)
Threads.@threads for j in eachindex(vEM)
    qGuessEM[j][4:5] = qGuessEM[j][4:5] .* vEM[j] ./ norm(qGuessEM[j][4:5])
end
qSB1::Vector{Vector{Float64}} = rotating12ToRotating41(EMDynamicsModel, qGuessEM, zeros(length(qGuessEM)))[1]
flagsSB1::Vector{Int64} = Vector{Int64}(undef, length(qSB1))
Threads.@threads for q in eachindex(qSB1)
    if any(x -> isnan(x), qSB1[q])
        flagsSB1[q] = 1
    else
        (arc::MBD.BCR4BP41Arc, event) = propagateWithEvents(propagator, endEventsBCR4BP, qSB1[q], [0, 4.0*pi], SB1DynamicsModel, [SB1DynamicsModel, get41MassRatio(systemData), R_H/get41CharLength(systemData), systemData.primaryData[1].bodyRadius/get41CharLength(systemData), systemData.primaryData[2].bodyRadius/get41CharLength(systemData), Periapsis(0, [])])
        if event == :earth
            flagsSB1[q] = 2
        elseif event == :moon
            flagsSB1[q] = 2
        elseif event == :capture
            flagsSB1[q] = 3
        elseif event == :escape
            flagsSB1[q] = 4
        elseif event == :escape1
            flagsSB1[q] = 5
        elseif event == :escape2
            flagsSB1[q] = 6
        else
            flagsSB1[q] = 0
        end
    end
end

xEM::Vector{Float64} = range(-radius, radius, n) .* get41CharLength(systemData) ./ get12CharLength(systemData)
# xEM::Vector{Float64} = range(-radius, radius, n) .* get41CharLength(systemData) ./ get12CharLength(systemData) .+ (1-get12MassRatio(systemData))
yEM::Vector{Float64} = range(-radius, radius, n) .* get41CharLength(systemData) ./ get12CharLength(systemData)
deltar::Vector{Vector{Float64}} = [[x, y] for x in xEM for y in yEM]
r::Vector{Float64} = LinearAlgebra.norm.(deltar)
rhat::Vector{StaticArrays.SVector{2, Float64}} = deltar ./ r
that::Vector{StaticArrays.SVector{2, Float64}} = [StaticArrays.SVector(-v[2], v[1]) for v in rhat] # Prograde
# that::Vector{StaticArrays.SVector{2, Float64}} = [StaticArrays.SVector(v[2], -v[1]) for v in rhat] # Retrograde
Omega::Vector{Float64} = map(q -> getPseudopotential(EMCR3BPDynamicsModel, push!(copy(q), 0.0)), deltar)
v2::Vector{Float64} = 2 .* Omega .- JCEM
v2[v2 .< 0] .= NaN
v::Vector{Float64} = sqrt.(v2)
qEM::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(v))
Threads.@threads for j in eachindex(v)
    qEM[j] = append!(copy(deltar[j]), [0.0], v[j] .* that[j], [0.0])
end
flags::Vector{Int64} = Vector{Int64}(undef, length(qEM))
Threads.@threads for q in eachindex(qEM)
    if any(x -> isnan(x), qEM[q])
        flags[q] = 1
    else
        (arc::MBD.CR3BPArc, event) = propagateWithEvents(propagator, endEventsCR3BP, qEM[q], [0, 4*pi*get41CharTime(systemData)/get12CharTime(systemData)], EMCR3BPDynamicsModel, [EMCR3BPDynamicsModel, R_H/get12CharLength(systemData), systemData.primaryData[1].bodyRadius/get12CharLength(systemData), systemData.primaryData[2].bodyRadius/get12CharLength(systemData), Periapsis(0, [])])
        if event == :earth
            flags[q] = 2
        elseif event == :moon
            flags[q] = 2
        elseif event == :capture
            flags[q] = 3
        elseif event == :escape
            flags[q] = 4
        elseif event == :escape1
            flags[q] = 5
        elseif event == :escape2
            flags[q] = 6
        else
            flags[q] = 0
        end
    end
end

targeter = PlanarPerpJCTargeter(EMCR3BPDynamicsModel)
familyFile::String = "FamilyData/CR3BPEML2Lyapunovs.csv"
orbit::MBD.CR3BPPeriodicOrbit = interpOrbit(targeter, familyFile, "JC", JCEM; choiceIndex = 1)
manifold::MBD.CR3BPManifold = getManifoldByArclength(orbit, "Stable", "Negative", 25/get12CharLength(systemData), n*10)
nThreads::Int64 = Threads.nthreads()
qMan_thread::Vector{Vector{Vector{Float64}}} = [Vector{Vector{Float64}}() for _ in 1:nThreads]
Threads.@threads for a in eachindex(manifold.orbitTimes)
    tid = mod1(Threads.threadid(), nThreads)
    peris = Periapsis(0, [])
    arc::MBD.CR3BPArc = propagateWithEvent(propagator, manifoldEvent, real(manifold.initialConditions[a]), [0, -4*pi*get41CharTime(systemData)/get12CharTime(systemData)], EMCR3BPDynamicsModel, [peris])
    for p in eachindex(peris.states)
        push!(qMan_thread[tid], peris.states[p])
    end
end
qMan::Vector{Vector{Float64}} = reduce(vcat, qMan_thread)

# sample::Int64 = 21567
# (arc41::MBD.BCR4BP41Arc, event) = propagateWithEvents(propagator, endEventsBCR4BP, qSB1[sample], [0, 4.0*pi], SB1DynamicsModel, [SB1DynamicsModel, get41MassRatio(systemData), R_H/get41CharLength(systemData), systemData.primaryData[1].bodyRadius/get41CharLength(systemData), systemData.primaryData[2].bodyRadius/get41CharLength(systemData), Periapsis(0, [])])
# exportBCR4BP41Trajectory(qSB1[sample]..., getTimeByIndex(arc41, -1), SB1DynamicsModel, mf, :Traj41)
# q::Vector{Float64} = rotating41ToRotating12(SB1DynamicsModel, [qSB1[sample]], [0.0])[1][1]
# exportBCR4BP12Trajectory(q..., getTimeByIndex(arc41, -1)*get41CharTime(systemData)/get12CharTime(systemData), EMDynamicsModel, mf, :Traj12)

MATLAB.put_variable(mf, :pointsSB1, qSB1)
MATLAB.put_variable(mf, :pointsEM, qEM)
MATLAB.put_variable(mf, :manifold, qMan)
MATLAB.put_variable(mf, :flagsSB1, flagsSB1)
MATLAB.put_variable(mf, :flagsEM, flags)
MATLAB.put_variable(mf, :JC, JCEM)
MATLAB.put_variable(mf, :moonAngle, thetaM)
MATLAB.close(mf)

println()
# end
