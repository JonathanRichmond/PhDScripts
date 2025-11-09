"""
Script for computing periapsis maps

Author: Jonathan Richmond
C: 10/9/25
U: 11/8/25
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
    if time < 1E-10
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
    if time < 1E-10
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
    if abs(time) < 1E-10
        -1.0
    else
        LinearAlgebra.dot(state[1:3], state[4:6])
    end
end

function endAffectBCR4BP!(integrator, index)
    if index == 1
        integrator.p[8].count += 1
    else
        if index == 2
            integrator.p[2].flag = Symbol("escape", integrator.p[8].count)
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
    else
        if index == 2
            integrator.p[2].flag = Symbol("escape", integrator.p[7].count)
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

function matRotating12ToRotating41(dynamicsModel::MBD.BCR4BP12DynamicsModel, states12::Vector{Vector{Float64}}, times12::Vector{Float64})
    numStates::Int64 = length(states12)
    states12_mat::Matrix{Float64} = Matrix{Float64}(undef, 7, numStates)
    Threads.@threads for s::Int64 in 1:numStates
        @inbounds @simd for j in 1:7
            states12_mat[j,s] = states12[s][j]
        end
    end
    theta4::Float64 = states12[1][7]
    mu41::Float64 = get41MassRatio(dynamicsModel.systemData)
    m4::Float64 = get4Mass(dynamicsModel.systemData)
    a4::Float64 = get4Distance(dynamicsModel.systemData)
    theta4dot::Float64 = sqrt((m4+1)/(a4^3))-1
    C::Matrix{Float64} = [-cos(theta4) -sin(theta4) 0; sin(theta4) -cos(theta4) 0; 0 0 1]
    Cdot::Matrix{Float64} = theta4dot.*[sin(theta4) -cos(theta4) 0; cos(theta4) sin(theta4) 0; 0 0 0]
    N::Matrix{Float64} = [(1/a4).*C zeros(Float64, 3, 4); sqrt(a4/(m4+1)).*Cdot sqrt(a4/(m4+1)).*C zeros(Float64, 3, 1); zeros(Float64, 1, 6) -1]
    rho::Vector{Float64} = [1-mu41; zeros(Float64, 5); pi]
    states41_mat::Matrix{Float64} = N*states12_mat .+ rho
    states41::Vector{Vector{Float64}} = [@view states41_mat[:, c] for c in axes(states41_mat, 2)]
    times41::Vector{Float64} = (get12CharTime(dynamicsModel.systemData)/get41CharTime(dynamicsModel.systemData)).*times12

    return (states41, times41)
end

function matRotating41ToRotating12(dynamicsModel::MBD.BCR4BP41DynamicsModel, states41::Vector{Vector{Float64}}, times41::Vector{Float64})
    numStates::Int64 = length(states41)
    states41_mat::Matrix{Float64} = Matrix{Float64}(undef, 7, numStates)
    Threads.@threads for s::Int64 in 1:numStates
        @inbounds @simd for j in 1:7
            states41_mat[j,s] = states41[s][j]
        end
    end
    theta2::Float64 = states41[1][7]
    mu41::Float64 = get41MassRatio(dynamicsModel.systemData)
    m4::Float64 = get4Mass(dynamicsModel.systemData)
    a4::Float64 = get4Distance(dynamicsModel.systemData)
    theta2dot::Float64 = sqrt((a4^3)/(m4+1))-1
    C::Matrix{Float64} = [cos(theta2) sin(theta2) 0; -sin(theta2) cos(theta2) 0; 0 0 1]
    Cdot::Matrix{Float64} = theta2dot.*[-sin(theta2) cos(theta2) 0; -cos(theta2) -sin(theta2) 0; 0 0 0]
    N::Matrix{Float64} = [a4.*C zeros(Float64, 3, 4); sqrt((m4+1)/a4).*Cdot sqrt((m4+1)/a4).*C zeros(Float64, 3, 1); zeros(Float64, 1, 6) -1]
    rho::Vector{Float64} = [-1+mu41; zeros(Float64, 6)]
    psi::Vector{Float64} = [zeros(Float64, 6); pi]
    states12_mat::Matrix{Float64} = N*(states41_mat .+ rho) .+ psi
    states12::Vector{Vector{Float64}} = [@view states12_mat[:, c] for c in axes(states12_mat, 2)]
    times12::Vector{Float64} = (get41CharTime(dynamicsModel.systemData)/get12CharTime(dynamicsModel.systemData)).*times41

    return (states12, times12)
end

mf = MATLAB.MatFile("Output/PeriapsisMap.mat", "w")

EMCR3BPSystemData = MBD.CR3BPSystemData("Earth", "Moon")
SB1CR3BPSystemData = MBD.CR3BPSystemData("Sun", "Earth_Barycenter")
systemData = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
EMCR3BPDynamicsModel = MBD.CR3BPDynamicsModel(EMCR3BPSystemData)
SB1CR3BPDynamicsModel = MBD.CR3BPDynamicsModel(SB1CR3BPSystemData)
EMDynamicsModel = MBD.BCR4BP12DynamicsModel(systemData)
SB1DynamicsModel = MBD.BCR4BP41DynamicsModel(systemData)

mu12::Float64 = get12MassRatio(systemData)
mu41::Float64 = get41MassRatio(systemData)
lstar12::Float64 = get12CharLength(systemData)
lstar41::Float64 = get41CharLength(systemData)
tstar12::Float64 = get12CharTime(systemData)
tstar41::Float64 = get41CharTime(systemData)

propagator = MBD.Propagator()
endEventsBCR4BP = DifferentialEquations.VectorContinuousCallback(endConditionsBCR4BP, endAffectBCR4BP!, nothing, 4)
endEventsCR3BP = DifferentialEquations.VectorContinuousCallback(endConditionsCR3BP, endAffectCR3BP!, nothing, 4)
manifoldEvent = DifferentialEquations.ContinuousCallback(periapsisCondition, nothing, endAffectManifold!)

JCEM::Float64 = 3.0663
R_H::Float64 = lstar41*(systemData.primaryData[1].mass/(3*(systemData.primaryData[3].mass+systemData.primaryData[1].mass)))^(1/3)
r_H_EM::Float64 = R_H/lstar12
r_H::Float64 = R_H/lstar41
r_E_EM::Float64 = systemData.primaryData[1].bodyRadius/lstar12
r_E::Float64 = systemData.primaryData[1].bodyRadius/lstar41
r_M_EM::Float64 = systemData.primaryData[2].bodyRadius/lstar12
r_M::Float64 = systemData.primaryData[2].bodyRadius/lstar41

radius::Float64 = 0.0035
# radius::Float64 = 0.00075
n::Int64 = 300
numAngles::Int64 = 37

thetaM::Vector{Float64} = range(0, 360, numAngles)*pi/180
x_B1::Float64 = 1-mu41
x_Moon::Float64 = 1-mu41+(1-mu12)*lstar12/lstar41
xSB1::Vector{Float64} = range(x_B1-radius, x_B1+radius, n)
# xSB1::Vector{Float64} = range(x_Moon-radius, x_Moon+radius, n)
ySB1::Vector{Float64} = range(-radius, radius, n)
posSB1::Vector{StaticArrays.SVector{2, Float64}} = vec([StaticArrays.SVector(x, y) for x in xSB1, y in ySB1])
deltarSB1::Vector{StaticArrays.SVector{2, Float64}} = posSB1 .- Ref(StaticArrays.SVector{2, Float64}([1-mu41, 0]))
rSB1::Vector{Float64} = LinearAlgebra.norm.(deltarSB1)
rhatSB1::Vector{StaticArrays.SVector{2, Float64}} = deltarSB1 ./ rSB1
thatSB1::Vector{StaticArrays.SVector{2, Float64}} = [StaticArrays.SVector(-v[2], v[1]) for v in rhatSB1] # Prograde
# thatSB1::Vector{StaticArrays.SVector{2, Float64}} = [StaticArrays.SVector(v[2], -v[1]) for v in rhatSB1] # Retrograde
qSB1::Vector{Vector{Vector{Float64}}} = Vector{Vector{Vector{Float64}}}(undef, numAngles)
qProp::Vector{Vector{Float64}} = []
map_theta::Vector{Int64} = []
map_q::Vector{Int64} = []
for t::Int64 in eachindex(thetaM)
    qGuessSB1::Vector{Vector{Float64}} = [[posSB1[j]..., 0, thatSB1[j]..., 0, thetaM[t]] for j in eachindex(posSB1)]
    qGuessEM::Vector{Vector{Float64}} = matRotating41ToRotating12(SB1DynamicsModel, qGuessSB1, zeros(length(qGuessSB1)))[1]
    OmegaEM::Vector{Float64} = map(q -> getPseudopotential(EMCR3BPDynamicsModel, q[1:3]), qGuessEM)
    v2EM::Vector{Float64} = 2 .* OmegaEM .- JCEM
    v2EM[v2EM .< 0] .= NaN
    vEM::Vector{Float64} = sqrt.(v2EM)
    @inbounds for j::Int64 in eachindex(vEM)
        vel::Vector{Float64} = qGuessEM[j][4:5]
        qGuessEM[j][4:5] = vel .* vEM[j] ./ LinearAlgebra.norm(vel)
    end
    qSB1[t] = matRotating12ToRotating41(EMDynamicsModel, qGuessEM, zeros(length(qGuessEM)))[1]
    for q::Int64 in eachindex(qSB1[t])
        if any(isnan, qSB1[t][q])
            continue
        else
            push!(qProp, qSB1[t][q])
            push!(map_theta, t)
            push!(map_q, q)
        end
    end
end
m::Int64 = length(qProp)
flagsSB1_vec::Vector{Int64} = zeros(Int64, m)
println("Propagating $m BCR4BP trajectories with $(Threads.nthreads()) threads...")
Threads.@threads for q::Int64 in 1:m
    IC::Vector{Float64} = qProp[q]
    (arc::MBD.BCR4BP41Arc, event::Symbol) = propagateWithEvents(propagator, endEventsBCR4BP, IC, [0, pi/2], SB1DynamicsModel, [SB1DynamicsModel, mu41, r_H, r_E, r_M, Periapsis(0, [])])
    if (event == :earth) || (event == :moon)
        flagsSB1_vec[q] = 8
    else
        e = String(event)
        if occursin("escape", e)
            numPeris::RegexMatch{String} = match(r"\d+$", e)
            flagsSB1_vec[q] = numPeris === nothing ? 9 : parse(Int, numPeris.match)
        else
            flagsSB1_vec[q] = 9
        end
    end
end
flagsSB1::Vector{Vector{Int64}} = Vector{Vector{Int64}}(undef, numAngles)
Threads.@threads for t::Int64 in 1:numAngles
    flagsSB1[t] = fill(7, length(qSB1[t]))
end
Threads.@threads for j::Int64 in 1:m
    flagsSB1[map_theta[j]][map_q[j]] = flagsSB1_vec[j]
end

xEM::Vector{Float64} = range(-radius, radius, n) .* lstar41 ./ lstar12
# xEM::Vector{Float64} = range(-radius, radius, n) .* lstar41 ./ lstar12) .+ (1-mu12)
yEM::Vector{Float64} = range(-radius, radius, n) .* lstar41 ./ lstar12
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
Threads.@threads for j::Int64 in eachindex(v)
    qEM[j] = append!(copy(deltar[j]), [0.0], v[j] .* that[j], [0.0])
end
qPropEM::Vector{Vector{Float64}} = []
map_qEM::Vector{Int64} = []
for q::Int64 in eachindex(qEM)
    if any(isnan, qEM[q])
        continue
    else
        push!(qPropEM, qEM[q])
        push!(map_qEM, q)
    end
end
mEM::Int64 = length(qPropEM)
flagsEM_vec::Vector{Int64} = zeros(Int64, mEM)
println("Propagating $mEM CR3BP trajectories with $(Threads.nthreads()) threads...")
Threads.@threads for q::Int64 in 1:mEM
    IC::Vector{Float64} = qPropEM[q]
    (arc::MBD.CR3BPArc, event::Symbol) = propagateWithEvents(propagator, endEventsCR3BP, IC, [0, pi/2*tstar41/tstar12], EMCR3BPDynamicsModel, [EMCR3BPDynamicsModel, r_H_EM, r_E_EM, r_M_EM, Periapsis(0, [])])
    if (event == :earth) || (event == :moon)
        flagsEM_vec[q] = 8
    else
        e = String(event)
        if occursin("escape", e)
            numPeris::RegexMatch{String} = match(r"\d+$", e)
            flagsEM_vec[q] = numPeris === nothing ? 9 : parse(Int, numPeris.match)
        else
            flagsEM_vec[q] = 9
        end
    end
end
flagsEM::Vector{Int64} = fill(7, length(qEM))
Threads.@threads for j::Int64 in 1:mEM
    flagsEM[map_qEM[j]] = flagsEM_vec[j]
end

# targeter = PlanarPerpJCTargeter(EMCR3BPDynamicsModel)
# familyFile::String = "FamilyData/CR3BPEML2Lyapunovs.csv"
# orbit::MBD.CR3BPPeriodicOrbit = interpOrbit(targeter, familyFile, "JC", JCEM; choiceIndex = 1)
# manifold::MBD.CR3BPManifold = getManifoldByArclength(orbit, "Stable", "Negative", 25/get12CharLength(systemData), n*10)
# nThreads::Int64 = Threads.nthreads()
# qMan_thread::Vector{Vector{Vector{Float64}}} = [Vector{Vector{Float64}}() for _ in 1:nThreads]
# Threads.@threads for a in eachindex(manifold.orbitTimes)
#     tid = mod1(Threads.threadid(), nThreads)
#     peris = Periapsis(0, [])
#     arc::MBD.CR3BPArc = propagateWithEvent(propagator, manifoldEvent, real(manifold.initialConditions[a]), [0, -pi/2*get41CharTime(systemData)/get12CharTime(systemData)], EMCR3BPDynamicsModel, [peris])
#     for p in eachindex(peris.states)
#         push!(qMan_thread[tid], peris.states[p])
#     end
# end
# qMan::Vector{Vector{Float64}} = reduce(vcat, qMan_thread)

# sample::Int64 = 22296
# (arc41::MBD.BCR4BP41Arc, event) = propagateWithEvents(propagator, endEventsBCR4BP, qSB1[sample], [0, pi/2], SB1DynamicsModel, [SB1DynamicsModel, get41MassRatio(systemData), R_H/get41CharLength(systemData), systemData.primaryData[1].bodyRadius/get41CharLength(systemData), systemData.primaryData[2].bodyRadius/get41CharLength(systemData), Periapsis(0, [])])
# exportBCR4BP41Trajectory(qSB1[sample]..., getTimeByIndex(arc41, -1), SB1DynamicsModel, mf, :Traj41)
# println(getTimeByIndex(arc41, -1))
# q::Vector{Float64} = rotating41ToRotating12(SB1DynamicsModel, [qSB1[sample]], [0.0])[1][1]
# exportBCR4BP12Trajectory(q..., getTimeByIndex(arc41, -1)*get41CharTime(systemData)/get12CharTime(systemData), EMDynamicsModel, mf, :Traj12)

MATLAB.put_variable(mf, :pointsSB1, qSB1)
MATLAB.put_variable(mf, :pointsEM, qEM)
# MATLAB.put_variable(mf, :manifold, qMan)
MATLAB.put_variable(mf, :flagsSB1, flagsSB1)
MATLAB.put_variable(mf, :flagsEM, flagsEM)
MATLAB.put_variable(mf, :JC, JCEM)
MATLAB.put_variable(mf, :moonAngles, thetaM)
MATLAB.close(mf)

println()
# end
