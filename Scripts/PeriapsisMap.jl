"""
Script for computing periapsis maps

Author: Jonathan Richmond
C: 10/9/25
U: 11/12/25
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
        q_E::Vector{Float64} = getPrimaryState(integrator.p[3], 1, state[7])
        q_M::Vector{Float64} = getPrimaryState(integrator.p[3], 2, state[7])
    
        out[1] = LinearAlgebra.dot(state[1:3], state[4:6])
        out[2] = LinearAlgebra.norm(state[1:3])-integrator.p[4]
        out[3] = integrator.p[5]-LinearAlgebra.norm(state[1:3]-q_E[1:3])
        out[4] = integrator.p[6]-LinearAlgebra.norm(state[1:3]-q_M[1:3])
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

# function periapsisCondition(state::Vector{Float64}, time::Float64, integrator)
#     if abs(time) < 1E-10
#         -1.0
#     else
#         LinearAlgebra.dot(state[1:3], state[4:6])
#     end
# end

function endAffectBCR4BP!(integrator, index)
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

# function endAffectManifold!(integrator)
#     integrator.p[2].count += 1
#     push!(integrator.p[2].states, copy(integrator.u))
#     if integrator.p[2].count >= 6
#         DifferentialEquations.terminate!(integrator)
#     end
# end

mf = MATLAB.MatFile("Output/PeriapsisMap.mat", "w")

CR3BPSystemData = MBD.CR3BPSystemData("Earth", "Moon")
BCR4BPSystemData = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
CR3BPDynamicsModel = MBD.CR3BPDynamicsModel(CR3BPSystemData)
BCR4BPDynamicsModel = MBD.BCR4BP12DynamicsModel(BCR4BPSystemData)

mu12::Float64 = get12MassRatio(BCR4BPSystemData)
lstar12::Float64 = get12CharLength(BCR4BPSystemData)
lstar41::Float64 = get41CharLength(BCR4BPSystemData)
tstar12::Float64 = get12CharTime(BCR4BPSystemData)
r_E::Float64 = BCR4BPSystemData.primaryData[1].bodyRadius/lstar12
r_M::Float64 = BCR4BPSystemData.primaryData[2].bodyRadius/lstar12

propagator = MBD.Propagator()
endEventsBCR4BP = DifferentialEquations.VectorContinuousCallback(endConditionsBCR4BP, endAffectBCR4BP!, nothing, 4)
endEventsCR3BP = DifferentialEquations.VectorContinuousCallback(endConditionsCR3BP, endAffectCR3BP!, nothing, 4)
# manifoldEvent = DifferentialEquations.ContinuousCallback(periapsisCondition, nothing, endAffectManifold!)

JC::Float64 = 3.0663
R_H::Float64 = lstar41*(BCR4BPSystemData.primaryData[1].mass/(3*(BCR4BPSystemData.primaryData[3].mass+BCR4BPSystemData.primaryData[1].mass)))^(1/3)
r_H::Float64 = R_H/lstar12
radius::Float64 = 1.25
# radius::Float64 = 0.00075
n::Int64 = 500 #500
numAngles::Int64 = 121 #73
thetaS::Vector{Float64} = collect(range(0, 360, numAngles))*pi/180

x_M::Float64 = 1-mu12
xGrid::Vector{Float64} = collect(range(-radius, radius, n))
# xGrid::Vector{Float64} = collect(range(x_M-radius, x_M+radius, n))
yGrid::Vector{Float64} = collect(range(-radius, radius, n))
rGrid::Vector{Vector{Float64}} = [[x, y] for x in xGrid for y in yGrid]
dGrid::Vector{Float64} = LinearAlgebra.norm.(rGrid)
rhatGrid::Vector{Vector{Float64}} = [v ./ d for (v, d) in zip(rGrid, dGrid)]
thatGrid::Vector{Vector{Float64}} = [StaticArrays.SVector(-v[2], v[1]) for v in rhatGrid] # Prograde
# thatGrid::Vector{Vector{Float64}} = [StaticArrays.SVector(v[2], -v[1]) for v in rhatGrid] # Retrograde
OmegaGrid::Vector{Float64} = map(q -> getPseudopotential(CR3BPDynamicsModel, push!(copy(q), 0.0)), rGrid)
v2Grid::Vector{Float64} = 2 .* OmegaGrid .- JC
v2Grid[v2Grid .< 0] .= NaN
vGrid::Vector{Float64} = sqrt.(v2Grid)
qGrid::Vector{StaticArrays.SVector{6, Float64}} = Vector{StaticArrays.SVector{6, Float64}}(undef, n^2)
xPoint::Vector{Float64} = Vector{Float64}(undef, n^2)
yPoint::Vector{Float64} = Vector{Float64}(undef, n^2)
Threads.@threads for j::Int64 in 1:n^2
    qGrid[j] = StaticArrays.SVector{6, Float64}(rGrid[j]..., 0.0, vGrid[j] .* thatGrid[j]..., 0.0)
    xPoint[j] = rGrid[j][1]
    yPoint[j] = rGrid[j][2]
end
qProp::Vector{StaticArrays.SVector{6, Float64}} = StaticArrays.SVector{6, Float64}[]
qMap::Vector{Int64} = []
for q::Int64 in eachindex(qGrid)
    if any(isnan, qGrid[q])
        continue
    else
        push!(qProp, qGrid[q])
        push!(qMap, q)
    end
end
p::Int64 = length(qProp)

flagsBCR4BP_vec::Vector{Int64} = zeros(Int64, p*numAngles)
println("Propagating $(p*numAngles) BCR4BP trajectories with $(Threads.nthreads()) threads...")
for t::Int64 in 1:numAngles
    Threads.@threads for j::Int64 in 1:p
        IC::Vector{Float64} = push!(collect(qProp[j]), thetaS[t])
        (arc::MBD.BCR4BP12Arc, event::Symbol) = propagateWithEvents(propagator, endEventsBCR4BP, IC, [0, pi*6.0], BCR4BPDynamicsModel, [BCR4BPDynamicsModel, r_H, r_E, r_M, Periapsis(0, [])])
        index::Int64 = (t-1)*p + j
        if (event == :earth) || (event == :moon)
            flagsBCR4BP_vec[index] = 7
        else
            e = String(event)
            if occursin("escape", e)
                numPeris::RegexMatch{String} = match(r"\d+$", e)
                tempFlag::Int64 = numPeris === nothing ? 6 : parse(Int64, numPeris.match)
                flagsBCR4BP_vec[index] = (tempFlag > 6) ? 6 : tempFlag
            else
                flagsBCR4BP_vec[index] = 6
            end
        end
    end
end
flagsBCR4BP::Vector{Vector{Int64}} = Vector{Vector{Int64}}(undef, numAngles)
Threads.@threads for t::Int64 in 1:numAngles
    flagsBCR4BP[t] = fill(8, n^2)
end
Threads.@threads for j::Int64 in 1:(p*numAngles)
    flagsBCR4BP[div(j-1, p)+1][qMap[mod1(j, p)]] = flagsBCR4BP_vec[j]
end

flagsCR3BP_vec::Vector{Int64} = zeros(Int64, p)
println("Propagating $p CR3BP trajectories with $(Threads.nthreads()) threads...")
Threads.@threads for j::Int64 in 1:p
    IC::Vector{Float64} = collect(qProp[j])
    (arc::MBD.CR3BPArc, event::Symbol) = propagateWithEvents(propagator, endEventsCR3BP, IC, [0, pi*6.0], CR3BPDynamicsModel, [CR3BPDynamicsModel, r_H, r_E, r_M, Periapsis(0, [])])
    if (event == :earth) || (event == :moon)
        flagsCR3BP_vec[j] = 7
    else
        e = String(event)
        if occursin("escape", e)
            numPeris::RegexMatch{String} = match(r"\d+$", e)
            tempFlag::Int64 = numPeris === nothing ? 6 : parse(Int64, numPeris.match)
            flagsCR3BP_vec[j] = (tempFlag > 6) ? 6 : tempFlag
        else
            flagsCR3BP_vec[j] = 6
        end
    end
end
flagsCR3BP::Vector{Int64} = fill(8, n^2)
Threads.@threads for j::Int64 in 1:p
    flagsCR3BP[qMap[j]] = flagsCR3BP_vec[j]
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

MATLAB.put_variable(mf, :xPoints, xPoint)
MATLAB.put_variable(mf, :yPoints, yPoint)
MATLAB.put_variable(mf, :flagsBCR4BP, flagsBCR4BP)
MATLAB.put_variable(mf, :flagsCR3BP, flagsCR3BP)
# MATLAB.put_variable(mf, :manifold, qMan)
MATLAB.put_variable(mf, :JC, JC)
MATLAB.put_variable(mf, :thetaS, thetaS)
MATLAB.close(mf)

println()
# end
