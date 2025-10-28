"""
Script for computing periapsis maps

Author: Jonathan Richmond
C: 10/9/25
"""
# module PeriMap
println("Running PeriapsisMap.jl...\n")

using MBD, DifferentialEquations, LinearAlgebra, Logging, MATLAB, StaticArrays

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../Utilities/Export.jl")

mutable struct Periapsis
    count::Int64
end

function endConditions(out::SubArray{Float64}, state::Vector{Float64}, time::Float64, integrator)
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

function endAffect!(integrator, index)
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

mf = MATLAB.MatFile("Output/PeriapsisMap.mat", "w")

EMCR3BPSystemData = MBD.CR3BPSystemData("Earth", "Moon")
SB1CR3BPSystemData = MBD.CR3BPSystemData("Sun", "Earth_Barycenter")
systemData = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
EMCR3BPDynamicsModel = MBD.CR3BPDynamicsModel(EMCR3BPSystemData)
SB1CR3BPDynamicsModel = MBD.CR3BPDynamicsModel(SB1CR3BPSystemData)
EMDynamicsModel = MBD.BCR4BP12DynamicsModel(systemData)
SB1DynamicsModel = MBD.BCR4BP41DynamicsModel(systemData)

propagator = MBD.Propagator()
endEvents = DifferentialEquations.VectorContinuousCallback(endConditions, endAffect!, nothing, 4)

x_B1::Float64 = 1-get41MassRatio(systemData)
xSB1::Vector{Float64} = range(x_B1-0.004, x_B1+0.004, 500)
ySB1::Vector{Float64} = range(-0.004, 0.004, 500)
JCEM::Float64 = 3.06633
thetaM::Float64 = pi*0
r_H::Float64 = (systemData.primaryData[1].mass/(3*(systemData.primaryData[3].mass+systemData.primaryData[1].mass)))^(1/3)

posSB1::Vector{StaticArrays.SVector{2, Float64}} = vec([StaticArrays.SVector(x, y) for x in xSB1, y in ySB1])
deltar::Vector{StaticArrays.SVector{2, Float64}} = posSB1 .- Ref(StaticArrays.SVector{2, Float64}([1-get41MassRatio(systemData), 0]))
r::Vector{Float64} = LinearAlgebra.norm.(deltar)
rhat::Vector{StaticArrays.SVector{2, Float64}} = deltar ./ r
that::Vector{StaticArrays.SVector{2, Float64}} = [StaticArrays.SVector(-v[2], v[1]) for v in rhat]
qGuessSB1::Vector{Vector{Float64}} = [[posSB1[j]..., 0, that[j]..., 0, thetaM] for j in eachindex(posSB1)]
qGuessEM::Vector{Vector{Float64}} = rotating41ToRotating12(SB1DynamicsModel, qGuessSB1, zeros(length(qGuessSB1)))[1]
OmegaEM::Vector{Float64} = map(q -> getPseudopotential(EMCR3BPDynamicsModel, q[1:3]), qGuessEM)
v2EM::Vector{Float64} = 2 .* OmegaEM .- JCEM
v2EM[v2EM .< 0] .= NaN
vEM::Vector{Float64} = sqrt.(v2EM)
Threads.@threads for j in eachindex(vEM)
    qGuessEM[j][4:5] = qGuessEM[j][4:5] .* vEM[j] ./ norm(qGuessEM[j][4:5])
end
qSB1::Vector{Vector{Float64}} = rotating12ToRotating41(EMDynamicsModel, qGuessEM, zeros(length(qGuessEM)))[1]

flags::Vector{Int64} = Vector{Int64}(undef, length(qSB1))
Threads.@threads for q in eachindex(qSB1)
    if any(x -> isnan(x), qSB1[q])
        flags[q] = 1
    else
        (arc::MBD.BCR4BP41Arc, event) = propagateWithEvents(propagator, endEvents, qSB1[q], [0.0, 4*pi], SB1DynamicsModel, [SB1DynamicsModel, get41MassRatio(systemData), r_H, systemData.primaryData[1].bodyRadius/get41CharLength(systemData), systemData.primaryData[2].bodyRadius/get41CharLength(systemData), Periapsis(0)])
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
            println(event)
            flags[q] = 0
        end
    end
end

# sample::Int64 = 18285
# (arc41::MBD.BCR4BP41Arc, event) = propagateWithEvents(propagator, endEvents, qSB1[sample], [0.0, 10*pi], SB1DynamicsModel, [SB1DynamicsModel, get41MassRatio(systemData), r_H, systemData.primaryData[1].bodyRadius/get41CharLength(systemData), systemData.primaryData[2].bodyRadius/get41CharLength(systemData), Periapsis(0)])
# exportBCR4BP41Trajectory(qSB1[sample]..., getTimeByIndex(arc41, -1), SB1DynamicsModel, mf, :Traj41)
# qEM::Vector{Float64} = rotating41ToRotating12(SB1DynamicsModel, [qSB1[sample]], [0.0])[1][1]
# exportBCR4BP12Trajectory(qEM..., getTimeByIndex(arc41, -1)*get41CharTime(systemData)/get12CharTime(systemData), EMDynamicsModel, mf, :Traj12)

MATLAB.put_variable(mf, :points, qSB1)
MATLAB.put_variable(mf, :flags, flags)
MATLAB.put_variable(mf, :JC, JCEM)
MATLAB.put_variable(mf, :moonAngle, thetaM)
MATLAB.close(mf)

println()
# end
