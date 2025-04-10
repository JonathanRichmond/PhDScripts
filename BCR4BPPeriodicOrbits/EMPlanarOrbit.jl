"""
Script for BCR4BP Earth-Moon planar orbits

Author: Jonathan Richmond
C: 4/9/25
"""
module EMPlanar
println()

using MBD, MATLAB

include("../BCR4BPTargeters/PlanarPerpP.jl")
include("../Utilities/Export.jl")

function homotopy(systemData::MBD.BCR4BPSystemData, targeter::PlanarPerpP12Targeter, initialStateGuess::Vector{Float64}, targetP::Float64, tol::Float64 = 1E-11)
    SunGravParam::Float64 = copy(systemData.primaryData[3].gravParam)
    systemData.primaryData[3].gravParam = 0.0
    solution::MBD.BCR4BP12MultipleShooterProblem = correct(targeter, initialStateGuess, targetP, tol)
    println("Converged CR3BP Orbit:\n\tState:$(solution.nodes[1].state.data[1:7])\n\tPeriod: $(getPeriod(targeter, solution))")
    
    epsilon = 0.0
    while epsilon < 1
        epsilon += 0.001
        systemData.primaryData[3].gravParam = epsilon*SunGravParam
        initialStateGuess = solution.nodes[1].state.data[1:7]
        solution = correct(targeter, initialStateGuess, targetP, tol)
        # println("\nConverged Homotopy Orbit:\n\tState:$(solution.nodes[1].state.data[1:7])")
    end

    return solution
end

systemData = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
dynamicsModel = MBD.BCR4BP12DynamicsModel(systemData)
CR3BPSystemData = MBD.CR3BPSystemData("Earth", "Moon")
CR3BPDynamicsModel = MBD.CR3BPDynamicsModel(CR3BPSystemData)

Earth = systemData.primaryData[1]
Moon = systemData.primaryData[2]

targeter = PlanarPerpP12Targeter(dynamicsModel)
propagator = MBD.Propagator()

initialStateGuess::Vector{Float64} = append!(getEquilibriumPoint(CR3BPDynamicsModel, 1), [0.0, 0.0, 0.0, 0.0])
targetP::Float64 = getSynodicPeriod(dynamicsModel)
solution::MBD.BCR4BP12MultipleShooterProblem = homotopy(systemData, targeter, initialStateGuess, targetP, 1E-8)
refinedSolution::MBD.BCR4BP12MultipleShooterProblem = correct(targeter, solution.nodes[1].state.data[1:7], targetP, 1E-10)
println("\nConverged BCR4BP Orbit:\n\tState:$(refinedSolution.nodes[1].state.data[1:7])\n\tPeriod: $(getPeriod(targeter, refinedSolution))")

arc::MBD.BCR4BP12Arc = propagate(propagator, refinedSolution.nodes[1].state.data[1:7], collect(range(0, targetP, 1001)), dynamicsModel)
nStates::Int64 = getStateCount(arc)
x::Vector{Float64} = Vector{Float64}(undef, nStates)
y::Vector{Float64} = Vector{Float64}(undef, nStates)
z::Vector{Float64} = Vector{Float64}(undef, nStates)
xdot::Vector{Float64} = Vector{Float64}(undef, nStates)
ydot::Vector{Float64} = Vector{Float64}(undef, nStates)
zdot::Vector{Float64} = Vector{Float64}(undef, nStates)
thetaS::Vector{Float64} = Vector{Float64}(undef, nStates)
t::Vector{Float64} = Vector{Float64}(undef, nStates)
for s::Int64 in 1:nStates
    state::Vector{Float64} = getStateByIndex(arc, s)
    x[s] = state[1]
    y[s] = state[2]
    z[s] = state[3]
    xdot[s] = state[4]
    ydot[s] = state[5]
    zdot[s] = state[6]
    thetaS[s] = state[7]
    t[s] = getTimeByIndex(arc, s)
end

mf = MATLAB.MatFile("Output/EMPlanarOrbit.mat", "w")
exportBCR4BP12Trajectory(x, y, z, xdot, ydot, zdot, thetaS, t, mf, :L1Orbit)
MATLAB.close(mf)

println()
end
