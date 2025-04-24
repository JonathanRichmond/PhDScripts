"""
Script for BCR4BP Earth-Moon planar orbits

Author: Jonathan Richmond
C: 4/23/25
"""
module EMPlanar
println()

using MBD, LinearAlgebra, MATLAB

include("../BCR4BPTargeters/PlanarPerpP.jl")
include("../Utilities/Export.jl")

function homotopy(systemData::MBD.BCR4BPSystemData, targeter::PlanarPerpP12Targeter, initialStateGuess::Vector{Float64}, revs::Int64; Deltaeps::Float64 = 0.001, tol::Float64 = 1E-11, JTol::Float64 = 2E-3)
    SunGravParam::Float64 = copy(systemData.primaryData[3].gravParam)
    systemData.primaryData[3].gravParam = 0.0
    targetP::Float64 = revs*getSynodicPeriod(targeter.dynamicsModel)
    solution::MBD.BCR4BP12MultipleShooterProblem = correct(targeter, initialStateGuess, targetP, revs; tol = tol, JTol = JTol)
    println("Converged CR3BP Orbit:\n\tState:$(solution.nodes[1].state.data[1:7])\n\tPeriod: $(getPeriod(targeter, solution))")
    
    epsilon = 0.0
    while epsilon < 1
        epsilon += Deltaeps
        systemData.primaryData[3].gravParam = epsilon*SunGravParam
        targetP = revs*getSynodicPeriod(targeter.dynamicsModel)
        initialStateGuess = solution.nodes[1].state.data[1:7]
        solution = correct(targeter, initialStateGuess, targetP, revs; tol = tol, JTol = JTol)
        # println("\nConverged Homotopy Orbit (epsilon = $epsilon):\n\tState:$(solution.nodes[1].state.data[1:7])\n\tPeriod: $(getPeriod(targeter, solution))")
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
solution::MBD.BCR4BP12MultipleShooterProblem = homotopy(systemData, targeter, initialStateGuess, 1; tol = 1E-9)
refinedSolution::MBD.BCR4BP12MultipleShooterProblem = correct(targeter, solution.nodes[1].state.data[1:7], targetP, 1; tol = 1E-10)
# refinedSolution::MBD.BCR4BP12MultipleShooterProblem = correct(targeter, [1.005405631673800, 0, 0.181390966578526, 0, -0.090738307446220, 0, 0], 2*targetP, 2, 5E-8)
println("\nConverged BCR4BP Orbit:\n\tState:$(refinedSolution.nodes[1].state.data[1:7])\n\tPeriod: $(getPeriod(targeter, refinedSolution))")
M::Matrix{Float64} = getMonodromy(targeter, refinedSolution)
E::LinearAlgebra.Eigen = LinearAlgebra.eigen(M)
println("\nEigenvalues: $(E.values)")

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
