"""
Script for BCR4BP Earth-Moon planar orbits

Author: Jonathan Richmond
C: 4/23/25
U: 6/5/25
"""
module EMPlanar
println()

using MBD, LinearAlgebra, MATLAB

include("../BCR4BPTargeters/PlanarPerpP.jl")
include("../CR3BPTargeters/PlanarPerpJC.jl")
include("../Utilities/Export.jl")

function homotopy(systemData::MBD.BCR4BPSystemData, targeter::PlanarPerpP12Targeter, initialStateGuess::Vector{Float64}, synRevs::Int64; Deltaeps::Float64 = 0.001, tol::Float64 = 1E-11, JTol::Float64 = 2E-3)
    SunGravParam::Float64 = copy(systemData.primaryData[3].gravParam)
    systemData.primaryData[3].gravParam = 0.0
    targetP::Float64 = synRevs*getSynodicPeriod(targeter.dynamicsModel)
    solution::MBD.BCR4BP12MultipleShooterProblem = correct(targeter, initialStateGuess, targetP, synRevs; tol = tol, JTol = JTol)
    println("\nConverged CR3BP Orbit:\n\tState:$(solution.nodes[1].state.data[1:7])\n\tPeriod: $(getPeriod(targeter, solution))\n")
    
    epsilon = 0.0
    while epsilon < 1
        epsilon += Deltaeps
        systemData.primaryData[3].gravParam = epsilon*SunGravParam
        targetP = getSynodicPeriod(targeter.dynamicsModel)
        initialStateGuess = solution.nodes[1].state.data[1:7]
        solution = correct(targeter, initialStateGuess, targetP, synRevs; tol = tol, JTol = JTol)
        # println("Converged Homotopy Orbit (epsilon = $epsilon):\n\tState:$(solution.nodes[1].state.data[1:7])\n\tPeriod: $(getPeriod(targeter, solution))\n")
    end

    return solution
end

systemData = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
dynamicsModel = MBD.BCR4BP12DynamicsModel(systemData)
CR3BPSystemData = MBD.CR3BPSystemData("Earth", "Moon")
CR3BPDynamicsModel = MBD.CR3BPDynamicsModel(CR3BPSystemData)

Earth = systemData.primaryData[1]
Moon = systemData.primaryData[2]

propagator = MBD.Propagator()
targeter = PlanarPerpP12Targeter(dynamicsModel)
CR3BPTargeter = PlanarPerpJCTargeter(CR3BPDynamicsModel)

targetP::Float64 = getSynodicPeriod(dynamicsModel)
guessOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, "FamilyData/CR3BPEML1Lyapunovs.csv", "Period", 2*pi)
compOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, "FamilyData/CR3BPEML1Lyapunovs.csv", "Period", targetP)
initialStateGuess::Vector{Float64} = push!(guessOrbit.initialCondition, 0.0)
solution::MBD.BCR4BP12MultipleShooterProblem = homotopy(systemData, targeter, initialStateGuess, 1)
refinedSolution::MBD.BCR4BP12MultipleShooterProblem = correct(targeter, solution.nodes[1].state.data[1:7], targetP, 1)
println("Converged BCR4BP Orbit:\n\tState:$(refinedSolution.nodes[1].state.data[1:7])\n\tPeriod: $(getPeriod(targeter, refinedSolution))\n")
M::Matrix{Float64} = getMonodromy(targeter, refinedSolution)
E::LinearAlgebra.Eigen = LinearAlgebra.eigen(M)
println("Eigenvalues: $(E.values)\n")

guessOrbit2::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, "FamilyData/CR3BPEML1Lyapunovs.csv", "Period", 1*pi)
compOrbit2::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, "FamilyData/CR3BPEML1Lyapunovs.csv", "Period", targetP/2)
initialStateGuess2::Vector{Float64} = push!(guessOrbit2.initialCondition, 0.0)
solution2::MBD.BCR4BP12MultipleShooterProblem = homotopy(systemData, targeter, initialStateGuess2, 1)
refinedSolution2::MBD.BCR4BP12MultipleShooterProblem = correct(targeter, solution2.nodes[1].state.data[1:7], targetP, 1)
println("Converged BCR4BP Orbit:\n\tState:$(refinedSolution2.nodes[1].state.data[1:7])\n\tPeriod: $(getPeriod(targeter, refinedSolution2))\n")

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

arc2::MBD.BCR4BP12Arc = propagate(propagator, refinedSolution2.nodes[1].state.data[1:7], collect(range(0, targetP, 1001)), dynamicsModel)
nStates2::Int64 = getStateCount(arc2)
x2::Vector{Float64} = Vector{Float64}(undef, nStates2)
y2::Vector{Float64} = Vector{Float64}(undef, nStates2)
z2::Vector{Float64} = Vector{Float64}(undef, nStates2)
xdot2::Vector{Float64} = Vector{Float64}(undef, nStates2)
ydot2::Vector{Float64} = Vector{Float64}(undef, nStates2)
zdot2::Vector{Float64} = Vector{Float64}(undef, nStates2)
thetaS2::Vector{Float64} = Vector{Float64}(undef, nStates2)
t2::Vector{Float64} = Vector{Float64}(undef, nStates2)
for s::Int64 in 1:nStates2
    state::Vector{Float64} = getStateByIndex(arc2, s)
    x2[s] = state[1]
    y2[s] = state[2]
    z2[s] = state[3]
    xdot2[s] = state[4]
    ydot2[s] = state[5]
    zdot2[s] = state[6]
    thetaS2[s] = state[7]
    t2[s] = getTimeByIndex(arc2, s)
end

mf = MATLAB.MatFile("Output/EMPlanarOrbit.mat", "w")
exportBCR4BP12Trajectory(x, y, z, xdot, ydot, zdot, thetaS, t, mf, :BCR4BPOrbit)
exportCR3BPOrbit(guessOrbit, mf, :CR3BPGuessOrbit)
exportCR3BPOrbit(compOrbit, mf, :CR3BPCompOrbit)
exportBCR4BP12Trajectory(x2, y2, z2, xdot2, ydot2, zdot2, thetaS2, t2, mf, :BCR4BPOrbit2)
exportCR3BPOrbit(guessOrbit2, mf, :CR3BPGuessOrbit2)
exportCR3BPOrbit(compOrbit2, mf, :CR3BPCompOrbit2)
MATLAB.close(mf)

println()
end
