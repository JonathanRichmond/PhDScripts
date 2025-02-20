"""
Script for BCR4BP code development

Author: Jonathan Richmond
C: 2/19/25
"""
module BCR4BPDev
println()

using MBD, MATLAB

include("../Targeters/PlanarPerpJC.jl")
include("../Utilities/Export.jl")

EMSystemData = MBD.CR3BPSystemData("Earth", "Moon")
SB1SystemData = MBD.CR3BPSystemData("Sun", "Earth_Barycenter")
systemData = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
CR3BPDynamicsModel = MBD.CR3BPDynamicsModel(EMSystemData)
EMDynamicsModel = MBD.BCR4BP12DynamicsModel(systemData)
Earth = systemData.primaryData[1]
Moon = systemData.primaryData[2]
Sun = systemData.primaryData[3]
B1 = systemData.primaryData[4]

(get12CharLength(systemData) == getCharLength(EMSystemData)) && (get12CharMass(systemData) == getCharMass(EMSystemData)) && (get12CharTime(systemData) == getCharTime(EMSystemData)) && (get12MassRatio(systemData) == getMassRatio(EMSystemData)) || ArgumentError("Earth-Moon system does not match between models")
(get41CharLength(systemData) == getCharLength(SB1SystemData)) && (get41CharMass(systemData) == getCharMass(SB1SystemData)) && (get41CharTime(systemData) == getCharTime(SB1SystemData)) && (get41MassRatio(systemData) == getMassRatio(SB1SystemData)) || ArgumentError("Sun-B1 system does not match between models")

CR3BPTargeter = PlanarPerpJCTargeter(CR3BPDynamicsModel)
CR3BPOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, "FamilyData/EML1Lyapunovs.csv", "JC", 3.0)
println("\nCR3BP Orbit:\n\tState:$(CR3BPOrbit.initialCondition)\n\tPeriod: $(CR3BPOrbit.period)\n\tJC: $(getJacobiConstant(CR3BPOrbit))\n\tStability: $(getStabilityIndex(CR3BPOrbit))")

propagator = MBD.Propagator()
arcEM::MBD.BCR4BP12Arc = propagate(propagator, push!(copy(CR3BPOrbit.initialCondition), 0.0), collect(range(0, 1.5*CR3BPOrbit.period, 1001)), EMDynamicsModel)
nStates::Int64 = getStateCount(arcEM)
xEM::Vector{Float64} = zeros(Float64, nStates)
yEM::Vector{Float64} = zeros(Float64, nStates)
zEM::Vector{Float64} = zeros(Float64, nStates)
xdotEM::Vector{Float64} = zeros(Float64, nStates)
ydotEM::Vector{Float64} = zeros(Float64, nStates)
zdotEM::Vector{Float64} = zeros(Float64, nStates)
thetaSEM::Vector{Float64} = zeros(Float64, nStates)
tEM::Vector{Float64} = zeros(Float64, nStates)
for s::Int64 in 1:nStates
    state::Vector{Float64} = getStateByIndex(arcEM, s)
    xEM[s] = state[1]
    yEM[s] = state[2]
    zEM[s] = state[3]
    xdotEM[s] = state[4]
    ydotEM[s] = state[5]
    zdotEM[s] = state[6]
    thetaSEM[s] = state[7]
    tEM[s] = getTimeByIndex(arcEM, s)
end

(statesSB1::Vector{Vector{Float64}}, tSB1::Vector{Float64}) = rotating122Rotating41(EMDynamicsModel, arcEM.states, arcEM.times)
xSB1::Vector{Float64} = zeros(Float64, nStates)
ySB1::Vector{Float64} = zeros(Float64, nStates)
zSB1::Vector{Float64} = zeros(Float64, nStates)
xdotSB1::Vector{Float64} = zeros(Float64, nStates)
ydotSB1::Vector{Float64} = zeros(Float64, nStates)
zdotSB1::Vector{Float64} = zeros(Float64, nStates)
thetaMSB1::Vector{Float64} = zeros(Float64, nStates)
for s::Int64 in 1:nStates
    xSB1[s] = statesSB1[s][1]
    ySB1[s] = statesSB1[s][2]
    zSB1[s] = statesSB1[s][3]
    xdotSB1[s] = statesSB1[s][4]
    ydotSB1[s] = statesSB1[s][5]
    zdotSB1[s] = statesSB1[s][6]
    thetaMSB1[s] = statesSB1[s][7]
end

mf = MATLAB.MatFile("Output/BCR4BPDev.mat", "w")
exportCR3BPOrbit(CR3BPOrbit, CR3BPDynamicsModel, mf, :orbitCR3BP)
exportBCR4BP12Trajectory(xEM, yEM, zEM, xdotEM, ydotEM, zdotEM, thetaSEM, tEM, mf, :trajBCR4BP12)
exportBCR4BP41Trajectory(xSB1, ySB1, zSB1, xdotSB1, ydotSB1, zdotSB1, thetaMSB1, tSB1, mf, :trajBCR4BP41)
MATLAB.close(mf)

println()
end
