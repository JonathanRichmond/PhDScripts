"""
Script for BCR4BP code development

Author: Jonathan Richmond
C: 2/19/25
U: 4/9/25
"""
module BCR4BPDev
println()

using MBD, MATLAB, SPICE

include("../BCR4BPTargeters/PlanarPerpP.jl")
include("../CR3BPTargeters/PlanarPerpJC.jl")
include("../Utilities/Export.jl")

function homotopy(systemData::MBD.BCR4BPSystemData, targeter::PlanarPerpP12Targeter, initialStateGuess::Vector{Float64}, targetP::Float64, tol::Float64 = 1E-11)
    SunGravParam::Float64 = copy(systemData.primaryData[3].gravParam)
    systemData.primaryData[3].gravParam = 0.0
    solution::MBD.BCR4BP12MultipleShooterProblem = correct(targeter, initialStateGuess, targetP, tol)
    println("\nConverged CR3BP L1 Orbit:\n\tState:$(solution.nodes[1].state.data[1:7])\n\tPeriod: $(getPeriod(targeter, solution))")
    
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

SPICE.furnsh("SPICEKernels/naif0012.tls", "SPICEKernels/de430.bsp", "SPICEKernels/de440.bsp")

EMSystemData = MBD.CR3BPSystemData("Earth", "Moon")
SB1SystemData = MBD.CR3BPSystemData("Sun", "Earth_Barycenter")
systemData = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
CR3BPDynamicsModel = MBD.CR3BPDynamicsModel(EMSystemData)
EMDynamicsModel = MBD.BCR4BP12DynamicsModel(systemData)
SB1DynamicsModel = MBD.BCR4BP41DynamicsModel(systemData)
Earth = systemData.primaryData[1]
Moon = systemData.primaryData[2]
Sun = systemData.primaryData[3]
B1 = systemData.primaryData[4]

(get12CharLength(systemData) == getCharLength(EMSystemData)) && (get12CharMass(systemData) == getCharMass(EMSystemData)) && (get12CharTime(systemData) == getCharTime(EMSystemData)) && (get12MassRatio(systemData) == getMassRatio(EMSystemData)) || ArgumentError("Earth-Moon system does not match between models")
(get41CharLength(systemData) == getCharLength(SB1SystemData)) && (get41CharMass(systemData) == getCharMass(SB1SystemData)) && (get41CharTime(systemData) == getCharTime(SB1SystemData)) && (get41MassRatio(systemData) == getMassRatio(SB1SystemData)) || ArgumentError("Sun-B1 system does not match between models")

targeter = PlanarPerpP12Targeter(EMDynamicsModel)
propagator = MBD.Propagator()
CR3BPTargeter = PlanarPerpJCTargeter(CR3BPDynamicsModel)
CR3BPOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(CR3BPTargeter, "FamilyData/EML1Lyapunovs.csv", "JC", 3.175)
println("\nCR3BP Orbit:\n\tState:$(CR3BPOrbit.initialCondition)\n\tPeriod: $(CR3BPOrbit.period)\n\tJC: $(getJacobiConstant(CR3BPOrbit))\n\tStability: $(getStabilityIndex(CR3BPOrbit))")
orbitArc::MBD.CR3BPArc = propagate(propagator, CR3BPOrbit.initialCondition, collect(range(0, 1.5*CR3BPOrbit.period, 1001)), CR3BPDynamicsModel)

arcEM::MBD.BCR4BP12Arc = propagate(propagator, push!(copy(CR3BPOrbit.initialCondition), pi*0/2), collect(range(0, 6.791174129196576, 1001)), EMDynamicsModel)
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

(statesSB1::Vector{Vector{Float64}}, tSB1::Vector{Float64}) = rotating12ToRotating41(EMDynamicsModel, arcEM.states, arcEM.times)
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

arcSB1::MBD.BCR4BP41Arc = propagate(propagator, statesSB1[1], [0, 6.791174129196576*get12CharTime(systemData)/get41CharTime(systemData)], SB1DynamicsModel)
nStates_v::Int64 = getStateCount(arcSB1)
xSB1_v::Vector{Float64} = zeros(Float64, nStates_v)
ySB1_v::Vector{Float64} = zeros(Float64, nStates_v)
zSB1_v::Vector{Float64} = zeros(Float64, nStates_v)
xdotSB1_v::Vector{Float64} = zeros(Float64, nStates_v)
ydotSB1_v::Vector{Float64} = zeros(Float64, nStates_v)
zdotSB1_v::Vector{Float64} = zeros(Float64, nStates_v)
thetaMSB1_v::Vector{Float64} = zeros(Float64, nStates_v)
tSB1_v::Vector{Float64} = zeros(Float64, nStates_v)
for s::Int64 in 1:nStates_v
    state::Vector{Float64} = getStateByIndex(arcSB1, s)
    xSB1_v[s] = state[1]
    ySB1_v[s] = state[2]
    zSB1_v[s] = state[3]
    xdotSB1_v[s] = state[4]
    ydotSB1_v[s] = state[5]
    zdotSB1_v[s] = state[6]
    thetaMSB1_v[s] = state[7]
    tSB1_v[s] = getTimeByIndex(arcSB1, s)
end

(statesEM_v::Vector{Vector{Float64}}, tEM_v::Vector{Float64}) = rotating41ToRotating12(SB1DynamicsModel, arcSB1.states, arcSB1.times)
xEM_v::Vector{Float64} = zeros(Float64, nStates_v)
yEM_v::Vector{Float64} = zeros(Float64, nStates_v)
zEM_v::Vector{Float64} = zeros(Float64, nStates_v)
xdotEM_v::Vector{Float64} = zeros(Float64, nStates_v)
ydotEM_v::Vector{Float64} = zeros(Float64, nStates_v)
zdotEM_v::Vector{Float64} = zeros(Float64, nStates_v)
thetaSEM_v::Vector{Float64} = zeros(Float64, nStates_v)
for s::Int64 in 1:nStates_v
    xEM_v[s] = statesEM_v[s][1]
    yEM_v[s] = statesEM_v[s][2]
    zEM_v[s] = statesEM_v[s][3]
    xdotEM_v[s] = statesEM_v[s][4]
    ydotEM_v[s] = statesEM_v[s][5]
    zdotEM_v[s] = statesEM_v[s][6]
    thetaSEM_v[s] = statesEM_v[s][7]
end

(statesEEclipJ2000::Vector{Vector{Float64}}, epochTimes::Vector{Float64}) = rotating12ToPrimaryEcliptic(EMDynamicsModel, "ECLIPJ2000", 1, "Jan 1 2030", arcEM.states, arcEM.times)
initialEpoch = SPICE.et2utc(epochTimes[1], :C, 0)
xEEclipJ2000::Vector{Float64} = zeros(Float64, nStates)
yEEclipJ2000::Vector{Float64} = zeros(Float64, nStates)
zEEclipJ2000::Vector{Float64} = zeros(Float64, nStates)
xdotEEclipJ2000::Vector{Float64} = zeros(Float64, nStates)
ydotEEclipJ2000::Vector{Float64} = zeros(Float64, nStates)
zdotEEclipJ2000::Vector{Float64} = zeros(Float64, nStates)
for s::Int64 in 1:nStates
    xEEclipJ2000[s] = statesEEclipJ2000[s][1]
    yEEclipJ2000[s] = statesEEclipJ2000[s][2]
    zEEclipJ2000[s] = statesEEclipJ2000[s][3]
    xdotEEclipJ2000[s] = statesEEclipJ2000[s][4]
    ydotEEclipJ2000[s] = statesEEclipJ2000[s][5]
    zdotEEclipJ2000[s] = statesEEclipJ2000[s][6]
end
println("\nInitial epoch: $initialEpoch")

(statesEM_v2::Vector{Vector{Float64}}, timesEM_v2::Vector{Float64}) = primaryEclipticToRotating12(EMDynamicsModel, "ECLIPJ2000", 1, statesEEclipJ2000, epochTimes)
xEM_v2::Vector{Float64} = zeros(Float64, nStates)
yEM_v2::Vector{Float64} = zeros(Float64, nStates)
zEM_v2::Vector{Float64} = zeros(Float64, nStates)
xdotEM_v2::Vector{Float64} = zeros(Float64, nStates)
ydotEM_v2::Vector{Float64} = zeros(Float64, nStates)
zdotEM_v2::Vector{Float64} = zeros(Float64, nStates)
thetaSEM_v2::Vector{Float64} = zeros(Float64, nStates)
for s::Int64 in 1:nStates
    xEM_v2[s] = statesEM_v2[s][1]
    yEM_v2[s] = statesEM_v2[s][2]
    zEM_v2[s] = statesEM_v2[s][3]
    xdotEM_v2[s] = statesEM_v2[s][4]
    ydotEM_v2[s] = statesEM_v2[s][5]
    zdotEM_v2[s] = statesEM_v2[s][6]
    thetaSEM_v2[s] = statesEM_v2[s][7]
end
println("Initial Sun Angle: $(thetaSEM_v2[1])")

(statesSB1EEclipJ2000::Vector{Vector{Float64}}, epochTimesSB1::Vector{Float64}) = rotating41ToPrimaryEcliptic(SB1DynamicsModel, "ECLIPJ2000", 1, "Jan 1 2030", arcSB1.states, arcSB1.times)
initialEpochSB1 = SPICE.et2utc(epochTimesSB1[1], :C, 0)
xSB1EEclipJ2000::Vector{Float64} = zeros(Float64, nStates_v)
ySB1EEclipJ2000::Vector{Float64} = zeros(Float64, nStates_v)
zSB1EEclipJ2000::Vector{Float64} = zeros(Float64, nStates_v)
xdotSB1EEclipJ2000::Vector{Float64} = zeros(Float64, nStates_v)
ydotSB1EEclipJ2000::Vector{Float64} = zeros(Float64, nStates_v)
zdotSB1EEclipJ2000::Vector{Float64} = zeros(Float64, nStates_v)
for s::Int64 in 1:nStates_v
    xSB1EEclipJ2000[s] = statesSB1EEclipJ2000[s][1]
    ySB1EEclipJ2000[s] = statesSB1EEclipJ2000[s][2]
    zSB1EEclipJ2000[s] = statesSB1EEclipJ2000[s][3]
    xdotSB1EEclipJ2000[s] = statesSB1EEclipJ2000[s][4]
    ydotSB1EEclipJ2000[s] = statesSB1EEclipJ2000[s][5]
    zdotSB1EEclipJ2000[s] = statesSB1EEclipJ2000[s][6]
end
println("Initial epoch check 1: $initialEpochSB1")

(statesSB1_v2::Vector{Vector{Float64}}, timesSB1_v2::Vector{Float64}) = primaryEclipticToRotating41(SB1DynamicsModel, "ECLIPJ2000", 1, statesSB1EEclipJ2000, epochTimesSB1)
xSB1_v2::Vector{Float64} = zeros(Float64, nStates_v)
ySB1_v2::Vector{Float64} = zeros(Float64, nStates_v)
zSB1_v2::Vector{Float64} = zeros(Float64, nStates_v)
xdotSB1_v2::Vector{Float64} = zeros(Float64, nStates_v)
ydotSB1_v2::Vector{Float64} = zeros(Float64, nStates_v)
zdotSB1_v2::Vector{Float64} = zeros(Float64, nStates_v)
thetaMSB1_v2::Vector{Float64} = zeros(Float64, nStates_v)
for s::Int64 in 1:nStates_v
    xSB1_v2[s] = statesSB1_v2[s][1]
    ySB1_v2[s] = statesSB1_v2[s][2]
    zSB1_v2[s] = statesSB1_v2[s][3]
    xdotSB1_v2[s] = statesSB1_v2[s][4]
    ydotSB1_v2[s] = statesSB1_v2[s][5]
    zdotSB1_v2[s] = statesSB1_v2[s][6]
    thetaMSB1_v2[s] = statesSB1_v2[s][7]
end
println("Initial Moon Angle: $(thetaMSB1_v2[1])")

thetaM::Float64 = (pi-getStateByIndex(arcEM, 1)[7] >= 0) ? (pi-getStateByIndex(arcEM, 1)[7]) : (3*pi-getStateByIndex(arcEM, 1)[7])
initialEpochTimeSB1::Float64 = getEpochTime(SB1DynamicsModel, "ECLIPJ2000", "Jan 10 2030", thetaM)
println("Initial epoch check 2: $(SPICE.et2utc(initialEpochTimeSB1, :C, 0))")

MoonEclipJ2000::Vector{Float64} = getEphemerides(initialEpoch, [0.0], "Moon", "Earth_Barycenter", "ECLIPJ2000")[1][1]
MoonSPICEElements::Vector{Float64} = SPICE.oscltx(MoonEclipJ2000, epochTimes[1], B1.gravParam)
MoonSPICEElements[1] = get12CharLength(systemData)*(1-get12MassRatio(systemData))
MoonSPICEElements[3] = 0.0
MoonJ2000::Vector{Float64} = SPICE.conics(MoonSPICEElements[1:8], epochTimes[1])
SunEclipJ2000::Vector{Float64} = getEphemerides(initialEpoch, [0.0], "Sun", "Earth_Barycenter", "ECLIPJ2000")[1][1]
SunSPICEElements::Vector{Float64} = SPICE.oscltx(SunEclipJ2000, epochTimes[1], B1.gravParam)
SunSPICEElements[1] = get41CharLength(systemData)
SunSPICEElements[3] = 0.0
SunJ2000::Vector{Float64} = SPICE.conics(SunSPICEElements[1:8], epochTimes[1])

statesCR3BPEclipJ2000::Vector{Vector{Float64}} = rotatingToPrimaryEclipJ2000(CR3BPDynamicsModel, initialEpoch, orbitArc.states, orbitArc.times)
nStates_CR3BP = getStateCount(orbitArc)
xCR3BPEclipJ2000::Vector{Float64} = zeros(Float64, nStates_CR3BP)
yCR3BPEclipJ2000::Vector{Float64} = zeros(Float64, nStates_CR3BP)
zCR3BPEclipJ2000::Vector{Float64} = zeros(Float64, nStates_CR3BP)
xdotCR3BPEclipJ2000::Vector{Float64} = zeros(Float64, nStates_CR3BP)
ydotCR3BPEclipJ2000::Vector{Float64} = zeros(Float64, nStates_CR3BP)
zdotCR3BPEclipJ2000::Vector{Float64} = zeros(Float64, nStates_CR3BP)
for s::Int64 in 1:nStates_CR3BP
    xCR3BPEclipJ2000[s] = statesCR3BPEclipJ2000[s][1]
    yCR3BPEclipJ2000[s] = statesCR3BPEclipJ2000[s][2]
    zCR3BPEclipJ2000[s] = statesCR3BPEclipJ2000[s][3]
    xdotCR3BPEclipJ2000[s] = statesCR3BPEclipJ2000[s][4]
    ydotCR3BPEclipJ2000[s] = statesCR3BPEclipJ2000[s][5]
    zdotCR3BPEclipJ2000[s] = statesCR3BPEclipJ2000[s][6]
end

(statesSEclipJ2000::Vector{Vector{Float64}}) = rotating12ToPrimaryEcliptic(EMDynamicsModel, "ECLIPJ2000", 4, "Jan 1 2030", arcEM.states, arcEM.times)[1]
xSEclipJ2000::Vector{Float64} = zeros(Float64, nStates)
ySEclipJ2000::Vector{Float64} = zeros(Float64, nStates)
zSEclipJ2000::Vector{Float64} = zeros(Float64, nStates)
xdotSEclipJ2000::Vector{Float64} = zeros(Float64, nStates)
ydotSEclipJ2000::Vector{Float64} = zeros(Float64, nStates)
zdotSEclipJ2000::Vector{Float64} = zeros(Float64, nStates)
for s::Int64 in 1:nStates
    xSEclipJ2000[s] = statesSEclipJ2000[s][1]
    ySEclipJ2000[s] = statesSEclipJ2000[s][2]
    zSEclipJ2000[s] = statesSEclipJ2000[s][3]
    xdotSEclipJ2000[s] = statesSEclipJ2000[s][4]
    ydotSEclipJ2000[s] = statesSEclipJ2000[s][5]
    zdotSEclipJ2000[s] = statesSEclipJ2000[s][6]
end

xE1::Vector{Float64} = zeros(Float64, nStates)
yE1::Vector{Float64} = zeros(Float64, nStates)
zE1::Vector{Float64} = zeros(Float64, nStates)
xdotE1::Vector{Float64} = zeros(Float64, nStates)
ydotE1::Vector{Float64} = zeros(Float64, nStates)
zdotE1::Vector{Float64} = zeros(Float64, nStates)
xE2::Vector{Float64} = zeros(Float64, nStates)
yE2::Vector{Float64} = zeros(Float64, nStates)
zE2::Vector{Float64} = zeros(Float64, nStates)
xdotE2::Vector{Float64} = zeros(Float64, nStates)
ydotE2::Vector{Float64} = zeros(Float64, nStates)
zdotE2::Vector{Float64} = zeros(Float64, nStates)
for th::Int64 = 1:nStates
    E1::Vector{Float64} = getInstantaneousEquilibriumPoint(EMDynamicsModel, 1, thetaSEM[th])
    E2::Vector{Float64} = getInstantaneousEquilibriumPoint(EMDynamicsModel, 2, thetaSEM[th])
    xE1[th] = E1[1]
    yE1[th] = E1[2]
    xE2[th] = E2[1]
    yE2[th] = E2[2]
end

initialStateGuess::Vector{Float64} = append!(getEquilibriumPoint(CR3BPDynamicsModel, 1), [0.0, 0.0, 0.0, 0.0])
targetP::Float64 = getSynodicPeriod(EMDynamicsModel)
solution::MBD.BCR4BP12MultipleShooterProblem = homotopy(systemData, targeter, initialStateGuess, targetP, 1E-8)
refinedSolution::MBD.BCR4BP12MultipleShooterProblem = correct(targeter, solution.nodes[1].state.data[1:7], targetP, 1E-10)
println("\nConverged BCR4BP L1 Orbit:\n\tState:$(refinedSolution.nodes[1].state.data[1:7])\n\tPeriod: $(getPeriod(targeter, refinedSolution))")
arcL1::MBD.BCR4BP12Arc = propagate(propagator, refinedSolution.nodes[1].state.data[1:7], tEM, EMDynamicsModel)
nStatesL1::Int64 = getStateCount(arcL1)
xL1::Vector{Float64} = Vector{Float64}(undef, nStatesL1)
yL1::Vector{Float64} = Vector{Float64}(undef, nStatesL1)
zL1::Vector{Float64} = Vector{Float64}(undef, nStatesL1)
xdotL1::Vector{Float64} = Vector{Float64}(undef, nStatesL1)
ydotL1::Vector{Float64} = Vector{Float64}(undef, nStatesL1)
zdotL1::Vector{Float64} = Vector{Float64}(undef, nStatesL1)
for s::Int64 in 1:nStatesL1
    state::Vector{Float64} = getStateByIndex(arcL1, s)
    xL1[s] = state[1]
    yL1[s] = state[2]
    zL1[s] = state[3]
    xdotL1[s] = state[4]
    ydotL1[s] = state[5]
    zdotL1[s] = state[6]
end

mf = MATLAB.MatFile("Output/BCR4BPDev.mat", "w")
exportCR3BPOrbit(CR3BPOrbit, CR3BPDynamicsModel, mf, :orbitCR3BP)
exportCR3BPTrajectory(xCR3BPEclipJ2000, yCR3BPEclipJ2000, zCR3BPEclipJ2000, xdotCR3BPEclipJ2000, ydotCR3BPEclipJ2000, zdotCR3BPEclipJ2000, orbitArc.times, mf, :trajCR3BPEclipJ2000)
exportCR3BPTrajectory([MoonJ2000[1]], [MoonJ2000[2]], [MoonJ2000[3]], [MoonJ2000[4]], [MoonJ2000[5]], [MoonJ2000[6]], [epochTimes[1]], mf, :MoonInitialState)
exportCR3BPTrajectory([SunJ2000[1]], [SunJ2000[2]], [SunJ2000[3]], [SunJ2000[4]], [SunJ2000[5]], [SunJ2000[6]], [epochTimes[1]], mf, :SunInitialState)
exportBCR4BP12Trajectory(xEM, yEM, zEM, xdotEM, ydotEM, zdotEM, thetaSEM, tEM, mf, :trajBCR4BPEM)
exportBCR4BP41Trajectory(xSB1, ySB1, zSB1, xdotSB1, ydotSB1, zdotSB1, thetaMSB1, tSB1, mf, :trajBCR4BPSB1)
exportBCR4BP41Trajectory(xSB1_v, ySB1_v, zSB1_v, xdotSB1_v, ydotSB1_v, zdotSB1_v, thetaMSB1_v, tSB1_v, mf, :validBCR4BPSB1)
exportBCR4BP12Trajectory(xEM_v, yEM_v, zEM_v, xdotEM_v, ydotEM_v, zdotEM_v, thetaSEM_v, tEM_v, mf, :validBCR4BPEM)
exportBCR4BP12Trajectory(xEEclipJ2000, yEEclipJ2000, zEEclipJ2000, xdotEEclipJ2000, ydotEEclipJ2000, zdotEEclipJ2000, thetaSEM, epochTimes, mf, :trajBCR4BPEEclipJ2000)
exportBCR4BP12Trajectory(xEM_v2, yEM_v2, zEM_v2, xdotEM_v2, ydotEM_v2, zdotEM_v2, thetaSEM_v2, timesEM_v2, mf, :validBCR4BPEM2)
exportBCR4BP41Trajectory(xSB1EEclipJ2000, ySB1EEclipJ2000, zSB1EEclipJ2000, xdotSB1EEclipJ2000, ydotSB1EEclipJ2000, zdotSB1EEclipJ2000, thetaMSB1_v, epochTimesSB1, mf, :trajBCR4BPSB1EEclipJ2000)
exportBCR4BP12Trajectory(xSB1_v2, ySB1_v2, zSB1_v2, xdotSB1_v2, ydotSB1_v2, zdotSB1_v2, thetaMSB1_v2, timesSB1_v2, mf, :validBCR4BPSB12)
exportBCR4BP12Trajectory(xSEclipJ2000, ySEclipJ2000, zSEclipJ2000, xdotSEclipJ2000, ydotSEclipJ2000, zdotSEclipJ2000, thetaSEM, epochTimes, mf, :trajBCR4BPSEclipJ2000)
exportBCR4BP12Trajectory(xE1, yE1, zE1, xdotE1, ydotE1, zdotE1, thetaSEM, tEM, mf, :E1BCR4BPEM)
exportBCR4BP12Trajectory(xE2, yE2, zE2, xdotE2, ydotE2, zdotE2, thetaSEM, tEM, mf, :E2BCR4BPEM)
exportBCR4BP12Trajectory(xL1, yL1, zL1, xdotL1, ydotL1, zdotL1, thetaSEM, tEM, mf, :L1Orbit)
MATLAB.close(mf)

SPICE.kclear()

println()
end
