"""
Auxiliary script for code development

Author: Jonathan LeFevre Richmond
U: 2/9/26
"""
module Aux
println()

using MBD, Logging, MATLAB

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../CR3BPTargeters/PlanarPerpJC.jl")

systemData = MBD.CR3BPSystemData("Sun", "Venus")
SSystemData = MBD.KSystemData("Sun")
dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
SDynamicsModel = MBD.KDynamicsModel(SSystemData)

targeter = PlanarPerpJCTargeter(dynamicsModel)
propagator = MBD.Propagator()

orbit::MBD.CR3BPPeriodicOrbit = interpOrbit(targeter, "FamilyData/CR3BPSVL1Lyapunovs.csv", "JC", 3.0003)
println("Orbit:\n\tIC:\t$(orbit.initialCondition)\n\tP:\t$(orbit.period)\n\tJC:\t$(getJacobiConstant(orbit))\n")
manifoldArc::MBD.CR3BPManifoldArc = getManifoldArcByTime(orbit, "Unstable", "Positive", 10000/getCharLength(systemData), 0.3)
arc::MBD.CR3BPArc = propagate(propagator, real(manifoldArc.initialCondition), [0, 2*pi], dynamicsModel)
states_SI::Vector{Vector{Float64}} = rotatingToPrimaryInertial(dynamicsModel, 1, arc.states, arc.times)
a::Vector{Float64} = zeros(Float64, getStateCount(arc))
e::Vector{Float64} = zeros(Float64, getStateCount(arc))
r::Vector{Float64} = zeros(Float64, getStateCount(arc))
for s::Int64 in 1:getStateCount(arc)
    Q_SI::Vector{Float64} = append!(states_SI[s][1:3].*getCharLength(systemData), states_SI[s][4:6].*getCharLength(systemData)./getCharTime(systemData))
    oe::Vector{Float64} = getOrbitalElements(SDynamicsModel, Q_SI)
    a[s] = oe[1]
    e[s] = oe[2]
    r[s] = getExcursion(dynamicsModel, 2, getStateByIndex(arc, s))/getCharLength(systemData)
end
r_a::Vector{Float64} = -0.001:-0.001:-0.25
a_V::Vector{Float64} = systemData.primaryData[2].gravParam ./ ((r_a .* getCharLength(systemData)) .^ 2)
a_S::Vector{Float64} = systemData.primaryData[1].gravParam ./ (((1 .+ r_a) .* getCharLength(systemData)) .^ 2)
d::Vector{Float64} = a_V ./ a_S

mf = MATLAB.MatFile("Output/Test.mat", "w")
MATLAB.put_variable(mf, :a, a)
MATLAB.put_variable(mf, :e, e)
MATLAB.put_variable(mf, :r, r)
MATLAB.put_variable(mf, :ra, r_a)
MATLAB.put_variable(mf, :d, d)
MATLAB.close(mf)

println()
end
