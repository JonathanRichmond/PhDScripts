"""
Script for Earth-Moon CR3BP lunar free return prograde orbit family

Author: Jonathan Richmond
C: 7/1/25
"""
module EMLFRPro
println()

using MBD, GLMakie, Logging

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../CR3BPTargeters/PlanarLFRX.jl")
include("../Utilities/Export.jl")
include("../Utilities/Plot.jl")

systemData = MBD.CR3BPSystemData("Earth", "Moon")
dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
Earth::MBD.BodyData, Moon::MBD.BodyData = systemData.primaryData[1], systemData.primaryData[2]
LagrangePoints::Vector{Vector{Float64}} = [getEquilibriumPoint(dynamicsModel, l) for l = 1:5]

propagator = MBD.Propagator()

node = MBD.CR3BPNode(0.0, [0.8, 0, 0, 0, 0.5, 0], dynamicsModel)
h::Float64 = 185/getCharLength(systemData)
altitudeConstraint = MBD.CR3BPAltitudeConstraint(node, 1, h)
flightPathAngleConstraint = MBD.CR3BPFlightPathAngleConstraint(node, 1, 0.0)
targeter = PlanarLFRXTargeter(dynamicsModel)

# h::Float64 = 185
# v::Float64 = 
# initialStateGuess::Vector{Float64} = getPrimaryState(dynamicsModel, 1)-[(Earth.bodyRadius+h)/getCharLength(systemData), 0, 0, 0, 0, 0]
# tSpanGuess::Vector{Float64} = [0, 7.00247]
# solution1::MBD.CR3BPMultipleShooterProblem = correct(targeter, initialStateGuess, tSpanGuess, targetJC)
# println("Converged Orbit 1:\n\tIC:\t$(solution1.nodes[1].state.data[1:6])\n\tP:\t$(getPeriod(targeter, solution1))\n\tJC:\t$(getJacobiConstant(dynamicsModel, solution1.nodes[1].state.data[1:6]))\n")

println()
end
