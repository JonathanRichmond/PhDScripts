"""
Script for Keplerian code development

Author: Jonathan Richmond
C: 9/22/25
"""
module KepDev
println()

using MBD

SSystemData = MBD.KSystemData("Sun")
SESystemData = MBD.CR3BPSystemData("Sun", "Earth")
SDynamicsModel = MBD.KDynamicsModel(SSystemData)
SEDynamicsModel = MBD.CR3BPDynamicsModel(SESystemData)

propagator = MBD.Propagator()

q0_dim::Vector{Float64} = [149600000, 10000, 10000, -1.5, 29.78, 1.0]
o0::Vector{Float64} = getOrbitalElements(SDynamicsModel, q0_dim)
println(o0)
arc::MBD.KArc = propagate(propagator, q0_dim, [0.0, 3600*24*30], SDynamicsModel)
qf_dim::Vector{Float64} = getStateByIndex(arc, -1)
of::Vector{Float64} = getOrbitalElements(SDynamicsModel, qf_dim)
println(of)

o_new::Vector{Float64} = append!(o0[1:5], of[6])
qf_new::Vector{Float64} = getCartesianState(SDynamicsModel, o_new)
println(qf_dim)
println(qf_new)

println()
end
