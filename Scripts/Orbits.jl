"""
Script for exporting select orbits

Author: Jonathan LeFevre Richmond
C: 4/9/26
"""
# module Orbits
println("Running Orbits.jl...\n")

using MBD, GLMakie, Logging, MATLAB

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../CR3BPTargeters/PlanarPerpJC.jl")
include("../CR3BPTargeters/SpatialAxialJC.jl")
include("../CR3BPTargeters/SpatialPerpJC.jl")
include("../CR3BPTargeters/SpatialVerticalJC.jl")
include("../Utilities/Export.jl")
include("../Utilities/Plot.jl")

mf = MATLAB.MatFile("Output/Orbits.mat", "w")

systemData = MBD.CR3BPSystemData("Earth", "Moon")
dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
Earth::MBD.BodyData, Moon::MBD.BodyData = systemData.primaryData[1], systemData.primaryData[2]
LagrangePoints::Vector{Vector{Float64}} = [getEquilibriumPoint(dynamicsModel, l) for l = 1:5]

propagator = MBD.Propagator()
LyapunovTargeter = PlanarPerpJCTargeter(dynamicsModel)
AxialTargeter = SpatialAxialJCTargeter(dynamicsModel)
HaloTargeter = SpatialPerpJCTargeter(dynamicsModel)
VerticalTargeter = SpatialVerticalJCTargeter(dynamicsModel)

orbits::Vector{MBD.CR3BPPeriodicOrbit} = []

# Lyapunovs
JCLyap::Vector{Float64} = [2.97, 3.0, 3.03, 3.07, 3.1, 3.13]
for j in eachindex(JCLyap)
    orbit::MBD.CR3BPPeriodicOrbit = interpOrbit(LyapunovTargeter, "FamilyData/CR3BPEML1Lyapunovs.csv", "JC", JCLyap[j])
    println("Earth-Moon orbit:\n\tIC:\t$(orbit.initialCondition)\n\tP:\t$(orbit.period)\n\tJC:\t$(JCLyap[j])\n")
    exportCR3BPOrbit(orbit, mf, Symbol("L1Lyapunov", j))
    push!(orbits, orbit)
end
for j in eachindex(JCLyap)
    orbit::MBD.CR3BPPeriodicOrbit = interpOrbit(LyapunovTargeter, "FamilyData/CR3BPEML2Lyapunovs.csv", "JC", JCLyap[j])
    println("Earth-Moon orbit:\n\tIC:\t$(orbit.initialCondition)\n\tP:\t$(orbit.period)\n\tJC:\t$(JCLyap[j])\n")
    exportCR3BPOrbit(orbit, mf, Symbol("L2Lyapunov", j))
    push!(orbits, orbit)
end

# Southern halos
JCHalo::Vector{Float64} = [3.03, 3.07, 3.1, 3.13]
for j in eachindex(JCHalo)
    orbit::MBD.CR3BPPeriodicOrbit = interpOrbit(HaloTargeter, "FamilyData/CR3BPEML1Halos.csv", "JC", JCHalo[j])
    println("Earth-Moon orbit:\n\tIC:\t$(orbit.initialCondition)\n\tP:\t$(orbit.period)\n\tJC:\t$(JCHalo[j])\n")
    exportCR3BPOrbit(orbit, mf, Symbol("L1Halo", j))
    push!(orbits, orbit)
end
for j in eachindex(JCHalo)
    orbit::MBD.CR3BPPeriodicOrbit = interpOrbit(HaloTargeter, "FamilyData/CR3BPEML2Halos.csv", "JC", JCHalo[j])
    println("Earth-Moon orbit:\n\tIC:\t$(orbit.initialCondition)\n\tP:\t$(orbit.period)\n\tJC:\t$(JCHalo[j])\n")
    exportCR3BPOrbit(orbit, mf, Symbol("L2Halo", j))
    push!(orbits, orbit)
end

# Axials
JCAx::Vector{Float64} = [2.97, 3.0]
for j in eachindex([JCAx[2]]) # SW
    orbit::MBD.CR3BPPeriodicOrbit = interpOrbit(AxialTargeter, "FamilyData/CR3BPEML1Axials.csv", "JC", JCAx[2])
    println("Earth-Moon orbit:\n\tIC:\t$(orbit.initialCondition)\n\tP:\t$(orbit.period)\n\tJC:\t$(JCAx[2])\n")
    exportCR3BPOrbit(orbit, mf, Symbol("L1Axial", j))
    push!(orbits, orbit)
end
for j in eachindex(JCAx) # NW
    orbit::MBD.CR3BPPeriodicOrbit = interpOrbit(AxialTargeter, "FamilyData/CR3BPEML2Axials.csv", "JC", JCAx[j])
    println("Earth-Moon orbit:\n\tIC:\t$(orbit.initialCondition)\n\tP:\t$(orbit.period)\n\tJC:\t$(JCAx[j])\n")
    exportCR3BPOrbit(orbit, mf, Symbol("L2Axial", j))
    push!(orbits, orbit)
end

# Verticals
JCVert::Vector{Float64} = [2.97, 3.0, 3.03, 3.07, 3.1, 3.13]
for j in eachindex(JCVert)
    orbit::MBD.CR3BPPeriodicOrbit = interpOrbit(VerticalTargeter, "FamilyData/CR3BPEML1Verticals.csv", "JC", JCVert[j])
    println("Earth-Moon orbit:\n\tIC:\t$(orbit.initialCondition)\n\tP:\t$(orbit.period)\n\tJC:\t$(JCVert[j])\n")
    exportCR3BPOrbit(orbit, mf, Symbol("L1Vertical", j))
    push!(orbits, orbit)
end
for j in eachindex(JCVert)
    orbit::MBD.CR3BPPeriodicOrbit = interpOrbit(VerticalTargeter, "FamilyData/CR3BPEML2Verticals.csv", "JC", JCVert[j])
    println("Earth-Moon orbit:\n\tIC:\t$(orbit.initialCondition)\n\tP:\t$(orbit.period)\n\tJC:\t$(JCVert[j])\n")
    exportCR3BPOrbit(orbit, mf, Symbol("L2Vertical", j))
    push!(orbits, orbit)
end

# Southern butterflies
JCButt::Vector{Float64} = [2.97, 3.0, 3.03, 3.07]
for j in eachindex(JCButt)
    orbit::MBD.CR3BPPeriodicOrbit = interpOrbit(HaloTargeter, "FamilyData/CR3BPEML2Butterflys.csv", "JC", JCButt[j])
    println("Earth-Moon orbit:\n\tIC:\t$(orbit.initialCondition)\n\tP:\t$(orbit.period)\n\tJC:\t$(JCButt[j])\n")
    exportCR3BPOrbit(orbit, mf, Symbol("L2Butterfly", j))
    push!(orbits, orbit)
end

MATLAB.close(mf)

println("Plotting orbit...")
(figure, axis) = set3DPlotParameters(L"Earth-Moon Orbits", L"$x$ [ndim]", L"$y$ [ndim]", L"$z$ [ndim]")
for o in orbits
    orbitArc::MBD.CR3BPArc = propagate(propagator, o.initialCondition, [0, o.period], dynamicsModel)
    xData::Vector{Float64}, yData::Vector{Float64}, zData::Vector{Float64} = Vector{Float64}(undef, getStateCount(orbitArc)), Vector{Float64}(undef, getStateCount(orbitArc)), Vector{Float64}(undef, getStateCount(orbitArc))
    for s::Int64 in 1:getStateCount(orbitArc)
        xData[s], yData[s], zData[s] = orbitArc.states[s][1], orbitArc.states[s][2], orbitArc.states[s][3]
    end
    GLMakie.lines!(axis, xData, yData, zData, color = :white)
end
GLMakie.scatter!(axis, LagrangePoints[1][1], LagrangePoints[1][2], LagrangePoints[1][3], color = :red, markersize = 5, label = L"$L_{1}$" => (; markersize = 20))
GLMakie.scatter!(axis, LagrangePoints[2][1], LagrangePoints[2][2], LagrangePoints[2][3], color = :orange, markersize = 5, label = L"$L_{2}$" => (; markersize = 20))
GLMakie.scatter!(axis, getPrimaryState(dynamicsModel, 1)[1], getPrimaryState(dynamicsModel, 1)[2], getPrimaryState(dynamicsModel, 1)[3], color = :blue, markerspace = :data, markersize = 2*Earth.bodyRadius/getCharLength(systemData), label = L"\mathrm{Earth}" => (; markersize = 20))
GLMakie.scatter!(axis, getPrimaryState(dynamicsModel, 2)[1], getPrimaryState(dynamicsModel, 2)[2], getPrimaryState(dynamicsModel, 2)[3], color = :gray, markerspace = :data, markersize = 2*Moon.bodyRadius/getCharLength(systemData), label = L"\mathrm{Moon}" => (; markersize = 20))
GLMakie.Legend(figure[1,2], axis)

println()
# end
