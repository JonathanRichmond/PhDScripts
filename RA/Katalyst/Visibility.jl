"""
Script for computing visibility

Author: Jonathan Richmond
C: 7/7/25
U: 7/17/25
"""
module Vis
println()

using MBD, LinearAlgebra, Logging, MATLAB, Statistics
using Base.Threads

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../../CR3BPTargeters/PlanarLFRX.jl")
include("../../CR3BPTargeters/PlanarPerpJC.jl")
include("../../CR3BPTargeters/SpatialPerpJCMS.jl")
include("../../CR3BPTargeters/SpatialPerpVy.jl")
include("../../Utilities/Export.jl")
include("../../Utilities/Visibility.jl")

function sampleRegion(nPoints::Int64, axisLimits::Matrix{Float64}, r_E::Vector{Float64}, r_M::Vector{Float64}, R_E::Float64, R_M::Float64)
    points::Matrix{Float64} = zeros(Float64, (3,nPoints))
    pointCount::Int64 = 0
    while pointCount < nPoints
        point::Vector{Float64} = axisLimits[:,1] + rand(3) .* (axisLimits[:,2]-axisLimits[:,1])
        # radius::Float64 = 2*(R_E+(R_M-R_E)*(point[1]-axisLimits[1,1]))
        radius::Float64 = (point[1]-axisLimits[1,1])*tan(10*pi/180)
        if (LinearAlgebra.norm(point[2:3]) <= radius) && (LinearAlgebra.norm(point-r_E) > R_E) && (LinearAlgebra.norm(point-r_M) > R_M)
            @inbounds points[:,pointCount] = point
            pointCount += 1
        end
    end

    return points
end

systemData = MBD.CR3BPSystemData("Earth", "Moon")
dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
BCR4BPSystemData = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
BCR4BPDynamicsModel = MBD.BCR4BP12DynamicsModel(BCR4BPSystemData)
L2::Vector{Float64} = getEquilibriumPoint(dynamicsModel, 2)

# obsTargeter = SpatialPerpVyTargeter(dynamicsModel)
obsTargeter = PlanarPerpJCTargeter(dynamicsModel)
# obsTargeter = SpatialPerpJCMSTargeter(dynamicsModel)
targTargeter = PlanarLFRXTargeter(dynamicsModel)
propagator = MBD.Propagator()
tSyn::Float64 = getSynodicPeriod(BCR4BPDynamicsModel)

p::Int64, q::Int64 = 1, 1
# obsOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(obsTargeter, "FamilyData/CR3BPEMResonant3_2ProsSpatial.csv", "Period", tSyn*q/p)
# obsOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(obsTargeter, "FamilyData/CR3BPEML1Halos.csv", "Period", tSyn*q/p)
# obsOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(obsTargeter, "FamilyData/CR3BPEML2Halos.csv", "Period", tSyn*q/p)
# obsOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(obsTargeter, "FamilyData/CR3BPEMResonant2_1Rets.csv", "Period", tSyn*q/p)
obsOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(obsTargeter, "FamilyData/CR3BPEMResonant2_1Pros.csv", "Period", tSyn*q/p)
# obsOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(obsTargeter, "FamilyData/CR3BPEMResonant3_2Pros.csv", "Period", tSyn*q/p)
# obsOrbit::MBD.CR3BPMSPeriodicOrbit = interpOrbit(obsTargeter, "FamilyData/CR3BPEMCyclersSpatial.csv", "Period", tSyn*q/p, 3)
# obsOrbit::MBD.CR3BPMSPeriodicOrbit = interpOrbit(obsTargeter, "FamilyData/CR3BPEML1P7HO1s.csv", "Period", tSyn*q/p, 15)
println("Converged Observer Orbit:\n\tIC:\t$(obsOrbit.initialCondition)\n\tP:\t$(obsOrbit.period) ($(obsOrbit.period*getCharTime(systemData)/3600/24) d.)\n\tJC:\t$(getJacobiConstant(obsOrbit))\n\tStab.:\t$(getStabilityIndex(obsOrbit))\n")
targTraj::MBD.CR3BPMultipleShooterProblem = interpSolution(targTargeter, "FamilyData/CR3BPEMLFRPros.csv", "x", getPrimaryState(dynamicsModel, 2)[1]+100/getCharLength(systemData), 3, 185/getCharLength(systemData), 0.0)
println("Converged Target Trajectory:\n\tIC:\t$(targTraj.nodes[end].state.data[1:6])\n\tTOF:\t$(getTOF(targTargeter, targTraj)) ($(getTOF(targTargeter, targTraj)*getCharTime(systemData)/3600/24) d.)\n\tJC:\t$(getJacobiConstant(dynamicsModel, targTraj.nodes[1].state.data[1:6]))\n")

println("Sampling region...")
nPoints::Int64 = 20000
EarthPos::Vector{Float64} = getPrimaryState(dynamicsModel, 1)[1:3]
MoonPos::Vector{Float64} = getPrimaryState(dynamicsModel, 2)[1:3]
r_E::Float64 = systemData.primaryData[1].bodyRadius/getCharLength(systemData)
r_M::Float64 = systemData.primaryData[2].bodyRadius/getCharLength(systemData)
# maxValue::Float64 = 2*r_E
maxValue::Float64 = (L2[1]-EarthPos[1])*tan(10*pi/180)
region::Matrix{Float64} = sampleRegion(nPoints, [EarthPos[1] L2[1]; -maxValue maxValue; -maxValue maxValue], EarthPos, MoonPos, r_E, r_M)

targR::Float64 = 3.5E-3/getCharLength(systemData)
Cd::Float64 = 0.5
thetaS0::Vector{Float64} = collect(range(0, 2*pi, 91))
xi0::Float64 = 0.0
times::Vector{Float64} = collect(range(0, tSyn, 501))
obsICs::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, q)
obsArcs::Vector{MBD.CR3BPArc} = Vector{MBD.CR3BPArc}(undef, q)
for o::Int64 = 1:q
    xi::Float64 = xi0+2*pi*(o-1)/q
    orbitArc::MBD.CR3BPArc = propagate(propagator, obsOrbit.initialCondition, [0, xi*obsOrbit.period/(2*pi)], dynamicsModel)
    obsICs[o] = getStateByIndex(orbitArc, -1)[1:6]
    obsArcs[o] = propagate(propagator, obsICs[o], times, dynamicsModel)
end

println("Analyzing phasing and visibility...")
thetaSdot::Float64 = evaluateEquations(BCR4BPDynamicsModel, MBD.SIMPLE, 0.0, [0.8, 0, 0, 0, 0, 0, thetaS0[1]])[7]
V::Matrix{Float64} = zeros(Float64, (length(thetaS0),length(times)))
Threads.@threads for index::Int64 = 1:(length(thetaS0)*length(times))
    theta::Int64 = div(index-1, length(times))+1
    t::Int64 = mod(index-1, length(times))+1
    thetaS::Float64 = thetaS0[theta]+thetaSdot*(times[t]-times[1])
    r_O::Matrix{Float64} = hcat([getStateByIndex(obsArcs[o], findfirst(obsArcs[o].times .== times[t]))[1:3] for o = 1:q]...)
    brightness::Matrix{Float64} = visibility(r_O, region, EarthPos, MoonPos, getPrimaryState(BCR4BPDynamicsModel, 4, thetaS)[1:3], 10*pi/180, 10*pi/180, 30*pi/180, r_M, getCharLength(systemData), targR, Cd, get4Distance(BCR4BPDynamicsModel))
    vis::Vector{Bool} = vec(any(brightness .<= 18.0, dims = 1))
    V[theta,t] = 100.0*count(vis)/nPoints
end
avgV::Vector{Float64} = vec(Statistics.mean(V, dims = 2))
maxIndex::Int64 = argmax(avgV)

mf = MATLAB.MatFile("RA/Katalyst/Visibility.mat", "w")
exportCR3BPOrbit(obsOrbit, mf, :ObserverOrbit)
exportArrays(obsICs, mf, :ObserverICs)
# exportCR3BPTrajectory(tarTraj, mf, :TargetTrajectory)
MATLAB.put_variable(mf, :Region, region)
MATLAB.put_variable(mf, :SunAngle, thetaS0)
MATLAB.put_variable(mf, :AverageVisibility, avgV)
MATLAB.put_variable(mf, :Index, maxIndex)
MATLAB.put_variable(mf, :Times, times./tSyn)
MATLAB.put_variable(mf, :Visible, V[maxIndex,:])
MATLAB.close(mf)

println()
end
