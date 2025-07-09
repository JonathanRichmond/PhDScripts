"""
Script for computing orbit metrics

Author: Jonathan Richmond
C: 7/7/25
"""
module Vis
println()

using MBD, LinearAlgebra, Logging, MATLAB
using Base.Threads

global_logger(ConsoleLogger(stderr, Logging.Warn)) # Debug, Info, Warn, Error

include("../../CR3BPTargeters/PlanarLFRX.jl")
include("../../CR3BPTargeters/SpatialPerpVy.jl")
include("../../Utilities/Export.jl")

function checkExclusion(r_O::Vector{Float64}, r_T::Vector{Float64}, primaryPos::Vector{Float64}, sigma::Float64)
    r_OT::Vector{Float64} = r_T-r_O
    r_OP::Vector{Float64} = primaryPos-r_O
    kappa::Float64 = acos(LinearAlgebra.dot(r_OT, r_OP)/(LinearAlgebra.norm(r_OT)*LinearAlgebra.norm(r_OP)))

    return (abs(kappa) <= sigma)
end

function getLunarExclusionAngle(dynamicsModel::MBD.CR3BPDynamicsModel, r_O::Vector{Float64}, MoonPos::Vector{Float64}, SunPos::Vector{Float64}, sigmaMDeg::Float64)
    r_OM::Vector{Float64} = MoonPos-r_O
    r_OS::Vector{Float64} = SunPos-r_O
    Chi::Float64 = acos(LinearAlgebra.dot(r_OM, r_OS)/(LinearAlgebra.norm(r_OM)*LinearAlgebra.norm(r_OS)))

    return sigmaMDeg*Chi/180+atan(dynamicsModel.systemData.primaryData[2].bodyRadius/(getCharLength(dynamicsModel)*LinearAlgebra.norm(MoonPos-r_O)))
end

function visibility(dynamicsModel::MBD.CR3BPDynamicsModel, R4BPModel::MBD.BCR4BP12DynamicsModel, r_O::Vector{Float64}, r_T::Vector{Float64}, objR::Float64, Cd::Float64, thetaS::Float64)
    SunPos::Vector{Float64} = getPrimaryState(R4BPModel, 4, thetaS)[1:3]
    isEarthExcluded = checkExclusion(r_O, r_T, getPrimaryState(dynamicsModel, 1)[1:3], 10*pi/180)
    isMoonExcluded = checkExclusion(r_O, r_T, getPrimaryState(dynamicsModel, 2)[1:3], getLunarExclusionAngle(dynamicsModel, r_O, getPrimaryState(dynamicsModel, 2)[1:3], SunPos, 10.0))
    isSunExcluded = checkExclusion(r_O, r_T, SunPos, 30*pi/180)
    (isEarthExcluded || isMoonExcluded || isSunExcluded) && return 100.0

    r_TO::Vector{Float64} = r_O-r_T
    r_TS::Vector{Float64} = SunPos-r_T
    alpha::Float64 = acos(LinearAlgebra.dot(r_TO, r_TS)/(LinearAlgebra.norm(r_TO)*LinearAlgebra.norm(r_TS)))

    return -26.832-2.5*log10(LinearAlgebra.norm(r_T-SunPos)^2*2*objR^2*Cd*(sin(alpha)+(pi-alpha)*cos(alpha))/(get4Distance(R4BPModel)^2*2*pi*LinearAlgebra.norm(-r_TO)^2*3*pi))
end

systemData = MBD.CR3BPSystemData("Earth", "Moon")
dynamicsModel = MBD.CR3BPDynamicsModel(systemData)
BCR4BPSystemData = MBD.BCR4BPSystemData("Earth", "Moon", "Sun", "Earth_Barycenter")
BCR4BPDynamicsModel = MBD.BCR4BP12DynamicsModel(BCR4BPSystemData)
L2::Vector{Float64} = getEquilibriumPoint(dynamicsModel, 2)

obsTargeter = SpatialPerpVyTargeter(dynamicsModel)
tarTargeter = PlanarLFRXTargeter(dynamicsModel)
propagator = MBD.Propagator()
tSyn::Float64 = getSynodicPeriod(BCR4BPDynamicsModel)

p::Int64 = 3
q::Int64 = 5
obsOrbit::MBD.CR3BPPeriodicOrbit = interpOrbit(obsTargeter, "FamilyData/CR3BPEMResonant3_2ProsSpatial.csv", "Period", tSyn*q/p)
println("Converged Observer Orbit:\n\tIC:\t$(obsOrbit.initialCondition)\n\tP:\t$(obsOrbit.period) ($(obsOrbit.period*getCharTime(systemData)/3600/24) d.)\n\tJC:\t$(getJacobiConstant(obsOrbit))\n\tStab.:\t$(getStabilityIndex(obsOrbit))\n")
tarTraj::MBD.CR3BPMultipleShooterProblem = interpSolution(tarTargeter, "FamilyData/CR3BPEMLFRPros.csv", "x", getPrimaryState(dynamicsModel, 2)[1]+100/getCharLength(systemData), 3, 185/getCharLength(systemData), 0.0)
println("Converged Target Trajectory:\n\tIC:\t$(tarTraj.nodes[end].state.data[1:6])\n\tTOF:\t$(getTOF(tarTargeter, tarTraj)) ($(getTOF(tarTargeter, tarTraj)*getCharTime(systemData)/3600/24) d.)\n\tJC:\t$(getJacobiConstant(dynamicsModel, tarTraj.nodes[1].state.data[1:6]))\n")

r_E::Vector{Float64} = getPrimaryState(dynamicsModel, 1)[1:3]
r_M::Vector{Float64} = getPrimaryState(dynamicsModel, 2)[1:3]
R_E::Float64 = systemData.primaryData[1].bodyRadius/getCharLength(systemData)
R_M::Float64 = systemData.primaryData[2].bodyRadius/getCharLength(systemData)
points::Vector{Vector{Float64}} = []
while length(points) < 4000
    d::Float64 = rand()*(L2[1]-r_E[1])
    theta::Float64 = rand()*2.0*pi
    radius::Float64 = sqrt(rand())*2*(R_E+(R_M-R_E)*d)
    point::Vector{Float64} = [r_E[1]+d, radius*cos(theta), radius*sin(theta)]
    ((LinearAlgebra.norm(point-r_E) > R_E) && (LinearAlgebra.norm(point-r_M) > R_M)) && push!(points, point)
end

tarR::Float64 = 3.54E-3/getCharLength(systemData)
Cd::Float64 = 0.5
thetaS0::Float64 = 0.0
xi0::Float64 = 0.0
obsICs::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, q)
for o::Int64 = 1:q
    xi::Float64 = xi0+2*pi*(o-1)/q
    orbitArc::MBD.CR3BPArc = propagate(propagator, obsOrbit.initialCondition, [0, xi*obsOrbit.period/(2*pi)], dynamicsModel)
    obsICs[o] = getStateByIndex(orbitArc, -1)
end
times::Vector{Float64} = collect(range(0, tSyn, 501))
V::Vector{Float64} = Vector{Float64}(undef, 501)
for t::Int64 = 1:length(times)
    thetaS::Float64 = thetaS0+evaluateEquations(BCR4BPDynamicsModel, MBD.SIMPLE, times[t], [0.8, 0, 0, 0, 0, 0, thetaS0])[7]*times[t]
    totalVis::Vector{Bool} = Vector{Bool}(undef, length(points))
    Threads.@threads for p::Int64 = 1:length(points)
        r_T::Vector{Float64} = points[p]
        vis::Bool = false
        for o::Int64 = 1:q
            arc::MBD.CR3BPArc = propagate(propagator, obsICs[o], [0, times[t]], dynamicsModel)
            r_O::Vector{Float64} = getStateByIndex(arc, -1)[1:3]
            if visibility(dynamicsModel, BCR4BPDynamicsModel, r_O, r_T, tarR, Cd, thetaS) <= 18.0
                vis = true
                break
            end
        end
        totalVis[p] = vis
    end
    V[t] = 100.0*count(identity, totalVis)/length(totalVis)
end

mf = MATLAB.MatFile("RA/Katalyst/Visibility.mat", "w")
exportCR3BPOrbit(obsOrbit, mf, :ObserverOrbit)
exportArrays(obsICs, mf, :ObserverICs)
# exportCR3BPTrajectory(tarTraj, mf, :TargetTrajectory)
exportArrays(points, mf, :Region)
MATLAB.put_variable(mf, :Times, times./tSyn)
MATLAB.put_variable(mf, :Visible, V)
MATLAB.close(mf)

println()
end
