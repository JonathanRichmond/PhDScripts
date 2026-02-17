"""
Frame rotation utility functions

Author: Jonathan LeFevre Richmond
C: 2/16/26
"""

using MBD, StaticArrays
using ..EphemeridesLoader: Ephemerides, FrameTransformations

export rotatingToPrimaryEclipJ2000, rotatingToSunEclipJ2000, secondaryEclipJ2000ToRotating

"""
    rotatingToPrimaryEclipJ2000(dynamicsModel, frame, initialEpochTime, states, times)

Return primary-centered Ecliptic J2000 inertial frame states [ndim]

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `frame::FrameSystem`: Frame system object
- `initialEpochTime::Float64`: Initial epoch time [s]
- `states::Vector{Vector{Float64}}`: Rotating states [ndim]
- `times::Vector{Float64}`: Epochs [ndim]
"""
function rotatingToPrimaryEclipJ2000(dynamicsModel::MBD.CR3BPDynamicsModel, frame::FrameTransformations.FrameSystem, initialEpochTime::Float64, states::Vector{Vector{Float64}}, times::Vector{Float64})
    numTimes::Int16 = Int16(length(times))
    (Int16(length(states)) == numTimes) || throw(ArgumentError("Number of state vectors, $(length(states)), must match number of times, $(length(times))"))
    lstar::Float64 = getCharLength(dynamicsModel)
    tstar::Float64 = getCharTime(dynamicsModel)
    bodyInitialStateDim::Vector{Float64} = FrameTransformations.vector6(frame, Int64(dynamicsModel.systemData.primaryData[1].spiceID), Int64(dynamicsModel.systemData.primaryData[2].spiceID), 17, initialEpochTime)
    bodyOrbitalElements::StaticArrays.MVector{6, Float64} = StaticArrays.MVector{6, Float64}(getOrbitalElements(dynamicsModel.systemData.primaryData[1].gravParam, bodyInitialStateDim))
    (dynamicsModel.systemData.primaryNames[2] == "Earth") && (bodyOrbitalElements[3] = 0.0)
    timesDim::Vector{Float64} = times.*tstar
    thetadotDim::Float64 = 1/tstar
    states_primaryInertial::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, numTimes)
    for i in Int16(1):numTimes
        state_primary::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(states[i]-getPrimaryState(dynamicsModel, 1))
        state_primaryDim::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(append!(state_primary[1:3].*lstar, state_primary[4:6].*lstar./tstar))
        bodyElements::Vector{Float64} = append!([lstar, 0.0], bodyOrbitalElements[3:5], [bodyOrbitalElements[6]+timesDim[i]/tstar])
        bodyStateDim::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(getCartesianState(dynamicsModel.systemData.primaryData[1].gravParam, bodyElements))
        xhat::StaticArrays.SVector{3, Float64} = StaticArrays.SVector{3, Float64}(bodyStateDim[1:3]./lstar)
        zhat::StaticArrays.SVector{3, Float64} = StaticArrays.SVector{3, Float64}(LinearAlgebra.cross(bodyStateDim[1:3], bodyStateDim[4:6])./LinearAlgebra.norm(LinearAlgebra.cross(bodyStateDim[1:3], bodyStateDim[4:6])))
        yhat::StaticArrays.SVector{3, Float64} = StaticArrays.SVector{3, Float64}(LinearAlgebra.cross(zhat, xhat))
        C::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}([xhat yhat zhat])
        Cdot::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}([thetadotDim.*yhat -thetadotDim.*xhat zeros(Float64, 3)])
        N::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([C zeros(Float64, (3,3)); Cdot C])
        state_primaryInertialDim::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(N*state_primaryDim)
        states_primaryInertial[i] = append!(state_primaryInertialDim[1:3]./lstar, state_primaryInertialDim[4:6].*tstar./lstar)
    end

    return states_primaryInertial
end

"""
    rotatingToSunEclipJ2000(dynamicsModel, frame, initialEpochTime, states, times)

Return Sun-centered Ecliptic J2000 inertial frame states [ndim]

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: Sun-Planet CR3BP dynamics model object
- `frame:FrameSystem`: Frame system object
- `initialEpochTime::Float64`: Initial epoch time [s]
- `states::Vector{Vector{Float64}}`: Rotating states [ndim]
- `times::Vector{Float64}`: Epochs [ndim]
"""
function rotatingToSunEclipJ2000(dynamicsModel::MBD.CR3BPDynamicsModel, frame::FrameTransformations.FrameSystem, initialEpochTime::Float64, states::Vector{Vector{Float64}}, times::Vector{Float64})
    numTimes::Int16 = Int16(length(times))
    (Int16(length(states)) == numTimes) || throw(ArgumentError("Number of state vectors, $(length(states)), must match number of times, $(length(times))"))
    lstar::Float64 = getCharLength(dynamicsModel)
    tstar::Float64 = getCharTime(dynamicsModel)
    bodyInitialStateDim::Vector{Float64} = FrameTransformations.vector6(frame, Int64(dynamicsModel.systemData.primaryData[1].spiceID), Int64(dynamicsModel.systemData.primaryData[2].spiceID), 17, initialEpochTime)
    bodyOrbitalElements::StaticArrays.MVector{6, Float64} = StaticArrays.MVector{6, Float64}(getOrbitalElements(dynamicsModel.systemData.primaryData[1].gravParam, bodyInitialStateDim))
    (dynamicsModel.systemData.primaryNames[2] == "Earth") && (bodyOrbitalElements[3] = 0.0)
    timesDim::Vector{Float64} = times.*tstar
    thetadotDim::Float64 = 1/tstar
    states_primaryInertial::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(times))
    for i in Int16(1):numTimes
        state_primary::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(states[i]-getPrimaryState(dynamicsModel, 1))
        state_primaryDim::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(append!(state_primary[1:3].*lstar, state_primary[4:6].*lstar./tstar))
        bodyElements::Vector{Float64} = append!([lstar, 0.0], bodyOrbitalElements[3:5], [bodyOrbitalElements[6]+timesDim[i]/tstar])
        bodyStateDim::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(getCartesianState(dynamicsModel.systemData.primaryData[1].gravParam, bodyElements))
        xhat::StaticArrays.SVector{3, Float64} = StaticArrays.SVector{3, Float64}(bodyStateDim[1:3]./lstar)
        zhat::StaticArrays.SVector{3, Float64} = StaticArrays.SVector{3, Float64}(LinearAlgebra.cross(bodyStateDim[1:3], bodyStateDim[4:6])./LinearAlgebra.norm(LinearAlgebra.cross(bodyStateDim[1:3], bodyStateDim[4:6])))
        yhat::StaticArrays.SVector{3, Float64} = StaticArrays.SVector{3, Float64}(LinearAlgebra.cross(zhat, xhat))
        C::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}([xhat yhat zhat])
        Cdot::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}([thetadotDim.*yhat -thetadotDim.*xhat zeros(Float64, 3)])
        N::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([C zeros(Float64, (3,3)); Cdot C])
        state_primaryInertialDim::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(N*state_primaryDim)
        states_primaryInertial[i] = append!(state_primaryInertialDim[1:3]./lstar, state_primaryInertialDim[4:6].*tstar./lstar)
    end

    return states_primaryInertial
end

"""
    secondaryEclipJ2000ToRotating(dynamicsModel, frame, initialEpochTime, states, times)

Return rotating frame states

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `frame::FrameSystem`: Frame system object
- `initialEpochTime::Float64`: Initial epoch time [s]
- `states_secondaryInertial::Vector{Vector{Float64}}`: Secondary-centered inertial states [ndim]
- `times::Vector{Float64}`: Epochs [ndim]
"""
function secondaryEclipJ2000ToRotating(dynamicsModel::MBD.CR3BPDynamicsModel, frame::FrameTransformations.FrameSystem, initialEpochTime::Float64, states_secondaryInertial::Vector{Vector{Float64}}, times::Vector{Float64})
    numTimes::Int16 = Int16(length(times))
    (Int16(length(states_secondaryInertial)) == numTimes) || throw(ArgumentError("Number of state vectors, $(length(states_primaryInertial)), must match number of times, $(length(times))"))
    lstar::Float64 = getCharLength(dynamicsModel)
    tstar::Float64 = getCharTime(dynamicsModel)
    bodyInitialStateDim::Vector{Float64} = FrameTransformations.vector6(frame, Int64(dynamicsModel.systemData.primaryData[1].spiceID), Int64(dynamicsModel.systemData.primaryData[2].spiceID), 17, initialEpochTime)
    bodyOrbitalElements::StaticArrays.MVector{6, Float64} = StaticArrays.MVector{6, Float64}(getOrbitalElements(dynamicsModel.systemData.primaryData[1].gravParam, bodyInitialStateDim))
    (dynamicsModel.systemData.primaryNames[2] == "Earth") && (bodyOrbitalElements[3] = 0.0)
    timesDim::Vector{Float64} = times.*tstar
    thetadotDim::Float64 = 1/tstar
    states::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, length(times))
    for i in Int16(1):numTimes
        state_secondaryInertialDim::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(append!(states_secondaryInertial[i][1:3].*lstar, states_secondaryInertial[i][4:6].*lstar./tstar))
        bodyElements::Vector{Float64} = append!([lstar, 0.0], bodyOrbitalElements[3:5], bodyOrbitalElements[6]+timesDim[i]/tstar)
        bodyStateDim::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(getCartesianState(dynamicsModel.systemData.primaryData[1].gravParam, bodyElements))
        state_primaryInertialDim::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(bodyStateDim+state_secondaryInertialDim)
        xhat::StaticArrays.SVector{3, Float64} = StaticArrays.SVector{3, Float64}(bodyStateDim[1:3]./lstar)
        zhat::StaticArrays.SVector{3, Float64} = StaticArrays.SVector{3, Float64}(LinearAlgebra.cross(bodyStateDim[1:3], bodyStateDim[4:6])./LinearAlgebra.norm(LinearAlgebra.cross(bodyStateDim[1:3], bodyStateDim[4:6])))
        yhat::StaticArrays.SVector{3, Float64} = StaticArrays.SVector{3, Float64}(LinearAlgebra.cross(zhat, xhat))
        C::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}([xhat yhat zhat])
        Cdot::StaticArrays.SMatrix{3, 3, Float64} = StaticArrays.SMatrix{3, 3, Float64}([thetadotDim.*yhat -thetadotDim.*xhat zeros(Float64, 3)])
        N::StaticArrays.SMatrix{6, 6, Float64} = StaticArrays.SMatrix{6, 6, Float64}([C zeros(Float64, (3,3)); Cdot C])
        state_primaryDim::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(N\state_primaryInertialDim)
        state_primary::StaticArrays.SVector{6, Float64} = StaticArrays.SVector{6, Float64}(append!(state_primaryDim[1:3]./lstar, state_primaryDim[4:6].*tstar./lstar))
        states[i] = state_primary+getPrimaryState(dynamicsModel, 1)
    end

    return states
end