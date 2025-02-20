"""
Export utility functions

Author: Jonathan Richmond
C: 2/19/25
"""

using MBD, MATLAB

export exportBCR4BP12Trajectory, exportBCR4BP12Trajectory, exportCR3BPOrbit, exportCR3BPTrajectory

"""
    BCR4BP12Traj(x, y, z, xdot, ydot, zdot, theta4, t)

BCR4BP P1-P2 trajectory export object

# Arguments
- `x::Vector{Float64}`: x data [ndim]
- `y::Vector{Float64}`: y data [ndim]
- `z::Vector{Float64}`: z data [ndim]
- `xdot::Vector{Float64}`: xdot data [ndim]
- `ydot::Vector{Float64}`: ydot data [ndim]
- `zdot::Vector{Float64}`: zdot data [ndim]
- `theta4::Vector{Float64}`: theta4 data [ndim]
- `t::Vector{Float64}`: time data [ndim]
"""
struct BCR4BP12Traj
    x::Vector{Float64}                                                  # x data [ndim]
    y::Vector{Float64}                                                  # y data [ndim]
    z::Vector{Float64}                                                  # z data [ndim]
    xdot::Vector{Float64}                                               # xdot data [ndim]
    ydot::Vector{Float64}                                               # ydot data [ndim]
    zdot::Vector{Float64}                                               # zdot data [ndim]
    theta4::Vector{Float64}                                             # theta4 data [ndim]
    t::Vector{Float64}                                                  # time data [ndim]

    function BCR4BP12Traj(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, theta4::Vector{Float64}, t::Vector{Float64})
        this = new(x, y, z, xdot, ydot, zdot, theta4, t)

        return this
    end
end

"""
    BCR4BP41Traj(x, y, z, xdot, ydot, zdot, theta2, t)

BCR4BP P4-B1 trajectory export object

# Arguments
- `x::Vector{Float64}`: x data [ndim]
- `y::Vector{Float64}`: y data [ndim]
- `z::Vector{Float64}`: z data [ndim]
- `xdot::Vector{Float64}`: xdot data [ndim]
- `ydot::Vector{Float64}`: ydot data [ndim]
- `zdot::Vector{Float64}`: zdot data [ndim]
- `theta2::Vector{Float64}`: theta2 data [ndim]
- `t::Vector{Float64}`: time data [ndim]
"""
struct BCR4BP41Traj
    x::Vector{Float64}                                                  # x data [ndim]
    y::Vector{Float64}                                                  # y data [ndim]
    z::Vector{Float64}                                                  # z data [ndim]
    xdot::Vector{Float64}                                               # xdot data [ndim]
    ydot::Vector{Float64}                                               # ydot data [ndim]
    zdot::Vector{Float64}                                               # zdot data [ndim]
    theta2::Vector{Float64}                                             # theta2 data [ndim]
    t::Vector{Float64}                                                  # time data [ndim]

    function BCR4BP41Traj(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, theta2::Vector{Float64}, t::Vector{Float64})
        this = new(x, y, z, xdot, ydot, zdot, theta2, t)

        return this
    end
end

"""
    CR3BPTraj(x, y, z, xdot, ydot, zdot, t)

CR3BP trajectory export object

# Arguments
- `x::Vector{Float64}`: x data [ndim]
- `y::Vector{Float64}`: y data [ndim]
- `z::Vector{Float64}`: z data [ndim]
- `xdot::Vector{Float64}`: xdot data [ndim]
- `ydot::Vector{Float64}`: ydot data [ndim]
- `zdot::Vector{Float64}`: zdot data [ndim]
- `t::Vector{Float64}`: time data [ndim]
"""
struct CR3BPTraj
    x::Vector{Float64}                                                  # x data [ndim]
    y::Vector{Float64}                                                  # y data [ndim]
    z::Vector{Float64}                                                  # z data [ndim]
    xdot::Vector{Float64}                                               # xdot data [ndim]
    ydot::Vector{Float64}                                               # ydot data [ndim]
    zdot::Vector{Float64}                                               # zdot data [ndim]
    t::Vector{Float64}                                                  # time data [ndim]

    function CR3BPTraj(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, t::Vector{Float64})
        this = new(x, y, z, xdot, ydot, zdot, t)

        return this
    end
end

"""
    exportBCR4BP12Trajectory(x, y, z, xdot, ydot, zdot, t, theta4, file, name)

Export BCR4BP P1-P2 trajectory data to MAT file

# Arguments
- `x::Vector{Float64}`: x data [ndim]
- `y::Vector{Float64}`: y data [ndim]
- `z::Vector{Float64}`: z data [ndim]
- `xdot::Vector{Float64}`: xdot data [ndim]
- `ydot::Vector{Float64}`: ydot data [ndim]
- `zdot::Vector{Float64}`: zdot data [ndim]
- `theta4::Vector{Float64}`: theta4 data [ndim]
- `t::Vector{Float64}`: time data [ndim]
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportBCR4BP12Trajectory(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, theta4::Vector{Float64}, t::Vector{Float64}, file::MATLAB.MatFile, name::Symbol)
    traj = BCR4BP12Traj(x, y, z, xdot, ydot, zdot, theta4, t)
    MATLAB.put_variable(file, name, traj)
end

"""
    exportBCR4BP41Trajectory(x, y, z, xdot, ydot, zdot, t, theta2, file, name)

Export BCR4BP P4-B1 trajectory data to MAT file

# Arguments
- `x::Vector{Float64}`: x data [ndim]
- `y::Vector{Float64}`: y data [ndim]
- `z::Vector{Float64}`: z data [ndim]
- `xdot::Vector{Float64}`: xdot data [ndim]
- `ydot::Vector{Float64}`: ydot data [ndim]
- `zdot::Vector{Float64}`: zdot data [ndim]
- `theta2::Vector{Float64}`: theta2 data [ndim]
- `t::Vector{Float64}`: time data [ndim]
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportBCR4BP41Trajectory(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, theta2::Vector{Float64}, t::Vector{Float64}, file::MATLAB.MatFile, name::Symbol)
    traj = BCR4BP41Traj(x, y, z, xdot, ydot, zdot, theta2, t)
    MATLAB.put_variable(file, name, traj)
end

"""
    exportCR3BPOrbit(orbit, dynamicsModel, file, name)

Export CR3BP orbit data to MAT file

# Arguments
- `orbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportCR3BPOrbit(orbit::MBD.CR3BPPeriodicOrbit, dynamicsModel::MBD.CR3BPDynamicsModel, file::MATLAB.MatFile, name::Symbol)
    propagator = MBD.Propagator()
    orbitArc::MBD.CR3BPArc = propagate(propagator, orbit.initialCondition, [0, orbit.period], dynamicsModel)
    nStates::Int64 = getStateCount(orbitArc)
    x::Vector{Float64} = zeros(Float64, nStates)
    y::Vector{Float64} = zeros(Float64, nStates)
    z::Vector{Float64} = zeros(Float64, nStates)
    xdot::Vector{Float64} = zeros(Float64, nStates)
    ydot::Vector{Float64} = zeros(Float64, nStates)
    zdot::Vector{Float64} = zeros(Float64, nStates)
    t::Vector{Float64} = zeros(Float64, nStates)
    for s::Int64 in 1:nStates
        state::Vector{Float64} = getStateByIndex(orbitArc, s)
        x[s] = state[1]
        y[s] = state[2]
        z[s] = state[3]
        xdot[s] = state[4]
        ydot[s] = state[5]
        zdot[s] = state[6]
        t[s] = getTimeByIndex(orbitArc, s)
    end
    exportCR3BPTrajectory(x, y, z, xdot, ydot, zdot, t, file, name)
end

"""
    exportCR3BPTrajectory(x, y, z, xdot, ydot, zdot, t, file, name)

Export CR3BP trajectory data to MAT file

# Arguments
- `x::Vector{Float64}`: x data [ndim]
- `y::Vector{Float64}`: y data [ndim]
- `z::Vector{Float64}`: z data [ndim]
- `xdot::Vector{Float64}`: xdot data [ndim]
- `ydot::Vector{Float64}`: ydot data [ndim]
- `zdot::Vector{Float64}`: zdot data [ndim]
- `t::Vector{Float64}`: time data [ndim]
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportCR3BPTrajectory(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, t::Vector{Float64}, file::MATLAB.MatFile, name::Symbol)
    traj = CR3BPTraj(x, y, z, xdot, ydot, zdot, t)
    MATLAB.put_variable(file, name, traj)
end
