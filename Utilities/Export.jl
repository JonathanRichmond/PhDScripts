"""
Export utility functions

Author: Jonathan Richmond
C: 2/17/25
"""

using MBD, MATLAB

export exportCR3BPOrbit, exportCR3BPTrajectory

"""
    exportCR3BPOrbit(orbit, dynamicsModel, filename)

Export CR3BP orbit data to MAT file

# Arguments
- `orbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `filename::String`: MAT file name
"""
function exportCR3BPOrbit(orbit::MBD.CR3BPPeriodicOrbit, dynamicsModel::MBD.CR3BPDynamicsModel, filename::String)
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
    exportCR3BPTrajectory(x, y, z, xdot, ydot, zdot, t, filename)
end

"""
    exportCR3BPTrajectory(x, y, z, xdot, ydot, zdot, t, filename)

Export CR3BP trajectory data to MAT file

# Arguments
- `x::Vector{Float64}`: x data [ndim]
- `y::Vector{Float64}`: y data [ndim]
- `z::Vector{Float64}`: z data [ndim]
- `xdot::Vector{Float64}`: xdot data [ndim]
- `ydot::Vector{Float64}`: ydot data [ndim]
- `zdot::Vector{Float64}`: zdot data [ndim]
- `t::Vector{Float64}`: t data [ndim]
- `filename::String`: MAT file name
"""
function exportCR3BPTrajectory(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, t::Vector{Float64}, filename::String)
    MATLAB.write_matfile("Output/"*filename, x = x, y = y, z = z, xdot = xdot, ydot = ydot, zdot = zdot, t = t)
end
