"""
Export utility functions

Author: Jonathan Richmond
C: 2/19/25
U: 6/9/25
"""

using MBD, CSV, DataFrames, LinearAlgebra, MATLAB

export CSVExportCR3BPFamily, exportBCR4BP12Orbit, exportBCR4BP12Trajectory, exportCR3BPOrbit
export exportCR3BPTrajectory, fullExportCR3BPFamily, MATExportCR3BPFamily

"""
    BCR4BP12Orb(x, y, z, xdot, ydot, zdot, theta4, t, H, P, varsig)

BCR4BP P1-P2 orbit export object

# Arguments
- `x::Vector{Float64}`: x data [ndim]
- `y::Vector{Float64}`: y data [ndim]
- `z::Vector{Float64}`: z data [ndim]
- `xdot::Vector{Float64}`: xdot data [ndim]
- `ydot::Vector{Float64}`: ydot data [ndim]
- `zdot::Vector{Float64}`: zdot data [ndim]
- `theta4::Vector{Float64}`: theta4 data [ndim]
- `t::Vector{Float64}`: Time data [ndim]
- `H::Vector{Float64}`: Hamiltonian data [ndim]
- `P::Float64`: Period [ndim]
- `varsig::Float64`: Stability index
"""
struct BCR4BP12Orb
    H::Vector{Float64}                                                  # Hamiltonian data [ndim]
    P::Float64                                                          # Period [ndim]
    t::Vector{Float64}                                                  # Time data [ndim]
    theta4::Vector{Float64}                                             # theta4 data [ndim]
    varsig::Float64                                                     # Stability index
    x::Vector{Float64}                                                  # x data [ndim]
    xdot::Vector{Float64}                                               # xdot data [ndim]
    y::Vector{Float64}                                                  # y data [ndim]
    ydot::Vector{Float64}                                               # ydot data [ndim]
    z::Vector{Float64}                                                  # z data [ndim]
    zdot::Vector{Float64}                                               # zdot data [ndim]

    function BCR4BP12Orb(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, theta4::Vector{Float64}, t::Vector{Float64}, H::Vector{Float64}, P::Float64, varsig::Float64)
        this = new(H, P, t, theta4, varsig, x, xdot, y, ydot, z, zdot)

        return this
    end
end

"""
    BCR4BP12Traj(x, y, z, xdot, ydot, zdot, theta4, t, H)

BCR4BP P1-P2 trajectory export object

# Arguments
- `x::Vector{Float64}`: x data [ndim]
- `y::Vector{Float64}`: y data [ndim]
- `z::Vector{Float64}`: z data [ndim]
- `xdot::Vector{Float64}`: xdot data [ndim]
- `ydot::Vector{Float64}`: ydot data [ndim]
- `zdot::Vector{Float64}`: zdot data [ndim]
- `theta4::Vector{Float64}`: theta4 data [ndim]
- `t::Vector{Float64}`: Time data [ndim]
- `H::Vector{Float64}`: Hamiltonian data [ndim]
"""
struct BCR4BP12Traj
    H::Vector{Float64}                                                  # Hamiltonian data [ndim]
    t::Vector{Float64}                                                  # Time data [ndim]
    theta4::Vector{Float64}                                             # theta4 data [ndim]
    x::Vector{Float64}                                                  # x data [ndim]
    xdot::Vector{Float64}                                               # xdot data [ndim]
    y::Vector{Float64}                                                  # y data [ndim]
    ydot::Vector{Float64}                                               # ydot data [ndim]
    z::Vector{Float64}                                                  # z data [ndim]
    zdot::Vector{Float64}                                               # zdot data [ndim]

    function BCR4BP12Traj(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, theta4::Vector{Float64}, t::Vector{Float64})
        this = new(H, t, theta4, x, xdot, y, ydot, z, zdot)

        return this
    end
end

"""
    CR3BPOrb(x, y, z, xdot, ydot, zdot, t, P, JC, varsig)

CR3BP orbit export object

# Arguments
- `x::Vector{Float64}`: x data [ndim]
- `y::Vector{Float64}`: y data [ndim]
- `z::Vector{Float64}`: z data [ndim]
- `xdot::Vector{Float64}`: xdot data [ndim]
- `ydot::Vector{Float64}`: ydot data [ndim]
- `zdot::Vector{Float64}`: zdot data [ndim]
- `t::Vector{Float64}`: Time data [ndim]
- `P::Float64`: Period [ndim]
- `JC::Float64`: Jacobi constant
- `varsig::Float64`: Stability index
"""
struct CR3BPOrb
    JC::Float64                                                         # Jacobi constant
    P::Float64                                                          # Period [ndim]
    t::Vector{Float64}                                                  # Time data [ndim]
    varsig::Float64                                                     # Stability index
    x::Vector{Float64}                                                  # x data [ndim]
    xdot::Vector{Float64}                                               # xdot data [ndim]
    y::Vector{Float64}                                                  # y data [ndim]
    ydot::Vector{Float64}                                               # ydot data [ndim]
    z::Vector{Float64}                                                  # z data [ndim]
    zdot::Vector{Float64}                                               # zdot data [ndim]

    function CR3BPOrb(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, t::Vector{Float64}, P::Float64, JC::Float64, varsig::Float64)
        this = new(JC, P, t, varsig, x, xdot, y, ydot, z, zdot)

        return this
    end
end

"""
    CR3BPTraj(x, y, z, xdot, ydot, zdot, t, JC)

CR3BP trajectory export object

# Arguments
- `x::Vector{Float64}`: x data [ndim]
- `y::Vector{Float64}`: y data [ndim]
- `z::Vector{Float64}`: z data [ndim]
- `xdot::Vector{Float64}`: xdot data [ndim]
- `ydot::Vector{Float64}`: ydot data [ndim]
- `zdot::Vector{Float64}`: zdot data [ndim]
- `t::Vector{Float64}`: Time data [ndim]
- `JC::Float64`: Jacobi constant
"""
struct CR3BPTraj
    JC::Float64                                                         # Jacobi constant
    t::Vector{Float64}                                                  # Time data [ndim]
    x::Vector{Float64}                                                  # x data [ndim]
    xdot::Vector{Float64}                                               # xdot data [ndim]
    y::Vector{Float64}                                                  # y data [ndim]
    ydot::Vector{Float64}                                               # ydot data [ndim]
    z::Vector{Float64}                                                  # z data [ndim]
    zdot::Vector{Float64}                                               # zdot data [ndim]

    function CR3BPTraj(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, t::Vector{Float64}, JC::Float64)
        this = new(JC, t, x, xdot, y, ydot, z, zdot)

        return this
    end
end

"""
    CSVExportCR3BPFamily(family, file)

Export CR3BP family data to CSV file

# Arguments
- `family::CR3BPOrbitFamily`: CR3BP periodic orbit family object
- `file::String`: CSV file name
"""
function CSVExportCR3BPFamily(family::MBD.CR3BPOrbitFamily, file::String)
    nMem::Int64 = getNumMembers(family)
    x::Vector{Float64} = Vector{Float64}(undef, nMem)
    y::Vector{Float64} = Vector{Float64}(undef, nMem)
    z::Vector{Float64} = Vector{Float64}(undef, nMem)
    xdot::Vector{Float64} = Vector{Float64}(undef, nMem)
    ydot::Vector{Float64} = Vector{Float64}(undef, nMem)
    zdot::Vector{Float64} = Vector{Float64}(undef, nMem)
    P::Vector{Float64} = Vector{Float64}(undef, nMem)
    JC::Vector{Float64} = Vector{Float64}(undef, nMem)
    varsig::Vector{Float64} = Vector{Float64}(undef, nMem)
    for o::Int64 = 1:nMem
        orbit::MBD.CR3BPPeriodicOrbit = getMember(family, o)
        x[o] = orbit.initialCondition[1]
        y[o] = orbit.initialCondition[2]
        z[o] = orbit.initialCondition[3]
        xdot[o] = orbit.initialCondition[4]
        ydot[o] = orbit.initialCondition[5]
        zdot[o] = orbit.initialCondition[6]
        P[o] = orbit.period
        JC[o] = getJacobiConstant(orbit)
        varsig[o] = getStabilityIndex(orbit)
    end
    familyData::DataFrames.DataFrame = DataFrames.DataFrame("x" => x, "y" => y, "z" => z, "xdot" => xdot, "ydot" => ydot, "zdot" => zdot, "Period" => P, "JC" => JC, "Stability" => varsig)
    CSV.write(file, familyData)
end

"""
    exportBCR4BP12Orbit(orbit, file, name)

Export BCR4BP P1-P2 orbit data to MAT file

# Arguments
- `orbit::BCR4BP12PeriodicOrbit`: BCR4BP P1-P2 periodic orbit object
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportBCR4BP12Orbit(orbit::MBD.BCR4BP12PeriodicOrbit, file::MATLAB.MatFile, name::Symbol)
    propagator = MBD.Propagator()
    orbitArc::MBD.BCR4BP12Arc = propagate(propagator, orbit.initialCondition, [0, orbit.period], orbit.dynamicsModel)
    nStates::Int64 = getStateCount(orbitArc)
    x::Vector{Float64} = zeros(Float64, nStates)
    y::Vector{Float64} = zeros(Float64, nStates)
    z::Vector{Float64} = zeros(Float64, nStates)
    xdot::Vector{Float64} = zeros(Float64, nStates)
    ydot::Vector{Float64} = zeros(Float64, nStates)
    zdot::Vector{Float64} = zeros(Float64, nStates)
    t::Vector{Float64} = zeros(Float64, nStates)
    theta4::Vector{Float64} = zeros(Float64, nStates)
    H::Vector{Float64} = zeros(Float64, nStates)
    for s::Int64 in 1:nStates
        state::Vector{Float64} = getStateByIndex(orbitArc, s)
        x[s] = state[1]
        y[s] = state[2]
        z[s] = state[3]
        xdot[s] = state[4]
        ydot[s] = state[5]
        zdot[s] = state[6]
        theta4[s] = state[7]
        t[s] = getTimeByIndex(orbitArc, s)
        H[s] = getHamiltonian(orbit.dynamicsModel, state)
    end
    varsig::Float64 = getStabilityIndex(orbit)
    exportBCR4BP12Orbit(x, y, z, xdot, ydot, zdot, theta4, t, H, orbit.period, varsig, file, name)
end

"""
    exportBCR4BP12Orbit(x, y, z, xdot, ydot, zdot, theta4, t, H, P, varsig, file, name)

Export BCR4BP P1-P2 orbit data to MAT file

# Arguments
- `x::Vector{Float64}`: x data [ndim]
- `y::Vector{Float64}`: y data [ndim]
- `z::Vector{Float64}`: z data [ndim]
- `xdot::Vector{Float64}`: xdot data [ndim]
- `ydot::Vector{Float64}`: ydot data [ndim]
- `zdot::Vector{Float64}`: zdot data [ndim]
- `theta4::Vector{Float64}`: theta4 data [ndim]
- `t::Vector{Float64}`: Time data [ndim]
- `H::Vector{Float64}`: Hamiltonian data [ndim]
- `P::Float64`: Period [ndim]
- `varsig::Float64`: Stability index
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportBCR4BP12Orbit(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, theta4::Vector{Float64}, t::Vector{Float64}, H::Vector{Float64}, P::Float64, varsig::Float64, file::MATLAB.MatFile, name::Symbol)
    orbit = BCR4BP12Orb(x, y, z, xdot, ydot, zdot, theta4, t, H, P, varsig)
    MATLAB.put_variable(file, name, orbit)
end

"""
    exportBCR4BP12Trajectory(x, y, z, xdot, ydot, zdot, theta4, t, H, file, name)

Export BCR4BP P1-P2 trajectory data to MAT file

# Arguments
- `x::Vector{Float64}`: x data [ndim]
- `y::Vector{Float64}`: y data [ndim]
- `z::Vector{Float64}`: z data [ndim]
- `xdot::Vector{Float64}`: xdot data [ndim]
- `ydot::Vector{Float64}`: ydot data [ndim]
- `zdot::Vector{Float64}`: zdot data [ndim]
- `theta4::Vector{Float64}`: theta4 data [ndim]
- `t::Vector{Float64}`: Time data [ndim]
- `H::Vector{Float64}`: Hamiltonian data [ndim]
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportBCR4BP12Trajectory(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, theta4::Vector{Float64}, t::Vector{Float64}, H::Vector{Float64}, file::MATLAB.MatFile, name::Symbol)
    traj = BCR4BP12Traj(x, y, z, xdot, ydot, zdot, theta4, t, H)
    MATLAB.put_variable(file, name, traj)
end

"""
    exportCR3BPOrbit(orbit, file, name)

Export CR3BP orbit data to MAT file

# Arguments
- `orbit::CR3BPPeriodicOrbit`: CR3BP periodic orbit object
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportCR3BPOrbit(orbit::MBD.CR3BPPeriodicOrbit, file::MATLAB.MatFile, name::Symbol)
    propagator = MBD.Propagator()
    orbitArc::MBD.CR3BPArc = propagate(propagator, orbit.initialCondition, [0, orbit.period], orbit.dynamicsModel)
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
    JC::Float64 = MBD.getJacobiConstant(orbit)
    varsig::Float64 = getStabilityIndex(orbit)
    exportCR3BPOrbit(x, y, z, xdot, ydot, zdot, t, orbit.period, JC, varsig, file, name)
end

"""
    exportCR3BPOrbit(x, y, z, xdot, ydot, zdot, t, P, JC, varsig, file, name)

Export CR3BP orbit data to MAT file

# Arguments
- `x::Vector{Float64}`: x data [ndim]
- `y::Vector{Float64}`: y data [ndim]
- `z::Vector{Float64}`: z data [ndim]
- `xdot::Vector{Float64}`: xdot data [ndim]
- `ydot::Vector{Float64}`: ydot data [ndim]
- `zdot::Vector{Float64}`: zdot data [ndim]
- `t::Vector{Float64}`: Time data [ndim]
- `P::Float64`: Period [ndim]
- `JC::Float64`: Jacobi constant
- `varsig::Float64`: Stability index
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportCR3BPOrbit(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, t::Vector{Float64}, P::Float64, JC::Float64, varsig::Float64, file::MATLAB.MatFile, name::Symbol)
    orbit = CR3BPOrb(x, y, z, xdot, ydot, zdot, t, P, JC, varsig)
    MATLAB.put_variable(file, name, orbit)
end

"""
    exportCR3BPTrajectory(x, y, z, xdot, ydot, zdot, t, JC, file, name)

Export CR3BP trajectory data to MAT file

# Arguments
- `x::Vector{Float64}`: x data [ndim]
- `y::Vector{Float64}`: y data [ndim]
- `z::Vector{Float64}`: z data [ndim]
- `xdot::Vector{Float64}`: xdot data [ndim]
- `ydot::Vector{Float64}`: ydot data [ndim]
- `zdot::Vector{Float64}`: zdot data [ndim]
- `t::Vector{Float64}`: Time data [ndim]
- `JC::Float64`: Jacobi constant
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportCR3BPTrajectory(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, t::Vector{Float64}, JC::Float64, file::MATLAB.MatFile, name::Symbol)
    traj = CR3BPTraj(x, y, z, xdot, ydot, zdot, t, JC)
    MATLAB.put_variable(file, name, traj)
end

"""
    fullExportCR3BPFamily(family, MATFile, CSVFile)

Export CR3BP family data to MAT and CSV files

# Arguments
- `family::CR3BPOrbitFamily`: CR3BP periodic orbit family object
- `MATFile::String`: MAT file name
- `CSVFile::String`: CSV file name
"""
function fullExportCR3BPFamily(family::MBD.CR3BPOrbitFamily, MATFile::String, CSVFile::String)
    mf = MATLAB.MatFile(MATFile, "w")
    MATExportCR3BPFamily(family, mf)
    MATLAB.close(mf)
    CSVExportCR3BPFamily(family, CSVFile)
end

"""
    MATExportCR3BPFamily(family, file)

Export CR3BP family data to MAT file

# Arguments
- `family::CR3BPOrbitFamily`: CR3BP periodic orbit family object
- `file::MatFile`: MAT file
"""
function MATExportCR3BPFamily(family::MBD.CR3BPOrbitFamily, file::MATLAB.MatFile)
    for o::Int64 = 1:getNumMembers(family)
        orbit::MBD.CR3BPPeriodicOrbit = getMember(family, o)
        exportCR3BPOrbit(orbit, file, Symbol("orbit"*string(o)))
    end
end
