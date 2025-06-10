"""
Export utility functions

Author: Jonathan Richmond
C: 2/19/25
U: 6/10/25
"""

using MBD, CSV, DataFrames, LinearAlgebra, MATLAB

export CSVExportCR3BPFamily, exportBCR4BP12Manifold, exportBCR4BP12Orbit, exportBCR4BP12Trajectory
export exportCR3BPManifold, exportCR3BPOrbit, exportCR3BPTrajectory, exportPseudoManifold
export fullExportCR3BPFamily, MATExportCR3BPFamily

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

    function BCR4BP12Traj(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, theta4::Vector{Float64}, t::Vector{Float64}, H::Vector{Float64})
        this = new(H, t, theta4, x, xdot, y, ydot, z, zdot)

        return this
    end
end

"""
    BCR4BP12Man(orbit, arcs, TOF)

BCR4BP P1-P2 manifold export object

# Arguments
- `orbit::BCR4BP12Orb`: BCR4BP P1-P2 periodic orbit export object
- `arcs::Vector{BCR4BP12Traj}`: BCR4BP P1-P2 manifold arc export objects
- `TOF::Float64`: Time-of-flight [ndim]
"""
struct BCR4BP12Man
    arcs::Vector{BCR4BP12Traj}                                          # Manifold arc export objects
    n::Int64                                                            # Number of manifold arcs
    orbit::BCR4BP12Orb                                                  # Periodic orbit export object
    TOF::Float64                                                        # Time-of-flight [ndim]

    function BCR4BP12Man(orbit::BCR4BP12Orb, arcs::Vector{BCR4BP12Traj}, TOF::Float64)
        this = new(arcs, length(arcs), orbit, TOF)

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
    CR3BPMan(orbit, arcs, TOF)

CR3BP manifold export object

# Arguments
- `orbit::CR3BPOrb`: CR3BP periodic orbit export object
- `arcs::Vector{CR3BPTraj}`: CR3BP manifold arc export objects
- `TOF::Float64`: Time-of-flight [ndim]
"""
struct CR3BPMan
    arcs::Vector{CR3BPTraj}                                             # Manifold arc export objects
    n::Int64                                                            # Number of manifold arcs
    orbit::CR3BPOrb                                                     # Periodic orbit export object
    TOF::Float64                                                        # Time-of-flight [ndim]

    function CR3BPMan(orbit::CR3BPOrb, arcs::Vector{CR3BPTraj}, TOF::Float64)
        this = new(arcs, length(arcs), orbit, TOF)

        return this
    end
end

"""
    PseudoMan(orbit, arcs, TOF)

Pseudo-manifold export object

# Arguments
- `orbit::CR3BPOrb`: CR3BP periodic orbit export object
- `arcs::Vector{BCR4BP12Traj}`: BCR4BP P1-P2 pseudo-manifold arc export objects
- `TOF::Float64`: Time-of-flight [ndim]
"""
struct PseudoMan
    arcs::Vector{BCR4BP12Traj}                                          # Manifold arc export objects
    n::Int64                                                            # Number of manifold arcs
    orbit::CR3BPOrb                                                     # Periodic orbit export object
    theta40::Float64                                                    # Initial P4 angle [ndim]
    TOF::Float64                                                        # Time-of-flight [ndim]

    function PseudoMan(orbit::CR3BPOrb, arcs::Vector{BCR4BP12Traj}, TOF::Float64)
        this = new(arcs, length(arcs), orbit, arcs[1].theta4[1], TOF)

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
    exportBCR4BP12Manifold(manifold, file, name)

Export BCR4BP P1-P2 manifold data to MAT file

# Arguments
- `manifold::BCR4BP12Manifold`: BCR4BP P1-P2 manifold object
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportBCR4BP12Manifold(manifold::MBD.BCR4BP12Manifold, file::MATLAB.MatFile, name::Symbol)
    propagator = MBD.Propagator()
    orbitArc::MBD.BCR4BP12Arc = propagate(propagator, manifold.periodicOrbit.initialCondition, [0, manifold.periodicOrbit.period], manifold.periodicOrbit.dynamicsModel)
    nStates::Int64 = getStateCount(orbitArc)
    orbitx::Vector{Float64} = zeros(Float64, nStates)
    orbity::Vector{Float64} = zeros(Float64, nStates)
    orbitz::Vector{Float64} = zeros(Float64, nStates)
    orbitxdot::Vector{Float64} = zeros(Float64, nStates)
    orbitydot::Vector{Float64} = zeros(Float64, nStates)
    orbitzdot::Vector{Float64} = zeros(Float64, nStates)
    orbitt::Vector{Float64} = zeros(Float64, nStates)
    orbittheta4::Vector{Float64} = zeros(Float64, nStates)
    orbitH::Vector{Float64} = zeros(Float64, nStates)
    for s::Int64 in 1:nStates
        state::Vector{Float64} = getStateByIndex(orbitArc, s)
        orbitx[s] = state[1]
        orbity[s] = state[2]
        orbitz[s] = state[3]
        orbitxdot[s] = state[4]
        orbitydot[s] = state[5]
        orbitzdot[s] = state[6]
        orbittheta4[s] = state[7]
        orbitt[s] = getTimeByIndex(orbitArc, s)
        orbitH[s] = getHamiltonian(manifold.periodicOrbit.dynamicsModel, state)
    end
    orb = BCR4BP12Orb(orbitx, orbity, orbitz, orbitxdot, orbitydot, orbitzdot, orbittheta4, orbitt, orbitH, manifold.periodicOrbit.period, getStabilityIndex(manifold.periodicOrbit))
    numArcs::Int64 = length(manifold.initialConditions)
    arcs::Vector{BCR4BP12Traj} = Vector{BCR4BP12Traj}(undef, numArcs)
    for a = 1:numArcs
        arc::MBD.BCR4BP12Arc = propagate(propagator, real(manifold.initialConditions[a]), [0, manifold.TOF], manifold.periodicOrbit.dynamicsModel)
        numStates::Int64 = getStateCount(arc)
        x::Vector{Float64} = zeros(Float64, numStates)
        y::Vector{Float64} = zeros(Float64, numStates)
        z::Vector{Float64} = zeros(Float64, numStates)
        xdot::Vector{Float64} = zeros(Float64, numStates)
        ydot::Vector{Float64} = zeros(Float64, numStates)
        zdot::Vector{Float64} = zeros(Float64, numStates)
        theta4::Vector{Float64} = zeros(Float64, numStates)
        t::Vector{Float64} = zeros(Float64, numStates)
        H::Vector{Float64} = zeros(Float64, numStates)
        for s::Int64 in 1:numStates
            state::Vector{Float64} = getStateByIndex(arc, s)
            x[s] = state[1]
            y[s] = state[2]
            z[s] = state[3]
            xdot[s] = state[4]
            ydot[s] = state[5]
            zdot[s] = state[6]
            theta4[s] = state[7]
            t[s] = getTimeByIndex(arc, s)
            H[s] = getHamiltonian(manifold.periodicOrbit.dynamicsModel, state)
        end
        arcs[a] = BCR4BP12Traj(x, y, z, xdot, ydot, zdot, theta4, t, H)
    end
    man = BCR4BP12Man(orb, arcs, manifold.TOF)
    MATLAB.put_variable(file, name, man)
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
    exportCR3BPManifold(manifold, file, name)

Export CR3BP manifold data to MAT file

# Arguments
- `manifold::CR3BPManifold`: CR3BP manifold object
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportCR3BPManifold(manifold::MBD.CR3BPManifold, file::MATLAB.MatFile, name::Symbol)
    propagator = MBD.Propagator()
    orbitArc::MBD.CR3BPArc = propagate(propagator, manifold.periodicOrbit.initialCondition, [0, manifold.periodicOrbit.period], manifold.periodicOrbit.dynamicsModel)
    nStates::Int64 = getStateCount(orbitArc)
    orbitx::Vector{Float64} = zeros(Float64, nStates)
    orbity::Vector{Float64} = zeros(Float64, nStates)
    orbitz::Vector{Float64} = zeros(Float64, nStates)
    orbitxdot::Vector{Float64} = zeros(Float64, nStates)
    orbitydot::Vector{Float64} = zeros(Float64, nStates)
    orbitzdot::Vector{Float64} = zeros(Float64, nStates)
    orbitt::Vector{Float64} = zeros(Float64, nStates)
    for s::Int64 in 1:nStates
        state::Vector{Float64} = getStateByIndex(orbitArc, s)
        orbitx[s] = state[1]
        orbity[s] = state[2]
        orbitz[s] = state[3]
        orbitxdot[s] = state[4]
        orbitydot[s] = state[5]
        orbitzdot[s] = state[6]
        orbitt[s] = getTimeByIndex(orbitArc, s)
    end
    orb = CR3BPOrb(orbitx, orbity, orbitz, orbitxdot, orbitydot, orbitzdot, orbitt, manifold.periodicOrbit.period, getJacobiConstant(manifold.periodicOrbit), getStabilityIndex(manifold.periodicOrbit))
    numArcs::Int64 = length(manifold.initialConditions)
    arcs::Vector{CR3BPTraj} = Vector{CR3BPTraj}(undef, numArcs)
    for a = 1:numArcs
        arc::MBD.CR3BPArc = propagate(propagator, real(manifold.initialConditions[a]), [0, manifold.TOF], manifold.periodicOrbit.dynamicsModel)
        numStates::Int64 = getStateCount(arc)
        x::Vector{Float64} = zeros(Float64, numStates)
        y::Vector{Float64} = zeros(Float64, numStates)
        z::Vector{Float64} = zeros(Float64, numStates)
        xdot::Vector{Float64} = zeros(Float64, numStates)
        ydot::Vector{Float64} = zeros(Float64, numStates)
        zdot::Vector{Float64} = zeros(Float64, numStates)
        t::Vector{Float64} = zeros(Float64, numStates)
        for s::Int64 in 1:numStates
            state::Vector{Float64} = getStateByIndex(arc, s)
            x[s] = state[1]
            y[s] = state[2]
            z[s] = state[3]
            xdot[s] = state[4]
            ydot[s] = state[5]
            zdot[s] = state[6]
            t[s] = getTimeByIndex(arc, s)
        end
        arcs[a] = CR3BPTraj(x, y, z, xdot, ydot, zdot, t, getJacobiConstant(manifold.periodicOrbit.dynamicsModel, getStateByIndex(arc, 1)))
    end
    man = CR3BPMan(orb, arcs, manifold.TOF)
    MATLAB.put_variable(file, name, man)
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
    exportPseudoManifold(pseudoManifold, file, name)

Export pseudo-manifold data to MAT file

# Arguments
- `pseudoManifold::BCR4BPPseudoManifold`: Pseudo-manifold object
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportPseudoManifold(pseudoManifold::MBD.BCR4BPPseudoManifold, file::MATLAB.MatFile, name::Symbol)
    propagator = MBD.Propagator()
    orbitArc::MBD.CR3BPArc = propagate(propagator, pseudoManifold.periodicOrbit.initialCondition, [0, pseudoManifold.periodicOrbit.period], pseudoManifold.periodicOrbit.dynamicsModel)
    nStates::Int64 = getStateCount(orbitArc)
    orbitx::Vector{Float64} = zeros(Float64, nStates)
    orbity::Vector{Float64} = zeros(Float64, nStates)
    orbitz::Vector{Float64} = zeros(Float64, nStates)
    orbitxdot::Vector{Float64} = zeros(Float64, nStates)
    orbitydot::Vector{Float64} = zeros(Float64, nStates)
    orbitzdot::Vector{Float64} = zeros(Float64, nStates)
    orbitt::Vector{Float64} = zeros(Float64, nStates)
    for s::Int64 in 1:nStates
        state::Vector{Float64} = getStateByIndex(orbitArc, s)
        orbitx[s] = state[1]
        orbity[s] = state[2]
        orbitz[s] = state[3]
        orbitxdot[s] = state[4]
        orbitydot[s] = state[5]
        orbitzdot[s] = state[6]
        orbitt[s] = getTimeByIndex(orbitArc, s)
    end
    orb = CR3BPOrb(orbitx, orbity, orbitz, orbitxdot, orbitydot, orbitzdot, orbitt, pseudoManifold.periodicOrbit.period, getJacobiConstant(pseudoManifold.periodicOrbit), getStabilityIndex(pseudoManifold.periodicOrbit))
    numArcs::Int64 = length(pseudoManifold.initialConditions)
    arcs::Vector{BCR4BP12Traj} = Vector{BCR4BP12Traj}(undef, numArcs)
    for a = 1:numArcs
        arc::MBD.BCR4BP12Arc = propagate(propagator, real(pseudoManifold.initialConditions[a]), [0, pseudoManifold.TOF], pseudoManifold.dynamicsModel)
        numStates::Int64 = getStateCount(arc)
        x::Vector{Float64} = zeros(Float64, numStates)
        y::Vector{Float64} = zeros(Float64, numStates)
        z::Vector{Float64} = zeros(Float64, numStates)
        xdot::Vector{Float64} = zeros(Float64, numStates)
        ydot::Vector{Float64} = zeros(Float64, numStates)
        zdot::Vector{Float64} = zeros(Float64, numStates)
        theta4::Vector{Float64} = zeros(Float64, numStates)
        t::Vector{Float64} = zeros(Float64, numStates)
        H::Vector{Float64} = zeros(Float64, numStates)
        for s::Int64 in 1:numStates
            state::Vector{Float64} = getStateByIndex(arc, s)
            x[s] = state[1]
            y[s] = state[2]
            z[s] = state[3]
            xdot[s] = state[4]
            ydot[s] = state[5]
            zdot[s] = state[6]
            theta4[s] = state[7]
            t[s] = getTimeByIndex(arc, s)
            H[s] = getHamiltonian(pseudoManifold.dynamicsModel, state)
        end
        arcs[a] = BCR4BP12Traj(x, y, z, xdot, ydot, zdot, theta4, t, H)
    end
    pseudoMan = PseudoMan(orb, arcs, pseudoManifold.TOF)
    MATLAB.put_variable(file, name, pseudoMan)
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
