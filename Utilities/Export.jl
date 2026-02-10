"""
Export utility functions

Author: Jonathan LeFevre Richmond
C: 2/19/25
U: 2/9/26
"""

using MBD, CSV, DataFrames, DifferentialEquations, LinearAlgebra, MATLAB

export CSVExportCR3BPFamily, exportArrays, exportBCR4BP12Manifold, exportBCR4BP12Orbit
export exportBCR4BP12Trajectory, exportBCR4BP41Trajectory, exportCR3BPManifold, exportCR3BPMMAT
export exportCR3BPOrbit, exportCR3BPTrajectory, exportInertialTrajectory, exportPseudoManifold
export fullExportCR3BPFamily, MATExportCR3BPOrbFamily, MATExportCR3BPTrajFamily

"""
    Arrays(vectors)

Arrays export object

# Arguments
- `vectors::Vector{Vector{Float64}}`: Vectors
"""
struct Arrays
    vectors::Vector{Vector{Float64}}                                    # Vectors

    function Arrays(vectors::Vector{Vector{Float64}})
        this = new(vectors)

        return this
    end
end

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
    BCR4BP12Traj(x, y, z, xdot, ydot, zdot, theta4, t, H, TOF)

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
- `TOF::Float64`: Time-of-flight [ndim]
"""
struct BCR4BP12Traj
    H::Vector{Float64}                                                  # Hamiltonian data [ndim]
    TOF::Float64                                                        # Time-of-flight [ndim]
    t::Vector{Float64}                                                  # Time data [ndim]
    theta4::Vector{Float64}                                             # theta4 data [ndim]
    x::Vector{Float64}                                                  # x data [ndim]
    xdot::Vector{Float64}                                               # xdot data [ndim]
    y::Vector{Float64}                                                  # y data [ndim]
    ydot::Vector{Float64}                                               # ydot data [ndim]
    z::Vector{Float64}                                                  # z data [ndim]
    zdot::Vector{Float64}                                               # zdot data [ndim]

    function BCR4BP12Traj(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, theta4::Vector{Float64}, t::Vector{Float64}, H::Vector{Float64}, TOF::Float64)
        this = new(H, TOF, t, theta4, x, xdot, y, ydot, z, zdot)

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
    BCR4BP41Traj(x, y, z, xdot, ydot, zdot, theta2, t, H, TOF)

BCR4BP P4-B1 trajectory export object

# Arguments
- `x::Vector{Float64}`: x data [ndim]
- `y::Vector{Float64}`: y data [ndim]
- `z::Vector{Float64}`: z data [ndim]
- `xdot::Vector{Float64}`: xdot data [ndim]
- `ydot::Vector{Float64}`: ydot data [ndim]
- `zdot::Vector{Float64}`: zdot data [ndim]
- `theta2::Vector{Float64}`: theta2 data [ndim]
- `t::Vector{Float64}`: Time data [ndim]
- `H::Vector{Float64}`: Hamiltonian data [ndim]
- `TOF::Float64`: Time-of-flight [ndim]
"""
struct BCR4BP41Traj
    H::Vector{Float64}                                                  # Hamiltonian data [ndim]
    TOF::Float64                                                        # Time-of-flight [ndim]
    t::Vector{Float64}                                                  # Time data [ndim]
    theta2::Vector{Float64}                                             # theta2 data [ndim]
    x::Vector{Float64}                                                  # x data [ndim]
    xdot::Vector{Float64}                                               # xdot data [ndim]
    y::Vector{Float64}                                                  # y data [ndim]
    ydot::Vector{Float64}                                               # ydot data [ndim]
    z::Vector{Float64}                                                  # z data [ndim]
    zdot::Vector{Float64}                                               # zdot data [ndim]

    function BCR4BP41Traj(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, theta2::Vector{Float64}, t::Vector{Float64}, H::Vector{Float64}, TOF::Float64)
        this = new(H, TOF, t, theta2, x, xdot, y, ydot, z, zdot)

        return this
    end
end

"""
    Conic(dynamicsModel, oe, TOF)

Keplerian conic export object

# Arguments
- `dynamicsModel::KDynamicsModel`: Keplerian dynamics model object
- `oe::Vector{Float64}`: Initial orbital elements [dim]
- `TOF::Float64`: Time-of-flight [s]
"""
struct Conic
    a::Float64                                                          # Semi-major axis [km]
    e::Float64                                                          # Eccentricity [ndim]
    i::Float64                                                          # Inclination [rad]
    Omega::Float64                                                      # RAAN [rad]
    omega::Float64                                                      # Argument of periapsis [rad]
    state::Vector{Float64}                                              # Cartesian state [dim]
    theta_0::Float64                                                    # Initial true anomaly [rad]
    TOF::Float64                                                        # Time-of-flight [s]

    function Conic(dynamicsModel::MBD.KDynamicsModel, oe::Vector{Float64}, TOF::Float64)
        this = new(oe[1], oe[2], oe[3], oe[4], oe[5], getCartesianState(dynamicsModel, oe), oe[6], TOF)

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
    CR3BPTraj(x, y, z, xdot, ydot, zdot, t, TOF, JC)

CR3BP trajectory export object

# Arguments
- `x::Vector{Float64}`: x data [ndim]
- `y::Vector{Float64}`: y data [ndim]
- `z::Vector{Float64}`: z data [ndim]
- `xdot::Vector{Float64}`: xdot data [ndim]
- `ydot::Vector{Float64}`: ydot data [ndim]
- `zdot::Vector{Float64}`: zdot data [ndim]
- `t::Vector{Float64}`: Time data [ndim]
- `TOF::Float64`: Time-of-flight [ndim]
- `JC::Float64`: Jacobi constant
"""
struct CR3BPTraj
    JC::Float64                                                         # Jacobi constant
    TOF::Float64                                                        # Time-of-flight [ndim]
    t::Vector{Float64}                                                  # Time data [ndim]
    x::Vector{Float64}                                                  # x data [ndim]
    xdot::Vector{Float64}                                               # xdot data [ndim]
    y::Vector{Float64}                                                  # y data [ndim]
    ydot::Vector{Float64}                                               # ydot data [ndim]
    z::Vector{Float64}                                                  # z data [ndim]
    zdot::Vector{Float64}                                               # zdot data [ndim]

    function CR3BPTraj(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, t::Vector{Float64}, TOF::Float64, JC::Float64)
        this = new(JC, TOF, t, x, xdot, y, ydot, z, zdot)

        return this
    end
end

"""
    CR3BPMMAT(initialEpoch, t_0, theta_dep_0, theta_arr_f, theta_arr_0, departureOrbit, departureManifoldArc1, departureManifoldArc2, departureConic, bridgeConic, arrivalConic, arrivalManifoldArc, arrivalOrbit, Deltav_1, Deltav_2, TOF)

CR3BP MMAT export object

# Arguments
- `initialEpoch::String`: Initial epoch
- `t_0::Float64`: Departure time from initial epoch [s]
- `theta_dep_0::Float64`: Initial departure body true anomaly [rad]
- `theta_arr_f::Float64`: Final arrival body true anomaly [rad]
- `theta_arr_0::Float64`: Initial arrival body true anomaly [rad]
- `departureOrbit::CR3BPOrb`: Departure CR3BP orbit object
- `departureManifoldArc1::CR3BPTraj`: Departure CR3BP manifold arc trajectory object
- `departureManifoldArc2::CR3BPTraj`: Intermediate CR3BP manifold arc trajectory object
- `departureConic::Conic`: Departure Keplerian conic object
- `bridgeConic::Conic`: Bridge Keplerian conic object
- `arrivalConic::Conic`: Arrival Keplerian conic object
- `arrivalManifoldArc::CR3BPTraj`: Arrival CR3BP manifold arc trajectory object
- `arrivalOrbit::CR3BPOrb`: Arrival CR3BP orbit object
- `Deltav_1::Float64`: First maneuver magnitude [km/s]
- `Deltav_2::Float64`: Second maneuver magnitude [km/s]
- `TOF::Float64`: Total transfer time-of-flight [s]
"""
struct CR3BPMMAT
    arrivalConic::Conic                                                 # Arrival conic object
    arrivalManifoldArc::CR3BPTraj                                       # Arrival manifold arc trajectory object
    arrivalOrbit::CR3BPOrb                                              # Arrival orbit object
    bridgeConic::Conic                                                  # Bridge conic object
    Deltav_1::Float64                                                   # First maneuver magnitude [km/s]
    Deltav_2::Float64                                                   # Second maneuver magnitude [km/s]
    departureConic::Conic                                               # Departure conic object
    departureManifoldArc1::CR3BPTraj                                    # Departure manifold arc trajectory object
    departureManifoldArc2::CR3BPTraj                                    # Intermediate manifold arc trajectory object
    departureOrbit::CR3BPOrb                                            # Departure orbit object
    initialEpoch::String                                                # Initial epoch
    theta_arr_f::Float64                                                # Final arrival body true anomaly [rad]
    theta_arr_0::Float64                                                # Initial arrival body true anomlay [rad]
    theta_dep_0::Float64                                                # Initial departure body true anomaly [rad]
    TOF::Float64                                                        # Transfer time-of-flight [s]
    t_0::Float64                                                        # Departure time from initial epoch [s]

    function CR3BPMMAT(initialEpoch::String, t_0::Float64, theta_dep_0::Float64, theta_arr_f::Float64, theta_arr_0::Float64, departureOrbit::CR3BPOrb, departureManifoldArc1::CR3BPTraj, departureManifoldArc2::CR3BPTraj, departureConic::Conic, bridgeConic::Conic, arrivalConic::Conic, arrivalManifoldArc::CR3BPTraj, arrivalOrbit::CR3BPOrb, Deltav_1::Float64, Deltav_2::Float64, TOF::Float64)
        this = new(arrivalConic, arrivalManifoldArc, arrivalOrbit, bridgeConic, Deltav_1, Deltav_2, departureConic, departureManifoldArc1, departureManifoldArc2, departureOrbit, initialEpoch, theta_arr_f, theta_arr_0, theta_dep_0, TOF, t_0)

        return this
    end
end

"""
    CR3BPOrbFam(orbits)

CR3BP obit family export object

# Arguments
- `orbits::Vector{CR3BPOrb}`: CR3BP periodic orbit export objects
"""
struct CR3BPOrbFam
    orbits::Vector{CR3BPOrb}                                            # Periodic orbit export objects

    function CR3BPOrbFam(orbits::Vector{CR3BPOrb})
        this = new(orbits)

        return this
    end
end

"""
    CR3BPTrajFam(trajs)

CR3BP trajectory family export object

# Arguments
- `trajs::Vector{CR3BPTraj}`: CR3BP trajectory export objects
"""
struct CR3BPTrajFam
    trajs::Vector{CR3BPTraj}                                            # Trajectory export objects

    function CR3BPTrajFam(trajs::Vector{CR3BPTraj})
        this = new(trajs)

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
    InertialTraj(x, y, z, xdot, ydot, zdot, t, TOF)

Inertial trajectory export object

# Arguments
- `x::Vector{Float64}`: x data [ndim]
- `y::Vector{Float64}`: y data [ndim]
- `z::Vector{Float64}`: z data [ndim]
- `xdot::Vector{Float64}`: xdot data [ndim]
- `ydot::Vector{Float64}`: ydot data [ndim]
- `zdot::Vector{Float64}`: zdot data [ndim]
- `t::Vector{Float64}`: Time data [ndim]
- `TOF::Float64`: Time-of-flight [ndim]
"""
struct InertialTraj
    TOF::Float64                                                        # Time-of-flight [ndim]
    t::Vector{Float64}                                                  # Time data [ndim]
    x::Vector{Float64}                                                  # x data [ndim]
    xdot::Vector{Float64}                                               # xdot data [ndim]
    y::Vector{Float64}                                                  # y data [ndim]
    ydot::Vector{Float64}                                               # ydot data [ndim]
    z::Vector{Float64}                                                  # z data [ndim]
    zdot::Vector{Float64}                                               # zdot data [ndim]

    function InertialTraj(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, t::Vector{Float64}, TOF::Float64)
        this = new(TOF, t, x, xdot, y, ydot, z, zdot)

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
    CSVExportCR3BPFamily(family, file)

Export CR3BP family data to CSV file

# Arguments
- `family::CR3BPMSOrbitFamily`: CR3BP multiple shooter periodic orbit family object
- `file::String`: CSV file name
"""
function CSVExportCR3BPFamily(family::MBD.CR3BPMSOrbitFamily, file::String)
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
        orbit::MBD.CR3BPMSPeriodicOrbit = getMember(family, o)
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
    CSVExportCR3BPFamily(family, file)

Export CR3BP family data to CSV file

# Arguments
- `family::CR3BPContinuationFamily`: CR3BP trajectory continuation family object
- `file::String`: CSV file name
"""
function CSVExportCR3BPFamily(family::MBD.CR3BPContinuationFamily, file::String)
    nMem::Int64 = getNumMembers(family)
    x::Vector{Float64} = Vector{Float64}(undef, nMem)
    y::Vector{Float64} = Vector{Float64}(undef, nMem)
    z::Vector{Float64} = Vector{Float64}(undef, nMem)
    xdot::Vector{Float64} = Vector{Float64}(undef, nMem)
    ydot::Vector{Float64} = Vector{Float64}(undef, nMem)
    zdot::Vector{Float64} = Vector{Float64}(undef, nMem)
    TOF::Vector{Float64} = Vector{Float64}(undef, nMem)
    JC::Vector{Float64} = Vector{Float64}(undef, nMem)
    for t::Int64 = 1:nMem
        solution = MBD.CR3BPMultipleShooterProblem()
        numNodes::Int64 = length(family.nodes[t])
        solution.nodes = [MBD.shallowClone(family.nodes[t][n]) for n = 1:numNodes]
        solution.segments = [MBD.shallowClone(family.segments[t][s]) for s = 1:numNodes-1]
        x[t] = solution.nodes[1].state.data[1]
        y[t] = solution.nodes[1].state.data[2]
        z[t] = solution.nodes[1].state.data[3]
        xdot[t] = solution.nodes[1].state.data[4]
        ydot[t] = solution.nodes[1].state.data[5]
        zdot[t] = solution.nodes[1].state.data[6]
        TOF[t] = 0.0
        for s::Int64 = 1:numNodes-1
            TOF[t] += solution.segments[s].TOF.data[1]
        end 
        JC[t] = getJacobiConstant(solution.nodes[1].dynamicsModel, solution.nodes[1].state.data[1:6])
    end
    familyData::DataFrames.DataFrame = DataFrames.DataFrame("x" => x, "y" => y, "z" => z, "xdot" => xdot, "ydot" => ydot, "zdot" => zdot, "TOF" => TOF, "JC" => JC)
    CSV.write(file, familyData)
end

"""
    exportArrays(vectors, MATFile, name)

Export vectors data to MAT file

# Arguments
- `vectors::Vector{Vector{Float64}}`: Vectors
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportArrays(vectors::Vector{Vector{Float64}}, file::MATLAB.MatFile, name::Symbol)
    arrays = Arrays(vectors)
    MATLAB.put_variable(file, name, arrays)
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
    exportBCR4BP12Manifold(manifold, manifoldArcs, file, name)

Export BCR4BP P1-P2 manifold data to MAT file

# Arguments
- `manifold::BCR4BP12Manifold`: BCR4BP P1-P2 manifold object
- `manifoldArcs::Vector{BCR4BP12ManifoldArc}`: BCR4BP P1-P2 manifold arcs
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportBCR4BP12Manifold(manifold::MBD.BCR4BP12Manifold, manifoldArcs::Vector{MBD.BCR4BP12ManifoldArc}, file::MATLAB.MatFile, name::Symbol)
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
    numArcs::Int64 = length(manifoldArcs)
    arcs::Vector{BCR4BP12Traj} = Vector{BCR4BP12Traj}(undef, numArcs)
    for a = 1:numArcs
        arc::MBD.BCR4BP12Arc = propagate(propagator, real(manifoldArcs[a].initialCondition), [0, manifoldArcs[a].TOF], manifold.periodicOrbit.dynamicsModel)
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
        arcs[a] = BCR4BP12Traj(x, y, z, xdot, ydot, zdot, theta4, t, H, manifoldArcs[a].TOF)
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
    orbitStates::Vector{Vector{Float64}} = []
    orbitEpochs::Vector{Float64} = []
    for n::Int64 in 1:length(orbit.nodeEpochs)-1
        arc::MBD.BCR4BP12Arc = propagate(propagator, orbit.nodeStates[n], [orbit.nodeEpochs[n], orbit.nodeEpochs[n+1]], orbit.dynamicsModel)
        append!(orbitStates, arc.states)
        append!(orbitEpochs, arc.times)
    end
    nStates::Int64 = length(orbitEpochs)
    x::Vector{Float64} = zeros(Float64, nStates)
    y::Vector{Float64} = zeros(Float64, nStates)
    z::Vector{Float64} = zeros(Float64, nStates)
    xdot::Vector{Float64} = zeros(Float64, nStates)
    ydot::Vector{Float64} = zeros(Float64, nStates)
    zdot::Vector{Float64} = zeros(Float64, nStates)
    theta4::Vector{Float64} = zeros(Float64, nStates)
    t::Vector{Float64} = zeros(Float64, nStates)
    H::Vector{Float64} = zeros(Float64, nStates)
    for s::Int64 in 1:nStates
        state::Vector{Float64} = orbitStates[s]
        x[s] = state[1]
        y[s] = state[2]
        z[s] = state[3]
        xdot[s] = state[4]
        ydot[s] = state[5]
        zdot[s] = state[6]
        theta4[s] = state[7]
        t[s] = orbitEpochs[s]
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
    exportBCR4BP12Trajectory(x0, y0, z0, xdot0, ydot0, zdot0, theta40, propTime, dynamicsModel, file, name)

Export BCR4BP P1-P2 trajectory data to MAT file

# Arguments
- `x0::Float64`: Initial x-position [ndim]
- `y0::Float64`: Initial y-position [ndim]
- `z0::Float64`: Initial z-position [ndim]
- `xdot0::Float64`: Initial x-velocity [ndim]
- `ydot0::Float64`: Initial y-velocity [ndim]
- `zdot0::Float64`: Initial z-velocity [ndim]
- `theta40::Float64`: Initial P4 angle [ndim]
- `propTime::Float64`: Propagation time [ndim]
- `dynamicsModel::BCR4BP12DynamicsModel`: BCR4BP P1-P2 dynamics model object
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportBCR4BP12Trajectory(x0::Float64, y0::Float64, z0::Float64, xdot0::Float64, ydot0::Float64, zdot0::Float64, theta40::Float64, propTime::Float64, dynamicsModel::MBD.BCR4BP12DynamicsModel, file::MATLAB.MatFile, name::Symbol)
    propagator = MBD.Propagator()
    arc::MBD.BCR4BP12Arc = propagate(propagator, [x0, y0, z0, xdot0, ydot0, zdot0, theta40], [0, propTime], dynamicsModel)
    nStates::Int64 = getStateCount(arc)
    x::Vector{Float64} = zeros(Float64, nStates)
    y::Vector{Float64} = zeros(Float64, nStates)
    z::Vector{Float64} = zeros(Float64, nStates)
    xdot::Vector{Float64} = zeros(Float64, nStates)
    ydot::Vector{Float64} = zeros(Float64, nStates)
    zdot::Vector{Float64} = zeros(Float64, nStates)
    theta4::Vector{Float64} = zeros(Float64, nStates)
    t::Vector{Float64} = zeros(Float64, nStates)
    H::Vector{Float64} = zeros(Float64, nStates)
    for s::Int64 in 1:nStates
        state::Vector{Float64} = getStateByIndex(arc, s)
        x[s] = state[1]
        y[s] = state[2]
        z[s] = state[3]
        xdot[s] = state[4]
        ydot[s] = state[5]
        zdot[s] = state[6]
        theta4[s] = state[7]
        t[s] = getTimeByIndex(arc, s)
        H[s] = getHamiltonian(dynamicsModel, state)
    end
    exportBCR4BP12Trajectory(x, y, z, xdot, ydot, zdot, theta4, t, H, propTime, file, name)
end

"""
    exportBCR4BP12Trajectory(x, y, z, xdot, ydot, zdot, theta4, t, H, TOF, file, name)

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
- `TOF::Float64`: Time-of-flight [ndim]
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportBCR4BP12Trajectory(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, theta4::Vector{Float64}, t::Vector{Float64}, H::Vector{Float64}, TOF::Float64, file::MATLAB.MatFile, name::Symbol)
    traj = BCR4BP12Traj(x, y, z, xdot, ydot, zdot, theta4, t, H, TOF)
    MATLAB.put_variable(file, name, traj)
end

"""
    exportBCR4BP41Trajectory(x0, y0, z0, xdot0, ydot0, zdot0, theta20, propTime, dynamicsModel, file, name)

Export BCR4BP P4-B1 trajectory data to MAT file

# Arguments
- `x0::Float64`: Initial x-position [ndim]
- `y0::Float64`: Initial y-position [ndim]
- `z0::Float64`: Initial z-position [ndim]
- `xdot0::Float64`: Initial x-velocity [ndim]
- `ydot0::Float64`: Initial y-velocity [ndim]
- `zdot0::Float64`: Initial z-velocity [ndim]
- `theta20::Float64`: Initial P2 angle [ndim]
- `propTime::Float64`: Propagation time [ndim]
- `dynamicsModel::BCR4BP41DynamicsModel`: BCR4BP P4-B1 dynamics model object
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportBCR4BP41Trajectory(x0::Float64, y0::Float64, z0::Float64, xdot0::Float64, ydot0::Float64, zdot0::Float64, theta20::Float64, propTime::Float64, dynamicsModel::MBD.BCR4BP41DynamicsModel, file::MATLAB.MatFile, name::Symbol)
    propagator = MBD.Propagator()
    arc::MBD.BCR4BP41Arc = propagate(propagator, [x0, y0, z0, xdot0, ydot0, zdot0, theta20], [0, propTime], dynamicsModel)
    nStates::Int64 = getStateCount(arc)
    x::Vector{Float64} = zeros(Float64, nStates)
    y::Vector{Float64} = zeros(Float64, nStates)
    z::Vector{Float64} = zeros(Float64, nStates)
    xdot::Vector{Float64} = zeros(Float64, nStates)
    ydot::Vector{Float64} = zeros(Float64, nStates)
    zdot::Vector{Float64} = zeros(Float64, nStates)
    theta2::Vector{Float64} = zeros(Float64, nStates)
    t::Vector{Float64} = zeros(Float64, nStates)
    H::Vector{Float64} = zeros(Float64, nStates)
    for s::Int64 in 1:nStates
        state::Vector{Float64} = getStateByIndex(arc, s)
        x[s] = state[1]
        y[s] = state[2]
        z[s] = state[3]
        xdot[s] = state[4]
        ydot[s] = state[5]
        zdot[s] = state[6]
        theta2[s] = state[7]
        t[s] = getTimeByIndex(arc, s)
        H[s] = getHamiltonian(dynamicsModel, state)
    end
    exportBCR4BP41Trajectory(x, y, z, xdot, ydot, zdot, theta2, t, H, t[end]-t[1], file, name)
end

"""
    exportBCR4BP41Trajectory(x, y, z, xdot, ydot, zdot, theta2, t, H, TOF, file, name)

Export BCR4BP P4-B1 trajectory data to MAT file

# Arguments
- `x::Vector{Float64}`: x data [ndim]
- `y::Vector{Float64}`: y data [ndim]
- `z::Vector{Float64}`: z data [ndim]
- `xdot::Vector{Float64}`: xdot data [ndim]
- `ydot::Vector{Float64}`: ydot data [ndim]
- `zdot::Vector{Float64}`: zdot data [ndim]
- `theta2::Vector{Float64}`: theta2 data [ndim]
- `t::Vector{Float64}`: Time data [ndim]
- `H::Vector{Float64}`: Hamiltonian data [ndim]
- `TOF::Float64`: Time-of-flight [ndim]
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportBCR4BP41Trajectory(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, theta2::Vector{Float64}, t::Vector{Float64}, H::Vector{Float64}, TOF::Float64, file::MATLAB.MatFile, name::Symbol)
    traj = BCR4BP41Traj(x, y, z, xdot, ydot, zdot, theta2, t, H, TOF)
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
        arcs[a] = CR3BPTraj(x, y, z, xdot, ydot, zdot, t, manifold.TOF, getJacobiConstant(manifold.periodicOrbit.dynamicsModel, getStateByIndex(arc, 1)))
    end
    man = CR3BPMan(orb, arcs, manifold.TOF)
    MATLAB.put_variable(file, name, man)
end

"""
    exportCR3BPManifold(manifold, manifoldArcs, file, name)

Export CR3BP manifold data to MAT file

# Arguments
- `manifold::CR3BPManifold`: CR3BP manifold object
- `manifoldArcs::Vector{CR3BPManifold}`: CR3BP manifold arcs
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportCR3BPManifold(manifold::MBD.CR3BPManifold, manifoldArcs::Vector{MBD.CR3BPManifoldArc}, file::MATLAB.MatFile, name::Symbol)
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
    numArcs::Int64 = length(manifoldArcs)
    arcs::Vector{CR3BPTraj} = Vector{CR3BPTraj}(undef, numArcs)
    for a = 1:numArcs
        arc::MBD.CR3BPArc = propagate(propagator, real(manifoldArcs[a].initialCondition), [0, manifoldArcs[a].TOF], manifold.periodicOrbit.dynamicsModel)
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
        arcs[a] = CR3BPTraj(x, y, z, xdot, ydot, zdot, t, manifoldArcs[a].TOF, getJacobiConstant(manifold.periodicOrbit.dynamicsModel, getStateByIndex(arc, 1)))
    end
    man = CR3BPMan(orb, arcs, manifold.TOF)
    MATLAB.put_variable(file, name, man)
end

"""
    exportCR3BPMMAT(env, t_0, intersect, theta_dep_0, theta_arr_f, theta_arr_0, departureManifoldArc, intermediateManifoldArc, oe_dep_SoI, t_depConic, oe_bridge_peri, arrivalManifoldArc, Deltav_1, TOF, file, name)

Export CR3BP MMAT data to MAT file

# Arguments
- `env::MMATEnv`: MMAT environment object
- `t_0::Float64`: Departure time from initial epoch [s]
- `intersect::Vector{Float64}`: Intersect data
- `theta_dep_0::Float64`: Initial departure body true anomaly [rad]
- `theta_arr_f::Float64`: Final arrival body true anomaly [rad]
- `theta_arr_0::Float64`: Initial arrival body true anomaly [rad]
- `departureManifoldArc::CR3BPManifoldArc`: Departure CR3BP manifold arc object
- `intermediateManifoldArc::CR3BPArc`: Intermediate CR3BP arc object
- `oe_dep_SoI::Vector{Float64}`: Dearture conic orbital elements at SoI [dim]
- `t_depConic::Float64`: Departure conic TOF [s]
- `oe_bridge_peri::Vector{Float64}`: Bridge conic orbital elements at periapsis [dim]
- `arrivalManifoldArc::CR3BPManifoldArc`: Arrival CR3BP manifold arc object
- `Deltav_1::Float64`: First maneuver magnitude [km/s]
- `TOF::Float64`: Total transfer time-fo-flight [s]
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportCR3BPMMAT(env, t_0::Float64, intersect::Vector{Float64}, theta_dep_0::Float64, theta_arr_f::Float64, theta_arr_0::Float64, departureManifoldArc::MBD.CR3BPManifoldArc, intermediateManifoldArc::MBD.CR3BPArc, oe_dep_SoI::Vector{Float64}, t_depConic::Float64, oe_bridge_peri::Vector{Float64}, arrivalManifoldArc::MBD.CR3BPManifoldArc, Deltav_1::Float64, TOF::Float64, file::MATLAB.MatFile, name::Symbol)
    departureOrbitArc::MBD.CR3BPArc = propagate(env.propagator, departureManifoldArc.periodicOrbit.initialCondition, [0, departureManifoldArc.periodicOrbit.period], departureManifoldArc.periodicOrbit.dynamicsModel)
    departureOrbitnStates::Int64 = getStateCount(departureOrbitArc)
    departureOrbitx::Vector{Float64} = zeros(Float64, departureOrbitnStates)
    departureOrbity::Vector{Float64} = zeros(Float64, departureOrbitnStates)
    departureOrbitz::Vector{Float64} = zeros(Float64, departureOrbitnStates)
    departureOrbitxdot::Vector{Float64} = zeros(Float64, departureOrbitnStates)
    departureOrbitydot::Vector{Float64} = zeros(Float64, departureOrbitnStates)
    departureOrbitzdot::Vector{Float64} = zeros(Float64, departureOrbitnStates)
    departureOrbitt::Vector{Float64} = zeros(Float64, departureOrbitnStates)
    for s::Int64 in 1:departureOrbitnStates
        state::Vector{Float64} = getStateByIndex(departureOrbitArc, s)
        departureOrbitx[s] = state[1]
        departureOrbity[s] = state[2]
        departureOrbitz[s] = state[3]
        departureOrbitxdot[s] = state[4]
        departureOrbitydot[s] = state[5]
        departureOrbitzdot[s] = state[6]
        departureOrbitt[s] = getTimeByIndex(departureOrbitArc, s)
    end
    departureJC::Float64 = getJacobiConstant(departureManifoldArc)
    departureOrbitvarsig::Float64 = getStabilityIndex(departureManifoldArc.periodicOrbit)
    departureOrb = CR3BPOrb(departureOrbitx, departureOrbity, departureOrbitz, departureOrbitxdot, departureOrbitydot, departureOrbitzdot, departureOrbitt, departureManifoldArc.periodicOrbit.period, departureJC, departureOrbitvarsig)

    departureTrajArc::MBD.CR3BPArc = propagate(env.propagator, real(departureManifoldArc.initialCondition), [0, departureManifoldArc.TOF], departureManifoldArc.periodicOrbit.dynamicsModel)
    departureTrajnStates::Int64 = getStateCount(departureTrajArc)
    departureTrajx::Vector{Float64} = zeros(Float64, departureTrajnStates)
    departureTrajy::Vector{Float64} = zeros(Float64, departureTrajnStates)
    departureTrajz::Vector{Float64} = zeros(Float64, departureTrajnStates)
    departureTrajxdot::Vector{Float64} = zeros(Float64, departureTrajnStates)
    departureTrajydot::Vector{Float64} = zeros(Float64, departureTrajnStates)
    departureTrajzdot::Vector{Float64} = zeros(Float64, departureTrajnStates)
    departureTrajt::Vector{Float64} = zeros(Float64, departureTrajnStates)
    for s::Int64 in 1:departureTrajnStates
        state::Vector{Float64} = getStateByIndex(departureTrajArc, s)
        departureTrajx[s] = state[1]
        departureTrajy[s] = state[2]
        departureTrajz[s] = state[3]
        departureTrajxdot[s] = state[4]
        departureTrajydot[s] = state[5]
        departureTrajzdot[s] = state[6]
        departureTrajt[s] = getTimeByIndex(departureTrajArc, s)
    end
    departureManifoldTraj1 = CR3BPTraj(departureTrajx, departureTrajy, departureTrajz, departureTrajxdot, departureTrajydot, departureTrajzdot, departureTrajt, departureManifoldArc.TOF, departureJC)

    intermediateTrajnStates::Int64 = getStateCount(intermediateManifoldArc)
    intermediateTrajx::Vector{Float64} = zeros(Float64, intermediateTrajnStates)
    intermediateTrajy::Vector{Float64} = zeros(Float64, intermediateTrajnStates)
    intermediateTrajz::Vector{Float64} = zeros(Float64, intermediateTrajnStates)
    intermediateTrajxdot::Vector{Float64} = zeros(Float64, intermediateTrajnStates)
    intermediateTrajydot::Vector{Float64} = zeros(Float64, intermediateTrajnStates)
    intermediateTrajzdot::Vector{Float64} = zeros(Float64, intermediateTrajnStates)
    intermediateTrajt::Vector{Float64} = zeros(Float64, intermediateTrajnStates)
    for s::Int64 in 1:intermediateTrajnStates
        state::Vector{Float64} = getStateByIndex(intermediateManifoldArc, s)
        intermediateTrajx[s] = state[1]
        intermediateTrajy[s] = state[2]
        intermediateTrajz[s] = state[3]
        intermediateTrajxdot[s] = state[4]
        intermediateTrajydot[s] = state[5]
        intermediateTrajzdot[s] = state[6]
        intermediateTrajt[s] = getTimeByIndex(intermediateManifoldArc, s)
    end
    intermediateTrajTOF::Float64 = getTimeByIndex(intermediateManifoldArc, -1)-getTimeByIndex(intermediateManifoldArc, 1)
    intermediateTrajJC::Float64 = getJacobiConstant(intermediateManifoldArc.dynamicsModel, getStateByIndex(intermediateManifoldArc, 1))
    departureManifoldTraj2 = CR3BPTraj(intermediateTrajx, intermediateTrajy, intermediateTrajz, intermediateTrajxdot, intermediateTrajydot, intermediateTrajzdot, intermediateTrajt, intermediateTrajTOF, intermediateTrajJC)

    departureConic = Conic(env.SDynamicsModel, oe_dep_SoI, t_depConic)
    bridgeConic = Conic(env.SDynamicsModel, oe_bridge_peri, intersect[1])
    arrivalConic = Conic(env.SDynamicsModel, intersect[3:8], intersect[9])

    arrivalOrbitArc::MBD.CR3BPArc = propagate(env.propagator, arrivalManifoldArc.periodicOrbit.initialCondition, [0, arrivalManifoldArc.periodicOrbit.period], arrivalManifoldArc.periodicOrbit.dynamicsModel)
    arrivalOrbitnStates::Int64 = getStateCount(arrivalOrbitArc)
    arrivalOrbitx::Vector{Float64} = zeros(Float64, arrivalOrbitnStates)
    arrivalOrbity::Vector{Float64} = zeros(Float64, arrivalOrbitnStates)
    arrivalOrbitz::Vector{Float64} = zeros(Float64, arrivalOrbitnStates)
    arrivalOrbitxdot::Vector{Float64} = zeros(Float64, arrivalOrbitnStates)
    arrivalOrbitydot::Vector{Float64} = zeros(Float64, arrivalOrbitnStates)
    arrivalOrbitzdot::Vector{Float64} = zeros(Float64, arrivalOrbitnStates)
    arrivalOrbitt::Vector{Float64} = zeros(Float64, arrivalOrbitnStates)
    for s::Int64 in 1:arrivalOrbitnStates
        state::Vector{Float64} = getStateByIndex(arrivalOrbitArc, s)
        arrivalOrbitx[s] = state[1]
        arrivalOrbity[s] = state[2]
        arrivalOrbitz[s] = state[3]
        arrivalOrbitxdot[s] = state[4]
        arrivalOrbitydot[s] = state[5]
        arrivalOrbitzdot[s] = state[6]
        arrivalOrbitt[s] = getTimeByIndex(arrivalOrbitArc, s)
    end
    arrivalJC::Float64 = getJacobiConstant(arrivalManifoldArc)
    arrivalOrbitvarsig::Float64 = getStabilityIndex(arrivalManifoldArc.periodicOrbit)
    arrivalOrb = CR3BPOrb(arrivalOrbitx, arrivalOrbity, arrivalOrbitz, arrivalOrbitxdot, arrivalOrbitydot, arrivalOrbitzdot, arrivalOrbitt, arrivalManifoldArc.periodicOrbit.period, arrivalJC, arrivalOrbitvarsig)

    arrivalTrajArc::MBD.CR3BPArc = propagate(env.propagator, real(arrivalManifoldArc.initialCondition), [0, arrivalManifoldArc.TOF], arrivalManifoldArc.periodicOrbit.dynamicsModel)
    arrivalTrajnStates::Int64 = getStateCount(arrivalTrajArc)
    arrivalTrajx::Vector{Float64} = zeros(Float64, arrivalTrajnStates)
    arrivalTrajy::Vector{Float64} = zeros(Float64, arrivalTrajnStates)
    arrivalTrajz::Vector{Float64} = zeros(Float64, arrivalTrajnStates)
    arrivalTrajxdot::Vector{Float64} = zeros(Float64, arrivalTrajnStates)
    arrivalTrajydot::Vector{Float64} = zeros(Float64, arrivalTrajnStates)
    arrivalTrajzdot::Vector{Float64} = zeros(Float64, arrivalTrajnStates)
    arrivalTrajt::Vector{Float64} = zeros(Float64, arrivalTrajnStates)
    for s::Int64 in 1:arrivalTrajnStates
        state::Vector{Float64} = getStateByIndex(arrivalTrajArc, s)
        arrivalTrajx[s] = state[1]
        arrivalTrajy[s] = state[2]
        arrivalTrajz[s] = state[3]
        arrivalTrajxdot[s] = state[4]
        arrivalTrajydot[s] = state[5]
        arrivalTrajzdot[s] = state[6]
        arrivalTrajt[s] = getTimeByIndex(arrivalTrajArc, s)
    end
    arrivalManifoldTraj = CR3BPTraj(arrivalTrajx, arrivalTrajy, arrivalTrajz, arrivalTrajxdot, arrivalTrajydot, arrivalTrajzdot, arrivalTrajt, arrivalManifoldArc.TOF, arrivalJC)

    transfer = CR3BPMMAT(env.initialEpoch, t_0, theta_dep_0, theta_arr_f, theta_arr_0, departureOrb, departureManifoldTraj1, departureManifoldTraj2, departureConic, bridgeConic, arrivalConic, arrivalManifoldTraj, arrivalOrb, Deltav_1, intersect[2], TOF)
    MATLAB.put_variable(file, name, transfer)
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
    exportCR3BPOrbit(orbit, file, name)

Export CR3BP multiple shooter orbit data to MAT file

# Arguments
- `orbit::CR3BPMSPeriodicOrbit`: CR3BP multiple shooter periodic orbit object
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportCR3BPOrbit(orbit::MBD.CR3BPMSPeriodicOrbit, file::MATLAB.MatFile, name::Symbol)
    propagator = MBD.Propagator()
    orbitStates::Vector{Vector{Float64}} = []
    orbitEpochs::Vector{Float64} = []
    for n::Int64 in 1:length(orbit.nodeEpochs)-1
        arc::MBD.CR3BPArc = propagate(propagator, orbit.nodeStates[n], [orbit.nodeEpochs[n], orbit.nodeEpochs[n+1]], orbit.dynamicsModel)
        append!(orbitStates, arc.states)
        append!(orbitEpochs, arc.times)
    end
    nStates::Int64 = length(orbitEpochs)
    x::Vector{Float64} = zeros(Float64, nStates)
    y::Vector{Float64} = zeros(Float64, nStates)
    z::Vector{Float64} = zeros(Float64, nStates)
    xdot::Vector{Float64} = zeros(Float64, nStates)
    ydot::Vector{Float64} = zeros(Float64, nStates)
    zdot::Vector{Float64} = zeros(Float64, nStates)
    t::Vector{Float64} = zeros(Float64, nStates)
    for s::Int64 in 1:nStates
        state::Vector{Float64} = orbitStates[s]
        x[s] = state[1]
        y[s] = state[2]
        z[s] = state[3]
        xdot[s] = state[4]
        ydot[s] = state[5]
        zdot[s] = state[6]
        t[s] = orbitEpochs[s]
    end
    JC::Float64 = getJacobiConstant(orbit)
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
    exportCR3BPTrajectory(solution, file, name)

Export CR3BP trajectory data to MAT file

# Arguments
- `solution::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportCR3BPTrajectory(solution::MBD.CR3BPMultipleShooterProblem, file::MATLAB.MatFile, name::Symbol)
    propagator = MBD.Propagator()
    states::Vector{Vector{Float64}} = []
    epochs::Vector{Float64} = []
    for n::Int64 in 1:length(solution.nodes)-1
        arc::MBD.CR3BPArc = propagate(propagator, solution.nodes[n].state.data[1:6], [solution.nodes[n].epoch.data[1], solution.nodes[n+1].epoch.data[1]], solution.nodes[1].dynamicsModel)
        append!(states, arc.states)
        append!(epochs, arc.times)
    end
    nStates::Int64 = length(epochs)
    x::Vector{Float64} = zeros(Float64, nStates)
    y::Vector{Float64} = zeros(Float64, nStates)
    z::Vector{Float64} = zeros(Float64, nStates)
    xdot::Vector{Float64} = zeros(Float64, nStates)
    ydot::Vector{Float64} = zeros(Float64, nStates)
    zdot::Vector{Float64} = zeros(Float64, nStates)
    t::Vector{Float64} = zeros(Float64, nStates)
    for s::Int64 in 1:nStates
        state::Vector{Float64} = states[s]
        x[s] = state[1]
        y[s] = state[2]
        z[s] = state[3]
        xdot[s] = state[4]
        ydot[s] = state[5]
        zdot[s] = state[6]
        t[s] = epochs[s]
    end
    JC::Float64 = getJacobiConstant(solution.nodes[1].dynamicsModel, solution.nodes[1].state.data[1:6])
    exportCR3BPTrajectory(x, y, z, xdot, ydot, zdot, t, t[end]-t[1], JC, file, name)
end

"""
    exportCR3BPTrajectory(x0, y0, z0, xdot0, ydot0, zdot0, propTime, dynamicsModel, file, name)

Export CR3BP trajectory data to MAT file

# Arguments
- `x0::Float64`: Initial x-position [ndim]
- `y0::Float64`: Initial y-position [ndim]
- `z0::Float64`: Initial z-position [ndim]
- `xdot0::Float64`: Initial x-velocity [ndim]
- `ydot0::Float64`: Initial y-velocity [ndim]
- `zdot0::Float64`: Initial z-velocity [ndim]
- `propTime::Float64`: Propagation time [ndim]
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportCR3BPTrajectory(x0::Float64, y0::Float64, z0::Float64, xdot0::Float64, ydot0::Float64, zdot0::Float64, propTime::Float64, dynamicsModel::MBD.CR3BPDynamicsModel, file::MATLAB.MatFile, name::Symbol)
    propagator = MBD.Propagator()
    arc::MBD.CR3BPArc = propagate(propagator, [x0, y0, z0, xdot0, ydot0, zdot0], [0, propTime], dynamicsModel)
    nStates::Int64 = getStateCount(arc)
    x::Vector{Float64} = zeros(Float64, nStates)
    y::Vector{Float64} = zeros(Float64, nStates)
    z::Vector{Float64} = zeros(Float64, nStates)
    xdot::Vector{Float64} = zeros(Float64, nStates)
    ydot::Vector{Float64} = zeros(Float64, nStates)
    zdot::Vector{Float64} = zeros(Float64, nStates)
    t::Vector{Float64} = zeros(Float64, nStates)
    for s::Int64 in 1:nStates
        state::Vector{Float64} = getStateByIndex(arc, s)
        x[s] = state[1]
        y[s] = state[2]
        z[s] = state[3]
        xdot[s] = state[4]
        ydot[s] = state[5]
        zdot[s] = state[6]
        t[s] = getTimeByIndex(arc, s)
    end
    JC::Float64 = getJacobiConstant(dynamicsModel, [x0, y0, z0, xdot0, ydot0, zdot0])
    exportCR3BPTrajectory(x, y, z, xdot, ydot, zdot, t, propTime, JC, file, name)
end

"""
    exportCR3BPTrajectory(x, y, z, xdot, ydot, zdot, t, TOF, JC, file, name)

Export CR3BP trajectory data to MAT file

# Arguments
- `x::Vector{Float64}`: x data [ndim]
- `y::Vector{Float64}`: y data [ndim]
- `z::Vector{Float64}`: z data [ndim]
- `xdot::Vector{Float64}`: xdot data [ndim]
- `ydot::Vector{Float64}`: ydot data [ndim]
- `zdot::Vector{Float64}`: zdot data [ndim]
- `t::Vector{Float64}`: Time data [ndim]
- `TOF::Float64`: Time-of-flight [ndim]
- `JC::Float64`: Jacobi constant
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportCR3BPTrajectory(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, t::Vector{Float64}, TOF::Float64, JC::Float64, file::MATLAB.MatFile, name::Symbol)
    traj = CR3BPTraj(x, y, z, xdot, ydot, zdot, t, TOF, JC)
    MATLAB.put_variable(file, name, traj)
end

"""
    exportInertalTrajectory(states, times, file, name)

Export inertial trajectory data to MAT file

# Arguments
- `states::Vector{Vector{Float64}}`: States [ndim]
- `times::Vector{Float64}`: Times [ndim]
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportInertialTrajectory(states::Vector{Vector{Float64}}, times::Vector{Float64}, file::MATLAB.MatFile, name::Symbol)
    nStates::Int64 = length(times)
    x::Vector{Float64} = zeros(Float64, nStates)
    y::Vector{Float64} = zeros(Float64, nStates)
    z::Vector{Float64} = zeros(Float64, nStates)
    xdot::Vector{Float64} = zeros(Float64, nStates)
    ydot::Vector{Float64} = zeros(Float64, nStates)
    zdot::Vector{Float64} = zeros(Float64, nStates)
    for s::Int64 in 1:nStates
        state::Vector{Float64} = states[s]
        x[s] = state[1]
        y[s] = state[2]
        z[s] = state[3]
        xdot[s] = state[4]
        ydot[s] = state[5]
        zdot[s] = state[6]
    end
    exportInertialTrajectory(x, y, z, xdot, ydot, zdot, times, times[end]-times[1], file, name)
end

"""
    exportInertialTrajectory(x, y, z, xdot, ydot, zdot, t, TOF, file, name)

Export inertial trajectory data to MAT file

# Arguments
- `x::Vector{Float64}`: x data [ndim]
- `y::Vector{Float64}`: y data [ndim]
- `z::Vector{Float64}`: z data [ndim]
- `xdot::Vector{Float64}`: xdot data [ndim]
- `ydot::Vector{Float64}`: ydot data [ndim]
- `zdot::Vector{Float64}`: zdot data [ndim]
- `t::Vector{Float64}`: Time data [ndim]
- `TOF::Float64`: Time-of-flight [ndim]
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportInertialTrajectory(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xdot::Vector{Float64}, ydot::Vector{Float64}, zdot::Vector{Float64}, t::Vector{Float64}, TOF::Float64, file::MATLAB.MatFile, name::Symbol)
    traj = InertialTraj(x, y, z, xdot, ydot, zdot, t, TOF)
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
    exportPseudoManifold(pseudoManifold, pseudoManifoldArcs, file, name)

Export pseudo-manifold data to MAT file

# Arguments
- `pseudoManifold::BCR4BPPseudoManifold`: BCR4BP pseudo-manifold object
- `pseudoManifoldArcs::Vector{BCR4BPPseudoManifoldArc}`: BCR4BP pseudo-manifold arcs
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function exportPseudoManifold(pseudoManifold::MBD.BCR4BPPseudoManifold, pseudoManifoldArcs::Vector{MBD.BCR4BPPseudoManifoldArc}, file::MATLAB.MatFile, name::Symbol)
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
    numArcs::Int64 = length(pseudoManifoldArcs)
    arcs::Vector{BCR4BP12Traj} = Vector{BCR4BP12Traj}(undef, numArcs)
    for a = 1:numArcs
        arc::MBD.BCR4BP12Arc = propagate(propagator, pseudoManifoldArcs[a].initialCondition, [0, pseudoManifoldArcs[a].TOF], pseudoManifold.dynamicsModel)
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
        arcs[a] = BCR4BP12Traj(x, y, z, xdot, ydot, zdot, theta4, t, H, pseudoManifoldArcs[a].TOF)
    end
    pseudoMan = PseudoMan(orb, arcs, pseudoManifold.TOF)
    MATLAB.put_variable(file, name, pseudoMan)
end

"""
    fullExportCR3BPFamily(family, MATFile, CSVFile, name)

Export CR3BP family data to MAT and CSV files

# Arguments
- `family::CR3BPOrbitFamily`: CR3BP periodic orbit family object
- `MATFile::String`: MAT file name
- `CSVFile::String`: CSV file name
- `name::Symbol`: Export object name
"""
function fullExportCR3BPFamily(family::MBD.CR3BPOrbitFamily, MATFile::String, CSVFile::String, name::Symbol)
    mf = MATLAB.MatFile(MATFile, "w")
    MATExportCR3BPOrbFamily(family, mf, name)
    MATLAB.close(mf)
    CSVExportCR3BPFamily(family, CSVFile)
end

"""
    fullExportCR3BPFamily(family, MATFile, CSVFile, name)

Export CR3BP multiple shooter family data to MAT and CSV files

# Arguments
- `family::CR3BPMSOrbitFamily`: CR3BP multiple shooter periodic orbit family object
- `MATFile::String`: MAT file name
- `CSVFile::String`: CSV file name
- `name::Symbol`: Export object name
"""
function fullExportCR3BPFamily(family::MBD.CR3BPMSOrbitFamily, MATFile::String, CSVFile::String, name::Symbol)
    mf = MATLAB.MatFile(MATFile, "w")
    MATExportCR3BPOrbFamily(family, mf, name)
    MATLAB.close(mf)
    CSVExportCR3BPFamily(family, CSVFile)
end

"""
    fullExportCR3BPFamily(family, MATFile, CSVFile, name)

Export CR3BP family data to MAT and CSV files

# Arguments
- `family::CR3BPContinuationFamily`: CR3BP trajectory continuation family object
- `MATFile::String`: MAT file name
- `CSVFile::String`: CSV file name
- `name::Symbol`: Export object name
"""
function fullExportCR3BPFamily(family::MBD.CR3BPContinuationFamily, MATFile::String, CSVFile::String, name::Symbol)
    mf = MATLAB.MatFile(MATFile, "w")
    MATExportCR3BPTrajFamily(family, mf, name)
    MATLAB.close(mf)
    CSVExportCR3BPFamily(family, CSVFile)
end

"""
    MATExportCR3BPOrbFamily(family, file, name)

Export CR3BP orbit family data to MAT file

# Arguments
- `family::CR3BPOrbitFamily`: CR3BP periodic orbit family object
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function MATExportCR3BPOrbFamily(family::MBD.CR3BPOrbitFamily, file::MATLAB.MatFile, name::Symbol)
    orbits::Vector{CR3BPOrb} = Vector{CR3BPOrb}(undef, getNumMembers(family))
    for o::Int64 = 1:getNumMembers(family)
        orbit::MBD.CR3BPPeriodicOrbit = getMember(family, o)
        propagator = MBD.Propagator()
        arc::MBD.CR3BPArc = propagate(propagator, orbit.initialCondition, [0, orbit.period], orbit.dynamicsModel)
        nStates::Int64 = length(arc.times)
        x::Vector{Float64} = zeros(Float64, nStates)
        y::Vector{Float64} = zeros(Float64, nStates)
        z::Vector{Float64} = zeros(Float64, nStates)
        xdot::Vector{Float64} = zeros(Float64, nStates)
        ydot::Vector{Float64} = zeros(Float64, nStates)
        zdot::Vector{Float64} = zeros(Float64, nStates)
        t::Vector{Float64} = zeros(Float64, nStates)
        for s::Int64 in 1:nStates
            state::Vector{Float64} = getStateByIndex(arc, s)
            x[s] = state[1]
            y[s] = state[2]
            z[s] = state[3]
            xdot[s] = state[4]
            ydot[s] = state[5]
            zdot[s] = state[6]
            t[s] = getTimeByIndex(arc, s)
        end
        JC::Float64 = getJacobiConstant(orbit)
        varsig::Float64 = getStabilityIndex(orbit)
        orbits[o] = CR3BPOrb(x, y, z, xdot, ydot, zdot, t, orbit.period, JC, varsig)
    end
    fam = CR3BPOrbFam(orbits)
    MATLAB.put_variable(file, name, fam)
end

"""
    MATExportCR3BPOrbFamily(family, file, name)

Export CR3BP multiple shooter orbit family data to MAT file

# Arguments
- `family::CR3BPMSOrbitFamily`: CR3BP multiple shooter periodic orbit family object
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function MATExportCR3BPOrbFamily(family::MBD.CR3BPMSOrbitFamily, file::MATLAB.MatFile, name::Symbol)
    orbits::Vector{CR3BPOrb} = Vector{CR3BPOrb}(undef, getNumMembers(family))
    for o::Int64 = 1:getNumMembers(family)
        orbit::MBD.CR3BPMSPeriodicOrbit = getMember(family, o)
        propagator = MBD.Propagator()
        orbitStates::Vector{Vector{Float64}} = []
        orbitEpochs::Vector{Float64} = []
        for n::Int64 in 1:length(orbit.nodeEpochs)-1
            arc::MBD.CR3BPArc = propagate(propagator, orbit.nodeStates[n], [orbit.nodeEpochs[n], orbit.nodeEpochs[n+1]], orbit.dynamicsModel)
            append!(orbitStates, arc.states)
            append!(orbitEpochs, arc.times)
        end
        nStates::Int64 = length(orbitEpochs)
        x::Vector{Float64} = zeros(Float64, nStates)
        y::Vector{Float64} = zeros(Float64, nStates)
        z::Vector{Float64} = zeros(Float64, nStates)
        xdot::Vector{Float64} = zeros(Float64, nStates)
        ydot::Vector{Float64} = zeros(Float64, nStates)
        zdot::Vector{Float64} = zeros(Float64, nStates)
        t::Vector{Float64} = zeros(Float64, nStates)
        for s::Int64 in 1:nStates
            state::Vector{Float64} = orbitStates[s]
            x[s] = state[1]
            y[s] = state[2]
            z[s] = state[3]
            xdot[s] = state[4]
            ydot[s] = state[5]
            zdot[s] = state[6]
            t[s] = orbitEpochs[s]
        end
        JC::Float64 = getJacobiConstant(orbit)
        varsig::Float64 = getStabilityIndex(orbit)
        orbits[o] = CR3BPOrb(x, y, z, xdot, ydot, zdot, t, orbit.period, JC, varsig)
    end
    fam = CR3BPOrbFam(orbits)
    MATLAB.put_variable(file, name, fam)
end

"""
    MATExportCR3BPTrajFamily(family, file, name)

Export CR3BP trajectory family data to MAT file

# Arguments
- `family::CR3BPContinuationFamily`: CR3BP trajectory continuation family object
- `file::MatFile`: MAT file
- `name::Symbol`: Export object name
"""
function MATExportCR3BPTrajFamily(family::MBD.CR3BPContinuationFamily, file::MATLAB.MatFile, name::Symbol)
    trajs::Vector{CR3BPTraj} = Vector{CR3BPTraj}(undef, getNumMembers(family))
    for m::Int64 = 1:getNumMembers(family)
        solution = MBD.CR3BPMultipleShooterProblem()
        numNodes::Int64 = length(family.nodes[m])
        solution.nodes = [MBD.shallowClone(family.nodes[m][n]) for n = 1:numNodes]
        solution.segments = [MBD.shallowClone(family.segments[m][s]) for s = 1:numNodes-1]
        propagator = MBD.Propagator()
        states::Vector{Vector{Float64}} = []
        epochs::Vector{Float64} = []
        for n::Int64 in 1:length(solution.nodes)-1
            arc::MBD.CR3BPArc = propagate(propagator, solution.nodes[n].state.data[1:6], [solution.nodes[n].epoch.data[1], solution.nodes[n+1].epoch.data[1]], solution.nodes[1].dynamicsModel)
            append!(states, arc.states)
            append!(epochs, arc.times)
        end
        nStates::Int64 = length(epochs)
        x::Vector{Float64} = zeros(Float64, nStates)
        y::Vector{Float64} = zeros(Float64, nStates)
        z::Vector{Float64} = zeros(Float64, nStates)
        xdot::Vector{Float64} = zeros(Float64, nStates)
        ydot::Vector{Float64} = zeros(Float64, nStates)
        zdot::Vector{Float64} = zeros(Float64, nStates)
        t::Vector{Float64} = zeros(Float64, nStates)
        for s::Int64 in 1:nStates
            state::Vector{Float64} = states[s]
            x[s] = state[1]
            y[s] = state[2]
            z[s] = state[3]
            xdot[s] = state[4]
            ydot[s] = state[5]
            zdot[s] = state[6]
            t[s] = epochs[s]
        end
        TOF::Float64 = 0.0
        for s::Int64 = 1:length(solution.segments)
            TOF += solution.segments[s].TOF.data[1]
        end
        JC::Float64 = getJacobiConstant(solution.nodes[1].dynamicsModel, solution.nodes[1].state.data[1:6])
        trajs[m] = CR3BPTraj(x, y, z, xdot, ydot, zdot, t, TOF, JC)
    end
    fam = CR3BPTrajFam(trajs)
    MATLAB.put_variable(file, name, fam)
end
