"""
Initial x-position targeter for CR3BP planar lunar free returns

Author: Jonathan Richmond
C: 7/1/25
"""

using MBD, CSV, DataFrames, DifferentialEquations, LinearAlgebra, StaticArrays

export PlanarPerpJCTargeter
export correct, getIndividualTrajectory, getTOF, interpTrajectory, propagateState

"""
    PlanarLFRXTargeter(dynamicsModel)

CR3BP planar lunar free return x-position targeter object

# Arguments
- `dynamicsModel::CR3BPDynamicsModel`: CR3BP dynamics model object
"""
struct PlanarLFRXTargeter
    dynamicsModel::MBD.CR3BPDynamicsModel                               # CR3BP dynamics model object

    function PlanarLFRXTargeter(dynamicsModel::MBD.CR3BPDynamicsModel)
        this = new(dynamicsModel)

        return this
    end
end
