module EphemeridesLoader

import Pkg

const _ephem_env = "../ephemerides_env"

function __init__()
    Pkg.activate(_ephem_env; shared = false)
    @eval using Ephemerides
    @eval using FrameTransformations
    Pkg.activate(".")
end

export Ephemerides, FrameTransformations

end
