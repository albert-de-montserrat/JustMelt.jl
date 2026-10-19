module JustMelt

using DiffResults: DiffResults
using ForwardDiff: ForwardDiff

include("types.jl")
include("melt.jl")
include("derivatives.jl")

export DryMeltingLaw, WetMeltingLaw, melt_fraction,
    melt_fraction_and_T_derivative, melt_fraction_and_P_derivative

end # module JustMelt
