"""
    melt_fraction_and_T_derivative(P, T, M_cpx, Xbulk, law::WetMeltingLaw)

Return `(F, ∂F/∂T)` for the wet equilibrium melt fraction, computed in one
forward-mode AD pass. Units follow [`melt_fraction`](@ref). Uses
`check=Val(false)`: invalid input or nonconvergence gives NaN.
"""
function melt_fraction_and_T_derivative(P, T, M_cpx, Xbulk, law::WetMeltingLaw)
    result = DiffResults.ImmutableDiffResult(T, (zero(T),))
    # ImmutableDiffResult is not updated in place; derivative! returns a new one.
    result = ForwardDiff.derivative!(result, T -> melt_fraction(P, T, M_cpx, Xbulk, law; check = Val(false)), T)
    return DiffResults.value(result), DiffResults.derivative(result)
end

"""
    melt_fraction_and_P_derivative(P, T, M_cpx, Xbulk, law::WetMeltingLaw)

Return `(F, ∂F/∂P)` for the wet equilibrium melt fraction, computed in one
forward-mode AD pass. Units follow [`melt_fraction`](@ref). Uses
`check=Val(false)`: invalid input or nonconvergence gives NaN.
"""
function melt_fraction_and_P_derivative(P, T, M_cpx, Xbulk, law::WetMeltingLaw)
    result = DiffResults.ImmutableDiffResult(P, (zero(P),))
    # ImmutableDiffResult is not updated in place; derivative! returns a new one.
    result = ForwardDiff.derivative!(result, P -> melt_fraction(P, T, M_cpx, Xbulk, law; check = Val(false)), P)
    return DiffResults.value(result), DiffResults.derivative(result)
end
