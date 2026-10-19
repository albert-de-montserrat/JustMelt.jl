# Melting API

Import names explicitly from `JustMelt`. Implementations are in `src/types.jl`, `src/melt.jl` and `src/derivatives.jl`.

## Laws and phase boundaries

`DryMeltingLaw()` and `WetMeltingLaw()` construct the Table 2 parameterization using Float64. `Law(Float32)` selects Float32 coefficients. Either law accepts `Law(solidus, lherzolite_liquidus, liquidus)` or `Law(solidus, lherzolite_liquidus, liquidus, Float64)`, where all three coefficient tuples have the same length. The three-argument custom constructor preserves floating-point coefficient precision (and converts integer coefficients to a floating-point type). Construct laws on the host before passing them into kernels.

Coefficients are ordered by ascending powers of pressure. `melting_curves(P, coefficients)` evaluates their polynomial. Default solidus, lherzolite liquidus and liquidus tuples are `(1085.7, 132.9, -5.1)`, `(1475, 80, -3.2)` and `(1780, 45, -2)`, respectively. The result is °C for pressure in GPa.

## Dry melting

```julia
melt_fraction(P, T, M_cpx, law)
```

Evaluates the dry melting function (even when `law` is a `WetMeltingLaw`). Below the solidus it returns zero; at or above the liquidus it returns one. The cpx-bearing branch uses exponent 1.5. After cpx exhaustion, a second power-law branch joins continuously to the liquidus.

`exhaustion_cpx(P)` returns the reaction coefficient `0.5 + 0.08P`. `melt_fraction_exhausted_cpx(P, M_cpx)` returns `M_cpx / (0.5 + 0.08P)`. `T_exhausted_cpx(P, M_cpx, β, Ts, Tlh, Tl)` returns the corresponding dry transition temperature; `Tl` is accepted but unused.

The dry kernel expects compatible floating-point argument types. Callers must supply physically meaningful modal cpx and pressures for which the phase boundaries are ordered; these constraints are not generally checked.

## Water partitioning

```julia
water_partitioning(P, Xbulk, D, F)
compute_ΔT(P, Xbulk, F)
```

The first function returns dissolved water in wt%:

```text
min(Xbulk / (D + F*(1-D)), 12*P^0.6 + P)
```

`compute_ΔT` fixes `D=0.01` and returns `43 * Xliquid^0.75` in °C. At `F=0`, this describes solidus depression. `Xbulk=Inf` may be used to evaluate the saturation limit. Do not pass water mass fractions without converting to wt%.

## Wet equilibrium

```julia
melt_fraction(P, T, M_cpx, Xbulk, law::WetMeltingLaw;
              it_max=100, atol=max(S(1e-12), 4eps(S)), check=Val(true))
```

Solves `F = Fdry(P, T + ΔT(P, Xbulk, F), M_cpx)` by Newton–Raphson, with iterates clamped to `[0,1]` and the step halved until the residual decreases. Shifting temperature is equivalent to lowering all three phase boundaries equally, as in equations (11)–(13). This computes closed-system equilibrium batch melting; it does not model melt extraction or conserve enthalpy along a trajectory.

The solver returns exact endpoints when appropriate and otherwise stops when either the absolute residual or the Newton step is at most `atol`. `S` is the number type of pressure; the endpoints and the default tolerance preserve it. `S` may be a dual number (any `Real`), so ForwardDiff can differentiate through the solver. Float32 defaults to `4eps(Float32)`; Float64 retains `1e-12`. Negative pressure and bulk water throw `ArgumentError`; exhausting `it_max` throws an error. Supply finite physical pressure, temperature and nonnegative modal cpx, positive finite tolerance and sufficient iterations. Infinite bulk water is permitted for the saturation limit. Phase-boundary ordering and the physical upper limit on modal cpx are not comprehensively validated.

`melt_fraction_residual(F, P, T, M_cpx, Xbulk, law)` returns `F-Fdry(...)`. `melt_fraction_residual_and_derivative` returns `(residual, derivative_wrt_F)` analytically, without mutable derivative buffers or an AD runtime. At saturation and phase transitions the model has kinks, so derivatives should be interpreted within the active branch. The wet solver uses this analytic derivative for its Newton steps.

## Derivatives of the wet melt fraction

```julia
melt_fraction_and_T_derivative(P, T, M_cpx, Xbulk, law::WetMeltingLaw)
melt_fraction_and_P_derivative(P, T, M_cpx, Xbulk, law::WetMeltingLaw)
```

Return `(F, ∂F/∂T)` and `(F, ∂F/∂P)` from one ForwardDiff pass through the wet solver; the value and derivative come from the same evaluation. Both call `melt_fraction` with `check=Val(false)`, so invalid input or nonconvergence returns NaN rather than throwing. Derivatives at phase and saturation kinks follow the branch the solver lands on.

Use matching Float32 inputs and coefficients for Metal, and `check=Val(false)` inside kernels to return typed NaNs on invalid input or nonconvergence instead of constructing exceptions. See [Metal usage and tests](metal.md).

## Figure workflow

`figures/Katz_figs.jl` calls the same package functions. It records figure handles in `FIGURES` and saves PNGs in the repository's `figures/` directory, overwriting files with the same names. See [coverage and assumptions](../figures/figures-notes.md). No experimental observations are synthesized.
