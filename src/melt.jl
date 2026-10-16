# Base's Float32 ^ Float32 can use Float64 intermediates on the CPU.
# The model only needs positive exponents of nonnegative physical quantities.
# Keep Float32 transcendental operations Float32, including the x=0 limit.
@inline _melting_power(x, y) = x^y
@inline _melting_power(x::Float32, y::Float32) = iszero(x) ? zero(x) : exp(y * log(x))

"""Evaluate a phase-boundary polynomial at pressure `P` (GPa), returning °C."""
@inline melting_curves(P, coefs) = evalpoly(P, coefs)

@inline exhaustion_cpx(P::T) where {T<:Real} = T(0.5) + T(0.08) * P
@inline melt_fraction_exhausted_cpx(P, M_cpx) = M_cpx / exhaustion_cpx(P)

@inline function T_exhausted_cpx(P, M_cpx, β, Ts, Tlh, Tl)
    return Ts + (Tlh - Ts) * _melting_power(melt_fraction_exhausted_cpx(P, M_cpx), inv(β))
end

# x^(3/2), for the nonnegative rescaled temperatures in the physical branches.
# sqrt avoids a general-purpose power call in these frequently evaluated branches.
@inline _power_three_halves(x) = x * sqrt(x)

"""
    melt_fraction(P, T, M_cpx, law)

Dry equilibrium melt mass fraction. Pressure is in GPa, temperature in °C,
and initial modal cpx is a mass fraction. Use Float32 inputs and a Float32 law
for Float32 device arithmetic.
"""
@inline function melt_fraction(P::S, T, M_cpx, law::AbstractMeltingLaw) where {S<:Real}
    Ts = law.T(P, law.coeffs_s)
    Tl = law.T(P, law.coeffs_l)
    Tlh = law.T(P, law.coeffs_lh)
    Tc = T_exhausted_cpx(P, M_cpx, S(1.5), Ts, Tlh, Tl)
    if T <= Ts
        return zero(Ts)
    elseif T >= Tl
        return one(Ts)
    elseif T <= Tc
        return _power_three_halves((T - Ts) / (Tlh - Ts))
    else
        Fc = melt_fraction_exhausted_cpx(P, M_cpx)
        return Fc + (one(Fc) - Fc) * _power_three_halves((T - Tc) / (Tl - Tc))
    end
end

"""Dissolved water (wt%) from bulk water (wt%), capped at saturation at `P` GPa."""
@inline function water_partitioning(P::S, Xbulk, D, F) where {S<:Real}
    return min(Xbulk / (D + F * (one(D) - D)), S(12) * _melting_power(P, S(0.6)) + P)
end

"""Melting-curve depression (°C), using bulk water in wt% and `D=0.01`."""
@inline function compute_ΔT(P::S, Xbulk, F) where {S<:Real}
    X = water_partitioning(P, Xbulk, S(0.01), F)
    return S(43) * _melting_power(X, S(0.75))
end

@inline melt_fraction_residual(F, P, T, M_cpx, Xbulk, law) =
    F - melt_fraction(P, T + compute_ΔT(P, Xbulk, F), M_cpx, law)

# Analytic branch derivatives keep the device path scalar and allocation-free.
@inline function _dry_productivity(P::S, T, M_cpx, law) where {S<:Real}
    Ts = law.T(P, law.coeffs_s)
    Tl = law.T(P, law.coeffs_l)
    Tlh = law.T(P, law.coeffs_lh)
    Tc = T_exhausted_cpx(P, M_cpx, S(1.5), Ts, Tlh, Tl)
    if T <= Ts || T >= Tl
        return zero(Ts)
    elseif T <= Tc
        return S(1.5) * sqrt((T - Ts) / (Tlh - Ts)) / (Tlh - Ts)
    else
        Fc = melt_fraction_exhausted_cpx(P, M_cpx)
        return (one(Fc) - Fc) * S(1.5) * sqrt((T - Tc) / (Tl - Tc)) / (Tl - Tc)
    end
end

"""
    melt_fraction_residual_and_derivative(F, P, T, M_cpx, Xbulk, law)

Return the wet equilibrium residual and its analytic derivative with respect
 to `F`. At phase and saturation kinks, the selected branch derivative is used.
"""
@inline function melt_fraction_residual_and_derivative(F, P::S, T, M_cpx, Xbulk, law) where {S<:Real}
    D = S(0.01)
    denominator = D + F * (one(D) - D)
    X = water_partitioning(P, Xbulk, D, F)
    shifted_T = T + S(43) * _melting_power(X, S(0.75))
    r = F - melt_fraction(P, shifted_T, M_cpx, law)
    # At zero water (or zero-pressure saturation), ΔT is constant in F.
    if iszero(X) || Xbulk / denominator >= S(12) * _melting_power(P, S(0.6)) + P
        return r, one(r)
    end
    dX = -Xbulk * (one(D) - D) / (denominator * denominator)
    dΔT = S(43) * S(0.75) * dX / sqrt(sqrt(X))
    dr = one(r) - _dry_productivity(P, shifted_T, M_cpx, law) * dΔT
    return r, dr
end

# Keep Float64's existing tolerance; Float32 needs a representable step tolerance.
@inline _melt_tolerance(::Type{S}) where {S<:AbstractFloat} = max(S(1e-12), S(4) * eps(S))

"""
    melt_fraction(P, T, M_cpx, Xbulk, law::WetMeltingLaw;
                  it_max=100, atol=max(S(1e-12), 4eps(S)), check=Val(true))

Solve hydrous batch-melting equilibrium by Newton–Raphson, with iterates
clamped to `[0,1]`. Pressure is
in GPa, temperature in °C, bulk water in wt%, and modal cpx a mass fraction.
`S` is the pressure's floating-point type; use matching scalar/law types.

Stops on residual or Newton step tolerance. The default checked mode throws
on invalid input or nonconvergence. In device kernels use `check=Val(false)`:
invalid input or nonconvergence returns a typed NaN without constructing errors.
"""
@inline function melt_fraction(P::S, T, M_cpx, Xbulk, law::WetMeltingLaw;
        it_max=100, atol=_melt_tolerance(S), check::Val{Check}=Val(true)) where {S<:AbstractFloat, Check}
    if !(P >= zero(P) && isfinite(P) && Xbulk >= zero(Xbulk) &&
         isfinite(T) && isfinite(M_cpx) && M_cpx >= zero(M_cpx) &&
         isfinite(atol) && atol > zero(atol))
        if Check
            throw(ArgumentError("invalid wet melting inputs or tolerance"))
        else
            return S(NaN)
        end
    end
    lo, hi = zero(S), one(S)
    melt_fraction_residual(lo, P, T, M_cpx, Xbulk, law) >= zero(S) && return lo
    melt_fraction_residual(hi, P, T, M_cpx, Xbulk, law) <= zero(S) && return hi
    F = (lo + hi) * S(0.5)
    tolerance = max(S(atol), eps(S))
    for _ in 1:it_max
        r, dr = melt_fraction_residual_and_derivative(F, P, T, M_cpx, Xbulk, law)
        abs(r) <= tolerance && return F
        next_F = clamp(F - r / dr, lo, hi)
        # Backtrack the Newton step across phase and saturation kinks.
        for _ in 1:32
            abs(melt_fraction_residual(next_F, P, T, M_cpx, Xbulk, law)) < abs(r) && break
            next_F = F + (next_F - F) * S(0.5)
        end
        # A representable Newton step can vanish at machine precision.
        abs(next_F - F) <= tolerance && return next_F
        F = next_F
    end
    if Check
        error("wet melting solver did not converge")
    else
        return S(NaN)
    end
end
