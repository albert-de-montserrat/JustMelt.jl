abstract type AbstractMeltingLaw end

## Anhydrous melting

struct DryMeltingLaw{N, _T, F} <: AbstractMeltingLaw
    coeffs_s::NTuple{N, _T}   # solidus
    coeffs_lh::NTuple{N, _T}  # auxiliar liquidos lherzolite
    coeffs_l::NTuple{N, _T}   # liquidus
    T::F                      # solidus/liquidus parameterization
end

@inline function DryMeltingLaw(coeffs_s::NTuple{N}, coeffs_lh::NTuple{N}, coeffs_l::NTuple{N}) where {N}
    DryMeltingLaw(coeffs_s, coeffs_lh, coeffs_l, float(promote_type(eltype(coeffs_s), eltype(coeffs_lh), eltype(coeffs_l))))
end

@inline DryMeltingLaw(coeffs_s, coeffs_lh, coeffs_l) = throw(ArgumentError("coeffs_s, coeffs_lh and coeffs_l must all be NTuple{N}"))

@inline function DryMeltingLaw(coeffs_s::NTuple{N, S}, coeffs_lh::NTuple{N, S}, coeffs_l::NTuple{N, S}, ::Type{T}) where {N, S, T <: Real}
    f         = melting_curves
    DryMeltingLaw{N, T, typeof(f)}(T.(coeffs_s), T.(coeffs_lh), T.(coeffs_l), f)
end

@inline function DryMeltingLaw(::Type{T}) where {T <: Real}
    coeffs_s  = T.((1085.7, 132.9, -5.1))
    coeffs_lh = T.((1475, 80, -3.2))
    coeffs_l  = T.((1780, 45, -2 ))
    f         = melting_curves
    DryMeltingLaw{3, T, typeof(f)}(coeffs_s, coeffs_lh, coeffs_l, f)
end

@inline DryMeltingLaw() = DryMeltingLaw(Float64)

## Hydrous melting

struct WetMeltingLaw{N, _T, F} <: AbstractMeltingLaw
    coeffs_s::NTuple{N, _T}   # solidus
    coeffs_lh::NTuple{N, _T}  # auxiliar liquidos lherzolite
    coeffs_l::NTuple{N, _T}   # liquidus
    T::F                      # solidus/liquidus parameterization
end


@inline function WetMeltingLaw(coeffs_s::NTuple{N}, coeffs_lh::NTuple{N}, coeffs_l::NTuple{N}) where {N}
    WetMeltingLaw(coeffs_s, coeffs_lh, coeffs_l, float(promote_type(eltype(coeffs_s), eltype(coeffs_lh), eltype(coeffs_l))))
end

@inline WetMeltingLaw(coeffs_s, coeffs_lh, coeffs_l) = throw(ArgumentError("coeffs_s, coeffs_lh and coeffs_l must all be NTuple{N}"))

@inline function WetMeltingLaw(coeffs_s::NTuple{N, S}, coeffs_lh::NTuple{N, S}, coeffs_l::NTuple{N, S}, ::Type{T}) where {N, S, T <: Real}
    f         = melting_curves
    WetMeltingLaw{N, T, typeof(f)}(T.(coeffs_s), T.(coeffs_lh), T.(coeffs_l), f)
end

@inline function WetMeltingLaw(::Type{T}) where {T <: Real}
    coeffs_s  = T.((1085.7, 132.9, -5.1))
    coeffs_lh = T.((1475, 80, -3.2))
    coeffs_l  = T.((1780, 45, -2 ))
    f         = melting_curves
    WetMeltingLaw{3, T, typeof(f)}(coeffs_s, coeffs_lh, coeffs_l, f)
end

@inline WetMeltingLaw() = WetMeltingLaw(Float64)