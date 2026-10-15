# Float32 and Metal.jl

The scalar calculations preserve Float32 when **all inputs and law coefficients are Float32**. No Metal dependency is needed for CPU calculations. Laws store coefficient tuples and a statically known function, so default Float32 laws are isbits values that can be passed as kernel arguments. Construct them on the host.

```julia
using JustMelt
law = WetMeltingLaw(Float32)
F = melt_fraction(1f0, 1200f0, 0.17f0, 0.1f0, law)
@assert F isa Float32
```

The zero-argument constructors still default to Float64. Transferring an array to Metal does not convert a separately constructed law or scalar arguments: explicitly use `Float32` laws and `f0` literals. Mixed Float64 arguments are not silently narrowed.

## Inside a Metal kernel

Use `check=Val(false)` for a device path that does not construct exceptions for invalid inputs or nonconvergence. It returns a Float32 NaN in those cases; inspect output for failures on the host. This mode does not make nonphysical phase boundaries or negative inputs to the dry helper formulas meaningful.

```julia
using Metal, JustMelt

function melting_kernel!(out, temperatures, law)
    i = Metal.thread_position_in_grid().x
    if i <= length(temperatures)
        @inbounds out[i] = melt_fraction(
            1f0, temperatures[i], 0.17f0, 0.1f0, law; check=Val(false))
    end
    return nothing
end

temperatures = Metal.MtlArray(collect(900f0:1f0:1400f0))
out = similar(temperatures)
law = WetMeltingLaw(Float32)
Metal.@sync begin
    @metal threads=128 groups=cld(length(out), 128) melting_kernel!(out, temperatures, law)
end
fractions = Array(out)
```

Indexing and synchronization follow the [Metal kernel programming documentation](https://metal.juliagpu.org/stable/usage/kernel/). The example uses the indexing API introduced in Metal 1.9. The explicit no-exception kernel path above is preferred when targeting runtimes with limited exception support.

## Precision and convergence

The default absolute tolerance is `max(S(1e-12), S(4)*eps(S))`, where `S` is the pressure type. This is about `4.77e-7` for Float32 and `1e-12` for Float64. The wet solver stops on residual tolerance, bracket width, or a bracket that cannot shrink further. A small bracket bounds the error in F; it does not guarantee the same tolerance on the residual near a steep branch. Asking for a sub-epsilon Float32 tolerance cannot recover Float64 accuracy.

The 3/2 melting powers use `x*sqrt(x)`. Other nonnegative Float32 powers use Float32 `exp` and `log`, with an explicit zero-base limit, because Base's general Float32 power implementation can internally promote to Float64. This introduces normal Float32 rounding differences from the Float64 model. No global fast-math mode is enabled.

Residual derivatives are analytic and return a scalar tuple; the implementation does not allocate a mutable DiffResults buffer or invoke ForwardDiff. At saturation and cpx-out kinks, derivatives follow the selected branch and are not globally smooth.

## Verification

The regular suite checks Float32 inference and return types, phase branches, saturation, convergence, agreement with Float64, analytic derivatives against CPU ForwardDiff, and absence of `double` arithmetic in inspected CPU LLVM entry points. CPU LLVM checks do not prove Metal compilation.

Metal compilation and execution have not been verified on the Linux development host. The current repository contains CPU tests, not an automated Metal hardware suite. Before claiming device validation, run the kernel example above on macOS with a functional Metal device and compare its output with CPU Float32 results across dry/wet branches, saturation and cpx exhaustion. Also check failure signaling with invalid inputs or exhausted iteration limits.
