using InteractiveUtils: code_llvm
using ForwardDiff

# CPU regression checks exercise every physical branch before device execution.
@testset "Float32 throughout" begin
    for Law in (DryMeltingLaw, WetMeltingLaw)
        law = Law(Float32)
        @test isbitstype(typeof(law))
        @test Law(law.coeffs_s, law.coeffs_lh, law.coeffs_l).coeffs_s isa NTuple{3,Float32}
        @test (@inferred melting_curves(1f0, law.coeffs_s)) isa Float32
        @test (@inferred exhaustion_cpx(1f0)) isa Float32
        @test (@inferred melt_fraction_exhausted_cpx(1f0, 0.17f0)) isa Float32
        @test (@inferred T_exhausted_cpx(1f0, 0.17f0, 1.5f0, 1213.5f0, 1551.8f0, 1823f0)) isa Float32
        for t in (700f0, 1213.5f0, 1300f0, 1400f0, 1800f0, 2000f0)
            @test (@inferred melt_fraction(1f0,t,0.17f0,law)) isa Float32
            @test melt_fraction(1f0,t,0.17f0,law) ≈ melt_fraction(1.0,Float64(t),0.17,Law()) atol=2e-6
        end
    end
    @test (@inferred water_partitioning(1f0,0.1f0,0.01f0,0.1f0)) isa Float32
    @test (@inferred compute_ΔT(1f0,0.1f0,0.1f0)) isa Float32
    law = WetMeltingLaw(Float32)
    reference = WetMeltingLaw()
    for p in (0f0,1f0,1.5f0,3f0), w in (0f0,0.02f0,0.1f0,0.3f0,0.5f0)
        fs = Float32[]
        for t in 700f0:25f0:2000f0
            f = @inferred melt_fraction(p,t,0.17f0,w,law)
            device_f = @inferred melt_fraction(p,t,0.17f0,w,law;check=Val(false))
            @test f isa Float32
            @test f == device_f
            @test 0 <= f <= 1
            @test f ≈ melt_fraction(Float64(p),Float64(t),0.17,Float64(w),reference) atol=3e-6
            @test abs(melt_fraction_residual(f,p,t,0.17f0,w,law)) < 1e-4
            push!(fs,f)
        end
        @test all(diff(fs) .>= -2f-6)
    end
    @test melt_fraction(1f0,1200f0,0.17f0,0.1f0,law;atol=1f-15) isa Float32
    @test isnan(melt_fraction(1f0,1200f0,0.17f0,0.1f0,law;it_max=0,check=Val(false)))
    @test isnan(melt_fraction(-1f0,1200f0,0.17f0,0.1f0,law;check=Val(false)))
    @test isnan(melt_fraction(1f0,1200f0,0.17f0,-0.1f0,law;check=Val(false)))
    @test isnan(melt_fraction(1f0,1200f0,0.17f0,0.1f0,law;atol=0f0,check=Val(false)))
    @test_throws ArgumentError melt_fraction(1f0,1200f0,0.17f0,0.1f0,law;atol=0f0)

    # Analytic derivatives match independent CPU AD in both dry branches,
    # the no-melt branch and the saturation plateau (away from kinks).
    for S in (Float32, Float64), (f,p,t,w) in (
            (0.1,1.0,1200.0,0.1), (0.4,1.0,1500.0,0.1),
            (0.1,1.0,700.0,0.1), (0.01,1.0,1200.0,1.0),
            (0.1,1.0,1200.0,0.0), (0.1,0.0,1200.0,0.1))
        l = WetMeltingLaw(S)
        f,p,t,w = S.((f,p,t,w))
        r,dr = @inferred melt_fraction_residual_and_derivative(f,p,t,S(0.17),w,l)
        @test r isa S && dr isa S
        ad = ForwardDiff.derivative(x -> melt_fraction_residual(x,p,t,S(0.17),w,l),f)
        @test dr ≈ ad rtol=1e-5 atol=1e-6
    end
end

# Inspect the reachable arithmetic in the scalar entry points. CPU LLVM is
# a promotion regression check, not a substitute for Metal device compilation.
fp32_device(p,t,m,w,law) = melt_fraction(p,t,m,w,law;check=Val(false))
@testset "Float32 LLVM entry points" begin
    L = typeof(WetMeltingLaw(Float32))
    for (fn, signature) in (
            (exhaustion_cpx, Tuple{Float32}),
            (water_partitioning, Tuple{Float32,Float32,Float32,Float32}),
            (compute_ΔT, Tuple{Float32,Float32,Float32}),
            (melt_fraction, Tuple{Float32,Float32,Float32,L}),
            (fp32_device, Tuple{Float32,Float32,Float32,Float32,L}),
            (melt_fraction_residual_and_derivative, Tuple{Float32,Float32,Float32,Float32,Float32,L}))
        ir = sprint(io -> code_llvm(io,fn,signature;debuginfo=:none))
        @test !occursin(r"\bdouble\b",ir)
        @test !occursin("jl_apply_generic",ir)
    end
end

function fp32_allocations(law::L) where {L}
    fp32_device(1f0,1200f0,0.17f0,0.1f0,law)
    melt_fraction_residual_and_derivative(0.1f0,1f0,1200f0,0.17f0,0.1f0,law)
    solver_bytes = @allocated fp32_device(1f0,1200f0,0.17f0,0.1f0,law)
    derivative_bytes = @allocated melt_fraction_residual_and_derivative(0.1f0,1f0,1200f0,0.17f0,0.1f0,law)
    return solver_bytes,derivative_bytes
end
@testset "Allocation-free Float32 device entry points" begin
    @test fp32_allocations(WetMeltingLaw(Float32)) == (0,0)
end
