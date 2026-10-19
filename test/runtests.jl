using Test, ForwardDiff
using JustMelt: DryMeltingLaw, WetMeltingLaw, melting_curves, melt_fraction,
    exhaustion_cpx, melt_fraction_exhausted_cpx, T_exhausted_cpx,
    water_partitioning, compute_ΔT, melt_fraction_residual,
    melt_fraction_residual_and_derivative,
    melt_fraction_and_T_derivative, melt_fraction_and_P_derivative

@testset "JustMelt" begin
    @testset "Constructors and Table 2 phase boundaries" begin
        for Law in (DryMeltingLaw, WetMeltingLaw)
            law = Law()
            @test law.coeffs_s == (1085.7, 132.9, -5.1)
            @test melting_curves(1.0, law.coeffs_s) ≈ 1213.5
            @test melting_curves(1.0, law.coeffs_lh) ≈ 1551.8
            @test melting_curves(1.0, law.coeffs_l) ≈ 1823.0
            @test Law(law.coeffs_s, law.coeffs_lh, law.coeffs_l).coeffs_s == law.coeffs_s
            typed = Law(law.coeffs_s, law.coeffs_lh, law.coeffs_l, Float32)
            @test typed.coeffs_s isa NTuple{3, Float32}
            @test Law(Float32).coeffs_s == typed.coeffs_s
            @test_throws ArgumentError Law([1,2], [3,4], [5,6])
            @test_throws ArgumentError Law((1,2), (3,4,5), (6,7))
        end
    end
    @testset "Dry branches and cpx-out continuity (equations 2–10)" begin
        law = DryMeltingLaw()
        for p in (0.0, 1.0, 1.5, 3.0)
            ts, th, tl = (melting_curves(p, cs) for cs in (law.coeffs_s, law.coeffs_lh, law.coeffs_l))
            fc = 0.15 / (0.5 + 0.08*p)
            tc = ts + (th-ts)*fc^(1/1.5)
            @test exhaustion_cpx(p) ≈ 0.5 + 0.08*p
            @test melt_fraction_exhausted_cpx(p, 0.15) ≈ fc
            @test T_exhausted_cpx(p, 0.15, 1.5, ts, th, tl) ≈ tc
            @test melt_fraction(p, ts-1, 0.15, law) == 0
            @test melt_fraction(p, ts, 0.15, law) == 0
            @test melt_fraction(p, tl, 0.15, law) == 1
            @test melt_fraction(p, tl+1, 0.15, law) == 1
            @test melt_fraction(p, tc, 0.15, law) ≈ fc
            @test melt_fraction(p, (ts+tc)/2, 0.15, law) ≈ fc/2^1.5
            @test melt_fraction(p, (tc+tl)/2, 0.15, law) ≈ fc+(1-fc)/2^1.5
            @test abs(melt_fraction(p, tc+1e-7, 0.15, law)-melt_fraction(p, tc-1e-7, 0.15, law)) < 1e-8
        end
    end
    @testset "Water mass balance and saturation (equations 16–18)" begin
        for f in (0.0, 0.1, 0.5, 1.0)
            x = water_partitioning(1.0, 0.05, 0.01, f)
            @test x*(0.01+f*0.99) ≈ 0.05
        end
        @test water_partitioning(1.0, Inf, 0.01, 0.0) == 13.0
        @test water_partitioning(0.0, 0.1, 0.01, 0.0) == 0.0
        @test water_partitioning(2.0, Inf, 0.01, 0.0) ≈ 12*2^0.6+2
        @test compute_ΔT(1.0, 0.0, 0.0) == 0.0
        @test compute_ΔT(1.0, 0.01, 0.0) ≈ 43.0
        @test compute_ΔT(1.0, Inf, 0.0) ≈ 43*13^0.75
    end
    @testset "Wet equilibrium and figure curves" begin
        law = WetMeltingLaw()
        for p in (0.0, 1.0, 1.5, 3.0), w in (0.0, 0.02, 0.05, 0.1, 0.3, 0.5)
            ts = 700.0:10.0:2000.0
            fs = [melt_fraction(p,t,0.17,w,law) for t in ts]
            @test all(f -> 0 <= f <= 1, fs)
            @test all(diff(fs) .>= -1e-10)
            @test maximum(abs(melt_fraction_residual(f,p,t,0.17,w,law)) for (f,t) in zip(fs,ts)) < 1e-10
            if iszero(w)
                @test fs ≈ [melt_fraction(p,t,0.17,law) for t in ts] atol=1e-10
            end
        end
        for t in (950.0, 1200.0, 1350.0, 1400.0)
            fs = [melt_fraction(1.0,t,0.17,w,law) for w in (0.0,0.02,0.05,0.1,0.3)]
            @test all(diff(fs) .>= -1e-10)
        end
        # Saturated solidus onset from equations 11, 16, 17.
        onset = 1213.5 - 43*13^0.75
        @test melt_fraction(1.0,onset-1e-6,0.17,0.3,law) == 0.0
        @test melt_fraction(1.0,onset+0.1,0.17,0.3,law) > 0.0
        @test_throws ArgumentError melt_fraction(-1.0,1200.0,0.17,0.1,law)
        @test_throws ArgumentError melt_fraction(1.0,1200.0,0.17,-0.1,law)
        @test_throws ErrorException melt_fraction(1.0,1200.0,0.17,0.1,law; it_max=0)
        # Derivative checked away from branch boundaries against finite differences.
        f, p, t, w = 0.1, 1.0, 1200.0, 0.1
        r, dr = melt_fraction_residual_and_derivative(f,p,t,0.17,w,law)
        residual(x) = melt_fraction_residual(x,p,t,0.17,w,law)
        @test r == residual(f)
        @test dr ≈ (residual(f+1e-6)-residual(f-1e-6))/2e-6 rtol=1e-6
    end

    @testset "Wet equilibrium accepts dual numbers" begin
        law = WetMeltingLaw()
        for (p, t, w) in ((0.5, 1200.0, 0.01), (1.0, 1200.0, 0.01), (1.0, 1300.0, 0.1))
            fP(x) = melt_fraction(x,t,0.17,w,law)
            fT(x) = melt_fraction(p,x,0.17,w,law)
            @test ForwardDiff.derivative(fP, p) ≈ (fP(p+1e-6)-fP(p-1e-6))/2e-6 rtol=1e-5
            @test ForwardDiff.derivative(fT, t) ≈ (fT(t+1e-4)-fT(t-1e-4))/2e-4 rtol=1e-5
            @test all(melt_fraction_and_P_derivative(p,t,0.17,w,law) .≈ (fP(p), ForwardDiff.derivative(fP, p)))
            @test all(melt_fraction_and_T_derivative(p,t,0.17,w,law) .≈ (fT(t), ForwardDiff.derivative(fT, t)))
        end
    end
end

include("float32.jl")
