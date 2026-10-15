# Katz, Spiegelman & Langmuir (2003), doi:10.1029/2002GC000433.
# Run: julia --project=. figures/Katz_figs.jl [--headless]
# Units: pressure GPa; temperature °C; water wt%; F and modal cpx mass fractions.
if "--headless" in ARGS
    @eval using CairoMakie
else
    @eval using GLMakie
end
using JustMelt
using JustMelt:  melting_curves, compute_ΔT, water_partitioning

const LAW = WetMeltingLaw()
const COLORS = (:royalblue, :seagreen, :orangered, :darkturquoise, :mediumorchid)
const OUTPUT = @__DIR__
mkpath(OUTPUT)
const FIGURES = Dict{Int, Figure}()

dry(P, T; cpx=0.15) = melt_fraction(Float64(P), T, cpx, LAW)
wet(P, T, water; cpx=0.17) = melt_fraction(Float64(P), T, cpx, water, LAW)
solidus(P) = melting_curves(P, LAW.coeffs_s)

function panel(fig, row, col; title="", xlabel="Temperature, °C", ylabel="Melt fraction, F", kwargs...)
    Axis(fig[row, col]; title, xlabel, ylabel,
         xgridstyle=:dot, ygridstyle=:dot, kwargs...)
end
function curve!(ax, x, y; label, color=:royalblue, kwargs...)
    @assert all(isfinite, x) && all(isfinite, y)
    lines!(ax, x, y; label, color, linewidth=2.5, kwargs...)
end
function finish(number, fig)
    FIGURES[number] = fig
    save(joinpath(OUTPUT, "fig$(number).png"), fig)
    println("Saved Figure $number")
    "--headless" in ARGS || display(fig)
end

# Figure 1: this model's phase boundaries; comparison models not implemented.
let
    fig = Figure(size=(760, 560))
    ax = panel(fig, 1, 1; title="Figure 1 — model phase boundaries (partial)",
               ylabel="Pressure, GPa", yreversed=true)
    ps = 0.0:0.01:8.0
    for (coeffs, label, style) in ((LAW.coeffs_s, "Solidus", :solid),
            (LAW.coeffs_lh, "Lherzolite liquidus", :dashdot),
            (LAW.coeffs_l, "Liquidus", :dash))
        curve!(ax, melting_curves.(ps, Ref(coeffs)), ps; label, color=:black, linestyle=style)
    end
    xlims!(ax, 1080, 2000); ylims!(ax, 0, 8); ax.yreversed = true
    axislegend(ax; position=:lb)
    finish(1, fig)
end

# Figure 2: anhydrous isobars, 15 wt% initial cpx.
let
    fig = Figure(size=(760, 560))
    ax = panel(fig, 1, 1; title="Figure 2 — isobaric anhydrous melting")
    ts = 1000.0:1.0:2000.0
    for (p, color) in zip((0.0, 1.0, 2.0, 3.0), COLORS)
        curve!(ax, ts, dry.(p, ts); label="$p GPa", color)
    end
    xlims!(ax, 1000, 2000); ylims!(ax, 0, 1)
    axislegend(ax; position=:lt)
    finish(2, fig)
end

# Figure 3: at the solidus F=0, so dissolved water is bulk water / D,
# capped at saturation. Experimental points are not available in this repo.
let
    fig = Figure(size=(760, 560))
    ax = panel(fig, 1, 1; title="Figure 3 — hydrous solidi (model curves only)",
               ylabel="Pressure, GPa", yreversed=true)
    ps = 0.0:0.005:8.0
    for (water, color) in zip((0.0, 0.05, 0.1, 0.3, Inf), COLORS)
        ts = [solidus(p) - compute_ΔT(p, water, 0.0) for p in ps]
        label = isinf(water) ? "Water saturated" : "$water bulk wt%"
        curve!(ax, ts, ps; label, color)
    end
    xlims!(ax, 850, 1900); ylims!(ax, 0, 8); ax.yreversed = true
    axislegend(ax; position=:rt)
    finish(3, fig)
end

# Figure 4: hydrous isobars, 1 GPa, 17 wt% initial cpx.
let
    fig = Figure(size=(760, 560))
    ax = panel(fig, 1, 1; title="Figure 4 — isobaric hydrous melting (1 GPa)")
    ts = 900.0:0.5:1400.0
    for (water, color) in zip((0.0, 0.02, 0.05, 0.1, 0.3), COLORS)
        curve!(ax, ts, wet.(1.0, ts, water); label="$water bulk wt%", color)
    end
    xlims!(ax, 900, 1400); ylims!(ax, 0, 0.4)
    axislegend(ax; position=:lt)
    finish(4, fig)
end

# Figure 5b specifies dissolved water directly: no implicit bulk-water solve.
let
    fig = Figure(size=(1100, 480))
    a = panel(fig, 1, 1; title="Figure 5a — 1.5 GPa", xlabel="Bulk water, wt%")
    b = panel(fig, 1, 2; title="Figure 5b — 1.5 GPa", xlabel="Water in the melt, wt%")
    bulk = 0.0:0.0005:0.25
    dissolved = 0.0:0.01:5.0
    for (t, color) in zip((1200.0, 1250.0, 1300.0, 1350.0), COLORS)
        curve!(a, bulk, wet.(1.5, t, bulk); label="$t °C", color)
        fs = [dry(1.5, t + 43*x^0.75; cpx=0.17) for x in dissolved]
        curve!(b, dissolved, fs; label="$t °C", color)
    end
    xlims!(a, 0, 0.25); ylims!(a, 0, 0.25)
    xlims!(b, 0, 5); ylims!(b, 0, 0.2)
    axislegend(a; position=:lt); axislegend(b; position=:rt)
    finish(5, fig)
end

# Figure 6: model curves only, without the experimental database.
let
    fig = Figure(size=(1050, 800))
    ts = 1000.0:1.0:2000.0
    for (i, p) in enumerate((0.0, 1.0, 1.5, 3.0))
        ax = panel(fig, (i-1)÷2+1, (i-1)%2+1;
                   title="Figure 6$(('a','b','c','d')[i]) — $p GPa (model only)")
        curve!(ax, ts, dry.(p, ts); label="Katz et al.")
        xlims!(ax, 1000, 2000); ylims!(ax, 0, 1)
    end
    finish(6, fig)
end

# Figure 8: calibrated functions only; experimental points are unavailable.
let
    fig = Figure(size=(1100, 480))
    a = panel(fig, 1, 1; title="Figure 8a — model calibration curve",
              xlabel="Water in the melt, wt%", ylabel="ΔT, °C")
    b = panel(fig, 1, 2; title="Figure 8b — model saturation curve",
              xlabel="Pressure, GPa", ylabel="Water at saturation, wt%")
    xs = 0.0:0.1:60.0
    # At 16 GPa the saturation cap exceeds 60 wt%; F=1 gives X_liquid=X_bulk.
    curve!(a, xs, [compute_ΔT(16.0, x, 1.0) for x in xs]; label="Parameterization", color=:black)
    ps = 0.0:0.01:8.0
    curve!(b, ps, [water_partitioning(p, Inf, 0.01, 0.0) for p in ps]; label="Parameterization", color=:black)
    xlims!(a, 0, 60); ylims!(a, 0, 900)
    xlims!(b, 0, 8); ylims!(b, 0, 50)
    finish(8, fig)
end

# Figures 9–10: only the Katz curves; other parameterizations and pMELTS
# are not implemented. Their 1 GPa cpx-out kink indicates 17 wt% cpx.
let
    fig = Figure(size=(760, 560))
    ax = panel(fig, 1, 1; title="Figure 9 — Katz model only (partial comparison)")
    ts = 1200.0:0.5:1800.0
    for (p, color) in zip((1.0, 3.0), COLORS)
        curve!(ax, ts, [dry(p, t; cpx=0.17) for t in ts]; label="$p GPa", color)
    end
    xlims!(ax, 1200, 1800); ylims!(ax, 0, 0.5)
    axislegend(ax; position=:lt)
    finish(9, fig)
end
let
    fig = Figure(size=(1100, 480))
    ts = 700.0:0.5:1800.0
    for (i, water) in enumerate((0.1, 0.5))
        ax = panel(fig, 1, i; title="Figure 10$(('a','b')[i]) — $water bulk wt% (Katz only)")
        for (p, color) in zip((1.0, 3.0), COLORS)
            curve!(ax, ts, wet.(p, ts, water); label="$p GPa", color)
        end
        xlims!(ax, 700, 1800); ylims!(ax, 0, 0.5)
        axislegend(ax; position=:lt)
    end
    finish(10, fig)
end
println("Figures written to $OUTPUT")
