# Katz et al. (2003) figure reproductions

Run `julia --project=. figures/Katz_figs.jl --headless` to save PNGs in `figures/` using CairoMakie. Omit `--headless` to use GLMakie and display figures interactively (requires GLMakie and a display server). `FIGURES` retains figure handles when included interactively.

All temperatures are in °C, pressures in GPa and water concentrations in wt%. Melt fraction and modal cpx are mass fractions. The script calls the existing functions in `src/melt.jl`; it does not digitize curves or invent experimental points.

| Figure | Coverage |
| --- | --- |
| 1 | Three Katz phase boundaries; other authors' solidi omitted. |
| 2 | Four dry isobars, initial modal cpx 0.15. |
| 3 | Five solidus curves including saturation; experimental points omitted. |
| 4 | Five hydrous isobars at 1 GPa, initial modal cpx 0.17. |
| 5 | Both panels at 1.5 GPa and 1200, 1250, 1300, 1350 °C; initial modal cpx 0.17. Panel b uses specified dissolved water directly. |
| 6 | Four model isobars; experimental points omitted. |
| 7 | Omitted: experimental P–T–F database is absent. Solidus already appears in Figure 1. |
| 8 | Both calibrated model functions; experimental points omitted. |
| 9 | Katz dry curves only; comparison models, pMELTS and experimental data omitted. |
| 10 | Katz hydrous curves in both panels; comparison models omitted. |
| 11 | Omitted: requires adding adiabatic integration of equations (20) and (23). |
| 12 | Omitted: requires the subduction-zone thermal field and enthalpy coupling. |

For Figures 4, 9 and 10 the initial modal cpx of 0.17 is inferred from the plotted cpx-out kink (F ≈ 0.293 at 1 GPa); their captions do not explicitly specify it. Figures 2, 5 and 6 specify their modal cpx in the captions.

The script produces model reproductions, not pixel-exact copies of the published plots. Partial plots are explicitly labelled.
