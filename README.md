# TIC 259541182.01 — a 1.15 Earth-radius planet candidate on a 1.00-day orbit

Discovered by **Advait Paliwal** (independent researcher, advaitpaliwal@gmail.com) on 2026-10-09 in TESS full-frame-image photometry.

| | |
|---|---|
| Host | TIC 259541182 = LEHPM 3836 = 2MASS J04470045−4836020, M dwarf, 45.6 pc, T = 11.95 |
| Period | 1.0044652 ± 0.0000009 d |
| Mid-transit | 2460968.5979 ± 0.0011 BJD_TDB |
| Depth / duration | 787 ± 69 ppm / 0.79 h |
| Radius | 1.15 ± 0.07 R⊕ |
| Equilibrium temperature | 875 ± 40 K |
| TESS sectors | 4, 5, 31, 32, 97, 98 (FFI only) |
| False-positive probability | 0.038 ± 0.001 (TRICERATOPS + MOLUSC Gaia limits, no neighbours removed, fitted depth; see note below); NFPP 2e-8 |

Status: **planet candidate**. Not a TOI or CTOI as of 2026-10-09. The remaining false-positive probability is dominated by an unresolved bound companion; high-resolution imaging is needed.

**Note (2026-10-16): corrected false-positive probability.** The TRICERATOPS runs behind the 0.034 ± 0.002 value were given the search depth of 586 ppm, not the fitted 787 ppm, and removed three neighbours. Rerun with all 44 sources kept, the fitted depth of 787 ppm and the wotan-detrended fit light curve folded on the fitted ephemeris (TRICERATOPS 1.1.0, ten runs of 10⁶ draws, seeds 1000–1009, same MOLUSC companion limits), it gives **FPP = 0.0376 ± 0.0014** (range 0.0351–0.0400; none of the ten runs is at or below 0.035) and NFPP = 1.9 × 10⁻⁸. Almost all of the FPP (0.036) is a planet orbiting an unresolved bound companion. This replaces the value in the table above and in `note/`; it is above both 0.035 and the 0.015 validation threshold. Source: the author's rerun of the TRICERATOPS script in `code/fpp/` with the depth and light curve changed as stated.

Contents: `note/` research note and figure; `code/` search, vetting, localization, false-positive and fit scripts; `data/` results (fit posterior summary, per-sector and per-pipeline depths, injections, systematics controls, difference-image localization, TRICERATOPS/MOLUSC summaries).

Data: TESS (NASA, via MAST: TESS-SPOC doi:10.17909/t9-wpz1-8s54, QLP doi:10.17909/t9-r086-e880) and ESA Gaia.

Licences: code MIT; text, figures and data CC BY 4.0.
