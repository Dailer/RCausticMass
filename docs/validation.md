# Validation against Tempel et al. (2017)

All comparisons use the Tempel et al. (2017) group catalogue, restricted to 0.02 < z < 0.1 (2139 groups), as the source of both candidate galaxies and the "real" R200/M200 values plotted on the x-axis. Two candidate-extraction conventions were tested, since (as found throughout development) the method's precision depends substantially on how the input sample is built -- not only on the caustic-fitting method itself.

Each panel below uses a random sample of 250 clusters for which `run_caustic()` converged in that specific mode (blind or informed) -- not necessarily the same 250 clusters across panels, since convergence depends on the mode and extraction. Error bars: `M200_err` is native (D99/Serra et al. 2011 formula); `R200_err` is derived by propagation from `M200_err` via M200 ∝ R200³, since R200 has no native uncertainty outside bayesian mode.

## Extraction conventions

**Constrained**: candidates within 3×R200 (R200 from Tempel's own catalogue), with the velocity window (|vlos|) constrained to the maximum |vlos| among Tempel's own confirmed members for that cluster. This ties the extraction directly to real cluster membership.

**Unconstrained**: candidates within a fixed 10 Mpc projected radius, with a fixed ±4000 km/s velocity window, regardless of cluster size. A more common convention when candidates are pulled from an external catalogue without prior knowledge of true membership.

`F_β(r)` was independently calibrated for each extraction (0.72 for Constrained, 0.50 for Unconstrained) -- see the main repository README for why a single default does not transfer between extraction geometries.

## Blind mode (no external prior)

R200 and M200 are both estimated from the data alone, with no R200/velocity-dispersion prior supplied.

### R200

| Constrained (N=250) | Unconstrained (N=250) |
|---|---|
| ![R200 blind constrained](images/R200_blind_original.png) | ![R200 blind unconstrained](images/R200_blind_extended.png) |

The unconstrained extraction shows a clear tendency to **overestimate R200 for intrinsically small clusters** (true R200 ~0.3-1 Mpc), by a factor of 2-3x in a substantial number of cases. This is consistent with a wider, fixed-radius field diluting the phase-space signal more severely for smaller systems -- see the main README's discussion of extraction-geometry effects.

### M200

| Constrained (N=250) | Unconstrained (N=250) |
|---|---|
| ![M200 blind constrained](images/M200_blind_original.png) | ![M200 blind unconstrained](images/M200_blind_extended.png) |

The R200 bias above propagates directly into M200 (mass is integrated up to the estimated R200), compounding the unconstrained extraction's lower precision.

## Informed mode (R200 fixed)

R200 is held fixed at Tempel's own catalogue value; only M200 (and concentration) are estimated from the data.

| Constrained (N=250) | Unconstrained (N=250) |
|---|---|
| ![M200 informed constrained](images/M200_informed_original.png) | ![M200 informed unconstrained](images/M200_informed_extended.png) |

Fixing R200 removes a major source of noise in both cases -- compare against the blind-mode M200 panels above. The constrained extraction still outperforms the unconstrained one even here, though the gap narrows substantially compared to blind mode.

## Summary

| Extraction | Mode | N |
|---|---|---|
| Constrained | Blind | 250 |
| Constrained | Informed | 250 |
| Unconstrained | Blind | 250 |
| Unconstrained | Informed | 250 |

**Takeaway**: both the extraction geometry (proportional vs. fixed-radius, membership-constrained vs. fixed velocity window) and the availability of an external R200 prior have a first-order effect on precision -- comparable to or larger than most algorithmic refinements to the fitting method itself. See the main README's "Known limitations" section, and [docs/parameter_calibration.md](parameter_calibration.md) for the detailed calibration process behind the current defaults.
