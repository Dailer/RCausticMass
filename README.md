# Caustic Mass Estimator for Galaxy Clusters

This is an R implementation of the Caustic Mass Estimator for Galaxy Clusters developed in Python by [Gifford et al.](https://github.com/giffordw/CausticMass) A pre-cleaning of interlopers is recommended, using e.g., the shifting-gapper technique (included in this file as `shifting_gapper()`).

The caustic technique is a powerful method to infer cluster mass profiles to clustrocentric distances well beyond the virial radius. It relies on the measure of the escape velocity of the system using only galaxy redshift information. This method was introduced by [Diaferio & Geller (1997)](https://arxiv.org/abs/astro-ph/9701034) and [Diaferio (1999)](http://arxiv.org/abs/astro-ph/9906331), with the NFW-fitting refinement of [Serra et al. (2011)](https://arxiv.org/abs/1109.5178).

## Installation

```r
install.packages(c("imager", "magicaxis"))
```

`magicaxis` is only needed for the optional diagnostic plot (`plot=TRUE`); everything else runs without it.

## Basic usage

Run the code, here using the sample data provided -- candidates around the Coma cluster (Abell 1656) -- with unknown values of R200 and cluster velocity dispersion (blind mode):

```r
source("RCausticMass.R")
data = read.table("coma_sample.txt", header = TRUE)
data = subset(data, CID == 1)
r = run_caustic(data$dproj, data$vlos, data$zclus[1], r200 = NA, clus_vdisp = NA, fbr = 0.50)
```

If R200 and the velocity dispersion are known from an external catalogue, pass them and set `fix_r200 = TRUE` for a more precise, informed run:

```r
r = run_caustic(data$dproj, data$vlos, data$zclus[1], r200 = 2.23, clus_vdisp = 947, fix_r200 = TRUE)
```

Setting `plot=TRUE` gives a full diagnostic figure -- the escape surface, the fitted NFW curve, its uncertainty band (shaded, from the same per-radius D99/Serra et al. 2011 error estimate used for the M200 error bar), R200, and the member/outlier split:

```r
r = run_caustic(data$dproj, data$vlos, data$zclus[1], fbr = 0.50, plot = TRUE)
```

<p align="center">
  <img src="docs/images/coma_diagnostic_plot.png" width="70%">
</p>

The smoothed phase-space itself can also be plotted directly:

```r
image(r$x_range, r$y_range, r$img_tot, asp = NA, las = 1, xlab = expression(R[proj] ~ (Mpc)),
      ylab = expression(v[proj] ~ (km/s)))
```

**How well did it do?** Coma is well-studied enough to check: the blind-mode run above (1743 candidates) gives

| | This run | Literature (Sohn et al. 2017, same caustic technique) | Difference |
|---|---|---|---|
| R200 | 2.03 Mpc | 2.23 Mpc | -9% |
| M200 | 9.68×10¹⁴ M☉ | 1.29×10¹⁵ M☉ | -25% |
| Velocity dispersion | 897 km/s | 947 km/s | -5% |

close to what the method's own reported uncertainty would suggest for this richness. See [docs/parameter_calibration.md](docs/parameter_calibration.md) for the same test on a second cluster (Abell 2029) and in informed mode, including a case where fixing R200 to its known value made the estimate *worse* -- a useful reminder that informed mode is more precise *on average*, not a guaranteed improvement for every individual cluster.

## Main functions

| Function | Purpose |
|---|---|
| `run_caustic()` | The core estimator: R200, M200, velocity dispersion, membership, and uncertainties from a cluster's projected phase-space (`dproj`, `vlos`). Works blind (no prior) or informed (`fix_r200=TRUE`, given R200/velocity dispersion). With `plot=TRUE`, also shows the fitted NFW escape curve's uncertainty band (D99/Serra et al. 2011 per-radius error, the same one behind the M200 error bar). |
| `run_caustic_robust()` | Wrapper around `run_caustic()` that retries with a smaller `rlimit` if the normal call fails with `"the contours do not expand to the radial limit"`. Recovered results are less precise on average — see the function's own documentation before trusting them. |
| `run_caustic_bootstrap()` | Runs `run_caustic()` over many galaxy resamples to obtain an empirical uncertainty estimate (distinct from the analytic D99 error). |
| `extend_outliers_nfw()` | Extends outlier/member classification to the *full* range of `dproj`, extrapolating the fitted NFW escape-velocity curve beyond the radius `run_caustic()` actually analysed. Also available directly as `caustic_outliers_extended` in `run_caustic()`'s own output. |
| `shifting_gapper()` | Interloper removal via the shifting-gapper technique (Fadda et al. 1996). |

## Advanced tuning parameters

Beyond `fbr`, a few lower-level parameters of `run_caustic()` were validated against three independent samples (two Tempel et al. 2017 extractions plus CIRS, Rines & Diaferio 2006) and are worth knowing about:

| Parameter | Default | Notes |
|---|---|---|
| `blur_gaussian` | `FALSE` | Controls `isoblur()`'s internal filter: `FALSE` uses the Deriche filter, `TRUE` uses Young-van Vliet (closer to a true Gaussian, but validated as never better and sometimes substantially worse — up to ~25-28% more M200 scatter on noisier extractions). |
| `gradu` | `1.0` | How fast a candidate contour's amplitude is allowed to grow with radius while tracing the escape surface. Inherited unchanged from causticpy at `0.5` originally; found never worse and sometimes clearly better at `1.0`. |
| `gradd` | `2.0` | Same idea for how fast the amplitude may fall. No detectable effect over the range tested; left at its original value. |
| `q` | `10` | Compresses the velocity axis relative to radius before the isotropic 2D kernel is applied. **Not recommended to change from a single default** — like `fbr`, its optimum is extraction-geometry-dependent and swings in *opposite* directions between datasets tested (favoring `q≈4-6` for a wide fixed-radius extraction, `q≈10` for a narrow membership-constrained one). Recalibrate jointly with `fbr` if you touch it. |

`neumann` (boundary conditions in `isoblur()`) and the density-grid resolution (`grid_by`, `nlevels`) were also tested and found to have negligible or ambiguous effect — left at their original defaults.

For the full calibration process (all sweeps, what was tried and didn't work, and the raw numbers behind these choices), see [docs/parameter_calibration.md](docs/parameter_calibration.md).

## Validation

Tested against the Tempel et al. (2017) group catalogue (0.02 < z < 0.1, 2139 groups), comparing blind-mode M200 estimates (no external prior) against the catalogue's own masses, for two different candidate-extraction conventions:

| | Extraction | N clusters | F_β(r) |
|---|---|---|---|
| **Constrained** | 3×R200, velocity window constrained by Tempel's own confirmed members | 250 | 0.72 |
| **Unconstrained** | Fixed 10 Mpc radius, ±4000 km/s velocity window | 250 | 0.50 |

<p align="center">
  <img src="docs/images/M200_blind_original.png" width="47%">
  <img src="docs/images/M200_blind_extended.png" width="47%">
</p>

Both use the same method with an F_β(r) calibrated for their respective extraction geometry, and the difference is visible: **how candidates are extracted matters as much as the method itself.** The proportional extraction (left), with a velocity window tied to real cluster membership, recovers M200 with substantially less scatter than the fixed-radius, wide-window extraction (right) -- even though both use blind mode (no external R200 prior) and both had F_β(r) independently optimized for their own geometry.

See [docs/validation.md](docs/validation.md) for the full comparison, including R200 recovery and informed-mode (R200 fixed) results for both extractions.

**Known caveat**: in blind mode, the fixed-radius extraction shows a clear tendency to overestimate R200 specifically for low-R200 clusters (a factor of 2-3x for some clusters with true R200 ~0.3-1 Mpc) -- see the R200 panel in [docs/validation.md](docs/validation.md). This is consistent with the wider extraction geometry diluting the phase-space signal more severely for intrinsically smaller systems.

## Known limitations

- **Blind mode is substantially less precise than informed mode.** If you have R200 and velocity dispersion from an external catalogue, use them (`fix_r200=TRUE`) — the difference in accuracy is large, not marginal.
- **`fbr` (the β-factor used in the mass integral) is not universal.** Its optimal value depends on the extraction geometry (fixed-radius vs. proportional-to-R200, wide vs. narrow velocity window), not just the cluster sample itself. Recalibrate against known masses when possible rather than trusting a single default across very different datasets.
- **`q` behaves the same way as `fbr`** — extraction-geometry-dependent, with no single default that serves all cases well (see Advanced tuning parameters above).
- **Convergence depends strongly on richness.** Clusters with fewer than ~15-20 candidate galaxies often can't be fit reliably in blind mode; this is a data limitation, not a bug.
- **`run_caustic_robust()`'s recovered results can be confidently wrong.** Their own reported uncertainty (`M200_err_frac`) is not a reliable indicator of how far off they may be — treat any result where `rlimit_frac_used` is not `NA` as lower-confidence, regardless of its error bar.
- **The caustic technique itself does not model Fingers-of-God distortion** (the velocity elongation of a cluster's core caused by peculiar-velocity dispersion when viewed nearly face-on). This is a structural limitation of the method (Serra & Diaferio 2013), not something any parameter here can fix -- it can cause genuine core members to be missed regardless of tuning.

## License

MIT — see [LICENSE](LICENSE).

## Citation

See [CITATION.cff](CITATION.cff), or cite Diaferio & Geller (1997), Diaferio (1999), and Serra et al. (2011) for the underlying method.

---

Author: Dailer F. Morell
