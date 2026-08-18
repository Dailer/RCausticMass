# Caustic Mass Estimator for Galaxy Clusters

This is an R implementation of the Caustic Mass Estimator for Galaxy Clusters developed in Python by [Gifford et al.](https://github.com/giffordw/CausticMass) A pre-cleaning of interlopers is recommended, using e.g., the shifting-gapper technique (included in this file as `shifting_gapper()`).

The caustic technique is a powerful method to infer cluster mass profiles to clustrocentric distances well beyond the virial radius. It relies on the measure of the escape velocity of the system using only galaxy redshift information. This method was introduced by [Diaferio & Geller (1997)](https://arxiv.org/abs/astro-ph/9701034) and [Diaferio (1999)](http://arxiv.org/abs/astro-ph/9906331), with the NFW-fitting refinement of [Serra et al. (2011)](https://arxiv.org/abs/1109.5178).

## Installation

```r
install.packages(c("magicaxis", "gplots", "imager", "pracma"))
source("RCausticMass.R")
```

## Basic usage

Run the code, here using the sample data provided with unknown values of R200 and cluster velocity dispersion (blind mode):

```r
data = read.table("sample_data.txt")
data = subset(data, CID == 1)
r = run_caustic(data$dproj, data$vlos, data$zclus, r200 = NA, clus_vdisp = NA)
```

If R200 and the velocity dispersion are known from an external catalogue, pass them and set `fix_r200 = TRUE` for a more precise, informed run:

```r
r = run_caustic(data$dproj, data$vlos, data$zclus, r200 = 1.2, clus_vdisp = 650, fix_r200 = TRUE)
```

The smoothed phase-space of galaxies can be plotted as well:

```r
image(r$x_range, r$y_range, r$img_tot, asp = NA, las = 1, xlab = expression(R[proj] ~ (Mpc)),
      ylab = expression(v[proj] ~ (km/s)))
```

## Main functions

| Function | Purpose |
|---|---|
| `run_caustic()` | The core estimator: R200, M200, velocity dispersion, membership, and uncertainties from a cluster's projected phase-space (`dproj`, `vlos`). Works blind (no prior) or informed (`fix_r200=TRUE`, given R200/velocity dispersion). |
| `run_caustic_robust()` | Wrapper around `run_caustic()` that retries with a smaller `rlimit` if the normal call fails with `"the contours do not expand to the radial limit"`. Recovered results are less precise on average — see the function's own documentation before trusting them. |
| `run_caustic_bootstrap()` | Runs `run_caustic()` over many galaxy resamples to obtain an empirical uncertainty estimate (distinct from the analytic D99 error). |
| `extend_outliers_nfw()` | Extends outlier/member classification to the *full* range of `dproj`, extrapolating the fitted NFW escape-velocity curve beyond the radius `run_caustic()` actually analysed. Also available directly as `caustic_outliers_extended` in `run_caustic()`'s own output. |
| `shifting_gapper()` | Interloper removal via the shifting-gapper technique (Fadda et al. 1996). |

## Validation

Tested against the Tempel et al. (2017) group catalogue (0.02 < z < 0.1), comparing blind-mode M200 estimates (no external prior) against the catalogue's own masses, for two different candidate-extraction conventions:

| | Extraction | N clusters | F_β(r) |
|---|---|---|---|
| **Original** | 3×R200, velocity window constrained by Tempel's own confirmed members | 713 | 0.72 |
| **Extended** | Fixed 10 Mpc radius, ±4000 km/s velocity window | 1485 | 0.50 |

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
- **Convergence depends strongly on richness.** Clusters with fewer than ~15-20 candidate galaxies often can't be fit reliably in blind mode; this is a data limitation, not a bug.
- **`run_caustic_robust()`'s recovered results can be confidently wrong.** Their own reported uncertainty (`M200_err_frac`) is not a reliable indicator of how far off they may be — treat any result where `rlimit_frac_used` is not `NA` as lower-confidence, regardless of its error bar.

## License

MIT — see [LICENSE](LICENSE).

## Citation

See [CITATION.cff](CITATION.cff), or cite Diaferio & Geller (1997), Diaferio (1999), and Serra et al. (2011) for the underlying method.

---

Author: Dailer F. Morell
