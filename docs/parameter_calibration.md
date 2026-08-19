# Parameter calibration: full process and results

This document records the exploration behind the defaults documented in the main README -- what was tested, on which samples, and what the numbers actually showed. The README states conclusions; this file shows the work.

**Validation samples used throughout:**
- **Constrained** (Tempel et al. 2017, 100 clusters): 3×R200 extraction, velocity window bounded by Tempel's own confirmed members.
- **Unconstrained** (same 100 Tempel clusters, re-extracted): fixed 10 Mpc radius, ±4000 km/s window.
- **CIRS** (Rines & Diaferio 2006, 74 clusters): fixed-radius extraction, independent dataset and centering convention.

All tables report blind mode (no R200/velocity-dispersion prior) unless noted. `sd` is the standard deviation of `R200_est/R200_true - 1` (R200) or `logM200_est - logM200_true` (M200) across the converged clusters in that run; `n` is how many converged.

---

## `fbr` (the F_β mass-integral factor)

No single value works across extractions -- this was the first and most consistent finding of the whole calibration effort. Final sweep, with the current full configuration (Deriche filter, `gradu=1.0`):

**Constrained** (n=38 throughout):
| fbr | R200 sesgo | M200 sd |
|---|---|---|
| 0.60 | -0.083 | 0.201 |
| 0.68 | -0.029 | 0.195 |
| **0.72** | **-0.004** | **0.191** |
| 0.75 | +0.015 | 0.188 |
| 0.80 | +0.044 | 0.184 |

**Unconstrained** (n=60 throughout):
| fbr | R200 sesgo | M200 sd |
|---|---|---|
| 0.40 | -0.121 | 0.380 |
| 0.48 | -0.036 | 0.371 |
| **0.50** | **-0.019** | **0.369** |
| 0.52 | -0.001 | 0.367 |
| 0.60 | +0.066 | 0.362 |

**CIRS**: fbr=0.44 chosen the same way (bias crossing zero), n≈65-67 depending on exact configuration.

Note the zero-bias crossing for M200 sits slightly above the values chosen for R200 in the Unconstrained case (~0.52-0.55 vs 0.50) -- 0.50 was kept since the difference is within run-to-run noise and R200 sesgo is already near zero there.

---

## `q` (velocity-axis compression before the isotropic kernel)

Behaves like `fbr` -- extraction-dependent, and here the disagreement between samples is *directional*, not just a matter of degree.

**A note on where `q=25`/`q=35` comes from**: the literature convention citing those specific values (e.g. GalWeight, Pisani 1993/1996) refers to the two-dimensional adaptive kernel method (2DAKM) -- i.e. `kernel="adaptive"` in this code, not the fixed-bandwidth Gaussian kernel used as the default here. The sweep below was run with `kernel="gaussian"` (this project's default), so it is not a direct test of the literature's own recommended value -- it may partly explain why no single q worked well across extractions, if the "right" q genuinely differs by kernel type. Not re-tested with `kernel="adaptive"` given that kernel's own lower precision ceiling (see below).

**Unconstrained** (fbr=0.50 fixed) -- prefers **low** q:
| q | n | R200 sd | M200 sd |
|---|---|---|---|
| 4 | 82 | 0.326 | 0.379 |
| 5 | 80 | 0.332 | 0.373 |
| 6 | 77 | 0.339 | **0.367** |
| 8 | 68 | 0.356 | 0.441 |
| 10 (default) | 59 | 0.357 | 0.385 |
| 15 | 42 | 0.417 | 0.469 |
| 20 | 26 | 0.480 | 0.551 |

**Constrained** (fbr=0.72 fixed) -- prefers **high** q, opposite direction:
| q | n | R200 sd | M200 sd |
|---|---|---|---|
| 4 | 37 | 0.228 | 0.484 |
| 6 | 39 | 0.195 | 0.303 |
| 8 | 37 | 0.161 | 0.220 |
| **10 (default)** | 38 | **0.154** | **0.196** |
| 15 | 39 | 0.178 | 0.211 |

**CIRS** (fbr=0.44 fixed) -- roughly flat, mild preference near q=8:
| q | n | R200 sd | M200 sd |
|---|---|---|---|
| 4 | 71 | 0.379 | 0.366 |
| 6 | 69 | 0.244 | 0.294 |
| 8 | 69 | 0.243 | **0.289** |
| 10 (default) | 67 | 0.244 | 0.294 |
| 15 | 60 | 0.272 | 0.336 |

**Conclusion**: left at `q=10`. No compromise value serves all three; recalibrate jointly with `fbr` for a new extraction convention rather than trusting either default.

---

## `blur_gaussian` (isoblur's internal filter: Deriche vs. Young-van Vliet)

`isoblur()`'s default (`gaussian=TRUE`) uses the Young-van Vliet recursive Gaussian approximation; `gaussian=FALSE` uses a Deriche filter instead (CImg's "pseudo-Gaussian", less accurate as a Gaussian approximation but tested here as more robust to noisy input). Full C++ source for both was retrieved from the `imager` GitHub repository to confirm exactly what each does (see `CImg.h`, `blur()`/`vanvliet()`/`deriche()`).

| Sample | Young-van Vliet M200 sd | Deriche M200 sd | Change |
|---|---|---|---|
| Constrained (fbr=0.72) | 0.198 | 0.196 | ~1% |
| **Unconstrained** (fbr=0.40-0.60) | 0.516-0.533 | 0.379-0.395 | **~25-28%** |
| CIRS (fbr=0.35-0.55) | 0.302-0.306 | 0.292-0.297 | ~3%, +2 converged |

Never worse in any sample; substantially better on the noisiest extraction. Adopted as the new default (`blur_gaussian = FALSE`).

**A custom from-scratch reimplementation** of the Young-van Vliet filter (exact coefficients and Triggs boundary conditions, derived directly from `CImg.h`) was also built and validated to match `imager::isoblur()` to floating-point precision (max difference ~1e-14) -- confirming the algorithm was understood correctly, though ultimately **not adopted**: it is 60-130x slower than the compiled C++ original (pure-R recursive loop vs. OpenMP-parallelized C++), and the project chose to keep the `imager` dependency rather than trade a rarely-binding dependency for a large runtime cost.

---

## `gradu` / `gradd` (contour gradient restriction)

Restricts how fast a candidate contour's amplitude can grow (`gradu`) or fall (`gradd`) between radial steps while `findcontours()` traces the escape surface outward. Both were inherited unchanged from causticpy (0.5 / 2.0) without independent validation.

**Unconstrained** (fbr=0.50): sharp transition between gradu=0.8 and 0.9, then a plateau:
| gradu | M200 sd |
|---|---|
| 0.5 (old default) | 0.385 |
| 0.8 | 0.383 |
| 0.9 | 0.354 |
| 1.0-3.0 (plateau) | 0.366-0.369 |

**Constrained** and **CIRS**: small but consistent improvement, never worse:
| Sample | gradu=0.5 | gradu=1.0 |
|---|---|---|
| Constrained | M200 sd=0.196 | M200 sd=0.191 |
| CIRS | M200 sd=0.294 | M200 sd=0.291 |

`gradd` showed no detectable effect anywhere in the range 1.25-3.0 tested. Adopted: `gradu = 1.0`, `gradd` unchanged.

---

## Parameters tested and found to have negligible or ambiguous effect

- **`neumann`** (isoblur boundary condition): differences at the 3rd decimal place (e.g. 0.389 vs 0.390) -- this code already reflect-pads manually before calling `isoblur()` and crops the result, which likely neutralizes most of the internal boundary behavior anyway.
- **External padding width** (currently 4σ before calling `isoblur()`): 4σ and 8σ gave matching results (0.385 vs 0.386); a 6σ test showed a small regression traced to a single marginal cluster, not a real trend. Left unchanged.
- **`grid_by`** (density-grid resolution, default 0.05): `0.08` gave a marginally better sd (0.278 vs 0.292 R200 sd) but converged 2 fewer clusters -- a real trade-off, not a clean win, and not confirmed on other samples. Left unchanged.
- **`nlevels`** (κ contour levels, default 200): already near a local optimum -- 50 and 100 are clearly worse, 300 showed an unexplained regression (M200 sd=0.476, worse than nlevels=50), 500 partially recovers. Left unchanged.

---

## `mirror` -- confirmed as a required bug fix, not a tunable parameter

Retested with the full current configuration (Deriche, gradu=1.0) across all three samples to make sure it still held after everything else changed:

| Sample | mirror=TRUE | mirror=FALSE |
|---|---|---|
| Constrained | M200 sesgo=+0.021, sd=0.191, n=38 | M200 sesgo=-0.166, sd=0.289, n=37 |
| Unconstrained | M200 sesgo=-0.001, sd=0.369, n=60 | M200 sesgo=-0.315, **sd=0.618**, n=53 |
| CIRS | M200 sesgo=+0.013, sd=0.291, n=67 | M200 sesgo=-0.078, sd=0.371, n=63 |

Worse in every metric, every sample, without exception. This mirrors galaxies to negative r before density estimation (Serra et al. 2011, Section 4.3) to avoid artificially depleted density near r=0; disabling it is not a legitimate configuration choice.

---

## Re-centering methods -- no clear winner found

Four independent ways to re-derive a cluster's center from its own candidate galaxies were compared (energy-tree via `sigma_plateau()`, 2D density peak, luminosity-weighted centroid, brightest-galaxy/BCG), each followed by a fresh `run_caustic()` call at the new center. Tested on the subset of Tempel Unconstrained clusters with ≥8 preliminary members (n≈25-38 depending on the exact comparison).

| Method | Typical R200 sd (varied by comparison run) |
|---|---|
| No re-centering (original catalogue center) | 0.221-0.414 |
| Energy tree (`sigma_plateau`) | 0.240-0.263 |
| 2D density peak | 0.232-0.246 |
| Luminosity-weighted centroid | 0.323 |
| BCG (brightest galaxy) | 0.447 (worst in every test) |

**BCG was consistently the worst option** -- a single galaxy's position is too noisy, and BCGs can be genuinely offset from the dynamical center by mergers or sloshing. Beyond that, **the ranking among the other three (and "no re-centering" itself) reordered every time the exact comparison sample changed**, at the n~25-38 scale available. This matches the out-of-sample instability finding below -- not enough clusters qualify for re-centering (≥8 preliminary members) to distinguish these methods with confidence. No re-centering method is implemented as a default in `run_caustic()`; if you want to experiment, re-derive a center externally and pass it as `clus_ra`/`clus_dec` before calling `rv_proj()`/`run_caustic()`.

**Iterative re-centering** (repeating the recenter-and-refit step, each time using the newly-estimated R200 to redefine members) was also tried and found **worse** than a single pass (R200 sd 0.302 vs 0.240) -- likely reinforcing the method's own bias rather than converging toward the true center.

---

## Kernel choice: gaussian vs. adaptive

The adaptive-bandwidth kernel (`kernel="adaptive"`, Serra et al. 2011 eq. 18-20) was compared against the standard fixed-bandwidth Gaussian kernel. The adaptive kernel showed a clearly lower precision ceiling regardless of its own tuning constant (~50-60% of clusters within a factor ~2 of the true mass, vs. a higher fraction for the Gaussian kernel under comparable conditions) and needs a much smaller `fbr` (~0.12-0.25) to avoid a large systematic bias. Not recommended; `kernel="gaussian"` remains the default.

---

## Concentration fitting: bounded vs. Bayesian

The default (`conc_method="bounded"`) fits concentration within a hard `[c_min, c_max]` range and was found to pile up against that boundary for a meaningful fraction of clusters (~25-30% in testing) when concentration is the only free parameter -- an artifact of the hard truncation, not a real physical result. `conc_method="bayesian"` (a log-normal prior centered on the Duffy et al. 2008 c-M relation) avoids this failure mode. Not changed as the default, since it requires both `M200_prior` and `R200_prior` to activate and is primarily relevant for informed-mode runs with external mass priors.

---

## Out-of-sample validation: a caution about all of the above

A held-out validation (calibrating `fbr` on one random half of the Constrained sample, testing on the other) found that **a single degenerate fit can double the reported `sd`** for a run of only ~25-30 clusters (0.337 calibration vs. 0.678 validation, driven almost entirely by one cluster with a 2.9 dex error; excluding it brought validation sd to 0.394, a much more modest and expected gap). Median-based statistics (bias, fraction within a factor of X) held up much better than `sd` under this test.

**Practical implication**: every `sd` comparison in this document (and the ones used to choose `q`, `gradu`, `blur_gaussian`, etc.) should be read with this in mind -- differences smaller than roughly 10-15% at these sample sizes (n~40-80) are not fully trustworthy on their own. The larger, more consistent effects documented above (e.g. `blur_gaussian`'s ~25% improvement on Unconstrained, or `mirror`'s effect) are far larger than this noise floor and can be trusted; several of the smaller, single-sample effects noted as "not adopted" above were left alone partly for this reason.

---

## External validation: literature clusters with independently published R200/M200

All the calibration above uses Tempel et al. (2017) and CIRS (Rines & Diaferio 2006) as ground truth. As an independent check, `run_caustic()` was also tested against two well-studied clusters with R200/M200 published from the same caustic technique (Sohn et al. 2017), not from Tempel/CIRS. Redshifts used (z=0.0235 for Coma, z=0.0784 for A2029) are the exact values from Sohn et al.'s own Table 1:

| Cluster | N candidates | R200 real | R200 est. (blind) | M200 real | M200 est. (blind) | M200 est. (informed) |
|---|---|---|---|---|---|---|
| Coma (A1656) | 1743 | 2.23 Mpc | 2.03 Mpc (-9%) | 1.29e15 M☉ | 9.70e14 M☉ (-25%) | 1.32e15 M☉ (**+2%**) |
| Abell 2029 | 627 | 1.97 Mpc | 2.09 Mpc (+6%) | 0.94e15 M☉ | 1.12e15 M☉ (+18%) | 1.65e15 M☉ (**+75%**) |

Blind-mode M200 errors (-25%, +18%) are broadly consistent with what Tempel/CIRS calibration would predict for this richness. The informed-mode result is a useful caution, though: fixing R200/velocity dispersion to their real values **helped Coma a lot** (error dropped to +2%) but **hurt Abell 2029 badly** (error grew to +75%, exceeding even the method's own wide reported uncertainty of 45.5%). A2029 has a known X-ray "sloshing" structure indicating it isn't fully dynamically relaxed, which may explain why fixing R200 didn't help the NFW fit find a better-behaved solution here. Informed mode remains better *on average* (as the large-sample Tempel/CIRS comparisons show clearly), but this is a concrete reminder that it is not a strict, guaranteed improvement for every individual cluster.

**A different kind of validation we haven't done**: Sohn et al. (2017) and related papers validate the caustic technique's *membership* completeness (not just R200/M200 accuracy) using Serra & Diaferio (2013) -- N-body mock catalogs with known true membership, recovering ~95% of true members within 3R200 with ~8% interloper contamination (better still within R200 itself: ~96%/~2%). All of our own validation compares against measured real-cluster values (which carry their own uncertainty), never against a simulation with exactly-known ground truth membership. We have not attempted an equivalent mock-catalog test.

---

## Bug fixes found during calibration (implemented, not configurable)

These aren't tunable parameters, but are worth listing since they were found *during* this calibration process and materially affect how trustworthy earlier or third-party results might be:

- `biwScale()`'s tuning constant didn't match Beers et al. (1990) (effectively ~20 instead of 9), and its `sqrt(n)` prefactor used the post-clipping count instead of the full sample size.
- `biwScale()` crashed (rather than returning `NA`) on empty input -- common for poor clusters where a radius/velocity cut leaves zero galaxies.
- Radial mirroring (r → -r) was missing from both `gaussian_kernel()` and `adaptive_kernel_2d()` (Serra et al. 2011, Sec. 4.3 requires it) -- caused artificially depleted density near r=0.
- `vlimit` was hardcoded to 4000 km/s inside `findcontours()`, disconnected from the value the caller actually passed to `run_caustic()`.
- `xmax`/`rlimit` used fixed defaults (6 / 5.8 Mpc) regardless of the actual data extent, causing systematic convergence failures for extractions at a different scale.
- κ contour levels were linearly spaced instead of log-spaced, wasting resolution in the sparse outer region of the density map.
- The gradient restriction on contour tracing (see `gradu`/`gradd` above) was entirely missing before being ported from causticpy.
- `run_caustic_robust()`'s recovered results can report a small, reassuring `M200_err_frac` while being wrong by a large factor -- this error estimate has no way to detect that the fitting radius was forced down artificially.
