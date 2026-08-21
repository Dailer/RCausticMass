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
| fbr | R200 bias | M200 sd |
|---|---|---|
| 0.60 | -0.083 | 0.201 |
| 0.68 | -0.029 | 0.195 |
| **0.72** | **-0.004** | **0.191** |
| 0.75 | +0.015 | 0.188 |
| 0.80 | +0.044 | 0.184 |

**Unconstrained** (n=60 throughout):
| fbr | R200 bias | M200 sd |
|---|---|---|
| 0.40 | -0.121 | 0.380 |
| 0.48 | -0.036 | 0.371 |
| **0.50** | **-0.019** | **0.369** |
| 0.52 | -0.001 | 0.367 |
| 0.60 | +0.066 | 0.362 |

**CIRS**: fbr=0.44 chosen the same way (bias crossing zero), n≈65-67 depending on exact configuration.

Note the zero-bias crossing for M200 sits slightly above the values chosen for R200 in the Unconstrained case (~0.52-0.55 vs 0.50) -- 0.50 was kept since the difference is within run-to-run noise and R200 bias is already near zero there.

### A single `fbr` also hides a real mass-dependent bias -- and it's mass, not richness

Investigated after noticing that M200 residuals correlate with true cluster mass (Pearson r=-0.35 to -0.40 in both Tempel and CIRS, using the single dataset-wide `fbr`): massive clusters are systematically *underestimated*, less massive ones *overestimated*. Several alternative explanations were ruled out first -- the radial extraction limit (`rlimit`), extrapolating the escape curve beyond the data with `nfw_escape_curve()`, and `kernel="adaptive"` all left the bias unchanged or made it worse (see "Kernel choice" and the discussion of `rlimit` elsewhere in this document). The real cause: **the optimal `fbr` itself depends on cluster mass**, confirmed independently in both CIRS and Tempel by splitting each sample into mass terciles and finding the zero-bias `fbr` separately in each:

| Mass tercile | CIRS optimal fbr | Tempel (unconstrained) optimal fbr |
|---|---|---|
| Low | ≈0.36 | <0.30 |
| Mid | ≈0.35 | ≈0.50 |
| **High** | **≈0.60** | **≈0.63** |

Monotonic and consistent in direction across both independent samples: the most massive clusters need close to *double* the `fbr` of the least massive ones. A single dataset-wide `fbr` (0.44 for CIRS, 0.50 for Tempel unconstrained) sits close to the low/mid terciles' optimum and is therefore too low for the high-mass tercile -- exactly reproducing the systematic underestimation found for massive clusters.

**Richness (N200) does not show the same pattern** -- tested the identical way (N200 terciles instead of mass terciles) in both samples:

| N200 tercile | CIRS optimal fbr | Tempel optimal fbr |
|---|---|---|
| Low | ≈0.50 | ≈0.50 |
| Mid | ≈0.37 | ≈0.59 |
| High | ≈0.39 | ≈0.40 |

Not monotonic, and not even in a consistent order between the two datasets -- confirming this is a genuine mass effect, not richness in disguise (richer clusters tend to be more massive, so the two are correlated, but only mass shows a clean, reproducible relationship with the needed `fbr`). Physically plausible explanation: `fbr`/F_β(r) corrects for velocity-anisotropy, and orbital anisotropy in cluster galaxy populations may itself depend systematically on halo mass.

**Not yet implemented as an automatic correction** -- doing so would need a preliminary mass estimate (chicken-and-egg with the thing `fbr` is used to compute) and calibration on a larger sample than the ~20-25 clusters per tercile used here. Documented as a known, real source of bias for now: if you have independent reason to expect a cluster is unusually massive, consider a higher `fbr` than the dataset-wide default.

### An iterative self-correction was tried -- works on average, adds per-cluster noise

Tested a two-pass approach on CIRS: run once with the dataset-wide `fbr` (0.44), look up a mass-informed `fbr` from a simple interpolation across the tercile anchor points (logM200=13.911→0.36, 14.364→0.35, 14.684→0.60), then re-run with that value. Per-tercile bias improved in all three groups (e.g. the high tercile went from -0.187 to -0.100 dex), confirming the correction works in the intended direction. But overall scatter got *worse* (sd 0.291→0.371) -- the correction depends on each cluster's own noisy first-pass mass estimate, so errors in that first pass propagate into a wrong `fbr` choice for that specific cluster. A damped version (blending the suggested `fbr` partially with the default, e.g. 30% toward the suggestion) gives a tunable trade-off -- less scatter penalty (0.291→0.311 at 30% damping) for a smaller bias improvement (high tercile: -0.187→-0.141) -- but there is no damping level that improves both at once. Not implemented in the code (no free lunch here), but documented as a real, usable option for someone willing to trade some scatter for less mass-dependent bias.

**Velocity dispersion was tested as an alternative predictor to mass** for the same purpose, on the idea that it's measured more directly (no dependence on the caustic fit itself). Using each cluster's *real* catalog vdisp, binned into quintiles: shows the same monotonic increase with `fbr` as mass does, but more gradually (CIRS: fbr≈0.39 at vdisp=361 km/s up to ≈0.61 at vdisp=784 km/s, roughly linear, R²=0.78; Tempel confirms the same slope within a factor ~1.2, though noisier, R²=0.32). The problem: this requires knowing the *true* vdisp, which isn't available in blind mode. Substituting the *preliminary* vdisp already computed internally (`vdisp0`, from the fixed ≤3 Mpc window used before any contour fitting) breaks the relationship entirely -- not monotonic, because `vdisp0` inherits its own richness/window-size-dependent bias (documented earlier: it over-estimates badly for intrinsically small clusters, where the fixed window is too wide) that has nothing to do with the true dynamical state the `fbr` relationship is tracking. The internally-iterated version of the preliminary vdisp (the one that self-consistently narrows its radius, similar to the blind window-narrowing method built for Coma) doesn't fix this either (correlation with true vdisp only 0.40). **No single-pass, no-external-data way to predict the right `fbr` was found.**

### A caution about what "true" mass even means here

All of the above treats Tempel/CIRS/Sohn et al. masses as ground truth. They aren't, and investigating this further changed the picture substantially:

- **Tempel's masses come from the virial theorem** (group velocity dispersion and virial radius from its own FoF group finder) -- genuinely independent of the caustic technique, unlike CIRS and Sohn et al. (2017), which both *are* caustic-technique measurements. So the CIRS/Sohn-based mass-`fbr` relationship carries a real circularity concern (calibrating a caustic-technique parameter against another caustic-technique measurement); the Tempel-based one does not.
- Recalibrating `fbr` for Coma specifically against a **multi-technique literature consensus** (median of X-ray, virial, Falco et al. 2014, Jeans, neural-flow, and deep-learning mass estimates -- deliberately excluding the two caustic-technique ones) gives `fbr≈0.39` -- barely different from the dataset-wide default (0.44), and far below the ≈0.63 needed to match Sohn et al.'s own caustic mass for Coma.
- But the independent-consensus median (0.70e15 M☉) is itself pulled down by X-ray hydrostatic estimates, which are well known in the literature to under-estimate true mass by omitting non-thermal pressure support (turbulence, bulk gas motions) -- so it isn't obviously "more correct" either.
- A controlled simulation study (Cluster-EAGLE, Old et al.) found that **both virial and caustic methods systematically over-estimate relative to X-ray hydrostatic mass** in mock clusters with known true mass (median ratios ≈1.3 and ≈1.5 respectively) -- consistent in direction with what we see, but this specific study found *no* significant mass-dependence in that bias for line-of-sight/projected samples, in tension with what we found on real clusters.
- A direct comparison against **weak lensing** (arguably the most model-independent technique available, and the only one tested here that's neither dynamical nor hydrostatic) found no consistent direction at all: of two clean, uncontaminated overlapping clusters (Kubo et al. 2009, Sloan Nearby Cluster Weak Lensing Survey) both massive enough to fall in the "high" tercile, `run_caustic()`'s default (`fbr=0.44`) over-estimated one (A1767: +0.17 dex) and under-estimated the other (A1066: -0.29 dex) by a similar margin; raising `fbr` to 0.60 made the first case worse (+0.35 dex) and the second better (-0.12 dex).

**Bottom line**: the mass-dependent trend in `fbr` is real and statistically well-established across dozens of clusters in two datasets. But there is no single "true" mass to calibrate it against -- every technique (caustic, virial, X-ray, weak lensing) carries its own known systematic, and at the level of any *individual* cluster, ordinary scatter is large enough to swamp whatever the mass-based correction would predict. Treat the mass-`fbr` relationship as a genuine population-level effect worth knowing about, not a per-cluster correction to apply with confidence.

**Update from simulation data**: a later test against TNG300 (see "Validation against a cosmological simulation" below), using true rather than measured mass, found this trend mostly disappears -- the strong monotonic pattern seen against Tempel/CIRS references was likely driven substantially by those references' own mass-dependent systematics, not a real property of the caustic technique. This doesn't overturn the observation that a single `fbr` is imperfect across a mass range, but it does mean the specific correction shape (roughly double `fbr` for the most massive clusters) should not be trusted as a property of the method itself.

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

## `mirror` -- velocity mirroring: keep on by default, but not universally

**A distinction worth being precise about**: this parameter controls *velocity* mirroring (duplicating each galaxy at both `+vlos` and `-vlos` before density estimation) -- it is unrelated to the *radial* mirroring (`r -> -r`, Serra et al. 2011 Section 4.3) that fixes the r=0 density-depletion problem. That radial mirroring is unconditional, always on, and not controlled by any exposed parameter; disabling it isn't possible and wasn't tested here. `mirror` is specifically about whether the velocity distribution is *assumed* symmetric about the cluster's systemic velocity.

Retested with the full current configuration (Deriche, gradu=1.0) across all three main samples to make sure the general default still held after everything else changed:

| Sample | mirror=TRUE | mirror=FALSE |
|---|---|---|
| Constrained | M200 bias=+0.021, sd=0.191, n=38 | M200 bias=-0.166, sd=0.289, n=37 |
| Unconstrained | M200 bias=-0.001, sd=0.369, n=60 | M200 bias=-0.315, **sd=0.618**, n=53 |
| CIRS | M200 bias=+0.013, sd=0.291, n=67 | M200 bias=-0.078, sd=0.371, n=63 |

Worse in every metric, every sample, without exception -- confirming `mirror=TRUE` (the default) is the right general-purpose choice: the assumption of velocity symmetry is a reasonable one for a typical, unclassified cluster sample.

**But this reverses for clusters specifically known to be unrelaxed.** Re-tested on the same CIRS relaxed/unrelaxed split used for the [systematic dynamical-state comparison](#) below (6 relaxed, 10 unrelaxed clusters):

| | mirror=TRUE | mirror=FALSE |
|---|---|---|
| **Relaxed** | R200 bias=-20.7%, sd=0.271 \| M200 bias=-0.271, sd=0.347 | R200 bias=-22.8%, sd=0.343 \| M200 bias=-0.310, sd=**0.573** |
| **Unrelaxed** | R200 bias=+16.5%, sd=0.290 \| M200 bias=+0.236, sd=0.309 | **R200 bias=+12.6%, sd=0.193** \| **M200 bias=+0.186, sd=0.256** |

For relaxed clusters, removing the mirror makes things worse, consistent with the general-sample result above (as expected: the symmetry assumption is a good one for a genuinely relaxed system). For unrelaxed clusters specifically, removing it **improves every metric at once** -- less bias and less scatter, in both R200 and M200. This makes physical sense: in a merging or infalling system, an asymmetric velocity distribution can be real signal (relative motion between substructures), not sampling noise, and forcing symmetry onto it actively discards that information rather than cleaning it up.

**Practical implication**: `mirror=TRUE` remains the right default for an unclassified sample (as the larger, general-purpose validation above shows). But if a cluster's dynamical state is known in advance (from the literature or an independent diagnostic) to be unrelaxed, `mirror=FALSE` is worth trying -- n=10 here is modest, so treat the exact numbers loosely, but the direction is clear and physically motivated, not just a small-sample fluctuation in one metric.

**Independent external evidence for the same effect**: Monteiro-Oliveira, Morell et al. (2021, MNRAS, "Unveiling the internal structure of Hercules supercluster") used `v1.0` of this code (`fbr=0.65`, following the Gifford et al. 2013 recipe) to estimate substructure masses in the Hercules supercluster (A2147, A2151, A2152), then tested the caustic technique against a controlled cluster-merger simulation with known true mass (ZuHone & Kowalik 2016, "Galaxy Cluster Merger Catalog"). They found the caustic mass estimate is measurably boosted for a period around pericentric passage during a merger -- contradicting the caustic technique's usual claim of being independent of dynamical state, and matching the direction of the systematic overestimation for unrelaxed clusters found independently in this document (A85, the CIRS relaxed/unrelaxed split, and the `mirror` result above). This is a genuinely independent confirmation, from a different code version, different clusters, and a different validation method (a controlled merger simulation rather than a real-cluster comparison), that the effect is real and not an artifact of this specific implementation.

**Update from simulation data**: a later, more direct test against TNG300 (see "Validation against a cosmological simulation" below) complicates the specific `mirror=TRUE`/`FALSE` recommendation above, even though it doesn't contradict the general finding that unrelaxed systems bias the mass estimate. Against *true* (not literature) substructure, `mirror=FALSE` did not reliably help, and whether it helps a given merging cluster turned out to depend on the viewing angle relative to the merger axis -- something no real dataset can control for. Read the `mirror`-vs-dynamical-state result above as a real pattern found in one specific sample (CIRS), not as settled operational guidance.

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

## Validation against a cosmological simulation (TNG300): genuine ground truth, no reference-mass ambiguity

Everything above compares against *measured* masses -- real-cluster references that, as the previous section found the hard way, carry their own technique-dependent systematics. TNG300 (IllustrisTNG, snapshot 94, z=0.06) removes that ambiguity: 43 massive halos (logM200 = 14.48-15.16, spanning the same mass range investigated for the `fbr`-mass question) with exactly known M200, R200, and 3D position/velocity for every member subhalo -- no measurement uncertainty anywhere. Candidates were mock-observed by projecting the true 3D data along a line of sight (dproj = projected separation perpendicular to the LOS; vlos = peculiar velocity along the LOS plus the local Hubble flow, H(z=0.06)≈69.7 km/s/Mpc for this simulation's cosmology). This is exactly the kind of test flagged as missing throughout this document (see "External validation" above, and the original discussion of Serra & Diaferio 2013's mock-catalog membership test).

**`fbr` recalibrated cleanly.** With full richness (median ~2300 members per cluster -- far richer than any real spectroscopic survey, since every simulated subhalo is "observed"), convergence was 43/43 and scatter far tighter than anything seen with real data (M200 sd≈0.10 dex vs. the usual 0.2-0.4 dex). The zero-bias `fbr` came out to ≈0.42-0.43 -- reassuringly close to the CIRS-calibrated default (0.44), a genuine independent cross-check that the earlier real-data calibration wasn't far off.

**The mass-`fbr` dependence found on Tempel/CIRS mostly disappears against true mass.** Splitting into mass terciles the same way as before: optimal `fbr` came out to ≈0.46 (low), ≈0.40 (mid), ≈0.41 (high) -- flat, not the strong monotonic trend (≈0.35→0.60) found against Tempel/CIRS references. This is a significant update to the earlier finding: it now looks like most or all of that mass-dependent trend was an artifact of the *reference* masses' own mass-dependent systematics (Tempel's virial masses, CIRS/Sohn's caustic masses), not a real property of the caustic technique itself, exactly the ambiguity flagged as unresolved in the `fbr` section above.

**Line-of-sight projection angle alone contributes as much scatter as everything else combined.** This is only testable with simulation data, since real observations never offer more than one viewing angle per cluster. The same 15 physical clusters were mock-observed along all three Cartesian axes: the resulting spread in `log M200_est` for a *single* cluster, from viewing angle alone, averaged sd≈0.091 dex -- essentially identical to the ≈0.097-0.10 dex *total* scatter found across different clusters in the main calibration. In other words, most of what looks like "measurement scatter" in any cluster mass comparison may not be reducible by better methodology at all -- it's a geometric fact about which direction a given cluster happens to be observed from, and no amount of parameter tuning changes that. This reframes essentially everything else in this document: the modest scatter improvements found from `blur_gaussian`, `gradu`, etc. (typically 5-25%) are working against a noise floor that projection angle alone can dominate.

**Richness matters continuously, confirming (and extending) the N200<6 finding.** Downsampling each cluster's full membership to a fixed N before running (same physical clusters, same viewing angle, only how many "galaxies" are kept varies): bias and scatter both degrade smoothly as N drops (M200 bias -0.04 dex / sd 0.095 at N=1000, worsening to bias -0.09 dex / sd 0.197 at N=50). Confirms the direction and rough scale of the richness-precision relationship found on real Tempel/CIRS data (see "Richness threshold" above), now isolated from any confound (same clusters, same true mass, only sample size differs).

**The `mirror`/dynamical-state finding does not replicate against true substructure, and now has an explanation.** A genuine (not literature-based) substructure diagnostic was built directly from the 3D data: split each cluster's members into two spatial groups via k-means on 3D position, and measure the velocity offset between the groups' centroids, normalized by their internal dispersion (`delta_v_norm`) -- a real, physically meaningful merger signature ranging from ≈0.1 (relaxed) to ≈2.9 (clear ongoing merger) across the 43 clusters. Splitting into terciles by this true diagnostic and re-running the `mirror=TRUE` vs. `FALSE` comparison from earlier in this document: for the most substructured tercile, `mirror=FALSE` did *not* help -- bias moved from +0.021 to -0.028 dex and scatter got *worse* (0.083→0.120), the opposite of the CIRS-based finding. Re-projecting the same 8 most-substructured clusters along all three axes showed why: whether `mirror=FALSE` helps a given *individual* cluster varies with viewing angle (2/8, 5/8, and 3/8 of clusters improved along the z, y, and x axes respectively, for the identical physical clusters) -- the benefit isn't a fixed property of "is this cluster unrelaxed," it depends on how the merger axis happens to align with the specific line of sight. Combined with the projection-angle finding above, this suggests the earlier CIRS result, while a real and reproducible pattern in that specific sample, may not generalize as a reliable "unrelaxed → use mirror=FALSE" rule; the true relationship is more likely mediated by projection geometry that isn't observable in real data at all. Treat the `mirror`/dynamical-state connection as a documented, sample-specific finding rather than settled guidance, pending a larger simulation-based test that could properly average over viewing angles.

---

## External validation: literature clusters with independently published R200/M200

All the calibration above uses Tempel et al. (2017) and CIRS (Rines & Diaferio 2006) as ground truth. As an independent check, `run_caustic()` was also tested against well-studied clusters with R200/M200 published independently (not from Tempel/CIRS), using several different techniques for comparison. Candidate galaxies for these individual clusters (Coma, Abell 2029, Abell 2199, Abell 1795, Abell 85, Abell 119) were extracted from **SDSS DR20**, separately from the Tempel/CIRS bulk samples.

| Cluster | N candidates | Reference technique | M200 real | M200 est. (blind) | M200 est. (informed) |
|---|---|---|---|---|---|
| Coma (A1656) | 1743 | Caustic (Sohn et al. 2017) | 1.29e15 M☉ | 9.70e14 M☉ (-25%) | 1.32e15 M☉ (**+2%**) |
| Abell 2029 | 627 | Caustic (Sohn et al. 2017) | 0.94e15 M☉ | 1.11e15 M☉ (+18%) | 1.65e15 M☉ (**+75%**) |
| Abell 2199 | 923 | Caustic (CIRS, same technique) | 3.41e14 M☉ | 4.08e14 M☉ (+20%) | -- |
| Abell 1795 | 890 | X-ray (Vikhlinin et al. 2006) | 6.03e14 M☉ | 5.47e14 M☉ (-9%) | -- |
| Abell 85 | 346 | Weak lensing, 2-halo fit (2025) | 4.14e14 M☉ (sum of both halos) | 6.56e14 M☉ (+58%) | -- |
| Abell 119 | 427 | Dynamical, substructure removed (Way et al. 1997) | ~3.05e14 M☉ | 4.49e14 M☉ (+47%) | -- |

Redshifts used for Coma/A2029 (z=0.0235 / z=0.0784) are the exact values from Sohn et al.'s own Table 1.

Blind-mode M200 errors for Coma/A2029/A2199/A1795 (-25%, +18%, +20%, -9%) are broadly consistent with what Tempel/CIRS calibration would predict for this richness. The informed-mode result for Coma/A2029 is a useful caution, though: fixing R200/velocity dispersion to their real values **helped Coma a lot** (error dropped to +2%) but **hurt Abell 2029 badly** (error grew to +75%, exceeding even the method's own wide reported uncertainty of 45.5%). A2029 has a known X-ray "sloshing" structure indicating it isn't fully dynamically relaxed, which may explain why fixing R200 didn't help the NFW fit find a better-behaved solution here. Informed mode remains better *on average* (as the large-sample Tempel/CIRS comparisons show clearly), but this is a concrete reminder that it is not a strict, guaranteed improvement for every individual cluster.

**A85 and A119 overestimate mass for a specific, understood reason**: both are known to have substructure the comparison reference explicitly excluded or separated (A85: an active 2:1-mass-ratio merger, compared against the sum of two independently-fit weak-lensing halos; A119: a literature dynamical mass computed "after elimination of subgroups", with the paper itself noting probable line-of-sight substructure). `run_caustic()` treats the whole phase-space as one system and has no way to separate substructure -- consistent with the systematic overestimation found for merging/unrelaxed clusters elsewhere in this document (see `mirror` below).

**A different kind of validation we haven't done**: Sohn et al. (2017) and related papers validate the caustic technique's *membership* completeness (not just R200/M200 accuracy) using Serra & Diaferio (2013) -- N-body mock catalogs with known true membership, recovering ~95% of true members within 3R200 with ~8% interloper contamination (better still within R200 itself: ~96%/~2%). All of our own validation compares against measured real-cluster values (which carry their own uncertainty), never against a simulation with exactly-known ground truth membership. We have not attempted an equivalent mock-catalog test.

---

## Predicting dynamical state from the data itself: substructure tests

*Note: the sections below (through "Edge-detection cross-check") describe several parameters (`combine_branches`, `fit_r_upper_mult`, `fit_r_lower_div`, `blind_refit`, `mass_curve`, `mass_smooth_spar`, `mass_r_min`, `level_spacing`) that were implemented specifically to run these tests and then removed from `RCausticMass.R` afterward, since none improved on the defaults -- they will not be found in the current file. `compute_edge` is the one exception: it's an original, still-present feature.*



Given the `mirror` finding above (relaxed vs. unrelaxed clusters need opposite settings), an obvious next question is whether dynamical state can be predicted from the candidate data itself, rather than requiring an external literature classification. Tested on the same CIRS relaxed/unrelaxed split (6 relaxed, 12 unrelaxed clusters with independent literature classifications agreeing across multiple papers).

**Simple univariate tests on the member velocity distribution (Anderson-Darling normality, skewness, a v-vs-(-v) Kolmogorov-Smirnov symmetry test) did not work** -- tested twice, once using cluster membership derived from a `mirror=TRUE` run (circular, since that membership definition already assumes symmetry) and once using a blind, symmetry-agnostic window-narrowing method to avoid that circularity. Both times the direction was mostly *opposite* to what was expected (relaxed clusters showing more apparent asymmetry than unrelaxed ones by these metrics), not just noisy. Likely explanation: the literature relaxed/unrelaxed classification comes from diagnostics unrelated to line-of-sight velocity symmetry specifically (X-ray morphology, 2D positional substructure statistics), so there's no strong reason these 1D velocity-only tests should track it.

**The Dressler & Shectman (1988) test (implemented in `ds_test.R`, included alongside this project) worked well on CIRS**: unlike the univariate tests, it combines projected *position* with velocity (comparing each galaxy's local neighbourhood mean velocity/dispersion against the cluster-wide values), which is the standard, literature-validated approach for exactly this purpose. Median p-value 0.202 for relaxed clusters vs. 0.000 for unrelaxed ones -- a real, promising separation (though not perfect: one relaxed cluster, A2244, showed strong spurious substructure evidence; 7 of 11 unrelaxed clusters were correctly flagged at p<0.1).

**But using it as an automatic decision rule for `mirror` on a noisier extraction failed badly.** Tested end-to-end on 39 real Tempel et al. (2017) clusters (unconstrained extraction): for each cluster, the DS test decided `mirror=TRUE`/`FALSE`, then `run_caustic()` ran with that choice, compared against a fixed `mirror=TRUE` baseline.

| | R200 bias | R200 sd | M200 bias (dex) | M200 sd |
|---|---|---|---|---|
| Baseline (`mirror=TRUE` always) | -9.7% | 0.269 | -0.105 | 0.356 |
| DS-adaptive choice | -16.7% | 0.325 | -0.217 | 0.656 |

Worse in every metric. Specifically for the 18 clusters where the DS test picked `mirror=FALSE`: had they used `mirror=TRUE` instead, R200 bias would have been -2.8% (sd=0.228); with the DS-chosen `mirror=FALSE` it was -32.3% (sd=0.258) -- a clear wrong call, not a marginal one. Likely explanation: DS was validated on CIRS, a curated sample with good per-cluster spectroscopic coverage; on Tempel's wide, contamination-heavy "unconstrained" extraction, a significant DS detection more often reflects field interlopers than genuine merger substructure, so it doesn't carry the same meaning it did on the cleaner sample it was validated against.

**Takeaway**: the DS test is a useful *manual* diagnostic on a reasonably clean candidate sample (which is what it was validated for in the literature), not a plug-in automatic classifier for noisy, wide-window extractions.

---

## Combining the upper/lower velocity branches: `combine_branches`

When `mirror=FALSE` allows the density map to be genuinely asymmetric, `findcontours()` has to decide how to combine the separately-traced positive-velocity and negative-velocity branches into one final amplitude curve at each radius. `combine_branches='min'` (default, matches causticpy) takes the smaller of the two magnitudes at each radius. Two alternatives were tested on the 11-15 worst-performing Tempel et al. (2017) clusters found under the default `mirror=TRUE` configuration, re-run with `mirror=FALSE`:

- **`'mean'`** (point-by-point average instead of minimum): worse on both R200 and M200 medians (dR 0.233→0.274, dM 0.315→0.387) -- one specific cluster degraded sharply (R200 error 0.006→0.426) where a noisier branch dragged the average away from an otherwise well-behaved fit. `min`'s conservatism (never let the amplitude leak through a gap on whichever side extends further) is doing real work, not just following convention.
- **`'best_nfw'`** (fit each whole branch to a simple NFW shape independently, keep whichever fits better): essentially indistinguishable from `min` in aggregate (identical median errors) -- in most individual cases it picks the same branch `min` effectively favours anyway, with one case notably worse (a branch that fit the NFW *shape* better was not the more *accurate* one, likely because the true curve there isn't a clean single NFW form).

Neither improves on `min`; kept as options for further experimentation, not as better defaults.

---

## The mass integral doesn't use the NFW-fit curve -- a structural finding, not a bug

Investigated while testing whether narrowing/widening the NFW-fit radial range (`fitting_radii`, normally `[R200/3, R200]`) affects M200/R200 accuracy. It doesn't, for a specific and important structural reason: under the default `mass_method='integral'`, `mass_from_Ar()` integrates `Ar_finalD` -- the raw contour curve exactly as `findcontours()` selected it -- directly over the *entire* radial grid, not the NFW-fit curve (`caustic_fit`) and not restricted to `fitting_radii` at all. The NFW fit (and therefore `fitting_radii`, `fit_r_upper_mult`, `fit_r_lower_div`) only affects membership classification (`memflag`) and the fit concentration (when a mass prior is supplied) -- never the R200/M200 numbers themselves in this mode.

This was confirmed exhaustively: extending `fit_r_upper_mult` from 1.0 to 2.0 produced *bit-for-bit identical* `log_M200_est` on real Tempel clusters, not just similar values. A `blind_refit=TRUE` option was added (refits once more using the actual `r200_est` from the first pass, rather than the necessarily-approximate preliminary R200, to fix the circularity that made this untestable in blind mode at all) -- it made no difference either, for the same underlying reason.

**Where `fitting_radii` genuinely does matter**: `mass_method='bayesian_joint'`, which fits `rii`/`ArD` directly for a joint (M200, concentration) estimate. Tested there (60 real Tempel clusters): extending `fit_r_upper_mult` made things *worse*, monotonically --

| fit_r_upper_mult | M200 bias (dex) | M200 sd |
|---|---|---|
| 1.0 (default) | +0.202 | 0.381 |
| 1.1 | +0.210 | 0.383 |
| 1.5 | +0.256 | 0.423 |
| 2.0 | +0.248 | 0.433 |

So the answer to "does widening/narrowing the fit range help" is no either way: irrelevant under the default method, mildly harmful under `bayesian_joint`.

---

## What curve to integrate for mass: raw contour vs. NFW fit vs. smoothed

Given the above, a genuinely different question is whether integrating something *other than* the raw `Ar_finalD` curve -- which can be noisy, especially for poorly-sampled clusters -- would help, since (per the previous section) `mass_method='integral'` never gets that option by default. Added `mass_curve` (`'raw'` default, `'nfw_fit'`, `'smooth'`) to test this directly. Tested on 95 real Tempel clusters (both extractions):

| `mass_curve` | R200 bias | R200 sd | M200 bias (dex) | M200 sd |
|---|---|---|---|---|
| `raw` (default) | -1.8% | 0.255 | -0.000 | 0.332 |
| `nfw_fit` (integrate `caustic_fit` instead) | +4.5% | 0.303 | +0.087 | 0.322 |
| `smooth` (`smooth.spline()` on `Ar_finalD`, auto spar) | -1.8% | 0.255 | -0.000 | 0.332 (**identical to raw** -- GCV picked ~zero smoothing) |
| `smooth`, `mass_smooth_spar=0.9` (forced heavy smoothing) | -2.9% | 0.251 | -0.011 | 0.313 |

`nfw_fit` is a clear loser (worse bias in both R200 and M200, no compensating scatter win). `smooth` with automatic smoothing-parameter selection changes nothing at all -- `smooth.spline()`'s GCV criterion found no excess roughness worth penalizing in these curves. Forcing heavy smoothing (`spar=0.9`) gives a small, genuine improvement in scatter (M200 sd 0.332→0.313, ~6%) at a small bias cost -- real but modest, not worth changing the default over on its own.

**`mass_r_min` (zero out the integrand below a radius, to test whether the innermost, most mirror/grid-exposed radii add noise rather than signal)**: no improvement at `r_min=0.05` Mpc (bias got slightly worse, scatter barely moved), and `r_min=0.1` Mpc **broke the fit entirely for nearly every cluster** ("need at least two non-NA values to interpolate" -- too little of the profile survives to locate where it crosses 200×critical density for many clusters, especially smaller ones). No usable middle ground was found between "no effect" and "breaks".

---

## An unresolved observation: some escape curves fall near the cluster centre

While investigating the above, `caustic_profile` (`Ar_finalD`) was found to *decline* monotonically from a near-central peak (~0.1-0.15 Mpc) all the way out to at least R200/2, in about 7% of real Tempel clusters (4 of 60 checked) -- the opposite of the more typical rise-then-fall shape. Initial hypothesis (that this happens specifically when the contour-selection step, which matches a single summary statistic `tot_avg` against `4×σ²`, can't find any good match among the few candidate contours available) did not hold up under a larger check: the correlation between match quality and this "falling" shape was not statistically significant (Wilcoxon p=0.94) across the full 60-cluster sample -- some perfectly normal-shaped clusters had far *worse* matches (up to 3114% off target) than any of the falling cases. The falling clusters were also not particularly poor in candidate count (167-693 candidates). **No confirmed explanation currently exists for this pattern.** Averaging the main curve with the independent edge-detection curve (`Ar_finalE`, see below) does not fix or explain it either. Left as an open item; revisit if a better-motivated hypothesis emerges (e.g. inspecting the 2D density map shape directly for these specific clusters, rather than only the 1D extracted curve).

---

## Edge-detection cross-check (`compute_edge`): doesn't help, including when averaged

`compute_edge=TRUE` computes a second, independent amplitude curve (`Ar_finalE`) from the empirical extremes of raw galaxy velocities per radial bin, rather than the phase-space density map `Ar_finalD` is drawn from. As already noted in the code's own documentation: alone, it's less precise than the main curve (77% vs. 90% of 100 real Tempel clusters within a factor ~2 of the true mass). Re-tested here specifically for whether *averaging* the two helps, two ways, on 60 real Tempel clusters:

| | R200 bias | R200 sd | M200 bias (dex) | M200 sd |
|---|---|---|---|---|
| Main curve only | -1.9% | 0.292 | -0.001 | 0.369 |
| Average of final R200/M200 results | -7.0% | 0.301 | -0.073 | 0.386 |
| Average of the two curves, then integrate | -6.9% | 0.303 | -0.073 | 0.385 |

Both averaging approaches give essentially the same (worse) answer -- makes sense, since integration is close to linear in the curve. Confirms the main curve alone remains the best available choice; the two curves are not independent enough of each other (both drawn from the same limited set of candidate contours) for averaging to cancel out noise the way it would for two genuinely independent estimates.

---

## Convergence: recovering clusters with too few galaxies

The single most common cause of `run_caustic()` failing to converge in blind mode turned out **not** to be the radial-contour issue `run_caustic_robust()` already handled, but a distinct, more common one: too few galaxies within the starting radial window to even attempt a preliminary R200 estimate. Previously this failed silently -- the resulting `NA` propagated several steps downstream into `gaussian_kernel()`, which crashed with an opaque `"NA/NaN argument"` error carrying no hint of the real cause, and `run_caustic_robust()` had no way to recognise or retry it.

Fixed in two parts: (1) `run_caustic()` now stops immediately with a specific, greppable error (`"too few galaxies for a blind preliminary R200 estimate"`) as soon as this happens, instead of letting the `NA` propagate; (2) `min_n` (the galaxy-count threshold, default 15) is now an exposed parameter, and `run_caustic_robust()` recognises this failure mode and retries automatically with a smaller `min_n` (12, 10, 8, 6 in turn), the same way it already retried `rlimit`. The recovered result carries `min_n_used` (analogous to the existing `rlimit_frac_used`) and the same "treat as lower-confidence" caveat.

**Tested against the full Tempel et al. (2017) sample, both extractions combined (n=200 clusters)**:

| | Convergence |
|---|---|
| `run_caustic()` alone | 95/200 (47.5%) |
| `run_caustic_robust()`, with this fix | 130/200 (65.0%) |

Of the 105 clusters that failed to converge with a plain `run_caustic()` call, 61 (58%) failed specifically due to this too-few-galaxies issue -- more common than the radial-contour failure (33 clusters, 31%) that `run_caustic_robust()` already handled before this fix. Of those 61, the `min_n` retry recovered 24 (39%). Note that not every `min_n`-type failure is recoverable this way: for genuinely tiny samples (single digits after the initial radius/velocity cuts), lowering `min_n` can let the preliminary estimate itself succeed while a *different* failure still occurs further downstream (e.g. `findcontours()` itself needs enough points to trace a contour) -- there is a hard floor below which no amount of threshold-lowering helps.

### Recovered results are measurably less accurate, not just lower-confidence in principle

Tested on the Tempel extended sample (unconstrained extraction, n=100): comparing the 60 clusters that converged with a plain `run_caustic()` call against the 10 additional ones recovered by `run_caustic_robust()` (8 via `rlimit`, 2 via `min_n`):

| | n | R200 bias | R200 sd | M200 bias (dex) | M200 sd |
|---|---|---|---|---|---|
| Converged normally | 60 | -1.9% | 0.292 | -0.001 | 0.369 |
| **Recovered by `run_caustic_robust()`** | 10 | **+8.5%** | **0.375** | **+0.130** | **0.402** |

Worse in every metric -- confirming the function's own "treat as lower-confidence" warning with a concrete number, not just as a general caution. n=10 is small, so treat the exact magnitude loosely, but the direction is consistent with everything else found here.

### Richness threshold: N200 < 6 is where things break down

Using `N200` (galaxies inside the *estimated* R200 -- computed from `caustic_outliers` and `r200_est`, not returned directly by `run_caustic()` but easy to derive) as the richness measure, tested across the full Tempel sample (both extractions, n=130 converged, `run_caustic_robust()` included):

| N200 | n | R200 bias | R200 sd | M200 bias (dex) | M200 sd |
|---|---|---|---|---|---|
| **0-5** | 16 | **-29.9%** | 0.229 | **-44.3%** | **0.601** |
| 6-10 | 69 | -7.3% | 0.302 | -7.7% | 0.363 |
| 11-20 | 28 | -14.8% | 0.158 | -18.2% | 0.216 |
| 21-40 | 13 | -8.9% | 0.256 | -9.4% | 0.293 |
| >40 | 4 | -12.9% | 0.112 | -15.4% | 0.153 |

There's a clear break at `N200 = 6`: below it, both bias and scatter are substantially worse than any richer bin (M200 scatter in particular is roughly double the next-worst group). Above it, precision settles into a broadly stable range without a strong further trend -- crossing this floor matters much more than how far above it a cluster sits. `N_total` (input candidate count, before fitting) was tested too but showed a much less clean pattern, likely confounded with extraction geometry (very high-N_total clusters are almost all from the wider unconstrained extraction, which has its own, separate scatter penalty unrelated to richness); `N200` is the more direct and more useful diagnostic since it's a property of the fit itself. **Practical takeaway**: treat any `run_caustic()` result with `N200 < 6` as low-confidence regardless of convergence or its own reported error bar, the same way results recovered by `run_caustic_robust()` already are.

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

---

## v1.0 vs. current: what the accumulated changes actually bought

All of the above compares parameter choices within already-modified code. A different, more basic question: how does the true original (`v1.0`, github.com/Dailer/RCausticMass release tag, 281 lines, a single `run_caustic()` function with everything inline -- no `run_caustic_robust()`, `mirror=FALSE` default, hardcoded `rlimit=5.8`/`xmax=6`, no gradient restriction on contour tracing, linear (not log-spaced) κ levels, `fbr=0.65` default) compare to the current file, tested with each version's own defaults on the same 200 real Tempel et al. (2017) clusters (both extractions)?

| | Convergence (`run_caustic()` alone) | R200 bias | R200 sd | M200 bias | M200 sd |
|---|---|---|---|---|---|
| v1.0 (original) | **108/200 (54.0%)** | -4.5% | 0.266 | **-12.0%** | 0.355 |
| Current | 95/200 (47.5%) | -1.8% | 0.255 | **-0.0%** | 0.332 |

**The original converges more often but is measurably less accurate when it does.** This isn't a contradiction -- it's the direct, explainable consequence of the bug fixes and default changes documented throughout this file: with no gradient restriction on contour tracing, more candidate contours qualify as "valid" (nothing constrains how fast their amplitude can grow or fall), making it easier to find *some* contour that satisfies `findcontours()`'s span requirements -- at the cost of that contour being a worse fit. The fixed `rlimit=5.8` (vs. the current adaptive `0.8×xmax`) crops the wide Unconstrained extraction more aggressively, which also tends to help raw convergence by excluding more of the noisy outer field, again without the calibration behind the current adaptive choice. `mirror=FALSE`'s effect on convergence specifically is small (as found elsewhere in this document); its cost is almost entirely to precision.

Once `run_caustic_robust()` (which doesn't exist at all in v1.0) is added on the current version, the comparison isn't close: 130/200 (65.0%) vs. v1.0's 108/200 (54.0%), with better precision in every converged case besides. The lesson: convergence rate alone is not a sufficient metric to judge a change by -- several of the changes that *reduced* raw convergence (gradient restriction, adaptive rlimit, mirror) were net improvements once precision is weighed too, and had to be explicitly compensated for with `run_caustic_robust()`'s recovery logic rather than by reverting them.

---

## Conclusion

**What actually changed from where this document's calibration record starts** (already past several earlier fixes -- radial mirroring, adaptive `rlimit`/`xmax`, log-spaced κ, gradient restriction, `mirror=TRUE` as default -- all already in place by that point): `blur_gaussian=FALSE` (Deriche instead of Young-van Vliet), `gradu=1.0` (was 0.5), and `run_caustic_robust()` gained a second recovery path (`min_n` search) alongside its original `rlimit` search. That's the complete list of default changes made *during* the process this file documents. Everything else tested here -- dozens of parameters, several genuinely good ideas among them -- either made no difference or made things worse, and was removed from the code again after testing rather than left in as unused options (see the note at the top of "Predicting dynamical state..." above for the specific list). For the full picture back to the true original release (`v1.0`), including what all of those earlier fixes bought together, see "v1.0 vs. current" just above.

**What the calibration numbers mean in practice**: `fbr` is the single most consequential number in this pipeline, and it is not one number -- it depends on extraction geometry (proportional-to-R200 vs. fixed-radius), and, more subtly, on cluster mass itself. Neither dependency is fully resolved into an automatic correction; both are documented as things to be aware of and recalibrate for when the context changes, not solved problems.

**What this process confirmed works well**: the core method, once these two default changes and the convergence-recovery improvement are in place, recovers M200 to within the precision the D99/Serra et al. (2011) uncertainty formula itself predicts, for reasonably rich (N200 ≥ 6), reasonably relaxed clusters, blind or informed. That's a real, validated result across hundreds of real clusters in four independent samples (Tempel Constrained, Tempel Unconstrained, CIRS, plus individual literature clusters).

**What this process could not resolve, honestly**: (1) why a minority of clusters show a declining escape curve near the centre -- open; (2) whether the mass-dependent `fbr` trend reflects a real physical effect (mass-dependent orbital anisotropy) or is partly an artifact of comparing against reference masses that themselves have mass-dependent systematics -- likely some of both, impossible to fully disentangle without a reference technique known to be bias-free, which doesn't exist; (3) how to turn the confirmed dynamical-state sensitivity of `mirror` into something usable without an external relaxed/unrelaxed classification, since the one automatic substructure test tried (Dressler-Shectman) that worked on curated data failed on noisier real extractions.

**Update**: the first of these was since attempted -- a genuine ground-truth test against 43 halos from the TNG300 cosmological simulation (see "Validation against a cosmological simulation" above). It substantially revised two earlier conclusions (the mass-`fbr` trend mostly disappears against true mass; the `mirror`/dynamical-state connection doesn't reliably replicate against true substructure and turns out to be viewing-angle-dependent), and surfaced an entirely new, arguably more important finding that wasn't on this list at all: line-of-sight projection angle alone contributes as much scatter as everything else in this document combined. That result reframes how much further precision gain should even be expected from continued parameter tuning -- a real, informative answer, if a humbling one.

**If extending this work further**, the highest-value next steps based on everything above would be: a larger TNG-based sample (43 halos was enough to reveal these effects but not enough to fully re-derive `fbr` or the richness/precision relationship with confidence); systematically averaging multiple viewing angles per real cluster where possible (e.g. via independent X-ray or lensing-based orientation constraints) to see whether *known* projection effects can be corrected for after the fact; a larger-sample, better-powered version of the mass-`fbr` investigation on real data now that the simulation result suggests the earlier trend was partly a reference-mass artifact; and a closer look at the 2D density map shape itself (not just the 1D extracted curve) for the unresolved "falling near centre" cases, which the 1D diagnostics used throughout this document weren't able to explain.
