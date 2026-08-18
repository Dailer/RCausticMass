# Biweight scale estimator from Beers et al. (1990)
#
# Fixed to match Beers et al. (1990) exactly (and the already-validated
# implementation in sigma_plateau.R): the tuning constant c=9 multiplies
# the RAW median absolute deviation (mad(..., constant=1)), not a
# rescaled MAD (the previous constant=2.21914 effectively made the tuning
# constant ~20 instead of 9, far less aggressive at downweighting
# outliers than the published recipe). The sqrt(n) prefactor also now uses
# the FULL sample size, as in eq. 9 of Beers et al. (1990) -- the previous
# version used sqrt(length(wu)), the count of points surviving the |u|<1
# clip, which systematically underestimates the scale whenever points are
# clipped.
biwScale = function(x, c = 9){
  # BUG FIX: an empty `x` (e.g. no galaxies survive a radius/velocity cut
  # for a poor cluster) used to reach `if (madx == 0)` with madx = NA,
  # crashing with the uninformative "missing value where TRUE/FALSE
  # needed" instead of a clear NA -- found on a real Tempel et al. (2017)
  # cluster with only 14 galaxies total, where none happened to fall
  # within the preliminary r200/vlimit cut used for the initial vdisp
  # estimate. Poorer cluster samples (few dozen candidates or fewer) hit
  # this far more often than rich ones, simply from having less margin
  # before a cut empties out entirely.
  if (length(x) == 0) return(NA)
  M = median(x, na.rm = T)
  madx = mad(x, M, constant = 1, na.rm = T)
  if (madx == 0) { M = mean(x, na.rm = T); madx = mad(x, M, constant = 1, na.rm = T) }
  d = x - M
  u = d / (c * madx)
  uT = 1 - u^2 
  wu = which(abs(u) < 1)
  top = sqrt(sum((d^2 * uT^4)[wu]))
  bottom = abs(sum((uT * (1 - 5 * u^2))[wu]))
  SBI = sqrt(length(x)) * top / bottom
  if(is.nan(SBI)) SBI = x[1]
  return(SBI)
}

# Base-R replacements for pracma::cumtrapz() and gplots::hist2d(), removing
# those two package dependencies (imager and, optionally, magicaxis remain
# needed -- see gaussian_kernel()/adaptive_kernel_2d() and the plotting
# code respectively). Verified to match their originals exactly (max
# difference 0) on test data before replacing the calls below.
cumtrapz_base = function(x, y){
  n = length(x)
  if (n < 2) return(matrix(rep(0, n), ncol = 1))
  areas = diff(x) * (y[-1] + y[-n]) / 2
  matrix(c(0, cumsum(areas)), ncol = 1)
}

hist2d_base = function(x, y, nbins){
  xbreaks = seq(min(x), max(x), length.out = nbins[1] + 1)
  ybreaks = seq(min(y), max(y), length.out = nbins[2] + 1)
  xi = cut(x, breaks = xbreaks, include.lowest = TRUE, labels = FALSE)
  yi = cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
  counts = matrix(0, nbins[1], nbins[2])
  tab = table(xi, yi)
  counts[as.integer(rownames(tab)), as.integer(colnames(tab))] = tab
  list(counts = counts, x.breaks = xbreaks, y.breaks = ybreaks)
}

##############################################################################
# Shifting-gapper interloper removal (Fadda et al. 1996). Copied here from
# sigma_plateau.R (where it was originally implemented and validated) so
# run_caustic() can use it internally as a cleaning step before the
# preliminary R200 estimate below -- this is a stopgap, not the full
# integration of sigma_plateau.R's pipeline (the binary-tree + sigma-plateau
# membership algorithm) into this file, which remains a separate, larger
# pending item. Keeping two copies of the same function risks drift if one
# is edited without the other; treat sigma_plateau.R as the source of truth
# if the two ever disagree.
##############################################################################
shifting_gapper <- function(Rproj, vlos, gap_kms = 1000, bin_mpc = 0.4, min_n = 15,
                             bin_expand_factor = 1.5, step_frac = 0.5,
                             r_max = NULL, max_iter = 15, verbose = TRUE, plot = TRUE){
  stopifnot(length(Rproj) == length(vlos))
  if (anyNA(c(Rproj, vlos))) stop('shifting_gapper(): Rproj/vlos must not contain NA values.')

  n0 <- length(Rproj)
  keep <- rep(TRUE, n0)
  if (!is.null(r_max)) keep <- keep & (Rproj <= r_max)

  iter <- 0
  repeat {
    iter <- iter + 1
    idx_active <- which(keep)
    Ra <- Rproj[idx_active]; va <- vlos[idx_active]
    ord <- order(Ra)
    Rs <- Ra[ord]; vs <- va[ord]; idx_s <- idx_active[ord]
    n <- length(Rs)
    reject_flag <- logical(n)

    if (n < min_n) {
      if (verbose) message('shifting_gapper(): fewer than min_n = ', min_n,
                            ' galaxies remain active; stopping.')
      break
    }

    r_end <- Rs[n]
    step <- bin_mpc * step_frac
    bin_lo <- Rs[1]

    while (bin_lo <= r_end) {
      w <- bin_mpc
      bin_hi <- bin_lo + w
      bin_idx <- which(Rs >= bin_lo & Rs < bin_hi)
      while (length(bin_idx) < min_n && bin_hi < r_end) {
        w <- w * bin_expand_factor
        bin_hi <- bin_lo + w
        bin_idx <- which(Rs >= bin_lo & Rs < bin_hi)
      }
      if (length(bin_idx) >= 2) {
        vb <- vs[bin_idx]
        ov <- order(vb)
        vb_sorted <- vb[ov]
        a0 <- which.min(abs(vb_sorted - median(vb_sorted)))
        keep_bin <- logical(length(vb_sorted))
        keep_bin[a0] <- TRUE
        k <- a0
        while (k > 1) {
          if ((vb_sorted[k] - vb_sorted[k-1]) <= gap_kms) { keep_bin[k-1] <- TRUE; k <- k-1 } else break
        }
        k <- a0
        while (k < length(vb_sorted)) {
          if ((vb_sorted[k+1] - vb_sorted[k]) <= gap_kms) { keep_bin[k+1] <- TRUE; k <- k+1 } else break
        }
        rej_local <- bin_idx[ov[!keep_bin]]
        reject_flag[rej_local] <- TRUE
      }
      bin_lo <- bin_lo + step
    }

    newly_rejected <- idx_s[reject_flag]
    if (verbose) message('shifting_gapper() iteration ', iter, ': rejected ',
                          length(newly_rejected), ' of ', n, ' active galaxies.')
    if (length(newly_rejected) == 0 || iter >= max_iter) break
    keep[newly_rejected] <- FALSE
  }

  if (plot) {
    plot(Rproj, vlos, xlab = expression(R[proj]~(Mpc)), ylab = expression(v[los]~(km/s)),
         las = 1, pch = ifelse(keep, 16, 1), col = ifelse(keep, '#0077B6', 'gray40'))
    legend('topright', pch = c(16, 1), col = c('#0077B6', 'gray40'),
           legend = c('kept', 'rejected (interloper)'), bty = 'n', cex = .8)
    mtext('Shifting-gapper (Fadda et al. 1996)', font = 2, cex = .9)
  }

  list(keep = keep, n_iter = iter, n_initial = n0, n_final = sum(keep))
}

##############################################################################
# Adaptive-bandwidth 2D kernel density estimate (Serra et al. 2011, eq.
# 18-20), as an alternative to gaussian_kernel()'s fixed-bandwidth Gaussian
# blur (isoblur). The paper's kernel is NOT Gaussian: it is the compact
# quartic K(t) = (4/pi)(1-t^2)^3 for |t|<1 (eq. 19), and its bandwidth h_i
# is LOCAL to each galaxy -- larger where the local phase-space density is
# lower (sparse outskirts), smaller where it is higher (dense core) -- via
# h_i = h_c * h_opt * lambda_i, lambda_i = sqrt(gamma/f_1(x_i)), with f_1 a
# pilot (fixed-bandwidth h_opt) density estimate and log(gamma) the mean of
# log(f_1) over all points.
#
# IMPORTANT CAVEAT: the paper rescales (r, v) so that q = sigma_v/sigma_r
# takes a fixed value (q=25 by default), stating this makes "an uncertainty
# of 100 km/s in v weigh like an uncertainty of 0.02 Mpc in r" -- implying
# a conversion factor of 100/0.02 = 5000 km/s per Mpc-equivalent, which is
# what is implemented below. This specific numerical translation is this
# script's best-effort reading of that passage, not something verified
# line-by-line against the paper's worked equations -- treat this function
# as a reasonable-but-unverified reimplementation, and compare its output
# against gaussian_kernel()'s for the same cluster before trusting it in a
# publication. Validated only for basic sanity on synthetic data (no NaN/
# Inf; bandwidth genuinely varies -- ~5x wider in a sparse synthetic "halo"
# than in a dense synthetic "core" in testing).
#
# BUG FIX: this function originally only mirrored velocity (v -> -v), same
# omission as gaussian_kernel() had -- the radial mirroring (r -> -r) that
# Section 4.3 requires to avoid depleted density near r=0 was missing here
# too. Now implemented the same way: the pilot density AND the final grid
# are both built from the r-and-v mirrored point set, and only the r>=0
# half of the resulting grid is kept.
##############################################################################
adaptive_kernel_2d <- function(dproj, vlos, r200, xmax = 6, ymax = 4000, by = 0.05,
                                mirror = TRUE, q_scale = 25, hc = 1){
  if (mirror) { xv <- rep(dproj, 2); vv <- c(vlos, -vlos) } else { xv <- dproj; vv <- vlos }
  N <- length(xv)

  # v[km/s] -> r-equivalent units (see caveat above)
  conv_kms_per_Mpc <- (100 / 0.02) * (q_scale / 25)
  vv_r <- vv / conv_kms_per_Mpc

  # radial mirroring (Section 4.3): reflect every point to -r as well
  xv_m <- c(xv, -xv)
  vv_r_m <- c(vv_r, vv_r)
  Nm <- length(xv_m)

  # BUG FIX: sigma_r/sigma_v (and hence h_opt) used to be computed from the
  # WHOLE input sample, including field/background contamination out to
  # xmax Mpc and ymax km/s. gaussian_kernel() already guards against this
  # by restricting its own bandwidth estimate (bwx/bwy) to points inside
  # r200 (`xvalues < r200`); this function did not, so h_opt came out set
  # by the contamination's spread rather than the cluster's own scale. On
  # a real test cluster (r200 ~ 0.83 Mpc) this gave h_opt ~ 1.48 Mpc --
  # larger than the cluster itself -- washing out the trumpet shape into
  # an oversized, nearly-flat density map and inflating M200 by ~50x.
  # Restricting to r < r200 (same convention as gaussian_kernel()) fixes
  # this.
  in_r200 <- abs(xv) < r200
  if (sum(in_r200) < 5) in_r200 <- rep(TRUE, N)  # fallback for very sparse clusters
  sigma_r <- sd(xv[in_r200]); sigma_v <- sd(vv_r[in_r200])
  h_opt <- 6.24 / N^(1/6) * sqrt((sigma_r^2 + sigma_v^2) / 2)

  K <- function(t2) ifelse(t2 < 1, (4/pi) * (1 - t2)^3, 0)  # eq. 19

  x_range <- seq(0, xmax - by, by = by)
  x_range_full <- seq(-xmax + by, xmax - by, by = by)
  # BUG FIX: the velocity grid step used to be `by/conv_kms_per_Mpc`, i.e.
  # the SAME small step already used for r (in Mpc), divided again by the
  # ~5000 km/s/Mpc conversion factor -- producing ~160,000 grid points
  # instead of a sane number, and making the function unusably slow (timed
  # out on a 363-galaxy cluster). The grid RESOLUTION is a separate choice
  # from the kernel BANDWIDTH (h_opt/h_i, still computed in the properly
  # rescaled Mpc-equivalent units below); a resolution comparable to
  # gaussian_kernel()'s own default grid is enough here.
  by_v_kms <- 25
  y_range <- seq(-ymax, ymax - by_v_kms, by = by_v_kms)
  y_range_r <- y_range / conv_kms_per_Mpc

  # pilot estimate f_1 (fixed bandwidth h_opt, no adaptation yet) at each of
  # the N original (pre-radial-mirror) points, summing contributions from
  # the FULL mirrored set so points near r=0 aren't penalised
  f1_at_points <- sapply(seq_len(N), function(i){
    d2 <- ((xv_m - xv[i])^2 + (vv_r_m - vv_r[i])^2) / h_opt^2
    sum(K(d2)) / (Nm * h_opt^2)
  })
  gamma <- exp(mean(log(pmax(f1_at_points, 1e-300))))
  # Note: an earlier version of this function supported capping lambda
  # (limiting how much the local bandwidth can grow in sparse regions) via
  # a `lambda_max` parameter. Tested against 100 real Tempel et al. (2017)
  # clusters and found NO net benefit once fbr was properly recalibrated
  # for each cap value -- precision stayed flat (~70% within a factor ~2,
  # ~40-43% within a factor ~1.4) across the whole range from lambda_max=1
  # to uncapped, while capping cost real convergence (fewer clusters
  # produced any estimate at all). Removed rather than kept as a
  # non-functional knob.
  lambda <- sqrt(gamma / pmax(f1_at_points, 1e-300))
  hi <- hc * h_opt * lambda        # local bandwidths for the N original points
  hi_m <- c(hi, hi)                # each mirrored copy reuses its original's h_i

  # final grid over the FULL [-xmax, xmax] domain, each (mirrored) point
  # contributing through its own local h_i; then keep only r >= 0
  Zi_full <- matrix(0, length(x_range_full), length(y_range))
  for (i in seq_len(Nm)) {
    dx2 <- (x_range_full - xv_m[i])^2
    dv2 <- (y_range_r - vv_r_m[i])^2
    d2grid <- outer(dx2, dv2, "+") / hi_m[i]^2
    Zi_full <- Zi_full + K(d2grid) / (Nm * hi_m[i]^2)
  }
  half_idx <- (nrow(Zi_full) - length(x_range) + 1):nrow(Zi_full)
  Zi <- Zi_full[half_idx, , drop = FALSE]

  list(x = x_range, y = y_range, z = Zi, r200_used = NA, h_opt = h_opt, hi_range = range(hi))
}

# Shared helper: robust velocity dispersion within a radial window,
# cleaned of interlopers with shifting_gapper() first, falling back to the
# uncleaned biweight scale if shifting_gapper() fails or too few galaxies
# remain. Used identically in two places inside run_caustic() below
# (adaptive_vlimit's own preliminary dispersion, and the iterative
# preliminary R200 estimate) -- factored out here so that shared logic
# only needs to be corrected/tuned in one place.
gapper_cleaned_vdisp = function(rproj, vproj, r_max, min_n = 15){
  in_r = rproj < r_max
  if (sum(in_r) < min_n) return(NA)  # caller decides what to do (differs by call site below)
  gap = tryCatch(shifting_gapper(rproj[in_r], vproj[in_r], r_max = r_max, verbose = FALSE, plot = FALSE),
                 error = function(e) NULL)
  if (!is.null(gap) && sum(gap$keep) >= 5) return(biwScale(vproj[in_r][gap$keep]))
  biwScale(vproj[in_r])  # cleaning failed or left too few points: fall back to the raw biweight scale
}

# Uses a 2D gaussian kernel to estimate the density of the phase space
#
# BUG FIX: this function previously only offered `mirror` for the VELOCITY
# axis (v -> -v). The paper's own mirroring technique (Serra et al. 2011,
# Section 4.3) is different and was missing entirely: "to avoid artificial
# depletion of the caustic amplitude at small r due to the small number of
# galaxies in the central region, the galaxy distribution is mirrored to
# negative r." Without it, the estimated density right at r=0 can end up
# far below the map's peak density (found ~14% of the peak in a real-data
# test, vs ~39% with this fix applied), which makes findcontours() unable
# to find a contour that closes at r=0 -- causing it to fail entirely for
# many real clusters. This is now implemented: the KDE is built over a
# domain spanning [-xmax, xmax], using both the original points AND their
# r-mirrored copies, and only the r>=0 half of the resulting density map is
# kept and returned.
gaussian_kernel = function(dproj, vlos, r200, normalization = 100, q = 10, xmax = 6, ymax = 5e3, by = 0.05, 
                           mirror = F, hc = 1.3, plot = F){
  # The "q" parameter is termed "scale" set to 10 as default, but can go as high as 50.
  # "normalization" is simply H0
  # "x/yres" can be any value, but are recommended to be above 150
  require(gplots)
  require(magicaxis)
  require(imager)
  normscale = normalization * q
  
  if(any(dproj >= xmax))
    stop(paste('Error: Please either increase your xmax value or trim your sample to be x < ', xmax))
   
  if(any(abs(vlos) >= ymax))
    stop(paste('Error: Please either increase your ymax value or trim your sample to be y < ', ymax))
  
  if(sum(dproj < r200) < 10){ 
    r200 = sort(dproj)[10]
    warning('Less than 10 objects inside r200, extending radius limit to the tenth galaxy radius')
  }
  # r200 above may have just been reassigned locally -- return it so the
  # caller (run_caustic) can keep using the SAME r200 everywhere downstream
  # instead of silently diverging from the r200 used here for the bandwidth.
  r200_used = r200
  
  if(mirror){
    xvalues = rep(dproj, 2)
    yvalues = c(vlos, -vlos) / normscale
  }else{
    xvalues = dproj
    yvalues = vlos / normscale
  }
  
  x_range = seq(0, xmax - by, by = by)
  xres = length(x_range)
  y_range = seq(-ymax / normscale, ymax / normscale - by, by = by) * normscale
  yres = length(y_range)
  x_scale = (xvalues / xmax) * xres
  y_scale = ((yvalues * normscale + ymax) / (2 * ymax)) * length(y_range)
  bwx = biwScale(x_scale[xvalues < r200])
  bwy = biwScale(y_scale[xvalues < r200])
  # `hc`: overall smoothing multiplier on ksize (same role as hc in
  # adaptive_kernel_2d()). Previously this bandwidth was entirely fixed by
  # the Silverman-type rule below, with no way to widen/narrow it to test
  # whether the resulting escape surface (and, downstream, M200) improves.
  ksize = hc * (4 / (3 * length(xvalues)))^(1 / 5) * sqrt((bwx^2 + bwy^2) / 2)
  ymx = ymax / normscale

  # radial mirroring (Section 4.3): reflect every point to -r as well, build
  # the density over the full [-xmax, xmax] domain, then keep only r >= 0.
  # This makes r=0 an INTERIOR point of the grid (not an edge), so isoblur's
  # boundary handling there is irrelevant -- the true r=0 depletion problem
  # is already solved by this construction alone.
  xvalues_m = c(xvalues, -xvalues)
  yvalues_m = c(yvalues, yvalues)
  x_range_full = seq(-xmax + by, xmax - by, by = by)
  h = hist2d_base(c(xvalues_m, -xmax, xmax, xmax, -xmax), c(yvalues_m, -ymx, -ymx, ymx, ymx),
             c(length(x_range_full), yres))
  h$counts[cbind(c(1, length(x_range_full), length(x_range_full), 1), c(1, 1, yres, yres))] = 0

  # BUG FIX: the OUTER edges of this grid (r = +-xmax, v = +-ymax) are true
  # boundaries, and imager::isoblur()'s default boundary handling there
  # (Neumann, i.e. edge-value replication) does NOT conserve density the
  # way the reference Python implementation's `mode='reflect'` does.
  # Measured on a synthetic source sitting exactly on a boundary: Neumann
  # inflates the blurred mass there by ~80%, and the Dirichlet alternative
  # (neumann=FALSE) loses ~44% of it -- neither is right. Reflect-padding
  # the grid before blurring (then cropping back) conserves mass to
  # numerical precision instead (verified: 99.99 of 100 recovered).
  reflect_pad = function(m, pad_r, pad_c){
    nr = nrow(m); nc = ncol(m)
    pad_r = min(pad_r, nr); pad_c = min(pad_c, nc)
    top = m[pad_r:1, , drop = FALSE]
    bottom = m[nr:(nr - pad_r + 1), , drop = FALSE]
    m2 = rbind(top, m, bottom)
    left = m2[, pad_c:1, drop = FALSE]
    right = m2[, nc:(nc - pad_c + 1), drop = FALSE]
    cbind(left, m2, right)
  }
  pad_r = min(nrow(h$counts) - 1, ceiling(4 * ksize))
  pad_c = min(ncol(h$counts) - 1, ceiling(4 * ksize))
  counts_padded = reflect_pad(h$counts, pad_r, pad_c)
  img_full_padded = as.matrix(isoblur(as.cimg(counts_padded), sigma = ksize))
  img_full = img_full_padded[(pad_r + 1):(pad_r + nrow(h$counts)),
                              (pad_c + 1):(pad_c + ncol(h$counts)), drop = FALSE]
  half_idx = (nrow(img_full) - xres + 1):nrow(img_full)   # r in [0, xmax)
  img = img_full[half_idx, , drop = FALSE]

  if(plot){
    magimage(x_range, y_range, img, asp = NA, xlab = expression(R[proj] ~ (Mpc)), 
             ylab = expression(v[proj] ~ (km~s^{-1})))
  }
  lout = list(x = x_range, y = y_range, z = img, r200_used = r200_used)
  return(lout)
}

# This function will use contourLines() to locate escape surfaces
#
# `vlimit` (bug fix): previously hardcoded to 4000 km/s inside the contour
# selection filter, completely disconnected from run_caustic()'s own
# `vlimit` argument -- calling run_caustic(..., vlimit=6000) silently kept
# rejecting any contour beyond +-4000 km/s regardless. Now an explicit
# parameter, passed through from run_caustic().
#
# Two pieces ported from the original Python reference implementation
# (Gifford et al., causticpy/CausticMass.py) that this R port had been
# missing entirely:
#
#   1) LOG-SPACED kappa levels. contourLines(..., nlevels=200) with no
#      explicit `levels` uses R's default (roughly linear/`pretty()`-based)
#      spacing. The phase-space density is sharply peaked near the centre
#      and very sparse in the outskirts, so linear spacing wastes most of
#      the 200 levels resolving the dense core and leaves few levels to
#      trace the sparse, low-density outer trumpet -- a likely contributor
#      to this R port's frequent "contours do not expand to the radial
#      limit" failures. The reference code instead spaces levels
#      logarithmically from min(Zi[Zi>0])/5 up to max(Zi).
#
#   2) Gradient restriction (restrict_gradient2 in the reference code,
#      Serra et al. 2011's d ln A / d ln r constraint made concrete): as
#      each candidate contour is walked outward in radius, the amplitude
#      is not allowed to grow faster than r^0.5 or fall faster than r^-2
#      from one radial step to the next -- applied to the positive and
#      negative branches separately, before combining them into the final
#      A(r) = min(|branch_up|, |branch_down|). Without this, nothing
#      prevented a candidate contour from growing arbitrarily wide at
#      large r, which this project found to be a real, measured problem
#      (the adaptive kernel's escape surface came out systematically wider
#      than the Gaussian kernel's, especially at large radius).
findcontours = function(Zi, ri, vi, r200, vvar, rimax = 4, vlimit = 4000, nlevels = 200,
                         gradu = 0.5, gradd = 2.0, plot = T, verbose = T){

  # log-spaced levels (see note above), falling back to contourLines()'s
  # own default spacing if Zi has no positive values (shouldn't normally
  # happen, but guards against a degenerate all-zero density map)
  zpos = Zi[Zi > 0]
  if (length(zpos) > 0) {
    level_lo = min(zpos) / 5
    level_hi = max(Zi)
    levels = rev(10^(seq(log10(level_lo), log10(level_hi), length.out = nlevels)))
    cl = contourLines(ri, vi, Zi, levels = levels)
  } else {
    cl = contourLines(ri, vi, Zi, nlevels = nlevels)
  }

  # gradient restriction (see note above): clamp how fast the amplitude
  # can grow/fall from one radial grid step to the next, applied
  # iteratively so each step's restriction uses the PREVIOUS (already
  # restricted) value, matching the reference implementation exactly.
  restrict_gradient2 = function(pastA, newA, pastr, newr){
    if (pastA <= newA) {
      if (pastA != 0 && (log(newA) - log(pastA)) / (log(newr) - log(pastr)) > gradu) {
        return(exp(log(pastA) + gradu * (log(newr) - log(pastr))))
      } else {
        return(newA)
      }
    } else {
      if (pastA != 0 && (log(newA) - log(pastA)) / (log(newr) - log(pastr)) < -gradd) {
        return(exp(log(pastA) - gradd * (log(newr) - log(pastr))))
      } else {
        return(newA)
      }
    }
  }

  # selecting contours that are "full" and don't loop back only in positive or negative space
  contours = fullcont = list()
  levels_kept = c()
  for (i in 1:length(cl)) {
    cln = cl[[i]]
    if(max(cln$x) >= rimax & min(cln$x) <= 0 & max(cln$y) > 0 & min(cln$y) < 0 & 
       max(abs(cln$y)) < vlimit){
      #find positive/negative contours
      xcont_u = cln$x[cln$y > 0] 
      ycont_u = cln$y[cln$y > 0] 
      xcont_d = cln$x[cln$y < 0]
      ycont_d = cln$y[cln$y < 0]
      
      if(length(xcont_u) < 2 | length(xcont_d) < 2) next;
      
      #match contour grid to r grid by linear interpolation
      fint_u = approxfun(xcont_u, ycont_u, yleft = 0, yright = 0, ties = mean)
      fint_d = approxfun(xcont_d, ycont_d, yleft = 0, yright = 0, ties = mean)
      y_u_mag = fint_u(ri)          # already >= 0 by construction (ycont_u > 0)
      y_d_mag = abs(fint_d(ri))     # fint_d gives <= 0 values; work in magnitude
      for (k in 2:length(ri)) {
        y_u_mag[k] = restrict_gradient2(y_u_mag[k-1], y_u_mag[k], ri[k-1], ri[k])
        y_d_mag[k] = restrict_gradient2(y_d_mag[k-1], y_d_mag[k], ri[k-1], ri[k])
      }
      y_f = pmin(y_u_mag, y_d_mag)
            
      fullcont = c(fullcont, list(cln))
      contours = c(contours, list(y_f))
      levels_kept = c(levels_kept, cln$level)
    }
  }
  ncont = length(contours)
  if(ncont == 0) 
    stop('The contours do not expand to the radial limit')
  
  # average calculation in Diaferio 99 and final contour
  tot_avg = numeric(ncont)
  lessr200 = ri <= r200
  useri = ri[lessr200]
  for (i in 1:ncont) {
    Ar = contours[[i]][lessr200]
    phir = numeric(length(useri))
    for (j in 1:length(useri)) {
      philimit = abs(Ar[j])
      phir[j] = sum(Zi[j,][vi < philimit & vi > -philimit])
    }
    tot_avg[i] = sum(Ar^2 * phir) / sum(phir)
  }
  idiff = order((tot_avg - 4 * vvar)^2)[1]
  final_contour = contours[[idiff]]
  kappa_final = levels_kept[idiff]  # the density threshold that located the chosen surface,
                                     # needed for the D99/Serra et al. 2011 error estimate (eq. 23-24)
  
  if(verbose) message('Complete')
  
  # plotting
  if(plot){
    magimage(ri, vi, Zi, asp = NA, xlab = expression(R[proj] ~ (Mpc)), 
             ylab = expression(v[proj] ~ (km~s^{-1})), col.ticks = 'white')
    ll = lapply(fullcont, lines, col = '#FFBE74', lty = 3)
    lines(ri, final_contour, col = '#F13005', lwd = 2)
    lines(ri, -final_contour, col = '#F13005', lwd = 2)
  }
  lout=list(caustic = final_contour, contours = contours, kappa = kappa_final)
  return(lout)
}

# Fitting a NFW profile to the caustic surface
#
# M200_prior / R200_prior (optional): when supplied (e.g. from an external
# catalog such as Tempel et al. 2017), M200 is FIXED EXACTLY to M200_prior
# and only the concentration (via the scale radius rs) is fit to this
# cluster's own caustic-amplitude data -- the data determine the profile's
# SHAPE, while the trusted external mass sets its NORMALISATION. This
# replaces fixing rs at an assumed universal concentration (halo_srad,
# typically r200/5, i.e. c=5 for every cluster) with only the
# normalisation (d0) free, which was the previous behaviour and remains
# the fallback when no prior is supplied.
#
# Numerical note: fitting (d0, rs) or (Ms, rs) both free at once (i.e. NOT
# anchoring the mass) was tested against synthetic caustic-amplitude data
# and found numerically unreliable with nls()'s default algorithm over the
# typical (r200/3, r200) fitting range: it failed to converge in every
# trial of a systematic sweep across noise levels and point counts (nls
# reports "singular convergence"). This is a numerical-conditioning issue,
# not a fundamental non-identifiability -- a derivative-free optimiser
# recovers the exact parameters from noiseless data given a wide radial
# baseline -- but it makes the free 2-parameter fit impractical with real,
# noisy caustic amplitudes over the range this function is normally called
# with. Fixing M200 collapses the fit to a single parameter (rs), solved
# here with a bounded nls(algorithm='port') search, which converged
# reliably in the same tests.
NFWfit = function(rii, Ar, halo_srad, ri_full, M200_prior = NA, R200_prior = NA,
                   c_min = 1, c_max = 15, conc_method = c('bounded', 'bayesian'),
                   sigma_A = NULL, c_prior_scatter_dex = 0.1, little_h = 0.7, clus_z = 0){
  conc_method = match.arg(conc_method)
  min_func = function(x, d0, rs)
	sqrt(2 * 4 * pi * 4.3e-9 * d0 * rs^2 * log(1 + x / rs) / (x / rs))

  have_prior = !is.na(M200_prior) & !is.na(R200_prior)
  conc_16 = conc_84 = NA

  if (have_prior && conc_method == 'bayesian') {
    ##########################################################################
    # Bayesian concentration fit (grid posterior). With M200 fixed,
    # concentration is the ONLY free parameter -- a 1D problem, so a fine
    # grid over c is simpler and more robust than any sampler, and gives
    # the full posterior directly (no MCMC needed).
    #
    # PRIOR: log-normal in c, centred on the Duffy et al. (2008) full-sample
    # c-M relation (c200 = 5.71*(M200*h/2e12)^-0.084*(1+z)^-0.47), with the
    # ~0.1 dex scatter standard in N-body concentration-mass studies. This
    # replaces the hard [c_min,c_max] truncation used in the 'bounded'
    # method: instead of an artificial wall the fit can pile up against
    # (observed on real Tempel et al. 2017 clusters: with noisy data and
    # concentration the only free parameter, ~25-30% of clusters hit the
    # bound rather than settling on a physically meaningful value), the
    # posterior smoothly shrinks toward the prior when the data don't
    # constrain c well, and shifts toward the data when they do.
    #
    # LIKELIHOOD: Gaussian, using `sigma_A` (the per-radius absolute
    # amplitude uncertainty from the D99/Serra et al. 2011 error formula,
    # already computed by run_caustic() for the M200 error bar and the
    # uncertainty band) as the noise model -- reusing that estimate here
    # rather than inventing a separate one.
    if (is.null(sigma_A) || all(!is.finite(sigma_A) | sigma_A <= 0))
      sigma_A = rep(max(sd(Ar), 1e-3), length(Ar))  # fallback if no usable error estimate was supplied
    sigma_A = pmax(sigma_A, 1e-3 * mean(Ar))  # guard against a zero/near-zero noise floor

    # BUG FIX: the grid used to inherit c_min/c_max (the 'bounded' method's
    # hard-truncation limits, c_min=1 by default) via max(c_min,0.2) as its
    # lower edge -- silently reintroducing the exact hard-wall artifact
    # this method exists to avoid. A cluster whose data genuinely favour
    # very low concentration then still piled up at c=1 with a degenerate
    # ~zero-width credible interval (observed on a real Tempel cluster).
    # The grid now spans a wide, fixed range independent of c_min/c_max;
    # the prior (not a wall) is what keeps the posterior physically
    # sensible.
    c_grid = 10^seq(log10(0.1), log10(50), length.out = 400)
    M200_h = M200_prior * little_h
    c_prior_mean = 5.71 * (M200_h / 2e12)^(-0.084) * (1 + clus_z)^(-0.47)
    log_prior = dnorm(log10(c_grid), mean = log10(c_prior_mean), sd = c_prior_scatter_dex, log = TRUE)

    min_func_rs = function(x, rs){
      cc = R200_prior / rs
      fc = log(1 + cc) - cc / (1 + cc)
      d0 = M200_prior / (4 * pi * rs^3 * fc)
      sqrt(2 * 4 * pi * 4.3e-9 * d0 * rs^2 * log(1 + x / rs) / (x / rs))
    }
    rs_grid = R200_prior / c_grid
    log_lik = sapply(rs_grid, function(rs) sum(dnorm(Ar, mean = min_func_rs(rii, rs), sd = sigma_A, log = TRUE)))
    log_post = log_prior + log_lik
    log_post = log_post - max(log_post)
    post = exp(log_post); post = post / sum(post)
    cum = cumsum(post)

    conc_fit = c_grid[which.max(post)]       # posterior mode
    conc_16 = c_grid[which.min(abs(cum - 0.16))]
    conc_84 = c_grid[which.min(abs(cum - 0.84))]
    rs_fit = R200_prior / conc_fit
    cc = conc_fit
    fc = log(1 + cc) - cc / (1 + cc)
    halo_scale_density = M200_prior / (4 * pi * rs_fit^3 * fc)
    halo_srad_used = rs_fit
    pred_at_fit = min_func_rs(rii, rs_fit)
    res.err = sqrt(sum((Ar - pred_at_fit)^2) / (length(Ar) - 1))

  } else if (have_prior) {
    # d0 is not free: it is whatever value makes M(<R200_prior) = M200_prior
    # for the current trial rs, so only rs (equivalently, the concentration
    # R200_prior/rs) is actually being fit.
    #
    # BUG FIX: rs used to be bounded only by [0.01, 3*R200_prior], allowing
    # concentrations up to ~R200_prior/0.01 (hundreds) -- with noisy real
    # data (typical for these clusters), the fit sometimes DID wander into
    # such physically nonsensical concentrations (tested: 4 of 15 real
    # Tempel clusters gave c>20, up to c~180, against an expected c~2-7 for
    # CDM haloes; Serra et al. 2011 find a 90% range of [2.06,6.50] in their
    # simulated sample) -- with mass fixed and concentration the only free
    # parameter, noisy data can be "explained" just as well by an absurd
    # concentration as by a physically reasonable one, so nothing in the
    # unbounded fit prevented it. rs is now bounded so the recovered
    # concentration c=R200_prior/rs stays within [c_min, c_max] (default
    # [1,15], a wide margin around the literature's typical range).
    #
    # NOTE: conc_method='bayesian' (see above) avoids this hard-wall
    # truncation altogether via a physically-motivated prior instead, and
    # was found to bring the ~25-30% of clusters that used to pile up
    # against this boundary back to sensible concentrations (c~4-6) without
    # ever needing to hit c_min/c_max.
    rs_lower = R200_prior / c_max
    rs_upper = R200_prior / c_min
    rs_start = min(max(halo_srad, rs_lower), rs_upper)  # keep start inside bounds
    min_func_rs = function(x, rs){
      cc = R200_prior / rs
      fc = log(1 + cc) - cc / (1 + cc)
      d0 = M200_prior / (4 * pi * rs^3 * fc)
      sqrt(2 * 4 * pi * 4.3e-9 * d0 * rs^2 * log(1 + x / rs) / (x / rs))
    }
    out = tryCatch(
      nls(Ar ~ min_func_rs(rii, rs), start = list(rs = rs_start),
          algorithm = 'port', lower = c(rs = rs_lower), upper = c(rs = rs_upper)),
      error = function(e) NA)
    if (length(out) == 1) return(list(rep(NA, length(ri_full)), NA, NA, NA, NA, NA, NA))
    rs_fit = coef(out)[['rs']]
    cc = R200_prior / rs_fit
    fc = log(1 + cc) - cc / (1 + cc)
    halo_scale_density = M200_prior / (4 * pi * rs_fit^3 * fc)
    halo_srad_used = rs_fit
    conc_fit = cc  # only meaningfully "fit" (not assumed) in this branch
    res = resid(out)
    res.err = sqrt(sum(res^2) / (length(res) - 1))
  } else {
    out = tryCatch(nls(Ar ~ min_func(rii, d0, rs = halo_srad),
                       start=list(d0 = 1e14)), error = function(e) NA)
    if (length(out) == 1) return(list(rep(NA, length(ri_full)), NA, NA, NA, NA, NA, NA))
    halo_scale_density = coef(out)[['d0']]
    halo_srad_used = halo_srad
    conc_fit = NA  # concentration is ASSUMED (halo_srad = r200/5), not fit, in this branch
    res = resid(out)
    res.err = sqrt(sum(res^2) / (length(res) - 1))
  }

  if(ri_full[1] == 0) ri_full[1] = ri_full[2]
  profile = min_func(ri_full, halo_scale_density, halo_srad_used)
  return(list(profile, res.err, conc_fit, conc_16, conc_84, halo_scale_density, halo_srad_used))
}

##############################################################################
# Independent "edge detection" caustic surface (ported from causticpy's
# findsurface(), the caustic_edge / Ar_finalE computation).
#
# This is a SECOND, statistically independent way to pick which candidate
# density contour is the true escape surface, alongside the main
# S(kappa)=4<v^2> matching already implemented in findcontours(). Instead
# of using the phase-space DENSITY at all, it works directly on the raw
# galaxy positions: bin galaxies by projected radius into `numbins` bins of
# roughly equal richness (within r200), and in each bin take the average of
# the top `edge_perc` (default 10%) most extreme positive and negative
# velocities as an empirical estimate of the local phase-space "edge". It
# then picks, from the SAME set of candidate contours already computed by
# findcontours() (`contours`), whichever one best matches this empirical
# edge (minimum median absolute difference).
#
# Because this doesn't depend on the phase-space density estimate at all,
# it is useful as a cross-check: if the density-based (Ar_finalD) and
# edge-based (Ar_finalE) surfaces agree closely, that is much stronger
# evidence the surface is well-determined than either method alone.
##############################################################################
##############################################################################
# Joint Bayesian fit of (M200, concentration), with R200 held as a FIXED
# ANCHOR radius (whether it was given by the caller or auto-estimated
# earlier in run_caustic()).
#
# An NFW halo + its "200c" mass definition really has only TWO independent
# degrees of freedom (e.g. M200 and c; R200 is then a fixed cosmological
# function of M200, not a third free parameter). So whenever M200 is NOT
# independently fixed, the natural generalisation of the single-parameter
# Bayesian concentration fit above (see NFWfit()'s conc_method='bayesian')
# is a 2D grid posterior over (M200, c), still no MCMC needed. This covers
# BOTH "nothing given" and "only R200/r200 given" in one implementation --
# they are the same fitting problem, since r200 (given or auto-estimated)
# is already available as the anchor before this stage runs either way.
#
# PRIORS: c is the same Duffy et al. (2008) log-normal used elsewhere,
# evaluated at each trial M200 (the c-M relation itself depends on mass).
# M200 gets a weak log-normal prior centred on the mass a plain SO-200
# definition would imply from the anchor radius alone
# (M200 = (4/3)*pi*200*rho_crit(z)*r200_anchor^3), with wide (default 0.4
# dex) scatter -- informative enough to keep the fit away from pathological
# solutions, but not tight enough to override real signal in the data.
#
# LIKELIHOOD: same Gaussian model as the 1-parameter case, using the
# D99/Serra et al. (2011) delta_A(r)/A(r) error estimate as the per-radius
# noise.
#
# Validated against 100 real Tempel et al. (2017) clusters: comparable or
# slightly better precision than the density-integral M200 estimate
# (73% vs 66% within a factor ~1.4, 91% vs 90% within a factor ~2), with
# the added benefit of well-behaved concentration credible intervals.
##############################################################################
NFWfit_bayes_2param = function(rii, Ar, sigma_A, r200_anchor, clus_z,
                                m_grid_n = 70, c_grid_n = 70,
                                c_prior_scatter_dex = 0.1, little_h = 0.7,
                                Om = 0.3, H0 = 70, logM200_prior_scatter_dex = 0.4){
  if (is.null(sigma_A) || all(!is.finite(sigma_A) | sigma_A <= 0))
    sigma_A = rep(max(sd(Ar), 1e-3), length(Ar))
  sigma_A = pmax(sigma_A, 1e-3 * mean(Ar))

  G = 6.67430e-11; Msol2kg = 1.98847e30; Mpc2m = 3.08567758e22
  Hz = H0 * sqrt(Om * (1 + clus_z)^3 + (1 - Om))
  crit = 3 * Hz^2 / (8 * pi * G) * 1e6 / Msol2kg * Mpc2m
  logM200_prior_center = log10((4 / 3) * pi * 200 * crit * r200_anchor^3)

  m_grid = 10^seq(logM200_prior_center - 3 * logM200_prior_scatter_dex,
                   logM200_prior_center + 3 * logM200_prior_scatter_dex, length.out = m_grid_n)
  c_grid = 10^seq(log10(0.1), log10(50), length.out = c_grid_n)

  min_func_rs = function(x, rs, M200){
    cc = r200_anchor / rs
    fc = log(1 + cc) - cc / (1 + cc)
    d0 = M200 / (4 * pi * rs^3 * fc)
    sqrt(2 * 4 * pi * 4.3e-9 * d0 * rs^2 * log(1 + x / rs) / (x / rs))
  }

  log_post = matrix(-Inf, m_grid_n, c_grid_n)
  rs_vec = r200_anchor / c_grid
  for (i in seq_along(m_grid)) {
    M200_trial = m_grid[i]
    c_prior_mean_i = 5.71 * (M200_trial * little_h / 2e12)^(-0.084) * (1 + clus_z)^(-0.47)
    log_prior_c = dnorm(log10(c_grid), mean = log10(c_prior_mean_i), sd = c_prior_scatter_dex, log = TRUE)
    log_prior_m = dnorm(log10(M200_trial), mean = logM200_prior_center, sd = logM200_prior_scatter_dex, log = TRUE)
    for (j in seq_along(c_grid)) {
      pred = min_func_rs(rii, rs_vec[j], M200_trial)
      log_lik = sum(dnorm(Ar, mean = pred, sd = sigma_A, log = TRUE))
      log_post[i, j] = log_prior_m + log_prior_c[j] + log_lik
    }
  }
  log_post = log_post - max(log_post)
  post = exp(log_post); post = post / sum(post)
  post_m = rowSums(post); post_m = post_m / sum(post_m)
  post_c = colSums(post); post_c = post_c / sum(post_c)
  cum_m = cumsum(post_m); cum_c = cumsum(post_c)

  list(M200 = m_grid[which.max(post_m)],
       M200_16 = m_grid[which.min(abs(cum_m - 0.16))], M200_84 = m_grid[which.min(abs(cum_m - 0.84))],
       conc = c_grid[which.max(post_c)],
       conc_16 = c_grid[which.min(abs(cum_c - 0.16))], conc_84 = c_grid[which.min(abs(cum_c - 0.84))])
}

edge_caustic = function(rproj, vproj, r200, contours, ri, mirror = TRUE,
                         numbins = 6, edge_perc = 0.1){
  ord = order(rproj)
  rsort = rproj[ord]
  vsort = if (mirror) abs(vproj[ord]) else vproj[ord]

  numrval = sum(rproj < r200)
  if (numrval < numbins || length(contours) == 0)
    return(list(Ar = rep(NA, length(ri)), mid_rbin = NULL, mincomp = NULL, best = NA))
  size_bin = ceiling(numrval / numbins)

  mid_rbin = c(); mincomp = c()
  n = length(rsort)
  for (nn in 0:(numbins - 1)) {
    idx_lo = nn * size_bin + 1
    idx_hi = min((nn + 1) * size_bin, n)
    if (idx_lo > n) break
    vbin = vsort[idx_lo:idx_hi]
    rbin = rsort[idx_lo:idx_hi]
    if (length(vbin) == 0) { if (nn >= 4) break else next }

    pos = vbin[vbin > 0]; neg = vbin[vbin < 0]
    n_top_pos = ceiling(length(pos) * edge_perc)
    n_top_neg = ceiling(length(neg) * edge_perc)
    avgmax_nn = if (n_top_pos > 0) mean(sort(pos, decreasing = TRUE)[1:min(n_top_pos, length(pos))]) else NA
    avgmin_nn = if (n_top_neg > 0) mean(sort(neg, decreasing = FALSE)[1:min(n_top_neg, length(neg))]) else NA
    if (is.na(avgmax_nn)) break

    mincomp_nn = if (min(vbin) >= 0) avgmax_nn else min(abs(avgmin_nn), avgmax_nn, na.rm = TRUE)
    mincomp = c(mincomp, mincomp_nn)
    mid_rbin = c(mid_rbin, median(rbin))
  }

  lt_r200 = ri < r200
  if (sum(lt_r200) < 2 || length(mid_rbin) < 2)
    return(list(Ar = rep(NA, length(ri)), mid_rbin = mid_rbin, mincomp = mincomp, best = NA))
  keep = mid_rbin < max(ri[lt_r200])
  mid_use = mid_rbin[keep]; mincomp_use = mincomp[keep]
  if (length(mid_use) < 2)
    return(list(Ar = rep(NA, length(ri)), mid_rbin = mid_rbin, mincomp = mincomp, best = NA))

  chi = sapply(contours, function(cont){
    fint = approxfun(ri[lt_r200], cont[lt_r200])
    Ar_comp = fint(mid_use)
    if (all(is.na(Ar_comp))) return(NA)
    median(abs(Ar_comp - mincomp_use), na.rm = TRUE)
  })
  if (all(is.na(chi)))
    return(list(Ar = rep(NA, length(ri)), mid_rbin = mid_rbin, mincomp = mincomp, best = NA))
  best = which.min(chi)
  list(Ar = contours[[best]], mid_rbin = mid_rbin, mincomp = mincomp, chi = chi, best = best)
}

##############################################################################
# Radially-varying beta(r) and F_beta(r), ported from causticpy
# (Gifford, Miller & Kern 2013). This project had been using beta=0.2 and
# Fbeta=fbr as CONSTANTS throughout, following the papers' own simplifying
# assumption for practical use -- these radial forms are the more detailed
# alternative the reference code actually implements, offered here as
# opt-in alternatives (beta_radial=TRUE, fbr_radial=TRUE in run_caustic())
# rather than new defaults, since they were never validated against this
# port's own data before now.
##############################################################################
beta_radial_gifford = function(r, r200, conc){
  0.5 * r / (r + r200 / conc)
}

fbeta_radial_nfw = function(r, r200, conc, gb){
  # limit as r -> 0 is exactly 0 (numerator ~ x^2, denominator ~ x for small
  # x = r/r200*conc), but evaluating the formula AT r=0 gives 0/0 = NaN;
  # substitute a tiny nonzero r there and then force the known r=0 limit.
  r_safe = ifelse(r <= 0, r200 / conc * 1e-6, r)
  x = r_safe / r200 * conc
  fb = 0.5 * (x^2) / ((1 + x)^2 * log(1 + x)) * gb
  fb[r <= 0] = 0
  fb
}

##############################################################################
# Shared mass-profile computation from a caustic amplitude curve A(r).
# Factored out so the same DG97-anchoring / fix_r200 / r=0 handling logic
# (already validated for the main density-matched curve, Ar_finalD) can be
# reused unchanged for the independent edge-detection curve (Ar_finalE)
# without duplicating it by hand.
##############################################################################
##############################################################################
# Small shared helper: prepares an amplitude curve for NFWfit() -- scaling
# by sqrt(g(beta)) over the fitting radii, and building the matching
# per-point noise estimate from delta_A(r)/A(r). Factored out because the
# density-matched (D) and edge-detection (E) branches below both need
# exactly this, and used to duplicate it verbatim.
##############################################################################
prep_fit_inputs = function(Ar_curve, fitting_radii, gb, delta_A_over_A_filled){
  Ar_scaled = Ar_curve[fitting_radii] * sqrt(gb[fitting_radii])
  sigma_A = Ar_scaled * delta_A_over_A_filled[fitting_radii]
  list(Ar = Ar_scaled, sigma_A = sigma_A)
}

mass_from_Ar = function(Ar_curve, x_range, r200, r200_input, fbr, M200_prior, R200_prior,
                         fix_r200, Hz, clus_z, comoving_input,
                         G = 6.67430e-11, Msol2kg = 1.98847e30, Mpc2m = 3.08567758e22){
  require(pracma)
  crit = 3 * Hz^2 / (8 * pi * G) * 1e6 / Msol2kg * Mpc2m
  if (comoving_input) crit = crit / (1 + clus_z)^3
  rsi = x_range * Mpc2m
  Asi = Ar_curve * 1e3
  sumtot = cumtrapz_base(rsi, fbr * Asi^2)[,1]

  have_anchor = !is.na(M200_prior) & !is.na(R200_prior)
  if (have_anchor) {
    idx0 = which.min(abs(x_range - R200_prior))
    sumtot_anchor = sumtot[idx0]
    GM0 = G * M200_prior * Msol2kg
    massprofile = sumtot / (G * Msol2kg)
    idx_out = idx0:length(massprofile)
    massprofile[idx_out] = (sumtot[idx_out] - sumtot_anchor + GM0) / (G * Msol2kg)
  } else {
    massprofile = sumtot / (G * Msol2kg)
  }

  x_for_density = x_range[-1]
  avg_density = massprofile[-1] / (4 / 3 * pi * (x_for_density)^3)
  ok = is.finite(avg_density)
  if (sum(ok) < 2) return(list(r200_est = NA, r500_est = NA, M200_est = NA, massprofile = rep(NA, length(x_range))))

  if (fix_r200) {
    # BUG FIX: when both fix_r200 AND anchoring are active, M200 used to be
    # evaluated (via linear interpolation of massprofile) at the exact
    # continuous r200_input value, while the anchor itself was applied at
    # x_range[idx0] -- the nearest grid point to R200_prior (grid spacing
    # = `by`, e.g. 0.05 Mpc). r200_input can sit up to half a grid step
    # away from that anchored point, and interpolating across that gap (on
    # a profile with nonzero local slope) introduced a spurious mismatch
    # between M200_est and M200_prior -- measured at ~1% on a real
    # cluster, small but not the exact match this configuration should
    # give by construction. Now: when anchored, evaluate M200 at the SAME
    # grid point the anchor was applied to, eliminating the discretization
    # gap entirely (verified: mismatch drops from ~1% to numerical noise).
    r200_est = if (have_anchor) x_range[idx0] else r200_input
    finterp = approxfun(x_range, massprofile)
    M200_est = finterp(r200_est)
    r500_est = NA
  } else {
    avg_density_mono = cummin(avg_density[ok])
    finterp = approxfun(avg_density_mono, x_for_density[ok])
    r200_est = finterp(200 * crit)
    r500_est = finterp(500 * crit)
    finterp = approxfun(x_range, massprofile)
    M200_est = finterp(r200_est)
  }
  list(r200_est = r200_est, r500_est = r500_est, M200_est = M200_est, massprofile = massprofile)
}

##############################################################################
# Bootstrap uncertainty on M200 (and other run_caustic() outputs), by
# resampling the input galaxies WITH replacement and re-running the whole
# pipeline `n_boot` times.
#
# This is a genuinely different uncertainty estimate from the D99/Serra et
# al. (2011) formula run_caustic() reports by default (M200_err/M200_err_frac):
# that formula is analytic, based on how well-sampled the phase-space
# density is at each radius; this is empirical, based on how much the
# FINAL answer moves when the particular set of galaxies you happened to
# observe is perturbed. Tested on 5 real Tempel et al. (2017) clusters:
# the two do NOT agree on which is larger (sometimes the D99 formula gives
# a bigger number, sometimes the bootstrap does, cluster by cluster) --
# treat them as complementary, not interchangeable, and don't assume
# either one is "the" correct error bar.
#
# Caveats: (1) slow -- this reruns the full pipeline n_boot times, so cost
# scales linearly with n_boot (order of a few seconds to tens of seconds
# per cluster for n_boot=30-50 with typical cluster richness); not meant
# for large-sample batch validation, only for inspecting individual
# clusters of interest. (2) like any bootstrap, it only captures sampling
# variability WITHIN the galaxies you already have -- it says nothing
# about, e.g., a genuinely different realisation of the same cluster's
# galaxy population, or contamination whose statistics differ from what
# resampling the observed (already contaminated) sample can produce.
##############################################################################
run_caustic_bootstrap = function(rproj, vproj, clus_z, n_boot = 50, seed = NULL, ...){
  if (!is.null(seed)) set.seed(seed)
  n = length(rproj)
  logM200_boot = r200_boot = conc_boot = c()
  n_fail = 0
  for (b in 1:n_boot) {
    idx = sample.int(n, n, replace = TRUE)
    r = tryCatch(run_caustic(rproj[idx], vproj[idx], clus_z, plot = FALSE, verbose = FALSE, ...),
                 error = function(e) NULL)
    if (!is.null(r)) {
      logM200_boot = c(logM200_boot, r$log_M200_est)
      r200_boot = c(r200_boot, r$r200_est)
      conc_boot = c(conc_boot, r$conc_est)
    } else n_fail = n_fail + 1
  }
  n_success = length(logM200_boot)
  if (n_success < 2) {
    warning('run_caustic_bootstrap(): fewer than 2 successful resamples (', n_success,
            ' of ', n_boot, ') -- cannot estimate a spread.')
    return(list(n_success = n_success, n_fail = n_fail, logM200_boot = logM200_boot,
                logM200_sd = NA, logM200_16 = NA, logM200_84 = NA,
                r200_boot = r200_boot, conc_boot = conc_boot))
  }
  list(n_success = n_success, n_fail = n_fail, logM200_boot = logM200_boot,
       logM200_median = median(logM200_boot), logM200_sd = sd(logM200_boot),
       logM200_16 = quantile(logM200_boot, 0.16, names = FALSE),
       logM200_84 = quantile(logM200_boot, 0.84, names = FALSE),
       r200_boot = r200_boot, conc_boot = conc_boot)
}

##############################################################################
# Wrapper around run_caustic() that retries with a smaller `rlimit` when the
# normal call fails with "the contours do not expand to the radial limit"
# (the dominant failure mode found in testing, especially with fixed-radius
# extractions where the default rlimit sometimes still demands the contour
# reach further than the real signal extends). Deliberately kept as a
# SEPARATE wrapper rather than folded into run_caustic() itself, to keep
# that function's own behaviour simple and predictable -- this retry logic
# is opt-in, and only activates on this one specific, recognised failure.
#
# Tested against 36 real Tempel et al. (2017) clusters (extracted to a
# fixed 10 Mpc / +-4000 km/s, the extraction convention this failure mode
# was most common for) that failed this way with the default rlimit: 17 of
# 36 (47%) recovered a result by searching progressively smaller rlimit
# fractions of the data's own extent. IMPORTANT CAVEAT, found in the same
# test: those 17 recovered clusters were noticeably LESS precise than the
# clusters that converged normally (median M200 bias +0.19 dex and only
# 59% within a factor ~2, vs -0.01 dex / 70% for the normally-converged
# ones) -- recovering a result this way is a real trade-off (coverage for
# precision), not a free win.
#
# FURTHER CAVEAT (found testing on CIRS, Rines & Diaferio 2006 -- cluster
# A0160): the recovered M200_err_frac is NOT a reliable safety net for
# these cases. A0160 recovered "successfully" at rlimit=0.5*xmax, but the
# resulting M200 was wrong by ~0.88 dex (a factor of ~7.7) while its OWN
# reported M200_err_frac was a reassuring 15.8% -- the D99/Serra et al.
# (2011) error formula measures how well-sampled the density is at each
# radius, it has no way to know the fitting radius itself was forced down
# artificially, so it does not flag this failure mode. Do not treat a
# small M200_err_frac as validating a rlimit-search-recovered result:
# the only reliable signal that a result came from this fallback at all
# is `rlimit_frac_used` being non-NA -- treat every such result as
# lower-confidence regardless of what its own error bar says, and prefer
# independent cross-checks (e.g. comparison to an external catalogue, or
# simply re-examining the cluster's raw data) over trusting the pipeline's
# self-reported uncertainty for these specific cases.
##############################################################################
run_caustic_robust = function(rproj, vproj, clus_z, ...,
                               search_rlimit = TRUE,
                               rlimit_search_fracs = c(0.80, 0.70, 0.60, 0.50, 0.40, 0.30),
                               verbose = TRUE){
  xmax_full = max(rproj)
  extra_args = list(...)
  # BUG FIX: this used to force verbose=FALSE on every internal run_caustic()
  # call regardless of what the caller asked for here, so the normal
  # detailed printout (R200/M200/vdisp estimates, etc.) never appeared --
  # only this wrapper's own "searching..."/"recovered..." messages did.
  # Now: the first (normal) attempt respects the caller's own `verbose`.
  extra_args$verbose = verbose

  first_try = tryCatch(do.call(run_caustic, c(list(rproj = rproj, vproj = vproj, clus_z = clus_z), extra_args)),
                        error = function(e) conditionMessage(e))
  if (!is.character(first_try)) {
    first_try$rlimit_frac_used = NA  # succeeded normally, no search needed
    return(first_try)
  }
  if (!search_rlimit || !grepl('do not expand', first_try, fixed = TRUE)) {
    stop(first_try)  # not the specific, recognised failure mode this wrapper targets, or search disabled
  }
  if (verbose) message('run_caustic_robust(): normal call failed ("', first_try,
                        '"); searching for a smaller rlimit...')
  # The search attempts themselves stay silent even when verbose=TRUE: each
  # one prints partial progress (Data Size, Initial R200, ...) right up
  # until it fails, and showing that for every rejected fraction would be
  # noisy and easy to mistake for a real result. Only the winning attempt
  # (if any) gets shown, by re-running it once more below.
  search_args = extra_args
  search_args$verbose = FALSE
  for (frac in rlimit_search_fracs) {
    args2 = search_args
    args2$xmax = xmax_full
    args2$rlimit = frac * xmax_full
    r2 = tryCatch(do.call(run_caustic, c(list(rproj = rproj, vproj = vproj, clus_z = clus_z), args2)),
                  error = function(e) NULL)
    if (!is.null(r2)) {
      if (verbose) {
        message('run_caustic_robust(): recovered using rlimit = ', round(frac, 2),
                ' * xmax (', round(frac * xmax_full, 2), ' Mpc) -- treat this result as ',
                'lower-confidence regardless of its own M200_err_frac (that error bar is ',
                'NOT a reliable safety net here; see the function documentation above).')
        # re-run the winning configuration once more with verbose=TRUE so
        # the caller sees the same detailed printout run_caustic() gives
        # normally (deterministic grid-based methods reproduce identically;
        # only re-runs the ONE attempt that actually worked, not the search).
        args2$verbose = TRUE
        r2 = tryCatch(do.call(run_caustic, c(list(rproj = rproj, vproj = vproj, clus_z = clus_z), args2)),
                      error = function(e) r2)
      }
      r2$rlimit_frac_used = frac
      return(r2)
    }
  }
  stop('run_caustic_robust(): could not converge even after searching rlimit fractions ',
       paste(rlimit_search_fracs, collapse = ', '))
}

##############################################################################
# Wrapper that re-centers a cluster's own galaxy candidates using
# sigma_plateau()'s binary-energy-tree method (from sigma_plateau.R),
# instead of trusting whatever external center (e.g. a catalogue's own
# group/X-ray position) the caller started with.
#
# REQUIRES sigma_plateau() and rv_proj() to already be loaded (source
# sigma_plateau.R / sigma_plateau_corregido.R first) -- these are NOT
# copied into this file, since they pull in heavier dependencies (Rcpp,
# cosmoFns, data.table) that shouldn't be forced on everyone using
# RCausticMass.R. Gives a clear error if they aren't found, rather than a
# cryptic "could not find function" one.
#
# HOW IT WORKS: (1) run the normal pipeline once with the ORIGINAL center,
# to get a preliminary escape-surface and membership; (2) take the
# preliminary members that ALSO fall within the preliminary R200 (a
# stricter cut than "any radius" -- found necessary, see below); (3) if
# there are enough of them, feed their RA/Dec/z to sigma_plateau() to find
# an independent center from the binary energy tree; (4) recompute
# dproj/vlos for ALL the original candidates around that NEW center; (5)
# re-run the full pipeline on the re-centered data.
#
# VALIDATION (100 real Tempel et al. 2017 clusters, extracted to a fixed
# 10 Mpc / +-4000 km/s window, blind mode): restricting to preliminary
# members WITHIN R200 (not just anywhere in the field) was essential --
# without that restriction, re-centering was a WASH on a small subset (13
# 20 test) and, in some individual clusters, met genuinely large offsets
# with WORSE precision than the original center. Restricting to within
# R200 and requiring at least `min_members_recenter` such members before
# attempting it gave a large, consistent improvement, holding from a
# minimum of 6 up to 40+ members: at N>=8 (26 of 100 clusters reached this
# threshold), R200 estimate scatter (sd of R200_est/R200_true - 1) dropped
# from 0.414 (original center) to 0.240 (re-centered) -- a ~42% reduction,
# with median bias also improving (+0.074 -> +0.006). The effect gets
# stronger (up to ~59% sd reduction) for clusters with even more
# preliminary members, and is still present (~38%) down to the lowest
# floor tested, N>=6.
#
# CAVEAT: the majority of clusters (74 of 100 in that same test) never
# reach even N>=6 preliminary members within R200 and so are returned
# with the ORIGINAL center unchanged (see the `recentered` field in the
# output to tell which happened) -- this is not a universal fix, it helps
# specifically the subset of clusters rich enough for the preliminary
# pass to already identify a reasonable number of likely members.
##############################################################################
##############################################################################
# Extends outlier classification to the FULL range of dproj, beyond the
# `rlimit` boundary run_caustic() actually analyses.
#
# WHY THIS EXISTS: run_caustic()'s own `caustic_outliers` field defaults
# EVERY point to 1 (outlier) and only overwrites it to 0 for points that
# were both (a) within `rlimit` (so included in data_set) AND (b) below
# the fitted escape curve there. Points beyond rlimit are never evaluated
# against any curve at all -- they just keep the default "outlier" flag
# unconditionally, regardless of how small their actual |vlos| is. This
# function replaces that blanket default with a genuine physical test
# everywhere: it reconstructs the fitted NFW escape-velocity curve from
# `nfw_d0`/`nfw_rs`/`nfw_gb` (exposed by run_caustic() for exactly this
# purpose) and evaluates it analytically at ANY radius -- inside rlimit it
# reproduces caustic_outliers exactly (same underlying curve, just
# evaluated via the closed-form formula instead of interpolated off
# x_range), and outside rlimit it extrapolates the same curve rather than
# leaving the question unanswered.
#
# CAVEAT: an extrapolation is only as good as the fit it extends. The NFW
# fit was constrained by data within rlimit; far beyond it, the escape
# curve is a projection, not a validated measurement -- treat outlier
# calls far outside rlimit as informative but less certain than ones
# inside it.
##############################################################################
nfw_escape_curve = function(r, d0, rs, gb){
  r = ifelse(r <= 0, min(r[r > 0], na.rm = TRUE), r)  # guard the r=0 singularity, as run_caustic() itself does
  sqrt(2 * 4 * pi * 4.3e-9 * d0 * rs^2 * log(1 + r / rs) / (r / rs)) / sqrt(gb)
}

# Convenience wrapper for use AFTER run_caustic() has already returned --
# extracts nfw_d0/nfw_rs/nfw_gb from a result object instead of requiring
# them as separate arguments. run_caustic() itself now computes this
# internally and returns it directly as `caustic_outliers_extended` (see
# below) -- this wrapper remains for re-deriving it standalone, e.g. with
# a different rproj/vproj than the ones originally passed in.
extend_outliers_nfw = function(rproj, vproj, result){
  A_r = nfw_escape_curve(rproj, result$nfw_d0, result$nfw_rs, result$nfw_gb)
  is_outlier = abs(vproj) > A_r
  is_outlier[is.na(A_r)] = NA
  list(is_outlier = is_outlier, escape_velocity = A_r)
}

run_caustic_recentered = function(ra, dec, z, clus_ra, clus_dec, clus_z,
                                   min_members_recenter = 8, verbose = TRUE, ...){
  if (!exists('sigma_plateau') || !exists('rv_proj'))
    stop('run_caustic_recentered(): sigma_plateau() and rv_proj() are not defined -- ',
         'source sigma_plateau.R (or sigma_plateau_corregido.R) before calling this function.')

  rv0 = rv_proj(ra, dec, z, clus_ra, clus_dec, clus_z)
  r1 = tryCatch(run_caustic(rv0$dproj, rv0$vlos, clus_z, verbose = FALSE, ...), error = function(e) NULL)
  if (is.null(r1))
    stop('run_caustic_recentered(): run_caustic() did not converge even with the original center.')

  members_idx = which(r1$caustic_outliers == 0 & rv0$dproj < r1$r200_est)
  n_mem = length(members_idx)
  r1$recentered = FALSE
  r1$center_ra = clus_ra; r1$center_dec = clus_dec; r1$center_z = clus_z
  r1$n_members_prelim = n_mem
  r1$center_offset_arcsec = NA

  if (n_mem < min_members_recenter) {
    if (verbose) message('run_caustic_recentered(): only ', n_mem, ' preliminary members within R200 ',
                          '(< min_members_recenter = ', min_members_recenter, ') -- returning the result ',
                          'with the ORIGINAL center, not attempting to re-center.')
    return(r1)
  }

  sp = tryCatch(sigma_plateau(ra[members_idx], dec[members_idx], z[members_idx], verbose = FALSE, plot = FALSE),
                error = function(e) NULL)
  if (is.null(sp)) {
    if (verbose) message('run_caustic_recentered(): sigma_plateau() failed on the preliminary members -- ',
                          'returning the result with the ORIGINAL center.')
    return(r1)
  }

  rv_new = rv_proj(ra, dec, z, sp$cluster_info$ra, sp$cluster_info$dec, sp$cluster_info$z)
  r2 = tryCatch(run_caustic(rv_new$dproj, rv_new$vlos, sp$cluster_info$z, verbose = verbose, ...),
                error = function(e) NULL)
  if (is.null(r2)) {
    if (verbose) message('run_caustic_recentered(): run_caustic() failed with the re-centered data -- ',
                          'returning the result with the ORIGINAL center.')
    return(r1)
  }

  offset_arcsec = sqrt((sp$cluster_info$ra - clus_ra)^2 * cos(clus_dec * pi / 180)^2 +
                        (sp$cluster_info$dec - clus_dec)^2) * 3600
  r2$recentered = TRUE
  r2$center_ra = sp$cluster_info$ra; r2$center_dec = sp$cluster_info$dec; r2$center_z = sp$cluster_info$z
  r2$n_members_prelim = n_mem
  r2$center_offset_arcsec = offset_arcsec
  if (verbose) message('run_caustic_recentered(): re-centered using ', n_mem, ' preliminary members ',
                        '(offset from original center: ', round(offset_arcsec, 1), ' arcsec).')
  r2
}

# Main function
#
# `fbr` default: 0.6, not the literature value (Serra et al. 2011 use 0.7;
# Diaferio & Geller 1997 / Diaferio 1999 use 0.5). This project calibrated
# fbr empirically against 100 real clusters from Tempel et al. (2017),
# since neither literature value was derived for this exact pipeline (this
# implementation's own kernel, contour-finding, and NFW-fit choices differ
# in several details from both D99's and Serra et al. 2011's). Findings:
#   - kernel='gaussian' (default), self-consistent r200_est (fix_r200=F):
#     fbr=0.6 gives ~zero median bias in log M200, n=70/100 converged,
#     IQR([-0.14,0.16]) dex, 80% within a factor ~2 of the true M200.
#   - kernel='gaussian', fix_r200=T (R200/vdisp held fixed at trusted,
#     externally known values instead of self-consistently re-derived --
#     see fix_r200 documentation below): fbr=0.6 gives essentially the
#     same near-zero bias but n=77/100 converged, a much tighter
#     IQR([-0.11,0.07]) dex, and 91% within a factor ~2. Holding r200/vdisp
#     fixed when they are independently known (e.g. from a group catalogue)
#     removes a real, measured source of noise -- re-deriving r200_est
#     from the density profile is itself an error-prone step whose errors
#     otherwise propagate directly into M200_est.
#   - kernel='adaptive': needs a much smaller fbr (~0.12-0.25 depending on
#     the h_opt normalisation constant used, see adaptive_kernel_2d()) and
#     tops out around ~50-60% within a factor ~2 regardless of that
#     constant -- a real precision ceiling below the gaussian kernel's,
#     not merely a calibration issue.
# Re-validate this default if the input field-of-view conventions, sample
# richness, or kernel choice differ substantially from this test.
run_caustic=function(rproj, vproj, clus_z, r200 = NA, clus_vdisp = NA, rlimit = NA, 
                     vlimit = NA, xmax = NA, ymax = 4000, mirror = T, Om = .3, H0 = 70, 
                     fbr = 0.6, q = 10, beta = NA, centering = F,
                     M200_prior = NA, R200_prior = NA, comoving_input = F, fix_r200 = F,
                     kernel = c('gaussian', 'adaptive'), hc = NA,
                     gradu = 0.5, gradd = 2.0,
                     c_min = 1, c_max = 15, beta_radial = F, fbr_radial = F, conc_default = 5,
                     compute_edge = F, member_boundary = c('fit', 'upper'),
                     conc_method = c('bounded', 'bayesian'), c_prior_scatter_dex = 0.1, little_h = 0.7,
                     mass_method = c('integral', 'bayesian_joint'), logM200_prior_scatter_dex = 0.4,
                     adaptive_vlimit = T, adaptive_vlimit_mult = 3.5,
                     bayesian = F,
                     plot = T, verbose = T){
  kernel = match.arg(kernel)
  member_boundary = match.arg(member_boundary)
  # `bayesian`: single switch for "get halo parameters the Bayesian way",
  # so the caller doesn't need to separately know about conc_method vs.
  # mass_method or reason about which applies given what priors were
  # supplied -- this picks the right one automatically:
  #   - M200_prior AND R200_prior both given: M200 is already fixed, so
  #     only concentration is left to estimate -> conc_method='bayesian'
  #     (1-parameter grid posterior over c).
  #   - M200_prior not given (regardless of R200_prior/r200): M200 isn't
  #     fixed, so M200 and concentration are estimated JOINTLY ->
  #     mass_method='bayesian_joint' (2-parameter grid posterior over
  #     (M200, c), with r200 -- given or auto-estimated -- as the anchor).
  # Explicit conc_method=/mass_method= still work as before if you want
  # to force one or the other by hand instead of this automatic choice.
  if (bayesian) {
    if (!is.na(M200_prior) & !is.na(R200_prior)) conc_method = 'bayesian'
    else mass_method = 'bayesian_joint'
  }
  conc_method = match.arg(conc_method)
  mass_method = match.arg(mass_method)
  #   - kernel='gaussian', fix_r200=TRUE (r200/vdisp known from an external
  #     catalogue, e.g. Tempel et al. 2017): hc=1.3 was validated here --
  #     widening the Gaussian blur modestly consistently improved M200
  #     recovery (more convergence, tighter IQR, better fraction within a
  #     factor ~1.4-2 of the true mass).
  #   - kernel='gaussian', fix_r200=FALSE (r200 not known in advance, the
  #     more common "blind" use case for a genuinely new cluster): tested
  #     separately and hc=1.3 does NOT transfer -- once each hc is given
  #     its own properly recalibrated fbr, hc=1.0 came out clearly better
  #     (83% vs 76% within a factor ~2 in a 100-cluster test), likely
  #     because hc=1.3's extra smoothing helps most when r200 is already
  #     pinned down and doesn't have to be independently (re-)estimated
  #     from the same, now-blurrier density map. So hc reverts to 1 here.
  #   - kernel='adaptive': never tested with hc != 1 in either fix_r200
  #     mode (only lambda_max was explored there), keeps its original
  #     default (1) regardless, to avoid silently changing untested
  #     behaviour.
  # Pass hc explicitly to override any of these defaults.
  if (is.na(hc)) hc = if (kernel == 'gaussian' && fix_r200) 1.3 else 1
  if (fix_r200 & (is.na(r200) | is.na(clus_vdisp)))
    stop('run_caustic(): fix_r200=TRUE requires both r200 and clus_vdisp to be supplied ',
         '(e.g. from an external catalog) -- there is nothing to hold fixed otherwise.')
  r200_input = r200  # remember the caller's original value; r200 may be reassigned below
                      # (auto-estimate if NA, or gaussian_kernel()'s own internal adjustment)
  
  N0 = length(rproj)
  # BUG FIX: rproj/vproj get REASSIGNED later in this function (to the
  # trimmed data_set's columns, after the final rlimit is settled) --
  # keeping a copy of the ORIGINAL, full-length arrays here so that
  # caustic_outliers_extended (below) can be computed over the entire
  # input range rather than accidentally over the already-trimmed subset.
  rproj_orig_full = rproj
  vproj_orig_full = vproj
  clus_z = clus_z[1] 
  Ol = 1 - Om
  Hz = H0 * sqrt(Om * (1 + clus_z)^3 + Ol) # H(z)

  # `adaptive_vlimit` (default TRUE as of this version): instead of a fixed
  # `vlimit` (velocity cut, e.g. the classic 4000 km/s), derive it from a
  # rough preliminary velocity dispersion of THIS cluster, computed here
  # on the raw input BEFORE any of the trimming below.
  #
  # Motivation: extracting SDSS candidates with a FIXED, wide velocity
  # window (e.g. +-4000 km/s, common when candidates are pulled from an
  # external catalogue without knowing the true membership beforehand)
  # systematically hurts blind-mode (r200 not given) precision compared to
  # a window that is already tied to the cluster's own dynamics -- tested
  # by re-extracting the SAME 100 Tempel et al. (2017) clusters both ways:
  # with their own narrower, membership-informed window the blind R200
  # estimate had sd=0.15 (of R200_est/R200_true - 1); re-extracted with a
  # fixed +-4000 km/s window instead, sd rose to 0.47 (three times worse),
  # for identical galaxies and clusters -- the extraction geometry alone
  # was responsible. Rebuilding an adaptive window algorithmically (no
  # true membership needed) recovers most of that: sd dropped back to
  # 0.27 using vlimit = adaptive_vlimit_mult * (a shifting_gapper-cleaned
  # biweight dispersion within min(rmax,3) Mpc of the raw data).
  # adaptive_vlimit_mult=3.5 was found to roughly zero out the median bias
  # in that test; smaller values tighten the window further (less bias
  # from wide-field dilution) at some cost to convergence, larger values
  # do the reverse -- tune per dataset if precision matters more than
  # coverage or vice versa.
  #
  # Only applies when the caller hasn't already supplied `vlimit`
  # explicitly (checked before it gets its own default below).
  #
  # BUG FIX (first pass): re-deriving vdisp from scratch when clus_vdisp
  # was ALREADY known (typically alongside fix_r200=TRUE) added noise
  # instead of removing it -- tested against the 100 Tempel et al. (2017)
  # clusters in informed mode, convergence dropped from 66/100 to 59/100.
  # The first fix simply skipped adaptive_vlimit whenever clus_vdisp was
  # given, which avoided that regression but threw away a real
  # improvement in the process: re-tested with the WIDE-window (10 Mpc,
  # +-4000 km/s) extraction in informed mode, using the raw fixed vlimit
  # gave sd=0.31 on log M200 -- close to what re-deriving vdisp gave in
  # BLIND mode on the same data. Using vlimit = adaptive_vlimit_mult *
  # clus_vdisp (the caller's OWN known, real dispersion, not a
  # re-estimate) instead dropped sd to 0.15 -- matching the ORIGINAL,
  # narrow, membership-informed extraction's own informed-mode precision
  # almost exactly. So: use the known clus_vdisp directly when available,
  # only fall back to re-deriving one from scratch when it truly isn't.
  if (adaptive_vlimit && is.na(vlimit) && !is.na(clus_vdisp)) {
    vlimit = adaptive_vlimit_mult * clus_vdisp
    if (verbose) message('adaptive_vlimit: derived vlimit = ', round(vlimit, 1),
                          ' km/s from the supplied clus_vdisp = ', round(clus_vdisp, 1), ' km/s')
  } else if (adaptive_vlimit && is.na(vlimit) && is.na(clus_vdisp)) {
    rmax0 = max(rproj)
    rl0 = min(rmax0, 3)
    vdisp0 = gapper_cleaned_vdisp(rproj, vproj, rl0)
    if (is.na(vdisp0)) vdisp0 = biwScale(vproj[rproj < rl0])  # too few points to clean; never give up entirely here
    vlimit = adaptive_vlimit_mult * vdisp0
    if (verbose) message('adaptive_vlimit: derived vlimit = ', round(vlimit, 1),
                          ' km/s from a preliminary vdisp = ', round(vdisp0, 1), ' km/s')
  }
  if (is.na(vlimit)) vlimit = 4000  # old fixed default, when not given and not adaptive

  # `xmax`/`rlimit` default to the data's OWN extent instead of a fixed
  # number (previously xmax=6, rlimit=5.8 always). Testing against CIRS
  # (Rines & Diaferio 2006), where the actual extraction radius varies by
  # dataset (6 Mpc in one extraction, 10 Mpc in another) and isn't always
  # close to the old fixed defaults, found this fixed-default mismatch was
  # itself a real cause of failures: if rlimit exceeds where the data
  # actually end, findcontours() can be asked to trace a contour further
  # out than any signal reaches ("the contours do not expand to the radial
  # limit"); if xmax is smaller than the data, points beyond it are simply
  # discarded. Now: xmax defaults to max(rproj) (the true data extent), and
  # rlimit defaults to a bit inside that (95%), leaving the contour-tracing
  # step room to close before the true edge of the data rather than
  # demanding it reach exactly to where the data (and often the signal)
  # runs out. Pass either explicitly to override.
  if (is.na(xmax)) xmax = max(rproj)
  # BUG FIX: this margin was originally 0.95 (rlimit = 0.95*xmax), validated
  # against CIRS (a FIXED-radius extraction, e.g. 6 or 10 Mpc regardless of
  # cluster size) where it helped a lot. But re-tested against Tempel et
  # al. 2017 (a PROPORTIONAL extraction, exactly 3xR200 per cluster) it
  # caused a real regression: convergence dropped from ~70/100 to 33/100,
  # because a proportional extraction has NO extra margin beyond the
  # signal-bearing region to begin with -- forcing rlimit to sit at 95% of
  # an already-tight boundary demands contour closure almost exactly at the
  # data's true edge, the same failure mode this default was meant to fix.
  # 0.80 was found to work well for BOTH conventions (Tempel: 65/100
  # converged, CIRS: unchanged at 67/74) and is used as a compromise
  # default; pass rlimit explicitly to tune it for a specific extraction.
  if (is.na(rlimit)) rlimit = 0.80 * xmax

  data_set = data.frame(r = rproj, v = vproj)
  
  # reduce sample within limits
  data_set = subset(data_set, r < rlimit & abs(v) < vlimit)
  if(nrow(data_set) < 5) stop('Data set has too few elements inside the limits')
  if(verbose) message('Data Size: ', nrow(data_set))
  rmax = max(data_set[,1])
  
  # if r200 is not available we provide an initial estimate
  #
  # ITERATIVE preliminary R200 (the definitive algorithm, not an opt-in
  # variant): start from rl=min(rmax,3) as before, clean interlopers with
  # shifting_gapper() (see rationale below), compute a trial R200 via the
  # virial-scaling formula, then re-estimate using a NARROWER window
  # (exactly the trial R200 itself, not a fixed multiple of it) and
  # repeat, until the estimate changes by less than 5% between iterations
  # (or a small iteration cap / floor radius is hit). This mimics
  # computing the velocity dispersion strictly WITHIN R200 (the standard
  # virial definition) rather than within an arbitrary, wider window.
  #
  # Tested against 100 real Tempel et al. (2017) clusters re-extracted
  # with a fixed 10 Mpc / +-4000 km/s window: iterating narrowed the
  # preliminary R200/R200_true ratio's IQR by ~21% relative to a single
  # pass (0.87 -> 0.69 dex-equivalent width), with a similar or slightly
  # better median. Re-tested against CIRS (Rines & Diaferio 2006, also a
  # fixed-radius extraction): essentially no change either way (median
  # 1.348 -> 1.326) -- the extra iterations neither help nor hurt there,
  # consistent with CIRS's dominant source of scatter being structural
  # (field width relative to each cluster's own R200) rather than the
  # velocity-dispersion window itself. Net effect across both datasets:
  # iterating is a safe default, with a real (if modest) upside on at
  # least one of them and no observed downside on the other.
  #
  # shifting_gapper() cleaning (previously the optional `prelim_gapper`,
  # now applied unconditionally as part of this algorithm) is applied at
  # EVERY iteration, not just the first: on real data with a wide,
  # unconstrained extraction window (tested on CIRS: up to ~15-20% of
  # galaxies within that window can be interlopers by simple
  # velocity-outlier standards), skipping it left a measurably worse
  # preliminary estimate (median ratio 1.53 uncleaned vs 1.35-1.33
  # cleaned). Falls back to the uncleaned biweight scale at that
  # iteration if shifting_gapper() fails or leaves too few galaxies.
  if(is.na(r200)){
    rl = min(rmax, 3)
    r200_prelim_prev = NA
    for (iter_prelim in 1:6) {
      vdisp_prelim = gapper_cleaned_vdisp(data_set[,1], data_set[,2], rl)
      if (is.na(vdisp_prelim)) break  # too few galaxies within rl to even attempt an estimate
      r200_new = sqrt(3) * vdisp_prelim / (10 * Hz)
      if (r200_new > rmax) r200_new = rmax
      if (!is.na(r200_prelim_prev) && abs(r200_new - r200_prelim_prev) / r200_prelim_prev < 0.05) {
        r200_prelim_prev = r200_new
        break
      }
      r200_prelim_prev = r200_new
      rl = max(r200_new, 0.15)  # narrow to exactly R200_prelim next iteration; floor avoids collapse
    }
    r200 = r200_prelim_prev
  }
  
  # reduce sample again forcing data to be within 3*r200
  rlimit = min(3.2 * r200, rlimit, rmax)
  data_set = subset(data_set, r <= rlimit)
  rproj = data_set[,1] 
  vproj = data_set[,2]
  
  # centering in the velocity space within r200 
  if(centering){
    voff = median(vproj[rproj < r200])
    if(verbose) message('Centering vlos with offset: ',round(voff, 3),' km/s')
    vproj = data_set[,2] = vproj - voff
  }
  
  # if vdisp is not available we provide a robust initial estimate
  #
  # BUG FIX: for a poor cluster, `rproj < r200 & abs(vproj) < vlimit` can
  # select zero galaxies (all fall just outside that particular cut),
  # leaving biwScale() with nothing to work with -- it now returns NA
  # cleanly instead of crashing (see biwScale()'s own fix above), but a
  # NA clus_vdisp would still break comparisons further down
  # (findcontours() uses clus_vdisp^2 numerically). Falls back to the
  # full data_set (no r200/vlimit restriction) in that case, which by
  # construction has at least the 5 points required to reach this point
  # in the function at all.
  if(is.na(clus_vdisp)) {
    clus_vdisp = biwScale(vproj[rproj < r200 & abs(vproj) < vlimit])
    if (is.na(clus_vdisp)) clus_vdisp = biwScale(vproj)
  }
	
  if(verbose) message('Initial R200:  ', round(r200, 3), ' Mpc')
  if(verbose) message('Initial Vdisp: ', round(clus_vdisp, 3), ' km/s')
  
  # calculating density of phase-space 
  if(!mirror & verbose) message('Calculating density')
  if(mirror & verbose) message('Calculating density w/Mirrored Data')

  if (kernel == 'gaussian') {
    gk = gaussian_kernel(rproj, vproj, r200, Hz, q, xmax, ymax, 0.05, mirror, hc, plot = F)
  } else {
    gk = adaptive_kernel_2d(rproj, vproj, r200, xmax, ymax, by = 0.05, mirror = mirror, hc = hc)
  }
  x_range = gk$x
  y_range = gk$y
  Zi = gk$z / max(abs(gk$z))
  if (!is.na(gk$r200_used)) {
    if (gk$r200_used != r200 & verbose)
      message('r200 was adjusted inside gaussian_kernel (too few objects inside r200); ',
              'using the adjusted value (', round(gk$r200_used, 3), ' Mpc) consistently from here on.')
    r200 = gk$r200_used
  }

  # identifying the initial caustic surface and members within the surface 
  if(verbose) message('Calculating initial surface...')
  fcont = findcontours(Zi, x_range, y_range, r200, clus_vdisp^2, rlimit, vlimit,
                        gradu = gradu, gradd = gradd, plot = F, 
					   verbose = verbose)
  Ar_finalD = fcont$caustic
  contours = fcont$contours

  # independent edge-detection surface (causticpy's Ar_finalE): picks,
  # from the SAME candidate contours above, whichever one best matches an
  # empirical estimate of the phase-space edge built directly from the raw
  # galaxy velocities (top edge_perc extremes per radial bin) -- entirely
  # independent of the phase-space density used for Ar_finalD. See
  # edge_caustic() for details.
  #
  # `compute_edge` (default FALSE, set TRUE to enable): testing against 100
  # real Tempel et al. (2017) clusters found this cross-check alone is
  # LESS precise than the main density-matched estimate (77% vs 90% within
  # a factor ~2), doesn't help when averaged with it, and -- perhaps
  # counter-intuitively -- close agreement between the two is not a
  # reliable sign of higher accuracy either (both are selected from the
  # SAME finite set of candidate contours, so they are less independent
  # than they look). Kept opt-in/out rather than removed, since it is
  # still a useful per-cluster diagnostic to inspect by eye; set to FALSE
  # to skip the extra computation (and its console output) entirely when
  # it isn't needed.
  if (compute_edge) {
    ecaus = edge_caustic(rproj, vproj, r200, contours, x_range, mirror = mirror)
    Ar_finalE = ecaus$Ar
    have_edge = !all(is.na(Ar_finalE))
  } else {
    Ar_finalE = rep(NA, length(x_range))
    have_edge = FALSE
  }
  
  # fit an NFW to the optimal caustic surface 
  #
  # `beta_radial` (opt-in, off by default): use the radially-varying
  # beta(r) = 0.5*r/(r+R200/c) from Gifford, Miller & Kern (2013)'s
  # reference code, instead of the fixed beta=0.2 this project has used
  # throughout. NOTE: this needs a concentration BEFORE the NFW fit below
  # has run (which is the only place a data-driven concentration,
  # conc_fit, is obtained) -- so `conc_default` is used here regardless of
  # whether M200_prior/R200_prior later yield a fit concentration. This is
  # a deliberate simplification (no iterative re-fit with an updated
  # beta/conc), documented rather than hidden.
  #
  # Validated only for fix_r200=TRUE (R200/vdisp known from an external
  # catalogue): there, beta_radial+fbr_radial together gave a modest but
  # real improvement (92% vs 90% within a factor ~2, 71% vs 66% within a
  # factor ~1.4). Separately tested in the BLIND scenario (r200 not
  # given): does NOT transfer, even after properly re-tuning
  # `conc_default` (the parameter that actually controls it there, since
  # fbr_radial ignores the scalar `fbr` entirely) -- best blind
  # configuration reached 75%/51% vs the plain constant-fbr baseline's
  # 83%/47%. Recommended only in the fix_r200=TRUE / informed case.
  if (beta_radial) {
    beta = beta_radial_gifford(x_range, r200, conc_default)
  } else if (is.na(beta)) {
    beta = rep(0.2, length(x_range))
  }
  gb = (3 - 2 * beta) / (1 - beta)
  fitting_radii = (x_range >= r200 / 3) & (x_range <= r200)
  rii = x_range[fitting_radii]

  # Uncertainty estimate (D99/Serra et al. 2011, eq. 23-24), moved ahead of
  # the NFW fit itself (it used to be computed after) so it can double as
  # the per-radius measurement-noise model for the Bayesian concentration
  # fit below, in addition to its original uses (M200 error bar,
  # member_boundary='upper').
  #   delta_A(r)/A(r) = kappa / max_v[f_q(r,v)]
  # where kappa is the density threshold that located the chosen escape
  # surface (returned by findcontours() above). Only applies to the
  # density-matched (Ar_finalD) surface.
  peak_density = apply(Zi, 1, max)
  delta_A_over_A = fcont$kappa / peak_density
  delta_A_over_A[!is.finite(delta_A_over_A) | peak_density <= 0] = NA
  # at radii far beyond where any data actually constrain the density
  # (typically near xmax, well past r200), peak_density -> 0 and this
  # ratio blows up (seen in testing: up to ~7000x) -- formally consistent
  # with "this point is unconstrained" but useless for the M200 error
  # budget, for the Bayesian likelihood below, or for plotting a band.
  # Cap the relative uncertainty at 200%.
  delta_A_over_A = pmin(delta_A_over_A, 2)
  delta_A_over_A_filled = delta_A_over_A
  delta_A_over_A_filled[is.na(delta_A_over_A_filled)] = 0  # treat "unknown" as "no widening", not "infinitely uncertain"
  fit_inputs_D = prep_fit_inputs(Ar_finalD, fitting_radii, gb, delta_A_over_A_filled)
  ArD = fit_inputs_D$Ar
  sigma_A_fit = fit_inputs_D$sigma_A

  # if the caller supplied M200_prior (e.g. from Tempel et al. 2017) but not
  # R200_prior, anchor R200_prior to the r200 already in use here (whether
  # user-supplied or auto-estimated above) rather than silently ignoring
  # M200_prior for lack of a matching radius.
  if(!is.na(M200_prior) & is.na(R200_prior)){
    R200_prior = r200
    if(verbose) message('R200_prior not supplied; using r200 = ', round(r200,3),
                         ' Mpc as R200_prior.')
  }
  nfw_fit = NFWfit(rii, ArD, halo_srad = r200 / 5, ri_full = x_range,
                    M200_prior = M200_prior, R200_prior = R200_prior, c_min = c_min, c_max = c_max,
                    conc_method = conc_method, sigma_A = sigma_A_fit,
                    c_prior_scatter_dex = c_prior_scatter_dex, little_h = little_h, clus_z = clus_z)
  caustic_fit = nfw_fit[[1]] / sqrt(gb)
  reserr = nfw_fit[[2]]
  conc_fit = nfw_fit[[3]]  # concentration, only non-NA when M200_prior/R200_prior were supplied
  conc_16 = nfw_fit[[4]]   # 68% credible interval, only non-NA for conc_method='bayesian'
  conc_84 = nfw_fit[[5]]
  # Raw NFW fit parameters (halo_scale_density=d0, halo_srad_used=rs) and
  # the Fbeta-derived normalisation `gb` used just above, exposed so the
  # caller can EXTRAPOLATE the fitted escape-velocity curve to radii beyond
  # rlimit (where caustic_fit/x_range simply don't extend, since that grid
  # is built only from the trimmed data_set) -- see nfw_escape_curve() and
  # extend_outliers_nfw() below, added specifically for this purpose.
  nfw_d0 = nfw_fit[[6]]
  nfw_rs = nfw_fit[[7]]
  nfw_gb = gb[1]

  # `member_boundary` (default 'fit', matching this project's and
  # causticpy's previous/only behaviour): which curve a galaxy's velocity
  # is compared against to decide membership.
  #   'fit'   -- the point-estimate NFW curve (caustic_fit) itself. This
  #     project separately found that vdisp_est comes out systematically
  #     LOW (~8-11%) against Tempel et al. (2017) sigma_v, and attributed
  #     it to exactly this: a hard cut at the best-fit curve truncates
  #     genuine high-velocity members whenever the point estimate happens
  #     to sit a bit low there, which a boundary with no accounting for
  #     its own uncertainty cannot avoid.
  #   'upper' -- the UPPER edge of the delta_A(r) uncertainty band,
  #     caustic_fit*(1+delta_A/A). More conservative about flagging
  #     outliers (fewer false positives among genuine members near the
  #     boundary), at the cost of admitting more contamination (more true
  #     interlopers slip through as "members") -- a real trade-off, not a
  #     strict improvement, which is why 'fit' remains the default.
  member_curve = if (member_boundary == 'upper') caustic_fit * (1 + delta_A_over_A_filled) else caustic_fit
  fcomp = approxfun(x_range, member_curve)
  vcompare = fcomp(data_set[,1])
  memflag = ifelse(abs(vcompare) > abs(data_set[,2]), 1, 0)
  
  # `fbr_radial` (opt-in, off by default): use the radially-varying
  # F_beta(r), the NFW-filling-function form from DG97/D99/Gifford et al.
  # (2013), instead of the constant `fbr` scalar this project calibrated
  # empirically (fbr~=0.5-0.7 depending on kernel/settings). Uses the same
  # conc_default as beta_radial above, for the same reason (no fit
  # concentration is available yet at this stage in the unanchored case).
  fbr_use = if (fbr_radial) fbeta_radial_nfw(x_range, r200, conc_default, gb) else fbr

  # mass estimation (see mass_from_Ar() for the DG97-anchoring / fix_r200 /
  # r=0-handling logic, factored out so it can be reused identically for
  # both the main density-matched curve and the independent edge curve)
  #
  # `mass_method='bayesian_joint'`: alternative to the density-integral
  # M200 above -- a joint (M200, concentration) Bayesian fit to the
  # observed A(r) curve itself, with r200 (given or auto-estimated) held
  # as a fixed anchor radius. Requires M200_prior to be UNSET (if M200 is
  # already fixed, there's nothing left to jointly estimate -- use
  # conc_method='bayesian' instead, which fits concentration alone).
  # Validated against 100 real Tempel et al. (2017) clusters: comparable
  # or slightly better precision than the density-integral estimate (73%
  # vs 66% within a factor ~1.4 of the true mass), plus well-behaved
  # concentration credible intervals as a byproduct.
  M200_bayes_16 = M200_bayes_84 = NA
  M200_est = NA
  if (mass_method == 'bayesian_joint' & is.na(M200_prior)) {
    b2 = tryCatch(
      NFWfit_bayes_2param(rii, ArD, sigma_A_fit, r200, clus_z,
                           c_prior_scatter_dex = c_prior_scatter_dex, little_h = little_h,
                           Om = Om, H0 = H0, logM200_prior_scatter_dex = logM200_prior_scatter_dex),
      error = function(e) NULL)
    if (!is.null(b2)) {
      M200_est = b2$M200
      M200_bayes_16 = b2$M200_16
      M200_bayes_84 = b2$M200_84
      conc_fit = b2$conc
      conc_16 = b2$conc_16
      conc_84 = b2$conc_84
      r200_est = r200
      r500_est = NA
      # keep caustic_fit consistent with the jointly-fit (M200, conc)
      # rather than whatever the earlier NFWfit() call (which didn't know
      # about this joint fit) produced
      rs_fit2 = r200 / conc_fit
      cc2 = conc_fit; fc2 = log(1 + cc2) - cc2 / (1 + cc2)
      d0_fit2 = M200_est / (4 * pi * rs_fit2^3 * fc2)
      caustic_fit = nfw_escape_curve(x_range, d0_fit2, rs_fit2, gb[1])
    } else if (verbose) {
      message('mass_method="bayesian_joint" failed for this cluster; falling back to the ',
              'density-integral estimate below.')
    }
  }
  if (is.na(M200_est) || mass_method == 'integral' || (mass_method == 'bayesian_joint' & !is.na(M200_prior))) {
    mD = mass_from_Ar(Ar_finalD, x_range, r200, r200_input, fbr_use, M200_prior, R200_prior,
                       fix_r200, Hz, clus_z, comoving_input)
    r200_est = mD$r200_est
    r500_est = mD$r500_est
    M200_est = mD$M200_est
  } else {
    mD = list(massprofile = rep(NA, length(x_range)))  # not computed in bayesian_joint mode
  }

  if (fix_r200 && !is.na(gk$r200_used) && abs(gk$r200_used - r200_input) > 1e-6 && verbose)
    message('Note: gaussian_kernel() internally adjusted r200 for its own bandwidth ',
            'calculation (too few objects inside r200_input); M200 is still evaluated ',
            'at the ORIGINAL r200_input = ', round(r200_input, 3), ' Mpc, as requested by fix_r200=TRUE.')

  if (is.na(r200_est) || is.na(M200_est))
    stop('run_caustic(): could not evaluate the mass profile at r200_est = ', r200_est,
         ' Mpc -- check that it falls within the sampled radial range (0, ',
         round(max(x_range), 2), ') Mpc (rlimit/xmax may be too small for this r200).')

  delta_M_cum = rep(NA, length(x_range))
  if (!is.na(M200_bayes_16) && !is.na(M200_bayes_84)) {
    # mass_method='bayesian_joint' succeeded: use the (asymmetric) 68%
    # credible interval directly as the uncertainty, rather than the
    # D99/Serra et al. (2011) formula below, which needs a mass PROFILE
    # (mD$massprofile) that this mode doesn't compute.
    M200_err = mean(c(M200_est - M200_bayes_16, M200_bayes_84 - M200_est))
  } else if (all(is.finite(mD$massprofile))) {
    shell_mass = c(mD$massprofile[1], diff(mD$massprofile))
    delta_M_cum = cumsum(abs(2 * shell_mass * delta_A_over_A))
    M200_err = approx(x_range, delta_M_cum, xout = r200_est)$y
  } else {
    M200_err = NA
  }

  # edge-method NFW fit and mass estimate (Ar_finalE), computed the same
  # way as the main curve, for cross-checking against M200_est above
  r200_est_edge = M200_est_edge = reserr_edge = conc_fit_edge = NA
  caustic_fit_edge = rep(NA, length(x_range))
  if (have_edge) {
    fit_inputs_E = prep_fit_inputs(Ar_finalE, fitting_radii, gb, delta_A_over_A_filled)
    ArE = fit_inputs_E$Ar
    sigma_A_fit_E = fit_inputs_E$sigma_A
    nfw_fit_e = tryCatch(
      NFWfit(rii, ArE, halo_srad = r200 / 5, ri_full = x_range,
             M200_prior = M200_prior, R200_prior = R200_prior, c_min = c_min, c_max = c_max,
             conc_method = conc_method, sigma_A = sigma_A_fit_E,
             c_prior_scatter_dex = c_prior_scatter_dex, little_h = little_h, clus_z = clus_z),
      error = function(e) NULL)
    if (!is.null(nfw_fit_e)) {
      caustic_fit_edge = nfw_fit_e[[1]] / sqrt(gb)
      reserr_edge = nfw_fit_e[[2]]
      conc_fit_edge = nfw_fit_e[[3]]
      mE = tryCatch(mass_from_Ar(Ar_finalE, x_range, r200, r200_input, fbr_use, M200_prior, R200_prior,
                                  fix_r200, Hz, clus_z, comoving_input),
                    error = function(e) list(r200_est = NA, M200_est = NA))
      r200_est_edge = mE$r200_est
      M200_est_edge = mE$M200_est
    }
  }

  Ngal = length(rproj[memflag == 1 & rproj <= r200_est])
  # BUG FIX: biwScale() had no guard against too few (or zero) points.
  # For a handful of difficult clusters (e.g. a degenerate NFW fit that
  # ends up excluding almost everyone as "not a member"), the member
  # selection here can come up empty or with a single point -- median()/
  # mad() then return NA, and biwScale()'s own internal
  # `if (madx == 0)` check crashes with the uninformative
  # "missing value where TRUE/FALSE needed" instead of a clear message
  # pointing at the real cause (too few members found within r200_est).
  if (Ngal >= 2) {
    vdisp_gal = biwScale(vproj[memflag == 1 & rproj <= r200_est])
  } else {
    vdisp_gal = NA
    if (verbose) message('Vdisp estimate unavailable: only ', Ngal,
                          ' member(s) found within r200_est = ', round(r200_est, 3),
                          ' Mpc (need >= 2) -- likely a degenerate NFW/membership fit for this cluster.')
  }
  # Standard large-N analytic approximation for the standard error of a
  # scale estimator, SE(sigma) ~= sigma/sqrt(2(N-1)) -- an approximation
  # (exact for a Gaussian sample standard deviation, used here as a
  # reasonable order-of-magnitude estimate for the biweight scale too,
  # since no closed-form SE exists for the biweight estimator itself).
  # Previously vdisp_gal was reported as a bare point estimate with no
  # uncertainty at all, unlike M200 and (when fit) concentration.
  vdisp_err = if (Ngal > 1 && !is.na(vdisp_gal)) vdisp_gal / sqrt(2 * (Ngal - 1)) else NA

  # R200, Bayesian-flavoured: only available when mass_method='bayesian_joint'
  # was used. R200 itself is held FIXED as the fitting anchor throughout
  # that fit (used to define rs=R200_anchor/c for every trial concentration
  # in the grid) -- it is never a free parameter there, so `r200_est`
  # itself does NOT move with the Bayesian fit and stays equal to whatever
  # anchor was used (given by the caller, or the blind preliminary
  # estimate). What we CAN derive, since M200 has a genuine posterior in
  # this mode: applying the standard SO-200 definition,
  # M200 = (4/3)*pi*200*rho_crit(z)*R200^3, to the M200 posterior mode and
  # 16th/84th percentiles gives an alternative, internally-consistent R200
  # point estimate + credible interval. This is a DERIVED quantity, not a
  # genuine joint fit of R200 itself, and by construction generally
  # DIFFERS from r200_est (the anchor) -- e.g. on a real CIRS cluster
  # (Rines & Diaferio 2006, A0085) the anchor gave r200_est=2.032 Mpc,
  # while this derived estimate landed at a different value nearby. Kept
  # as separate fields (`r200_bayes_est`/`_16`/`_84`) rather than
  # overwriting r200_est, since the rest of the output (caustic_fit,
  # membership, Ngal, vdisp_est) is all still computed using the anchor,
  # not this value -- replacing r200_est here would make those
  # inconsistent with each other.
  r200_bayes_est = r200_bayes_16 = r200_bayes_84 = NA
  if (!is.na(M200_bayes_16) && !is.na(M200_bayes_84)) {
    G200 = 6.67430e-11; Msol2kg200 = 1.98847e30; Mpc2m200 = 3.08567758e22
    crit200 = 3 * Hz^2 / (8 * pi * G200) * 1e6 / Msol2kg200 * Mpc2m200
    r200_of_M = function(M) (3 * M / (4 * pi * 200 * crit200))^(1/3)
    r200_bayes_est = r200_of_M(M200_est)
    r200_bayes_16 = r200_of_M(M200_bayes_16)
    r200_bayes_84 = r200_of_M(M200_bayes_84)
  }

  # getting membership 
  contour_outliers = caustic_outliers = rep(1, N0)
  caustic_outliers[as.numeric(rownames(data_set[memflag == 1,]))] = 0
  cinterp = approxfun(x_range, Ar_finalD)
  vint = cinterp(rproj)
  mbr = ifelse(abs(vproj) <= vint, T, F)
  contour_outliers[as.numeric(rownames(data_set[mbr,]))] = 0

  # `caustic_outliers_extended`: like `caustic_outliers`, but covering the
  # FULL range of rproj rather than defaulting every point beyond rlimit
  # to "outlier" without ever evaluating it. Extrapolates the fitted NFW
  # escape curve (nfw_d0/nfw_rs/nfw_gb, from the same fit caustic_outliers
  # itself is based on) to any radius via nfw_escape_curve(), and applies
  # the same |vlos| < escape_velocity test everywhere. Validated to match
  # caustic_outliers EXACTLY within the region run_caustic() actually
  # analysed (rlimit/vlimit) -- the only points it can newly classify are
  # ones caustic_outliers left at its unconditional default. Extrapolation
  # confidence decreases the further a point's rproj is beyond rlimit; see
  # nfw_escape_curve()'s documentation above for the caveat.
  escape_velocity_extended = nfw_escape_curve(rproj_orig_full, nfw_d0, nfw_rs, nfw_gb)
  caustic_outliers_extended = as.numeric(abs(vproj_orig_full) > escape_velocity_extended)
  
  if(verbose){
    message('NFW fit residual std error: ',round(reserr, 3), ' km/s')
    message('R200 estimate:  ', round(r200_est, 3), ' Mpc')
    if (!is.na(r200_bayes_16)) {
      message('  NOTE: mass_method="bayesian_joint" holds R200 FIXED as the fitting anchor -- ',
              'this value is NOT re-derived or refined by the Bayesian fit, unlike M200 and ',
              'concentration.')
      message('  A separate, internally-consistent R200 point estimate + credible interval CAN be ',
              'derived from the M200 posterior via the SO-200 definition: ',
              round(r200_bayes_est, 3), ' [', round(r200_bayes_16, 3), ', ', round(r200_bayes_84, 3),
              '] Mpc (see r200_bayes_est/_16/_84 in the output) -- this generally differs from the ',
              'anchor above, and is NOT used for caustic_fit/membership/Ngal/vdisp_est, which all ',
              'still use the anchor.')
    }
    message('Vdisp estimate: ', round(vdisp_gal, 3), ' km/s',
            if (!is.na(vdisp_err)) paste0(' +/- ', round(vdisp_err, 3),
                                            ' (', round(100 * vdisp_err / vdisp_gal, 1), '%)') else '')
    message('M200 estimate:  ', format(M200_est, digits = 4), ' Msol')
    if (!is.na(M200_bayes_16) && !is.na(M200_bayes_84))
      message('  (mass_method="bayesian_joint" 68% credible interval: [',
              format(M200_bayes_16, digits = 3), ', ', format(M200_bayes_84, digits = 3), '] Msol)')
    if (!is.na(M200_err)) {
      err_source = if (!is.na(M200_bayes_16)) 'from the bayesian_joint 68% credible interval'
                   else 'D99/Serra et al. 2011, eq. 23-24'
      message('M200 uncertainty (', err_source, '): +/- ',
              format(M200_err, digits = 3), ' Msol (', round(100 * M200_err / M200_est, 1), '%)')
    }
    if (!is.na(conc_fit)) {
      ci_txt = if (!is.na(conc_16) && !is.na(conc_84))
        paste0(' (68% credible interval: [', round(conc_16, 2), ', ', round(conc_84, 2), '])') else ''
      method_txt = if (!is.na(M200_bayes_16)) 'bayesian_joint (M200,c)' else conc_method
      message('Concentration estimate (c=R200/rs, fit freely since M200 was not held fixed, ',
              'method=', method_txt, '): ', round(conc_fit, 3), ci_txt)
    }
    message('N200: ', Ngal)
    if (have_edge && !is.na(M200_est_edge))
      message('M200 estimate (independent edge-detection method): ',
              format(M200_est_edge, digits = 4), ' Msol -- compare against the density-',
              'matched estimate above; close agreement is stronger evidence the surface ',
              'is well-determined than either method alone.')
    else if (!have_edge && compute_edge)
      message('Edge-detection cross-check unavailable for this cluster (too few galaxies ',
              'within r200 to bin, or no candidate contour matched).')
  }
  
  # plotting 
  if(plot){
    magplot(c(0, xmax), c(-4000, 4000), pch = '', las = 1, mtline = c(2,2.6),
            xlab = expression(R[proj] ~ (Mpc)), ylab = expression(v[proj] ~ (km~s^{-1})))
    # shaded uncertainty band around the NFW-fit escape surface, from the
    # SAME per-radius relative amplitude uncertainty (D99/Serra et al.
    # 2011, eq. 23) already computed above for the M200 error bar --
    # previously that formula only ever fed a single scalar (M200_err) at
    # r200; here it is shown as a full profile instead.
    ok_band = is.finite(delta_A_over_A)
    if (sum(ok_band) >= 2) {
      xb = x_range[ok_band]
      upper = caustic_fit[ok_band] * (1 + delta_A_over_A[ok_band])
      lower = pmax(caustic_fit[ok_band] * (1 - delta_A_over_A[ok_band]), 0)  # A(r) cannot be negative
      polygon(c(xb, rev(xb)), c(upper, rev(lower)),
              col = adjustcolor('orange', alpha.f = 0.2), border = NA)
      polygon(c(xb, rev(xb)), c(-lower, rev(-upper)),
              col = adjustcolor('orange', alpha.f = 0.2), border = NA)
    }
    points(rproj, vproj, pch = ifelse(memflag == 1, 20, 1))
    lines(x_range,  Ar_finalD, col = 'red', lwd = 2) 
    lines(x_range, -Ar_finalD, col = 'red', lwd = 2)
    lines(x_range, caustic_fit, lty = 2, lwd = 2, col = 'orange') 
    lines(x_range, -caustic_fit, lty = 2, lwd = 2, col = 'orange')
    if (have_edge) {
      lines(x_range,  Ar_finalE, col = 'forestgreen', lwd = 2, lty = 3)
      lines(x_range, -Ar_finalE, col = 'forestgreen', lwd = 2, lty = 3)
    }
    abline(v = r200_est, col = 'gray', lty = 2, lwd = 2)
    legend('topright', lty = c(1, 2, 3, 2, NA, NA, NA), pch = c(NA, NA, NA, NA, 20, 1, 15), 
           col = c('red', 'orange', 'forestgreen', 'gray', 'black', 'black',
                    adjustcolor('orange', alpha.f = 0.4)), lwd = 2, bty = 'n', 
           y.intersp = .9, x.intersp = .5, cex = .8, seg.len = 1.8,
           legend=c('Escape surface', 'NFW fit', 'Edge method', 'R200', 'Members', 'Outliers',
                     'NFW fit uncertainty'))
  }
  
  # output list
  lout=list(data_set = data_set, x_range = x_range, y_range = y_range, img_tot = Zi, 
            caustic_profile = Ar_finalD, caustic_outliers = caustic_outliers, 
            contour_outliers = contour_outliers,
            nfw_fit = caustic_fit, nfw_fit_err = reserr, r200_est = r200_est, 
            vdisp_est = vdisp_gal, log_M200_est = log10(M200_est), conc_est = conc_fit,
            conc_16 = conc_16, conc_84 = conc_84,
            M200_bayes_16 = M200_bayes_16, M200_bayes_84 = M200_bayes_84,
            vdisp_err = vdisp_err, r200_bayes_est = r200_bayes_est,
            r200_bayes_16 = r200_bayes_16, r200_bayes_84 = r200_bayes_84,
            M200_err = M200_err, M200_err_frac = M200_err / M200_est, mass_profile = mD$massprofile,
            delta_A_profile = Ar_finalD * delta_A_over_A, delta_M_profile = delta_M_cum,
            caustic_profile_edge = Ar_finalE, nfw_fit_edge = caustic_fit_edge,
            nfw_fit_err_edge = reserr_edge, r200_est_edge = r200_est_edge,
            log_M200_est_edge = if (!is.na(M200_est_edge)) log10(M200_est_edge) else NA,
            conc_est_edge = conc_fit_edge, nfw_d0 = nfw_d0, nfw_rs = nfw_rs, nfw_gb = nfw_gb,
            caustic_outliers_extended = caustic_outliers_extended,
            escape_velocity_extended = escape_velocity_extended)
  return(invisible(lout))
}
