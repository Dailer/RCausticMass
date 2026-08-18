# Biweight scale estimator (Beers et al. 1990, eq. 9). Robust alternative
# to the standard deviation, downweighting outliers via a tuning constant
# c=9 applied to the raw median absolute deviation.
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

##############################################################################
# Shifting-gapper interloper removal (Fadda et al. 1996). Sorts galaxies
# into radial bins and rejects points separated from the local velocity
# distribution by a gap exceeding gap_kms, iterating until stable.
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
# 18-20), alternative to gaussian_kernel()'s fixed-bandwidth blur. Uses a
# compact quartic kernel with a LOCAL bandwidth per galaxy -- wider in
# sparse regions, narrower in dense ones -- via h_i = h_c*h_opt*lambda_i,
# lambda_i = sqrt(gamma/f_1(x_i)), with f_1 a pilot fixed-bandwidth
# density estimate. Rescales (r,v) so q=sigma_v/sigma_r takes a fixed
# value (q=25 default) before estimating the density, then rescales back.
# Both the pilot density and final grid mirror r and v (r -> -r, v -> -v)
# to avoid depleted density near r=0, matching gaussian_kernel().
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

# 2D Gaussian kernel density estimate of the phase space (dproj, vlos).
# Mirrors the galaxy distribution to negative r before estimating the
# density (Serra et al. 2011, Section 4.3), avoiding artificial depletion
# of the density near r=0 from the small number of galaxies there; only
# the r>=0 half of the resulting map is kept.
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
  h = hist2d(c(xvalues_m, -xmax, xmax, xmax, -xmax), c(yvalues_m, -ymx, -ymx, ymx, ymx),
             c(length(x_range_full), yres), show = F)
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

# Locates escape-velocity surfaces via contourLines() over the phase-space
# density map. Kappa levels are log-spaced (from min(Zi[Zi>0])/5 to
# max(Zi)) rather than linear, since the density is sharply peaked near
# the centre and sparse in the outskirts. As each candidate contour is
# walked outward, its amplitude is gradient-restricted (Serra et al. 2011):
# not allowed to grow faster than r^0.5 or fall faster than r^-2 between
# radial steps, applied to the positive/negative branches separately
# before combining as A(r) = min(|branch_up|, |branch_down|).
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

# Fits a NFW profile to the caustic amplitude curve.
#
# M200_prior / R200_prior (optional): when supplied (e.g. from an external
# catalog), M200 is fixed exactly to M200_prior and only the concentration
# (via the scale radius rs) is fit to this cluster's own data -- the data
# set the profile's shape, the trusted external mass sets its
# normalisation. Without a prior, an assumed universal concentration is
# used (halo_srad, typically r200/5) and only the normalisation (d0) is
# fit -- fitting (d0, rs) both free at once was found numerically
# unreliable with nls()'s default algorithm over typical fitting ranges.
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
    # concentration is the only free parameter -- a fine 1D grid gives the
    # full posterior directly, no MCMC needed.
    #
    # PRIOR: log-normal in c, centred on the Duffy et al. (2008) c-M
    # relation (c200 = 5.71*(M200*h/2e12)^-0.084*(1+z)^-0.47), ~0.1 dex
    # scatter. Replaces the hard [c_min,c_max] truncation of the 'bounded'
    # method, which the fit can pile up against with noisy data -- the
    # posterior instead shrinks toward the prior when data don't constrain
    # c well, and toward the data when they do.
    #
    # LIKELIHOOD: Gaussian, using `sigma_A` (the per-radius amplitude
    # uncertainty from the D99/Serra et al. 2011 error formula) as noise.
    if (is.null(sigma_A) || all(!is.finite(sigma_A) | sigma_A <= 0))
      sigma_A = rep(max(sd(Ar), 1e-3), length(Ar))  # fallback if no usable error estimate was supplied
    sigma_A = pmax(sigma_A, 1e-3 * mean(Ar))  # guard against a zero/near-zero noise floor

    # Grid spans a wide, fixed range independent of c_min/c_max, so the
    # prior (not a wall) is what keeps the posterior physically sensible.
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
# findsurface()). A SECOND, independent way to pick the escape surface,
# alongside the main S(kappa)=4<v^2> density matching: bins galaxies by
# projected radius, takes the average of the top `edge_perc` most extreme
# positive/negative velocities per bin as an empirical phase-space edge,
# and picks whichever candidate contour (from findcontours()) best
# matches it. Useful as a cross-check against the density-based surface.
##############################################################################
##############################################################################
# Joint Bayesian fit of (M200, concentration), with R200 held as a FIXED
# ANCHOR radius (given by the caller or auto-estimated earlier). An NFW
# halo has only two independent degrees of freedom (e.g. M200 and c; R200
# then follows from M200 via the SO-200 definition) -- this is a 2D grid
# posterior over (M200, c), generalising the single-parameter Bayesian
# concentration fit above. Priors: c uses the same Duffy et al. (2008)
# log-normal; M200 gets a weak log-normal prior centred on the SO-200 mass
# implied by the anchor radius (0.4 dex scatter default). Likelihood: same
# Gaussian model as the 1-parameter case, using the D99/Serra et al.
# (2011) delta_A(r)/A(r) error estimate as the per-radius noise.
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
  sumtot = cumtrapz(rsi, fbr * Asi^2)[,1]

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
# pipeline `n_boot` times. This is empirical (how much the answer moves
# under resampling), distinct from the analytic D99/Serra et al. (2011)
# error the main function reports by default -- the two need not agree,
# treat them as complementary rather than interchangeable. Slow: cost
# scales linearly with n_boot, meant for inspecting individual clusters
# rather than large-sample batch validation.
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
# (the dominant failure mode with fixed-radius extractions, where the
# default rlimit can still demand the contour reach further than the real
# signal extends). Kept as a separate wrapper rather than folded into
# run_caustic() itself, so that function's own behaviour stays simple and
# predictable. Recovered results are LESS PRECISE on average than normal
# convergence -- this is a coverage/precision trade-off, not a free win.
#
# IMPORTANT: a recovered result's own M200_err_frac is NOT a reliable
# safety net -- it can be small and reassuring even when the result is
# wrong by a large factor, since the D99/Serra et al. (2011) error formula
# has no way to know the fitting radius was forced down artificially. The
# only reliable signal that a result came from this fallback is
# `rlimit_frac_used` being non-NA -- treat every such result as
# lower-confidence regardless of its own error bar.
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
# Extends outlier classification to the FULL range of dproj, beyond the
# `rlimit` boundary run_caustic() actually analyses. `caustic_outliers`
# defaults every point to 1 (outlier) and only overwrites it for points
# both within `rlimit` AND below the fitted escape curve there -- points
# beyond rlimit keep the default regardless of their actual |vlos|. This
# reconstructs the fitted NFW escape-velocity curve from `nfw_d0`/`nfw_rs`/
# `nfw_gb` and evaluates it analytically at any radius: inside rlimit it
# reproduces caustic_outliers exactly, outside it extrapolates the same
# curve. CAVEAT: an extrapolation is only as good as the fit it extends --
# treat outlier calls far outside rlimit as less certain than ones inside.
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

# Main function
#
# `fbr` (the F_beta factor in the mass integral) has no universal default
# that works well across datasets -- its optimal value depends on the
# extraction geometry (fixed-radius vs. proportional-to-R200, wide vs.
# narrow velocity window) and on kernel choice, not just on the cluster
# sample itself. Recalibrate against known masses/R200 when possible
# rather than trusting a single default for a new dataset.
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
  # `bayesian`: single switch for "get halo parameters the Bayesian way" --
  # picks conc_method='bayesian' (1-parameter fit) if both M200_prior and
  # R200_prior are given, or mass_method='bayesian_joint' (2-parameter,
  # with r200 as anchor) otherwise. Explicit conc_method=/mass_method=
  # still work if you want to force one by hand.
  if (bayesian) {
    if (!is.na(M200_prior) & !is.na(R200_prior)) conc_method = 'bayesian'
    else mass_method = 'bayesian_joint'
  }
  conc_method = match.arg(conc_method)
  mass_method = match.arg(mass_method)
  # `hc` (Gaussian blur width): defaults to 1.3 when fix_r200=TRUE
  # (r200/vdisp externally known) and 1.0 otherwise -- validated
  # separately for each case; the two do not transfer. Pass explicitly to
  # override.
  if (is.na(hc)) hc = if (kernel == 'gaussian' && fix_r200) 1.3 else 1
  if (fix_r200 & (is.na(r200) | is.na(clus_vdisp)))
    stop('run_caustic(): fix_r200=TRUE requires both r200 and clus_vdisp to be supplied ',
         '(e.g. from an external catalog) -- there is nothing to hold fixed otherwise.')
  r200_input = r200  # remember the caller's original value; r200 may be reassigned below
                      # (auto-estimate if NA, or gaussian_kernel()'s own internal adjustment)
  
  N0 = length(rproj)
  # rproj/vproj get REASSIGNED later (to the trimmed data_set's columns) --
  # keep a copy of the original, full-length arrays for
  # caustic_outliers_extended below.
  rproj_orig_full = rproj
  vproj_orig_full = vproj
  clus_z = clus_z[1] 
  Ol = 1 - Om
  Hz = H0 * sqrt(Om * (1 + clus_z)^3 + Ol) # H(z)

  # `adaptive_vlimit` (default TRUE): instead of a fixed `vlimit` (e.g. the
  # classic 4000 km/s), derive it from the cluster's own velocity
  # dispersion -- the supplied `clus_vdisp` if known, otherwise a rough
  # preliminary estimate (shifting_gapper-cleaned biweight scale within
  # min(rmax,3) Mpc). A fixed, wide velocity window (common when
  # candidates are pulled from an external catalogue without knowing true
  # membership beforehand) systematically hurts precision compared to a
  # window tied to the cluster's own dynamics. `adaptive_vlimit_mult`
  # (default 3.5) sets the window as a multiple of that dispersion --
  # smaller values tighten the window (less bias, worse coverage), larger
  # values do the reverse. Only applies when `vlimit` isn't already given.
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

  # `xmax`/`rlimit` default to the data's OWN extent rather than a fixed
  # number: xmax = max(rproj) (the true data extent), rlimit = 0.80*xmax,
  # leaving room for the contour-tracing step to close before the true
  # edge of the data instead of demanding it reach exactly there ("the
  # contours do not expand to the radial limit" otherwise). Pass either
  # explicitly to override -- 0.80 is a compromise that works reasonably
  # for both fixed-radius and proportional-to-R200 extractions, but a
  # value tuned to your own extraction convention may do better.
  if (is.na(xmax)) xmax = max(rproj)
  if (is.na(rlimit)) rlimit = 0.80 * xmax

  data_set = data.frame(r = rproj, v = vproj)
  
  # reduce sample within limits
  data_set = subset(data_set, r < rlimit & abs(v) < vlimit)
  if(nrow(data_set) < 5) stop('Data set has too few elements inside the limits')
  if(verbose) message('Data Size: ', nrow(data_set))
  rmax = max(data_set[,1])
  
  # if r200 is not available we provide an initial estimate
  #
  # Iterative preliminary R200: start from rl=min(rmax,3), clean
  # interlopers with shifting_gapper(), compute a trial R200 via the
  # virial-scaling formula, then re-estimate using a NARROWER window
  # (exactly the trial R200, not a fixed multiple) and repeat until the
  # estimate changes by less than 5% between iterations (or a small
  # iteration cap / floor radius is hit) -- mimics computing the velocity
  # dispersion strictly within R200 rather than an arbitrary wider window.
  # shifting_gapper() cleaning is applied at every iteration, falling back
  # to the uncleaned biweight scale if it fails or too few galaxies remain.
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
  
  # if vdisp is not available we provide a robust initial estimate.
  # For a poor cluster the radius/velocity cut can select zero galaxies;
  # falls back to the full data_set in that case rather than propagating NA.
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
  # galaxy velocities -- independent of the density used for Ar_finalD.
  # `compute_edge` (default FALSE): this cross-check is less precise than
  # the main density-matched estimate and doesn't help when averaged with
  # it, but remains a useful per-cluster diagnostic; set TRUE to compute it.
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
  # beta(r) = 0.5*r/(r+R200/c) from Gifford, Miller & Kern (2013), instead
  # of the fixed beta=0.2 used by default. Needs a concentration BEFORE
  # the NFW fit runs, so `conc_default` is used regardless of whatever
  # concentration the fit itself later finds (no iterative re-fit).
  # Validated as a modest improvement only when fix_r200=TRUE (R200/vdisp
  # known externally); does not transfer to the blind case.
  if (beta_radial) {
    beta = beta_radial_gifford(x_range, r200, conc_default)
  } else if (is.na(beta)) {
    beta = rep(0.2, length(x_range))
  }
  gb = (3 - 2 * beta) / (1 - beta)
  fitting_radii = (x_range >= r200 / 3) & (x_range <= r200)
  rii = x_range[fitting_radii]

  # Uncertainty estimate (D99/Serra et al. 2011, eq. 23-24):
  #   delta_A(r)/A(r) = kappa / max_v[f_q(r,v)]
  # where kappa is the density threshold that located the chosen escape
  # surface. Doubles as the per-radius noise model for the Bayesian
  # concentration fit below, and for the M200 error bar / member_boundary
  # ='upper'. Only applies to the density-matched (Ar_finalD) surface.
  peak_density = apply(Zi, 1, max)
  delta_A_over_A = fcont$kappa / peak_density
  delta_A_over_A[!is.finite(delta_A_over_A) | peak_density <= 0] = NA
  # At radii far beyond where data constrain the density, this ratio can
  # blow up; cap the relative uncertainty at 200%.
  delta_A_over_A = pmin(delta_A_over_A, 2)
  delta_A_over_A_filled = delta_A_over_A
  delta_A_over_A_filled[is.na(delta_A_over_A_filled)] = 0  # treat "unknown" as "no widening", not "infinitely uncertain"
  fit_inputs_D = prep_fit_inputs(Ar_finalD, fitting_radii, gb, delta_A_over_A_filled)
  ArD = fit_inputs_D$Ar
  sigma_A_fit = fit_inputs_D$sigma_A

  # if the caller supplied M200_prior but not R200_prior, anchor
  # R200_prior to the r200 already in use here rather than silently
  # ignoring M200_prior for lack of a matching radius.
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
  # Raw NFW fit parameters (d0, rs) and the Fbeta normalisation `gb`,
  # exposed so the caller can extrapolate the fitted escape-velocity curve
  # beyond rlimit -- see nfw_escape_curve() / extend_outliers_nfw() above.
  nfw_d0 = nfw_fit[[6]]
  nfw_rs = nfw_fit[[7]]
  nfw_gb = gb[1]

  # `member_boundary`: which curve a galaxy's velocity is compared against
  # to decide membership. 'fit' (default) uses the point-estimate NFW
  # curve directly. 'upper' uses the upper edge of the delta_A(r)
  # uncertainty band -- more conservative about flagging outliers, at the
  # cost of admitting more contamination; a real trade-off, not a strict
  # improvement.
  member_curve = if (member_boundary == 'upper') caustic_fit * (1 + delta_A_over_A_filled) else caustic_fit
  fcomp = approxfun(x_range, member_curve)
  vcompare = fcomp(data_set[,1])
  memflag = ifelse(abs(vcompare) > abs(data_set[,2]), 1, 0)
  
  # `fbr_radial` (opt-in, off by default): use the radially-varying
  # F_beta(r) form (DG97/D99/Gifford et al. 2013), instead of the constant
  # `fbr` scalar. Uses conc_default, same reasoning as beta_radial above.
  fbr_use = if (fbr_radial) fbeta_radial_nfw(x_range, r200, conc_default, gb) else fbr

  # mass estimation (see mass_from_Ar() for the DG97-anchoring / fix_r200 /
  # r=0-handling logic, shared between the density-matched and edge curves)
  #
  # `mass_method='bayesian_joint'`: alternative to the density-integral
  # M200 -- a joint (M200, concentration) Bayesian fit to the observed
  # A(r) curve, with r200 held as a fixed anchor. Requires M200_prior to
  # be unset (if M200 is already fixed, use conc_method='bayesian' instead
  # to fit concentration alone).
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
    # bayesian_joint succeeded: use the credible interval directly.
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
  # biwScale() needs a guard against too few/zero points: for a degenerate
  # fit that excludes almost everyone as "not a member", the member
  # selection here can come up empty.
  if (Ngal >= 2) {
    vdisp_gal = biwScale(vproj[memflag == 1 & rproj <= r200_est])
  } else {
    vdisp_gal = NA
    if (verbose) message('Vdisp estimate unavailable: only ', Ngal,
                          ' member(s) found within r200_est = ', round(r200_est, 3),
                          ' Mpc (need >= 2) -- likely a degenerate NFW/membership fit for this cluster.')
  }
  # Standard large-N approximation for the standard error of a scale
  # estimator, SE(sigma) ~= sigma/sqrt(2(N-1)).
  vdisp_err = if (Ngal > 1 && !is.na(vdisp_gal)) vdisp_gal / sqrt(2 * (Ngal - 1)) else NA

  # R200, Bayesian-flavoured: only available when mass_method='bayesian_joint'
  # was used. R200 is held FIXED as the fitting anchor throughout that fit
  # (used to define rs=R200_anchor/c for every trial concentration), so
  # `r200_est` itself does NOT move with the fit. What CAN be derived,
  # since M200 has a genuine posterior in this mode: applying the SO-200
  # definition (M200 = (4/3)*pi*200*rho_crit(z)*R200^3) to the M200
  # posterior mode and 16th/84th percentiles gives an alternative,
  # internally-consistent R200 point estimate + credible interval. This is
  # a DERIVED quantity, not a genuine fit of R200, and generally DIFFERS
  # from r200_est (the anchor) -- kept as separate fields
  # (`r200_bayes_est`/`_16`/`_84`) since the rest of the output
  # (caustic_fit, membership, Ngal, vdisp_est) still uses the anchor.
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
