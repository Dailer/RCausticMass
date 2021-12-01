# Biweight scale estimator from Beers (1990)
biwScale = function(x, c = 9){
  M = median(x, na.rm = T) 
  d = x - M
  u = d / (c * mad(x, M, 2.21914, na.rm = T))
  uT = 1 - u^2 
  wu = which(abs(u) < 1)
  top = sqrt(sum((d^2 * uT^4)[wu]))
  bottom = abs(sum((uT * (1 - 5 * u^2))[wu]))
  SBI = sqrt(length(wu)) * top / bottom
  if(is.nan(SBI)) SBI = x[1]
  return(SBI)
}

# Uses a 2D gaussian kernel to estimate the density of the phase space
gaussian_kernel = function(dproj, vlos, r200, normalization = 100, q = 10, xmax = 6, ymax = 5e3, by = 0.05, 
                           mirror = F, plot = F){
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
  ksize = (4 / (3 * length(xvalues)))^(1 / 5) * sqrt((bwx^2 + bwy^2) / 2)
  ymx = ymax / normscale
  h = hist2d(c(xvalues, 0, xmax, xmax, 0), c(yvalues, -ymx, -ymx, ymx, ymx), c(xres, yres), show = F)
  h$counts[cbind(c(1, xres, xres, 1), c(1, 1, yres, yres))] = 0
  img = as.matrix(isoblur(as.cimg(h$counts), sigma = ksize))
     
  if(plot){
    magimage(x_range, y_range, img, asp = NA, xlab = expression(R[proj] ~ (Mpc)), 
             ylab = expression(v[proj] ~ (km~s^{-1})))
  }
  lout = list(x = x_range, y = y_range, z = img)
  return(lout)
}

# This function will use contourLines() to locate escape surfaces
findcontours = function(Zi, ri, vi, r200, vvar, rimax = 4, nlevels = 200, plot = T, verbose = T){

  cl = contourLines(ri, vi, Zi, nlevels = nlevels)
   
  # selecting contours that are "full" and don't loop back only in positive or negative space
  contours = fullcont = list()
  for (i in 1:length(cl)) {
    cln = cl[[i]]
    if(max(cln$x) >= rimax & min(cln$x) <= 0 & max(cln$y) > 0 & min(cln$y) < 0 & 
       max(abs(cln$y)) < 4000){
      #find positive/negative contours
      xcont_u = cln$x[cln$y > 0] 
      ycont_u = cln$y[cln$y > 0] 
      xcont_d = cln$x[cln$y < 0]
      ycont_d = cln$y[cln$y < 0]
      
      if(length(xcont_u) < 2 | length(xcont_d) < 2) next;
      
      #match contour grid to r grid by linear interpolation
      fint_u = approxfun(xcont_u, ycont_u, yleft = 0, yright = 0, ties = mean)
      fint_d = approxfun(xcont_d, ycont_d, yleft = 0, yright = 0, ties = mean)
      y_u = fint_u(ri)
      y_d = fint_d(ri)
      y_f = apply(cbind(y_u, abs(y_d)), 1, min)
            
      fullcont = c(fullcont, list(cln))
      contours = c(contours, list(y_f))
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
  
  if(verbose) message('Complete')
  
  # plotting
  if(plot){
    magimage(ri, vi, Zi, asp = NA, xlab = expression(R[proj] ~ (Mpc)), 
             ylab = expression(v[proj] ~ (km~s^{-1})), col.ticks = 'white')
    ll = lapply(fullcont, lines, col = '#FFBE74', lty = 3)
    lines(ri, final_contour, col = '#F13005', lwd = 2)
    lines(ri, -final_contour, col = '#F13005', lwd = 2)
  }
  lout=list(caustic = final_contour, contours = contours)
  return(lout)
}

# Fitting a NFW profile to the caustic surface
NFWfit = function(rii, Ar, halo_srad, ri_full){
  min_func = function(x, d0, rs)
	sqrt(2 * 4 * pi * 4.3e-9 * d0 * rs^2 * log(1 + x / rs) / (x / rs))
	
  out = tryCatch(nls(Ar ~ min_func(rii, d0, rs = halo_srad),
                     start=list(d0 = 1e14)), error = function(e) NA)  
  res = resid(out)
  res.err = sqrt(sum(res^2) / (length(res) - 1))
  halo_scale_density = coef(out)
  halo_scale_density_e = tryCatch(sqrt(halo_scale_density), error = function(e) 1e14)
  if(ri_full[1] == 0) ri_full[1] = ri_full[2]
  profile = min_func(ri_full, halo_scale_density, halo_srad)
  return(list(profile, res.err))
}

# Main function
run_caustic=function(rproj, vproj, clus_z, r200 = NA, clus_vdisp = NA, rlimit = 5.8, 
                     vlimit = 4000, xmax = 6, ymax = 4000, mirror = F, Om = .3, H0 = 70, 
                     fbr = 0.65, q = 10, beta = NA, centering = F, plot = T, verbose = T){
  
  N0 = length(rproj)
  clus_z = clus_z[1] 
  Ol = 1 - Om
  Hz = H0 * sqrt(Om * (1 + clus_z)^3 + Ol) # H(z)
  data_set = data.frame(r = rproj, v = vproj)
  
  # reduce sample within limits
  data_set = subset(data_set, r < rlimit & abs(v) < vlimit)
  if(nrow(data_set) < 5) stop('Data set has too few elements inside the limits')
  if(verbose) message('Data Size: ', nrow(data_set))
  rmax = max(data_set[,1])
  
  # if r200 is not available we provide an initial estimate
  if(is.na(r200)){
    rl = min(rmax, 3)
    vdisp_prelim = biwScale(data_set[,2][data_set[,1] < rl])
    r200 = sqrt(3) * vdisp_prelim / (10 * Hz)
    if(r200 > rmax) r200 = rmax
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
  if(is.na(clus_vdisp)) 
	clus_vdisp = biwScale(vproj[rproj < r200 & abs(vproj) < vlimit])
	
  if(verbose) message('Initial R200:  ', round(r200, 3), ' Mpc')
  if(verbose) message('Initial Vdisp: ', round(clus_vdisp, 3), ' km/s')
  
  # calculating density of phase-space 
  if(!mirror & verbose) message('Calculating density')
  if(mirror & verbose) message('Calculating density w/Mirrored Data')
    
  gk = gaussian_kernel(rproj, vproj, r200, Hz, q, xmax, ymax, 0.05, mirror, plot = F)
  x_range = gk$x
  y_range = gk$y
  Zi = gk$z / max(abs(gk$z))
  
  # identifying the initial caustic surface and members within the surface 
  if(verbose) message('Calculating initial surface...')
  fcont = findcontours(Zi, x_range, y_range, r200, clus_vdisp^2, rlimit, plot = F, 
					   verbose = verbose)
  Ar_finalD = fcont$caustic
  contours = fcont$contours
  
  # fit an NFW to the optimal caustic surface 
  if(is.na(beta)) beta = rep(0.2, length(x_range))
  gb = (3 - 2 * beta) / (1 - beta)
  fitting_radii = (x_range >= r200 / 3) & (x_range <= r200)
  rii = x_range[fitting_radii]
  ArD = Ar_finalD[fitting_radii] * sqrt(gb[fitting_radii])
  nfw_fit = NFWfit(rii, ArD, halo_srad = r200 / 5, ri_full = x_range)
  caustic_fit = nfw_fit[[1]] / sqrt(gb)
  reserr = nfw_fit[[2]]
  
  fcomp = approxfun(x_range, caustic_fit)
  vcompare = fcomp(data_set[,1])
  memflag = ifelse(abs(vcompare) > abs(data_set[,2]), 1, 0)
  
  # mass estimation 
  require(pracma)
  G = 6.67430e-11
  Msol2kg = 1.98847e30
  Mpc2m = 3.08567758e22
  # critical density in Msun/Mpc^3
  crit = 3 * Hz^2 / (8 * pi * G) * 1e6 / Msol2kg * Mpc2m / (1 + clus_z)^3 
  rsi = x_range * Mpc2m   # in meters
  Asi = Ar_finalD * 1e3   # in m/s
  sumtot = cumtrapz(rsi, fbr * Asi^2)[,1]
  massprofile = sumtot / (G * Msol2kg)
  avg_density = massprofile / (4 / 3 * pi * (x_range)^3)
  finterp = approxfun(avg_density, x_range)
  r200_est = finterp(200 * crit)
  r500_est = finterp(500 * crit)
  finterp = approxfun(x_range, massprofile)
  M200_est = finterp(r200_est)

  Ngal = length(rproj[memflag == 1 & rproj <= r200_est])
  vdisp_gal = biwScale(vproj[memflag == 1 & rproj <= r200_est])
  
  # getting membership 
  contour_outliers = caustic_outliers = rep(1, N0)
  caustic_outliers[as.numeric(rownames(data_set[memflag == 1,]))] = 0
  cinterp = approxfun(x_range, Ar_finalD)
  vint = cinterp(rproj)
  mbr = ifelse(abs(vproj) <= vint, T, F)
  contour_outliers[as.numeric(rownames(data_set[mbr,]))] = 0
  
  if(verbose){
    message('NFW fit residual std error: ',round(reserr, 3), ' km/s')
    message('R200 estimate:  ', round(r200_est, 3), ' Mpc')
    message('Vdisp estimate: ', round(vdisp_gal, 3), ' km/s')
    message('M200 estimate:  ', format(M200_est, digits = 4), ' Msol')
    message('N200: ', Ngal)
  }
  
  # plotting 
  if(plot){
    magplot(c(0, xmax), c(-4000, 4000), pch = '', las = 1, mtline = c(2,2.6),
            xlab = expression(R[proj] ~ (Mpc)), ylab = expression(v[proj] ~ (km~s^{-1})))
    points(rproj, vproj, pch = ifelse(memflag == 1, 20, 1))
    lines(x_range,  Ar_finalD, col = 'red', lwd = 2) 
    lines(x_range, -Ar_finalD, col = 'red', lwd = 2)
    lines(x_range, caustic_fit, lty = 2, lwd = 2, col = 'orange') 
    lines(x_range, -caustic_fit, lty = 2, lwd = 2, col = 'orange')
    abline(v = r200_est, col = 'gray', lty = 2, lwd = 2)
    legend('topright', lty = c(1, 2, 2, NA, NA), pch = c(NA, NA, NA, 20, 1), 
           col = c('red', 'orange', 'gray', 'black', 'black'), lwd = 2, bty = 'n', 
           y.intersp = .8, x.intersp = .5, cex = .8, seg.len = 1.8,
           legend=c('Escape surface', 'NFW fit', 'R200', 'Members', 'Outliers'))
  }
  
  # output list
  lout=list(data_set = data_set, x_range = x_range, y_range = y_range, img_tot = Zi, 
            caustic_profile = Ar_finalD, caustic_outliers = caustic_outliers, 
            nfw_fit = caustic_fit, nfw_fit_err = reserr, r200_est = r200_est, 
            vdisp_est = vdisp_gal, log_M200_est = log10(M200_est))
  return(invisible(lout))
}
