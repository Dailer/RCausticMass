# Caustic Mass Estimator for Galaxy Clusters 
This is an R implementation of the Caustic Mass Estimator for Galaxy Clusters developed in Python by [Gifford et al.](https://github.com/giffordw/CausticMass)
A pre-cleaning of interlopers is recommended, using e.g., the shifting-gapper technique.

The caustic technique is a powerful method to infer cluster mass profiles to clustrocentric distances well beyond the virial radius. It relies in the measure of the escape velocity of the system using only galaxy redshift information. This method was introduced by [Diaferio & Geller (1997)](https://arxiv.org/abs/astro-ph/9701034) and [Diaferio (1999)](http://arxiv.org/abs/astro-ph/9906331).

Code Usage
----------
Install the required libraries
```
install.packages(c("magicaxis","gplots","imager","pracma"))
```
Source the R file functions
```
source("RCausticMass.R")
```
Run the code, here using the sample data provided with unknown values of R200 and cluster velocity dispersion
```
data = read.table("sample_data.txt")
data = subset(data, CID == 1)
r = run_caustic(data$dproj, data$vlos, data$zclus, r200 = NA, clus_vdisp = NA)
```
The smoothed phase-space of galaxies can be plotted
```
image(r$x_range, r$y_range, r$img_tot, asp = NA, las = 1, xlab = expression(R[proj] ~ (Mpc)), 
      ylab = expression(v[proj] ~ (km/s)))
```
----------
Author: Dailer F. Morell

This code is public and may be used for research purposes.
