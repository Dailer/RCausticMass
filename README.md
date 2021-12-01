# Caustic Mass Estimator for Galaxy Clusters 
This is an R implementation of the Caustic Mass Estimator for Galaxy Clusters developed in Python by [Gifford et al.](https://github.com/giffordw/CausticMass)

The caustic technique is a powerful method to infer cluster mass profiles to clustrocentric distances well beyond the virial radius, by measuring the escape velocity of the sistem using only galaxy redshift information.

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
Run the code, here using the sample data
```
data = read.table("sample_data.txt")
data = subset(data, CID == 1)
run_caustic(data$dproj, data$vlos, data$zclus)
```
----------
Author: Dailer F. Morell

**Copyright 2021 the author**

This code is public and may be used for research purposes.
