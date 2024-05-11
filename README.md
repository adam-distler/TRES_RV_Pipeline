# TRES_RV_Pipeline
A pipeline for determining radial velocities of targets using either synthetic spectra or actual spectra from stars with known RV's. This was created by Adam Distler and Dr. Melinda Soares Furtado for the use of dating the Ursa Major moving group.
The pipeline has the following dependencies: numpy, math, astropy.io, matplotlib.pyplot, scipy.interpolate, specutils, barycorrpy, astroquery.simbad, PyAstronomy, specutils.manipulation, keplersplinev2, and warnings.
This code consists of two functions: abs_rv() and abs_rv_complete(), each with similar but distinct imputs. abs_rv() has the following parameters:
target - this is the target file, in a .fits format that can be read by SpectrumCollection.read() in Specutils.
