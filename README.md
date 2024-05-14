# TRES_RV_Pipeline
A pipeline for determining radial velocities of targets using either synthetic spectra or actual spectra from stars with known RV's. This was created by Adam Distler and Dr. Melinda Soares Furtado for the use of dating the Ursa Major moving group.
The pipeline has the following dependencies: numpy, math, astropy.io, matplotlib.pyplot, scipy.interpolate, specutils, barycorrpy, astroquery.simbad, PyAstronomy, specutils.manipulation, keplersplinev2, and warnings.
This code consists of two functions: abs_rv() and abs_rv_complete(), each with similar but distinct imputs. abs_rv() has the following parameters:
target - this is the target file, in a .fits format that can be read by SpectrumCollection.read() in Specutils. You will want to enter the path + file name in a string format, i.e target_file = 'my_file.fits'
template - the template file. Our team personally uses Pollux for syntehtic spectra, as they are easier to use than another star whose RV is known. If one uses a synthetic spectra, then it should be in a .fits format in a simple 2D array. For an actual star, we require that the star is in a .fits format with each axis is a specific order.

order - this is the parameter that one can change to change the specific order under consideration for the pipeline. We have found that different orders yield better or worse results, as different stars have more prominent features at different parts of the EM spectrum.
