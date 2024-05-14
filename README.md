# TRES_RV_Pipeline
A pipeline for determining radial velocities of targets using either synthetic spectra or actual spectra from stars with known RV's. This was created by Adam Distler and Dr. Melinda Soares Furtado for the use of dating the Ursa Major moving group.
The pipeline has the following dependencies: numpy, math, astropy.io, matplotlib.pyplot, scipy.interpolate, specutils, barycorrpy, astroquery.simbad, PyAstronomy, specutils.manipulation, keplersplinev2, and warnings.
This code consists of two functions: abs_rv() and abs_rv_complete(), each with similar but distinct imputs. abs_rv() and abs_rv_complete() have similar parameters, except those marked with a *. They have the following parameters:

target - this is the target file, in a .fits format that can be read by SpectrumCollection.read() in Specutils. You will want to enter the path + file name in a string format, i.e target_file = 'my_file.fits'

template - the template file. Our team personally uses Pollux for syntehtic spectra, as they are easier to use than another star whose RV is known. If one uses a synthetic spectra, then it should be in a .fits format in a simple 2D array. For an actual star, we require that the star is in a .fits format with each axis is a specific order.

* order - for only abs_rv(). This is the parameter that one can change to change the specific order under consideration for the pipeline. We have found that different orders yield better or worse results, as different stars have more prominent features at different parts of the EM spectrum.

* order_list - for abs_rv_complete(). This will be a list of order numbers that one wants to iterate over. We recommend using order_list = np.linspace(a,b,n). It is important to note that these should all be integers.

trim - the trim feature allows one to slice off the ends of each order to remove the drop-off in quantum efficiency. A trim of 0.1 would slice off the first and last 10% of the spectrum.

synthetic_template - either a 0 or 1, 1 indicates a synthetic spectrum. A 0 represents a real star.

instrument - the default is 'TRES', but more instruments that the team uses will be added as needed.

minrv, maxrv - the min and max rv's one will search. Input a simple number, i.e minrv=-50 

stepsize - the stepsize between RV's for the cross-correlation function.

smoothing - a binary input, as 0 indicates no smoothing, while a 1 indicates a smoothing of the synthetic spectrum. This has been found to make a sizeable difference, with a trend we have found is to match the depths of the features gives us a good result.

kw - a smoothing factor. 0 indicates no smoothing, while a higher number (i.e 20, 30) will show significant smoothing of the template.

fit - either 'spline', or 'gaussian'. this is used to counter the blaze corrections needed for the target.

printing - either a 0 or 1. 1 indicates the graphs and resulting text will be printed when calling the function. For calling abs_rv_complete(), there would be a lot of graphs printed, so we recommend setting this to 0 for this function.
