import numpy as np
import math
import astropy.io
from astropy.io import fits
import matplotlib.pyplot as plt
import scipy.interpolate as interp
from specutils import Spectrum1D
from barycorrpy import get_BC_vel , exposure_meter_BC_vel
from astropy.wcs import WCS
from specutils import Spectrum1D, SpectrumCollection, SpectralRegion, SpectrumList
from astroquery.simbad import Simbad
from astropy.modeling import models
from astropy import units as u
from PyAstronomy import pyasl
from astropy.convolution import Box1DKernel
from specutils.manipulation import convolution_smooth
from keplersplinev2 import *
import warnings
warnings.filterwarnings('ignore')
from specutils.fitting import fit_generic_continuum, fit_continuum


