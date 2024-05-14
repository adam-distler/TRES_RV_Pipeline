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

def abs_rv(target,template,order=32,trim=0.15,synth_template=0,instrument='TRES',
           minrv=-30.,maxrv=100,stepsize=0.05,smoothing=0,kw=3,fit='spline', printing=1):
    plt.rcParams.update({'font.size': 16, 'axes.labelsize': 16, 'axes.titlesize': 20})
    target= fits.open(target)
    template = fits.open(template)
    vc=299792.458# speed of light in km/s
    # Target --- read in our target file
    
    if instrument=='TRES':
        multispec = SpectrumCollection.read(target)
        #trim data to cut out some of the edges (set by "trim" value)
        cut_size = int(trim * len(multispec[order].wavelength))    
        trimmed_array = multispec[order][cut_size:-cut_size]
        #perform barycentric correction using header information
        target_head=target[0].header
        target_bcv=float(target_head['BCV'])
        if printing ==1:
            print('Performing a target barycentric correction of',target_bcv,'km/s')
        modified_wavelength=trimmed_array.wavelength*(1+((target_bcv)/vc)) #barycenter 
        target_finalwl = modified_wavelength
        target_finaladu=trimmed_array.flux
        #normalize spectra
        if printing ==1:
            print('You are normalizing your spectra in',fit,'mode.')
        if fit=='spline':
            s = keplersplinev2(trimmed_array.wavelength.value, trimmed_array.flux.value,bkspace = 3)
            target_finaladu=trimmed_array.flux.value/s
        if fit=='gaussian':
            model = models.Gaussian1D(max(trimmed_array.flux.value), np.nanmean(trimmed_array.wavelength).value, 30)
            g1_fit = fit_continuum(trimmed_array, model=model)
            g1_blaze = g1_fit(trimmed_array.wavelength)
            target_finaladu=trimmed_array.flux/g1_blaze
    minwl=min(target_finalwl.value)
    maxwl=max(target_finalwl.value)
    
    if synth_template==1:
        if printing ==1:
            print('You are incorporating a template file using synthetic data')
        template[1].data=template[1].data[(template[1].data.wavelength>minwl)&(template[1].data.wavelength<maxwl)]###leave off here...
        template[1].data.flux=template[1].data.flux/np.nanmedian(template[1].data.flux)
        template_finalwl=template[1].data.wavelength*u.AA
        template_finalflux=template[1].data.flux
        if fit=='spline':
            if printing ==1:
                print('Fitting a spline to your template too')
            s = keplersplinev2(template_finalwl.value, template_finalflux,bkspace = 3)
            template_finalflux=template_finalflux/s
        
        if smoothing==1:
            if printing ==1:
                print('Smoothing the template spectra with kw=',kw)
            box1d_kernel = Box1DKernel(width=kw,mode='linear_interp')
            spec1 = Spectrum1D(spectral_axis=template[1].data.wavelength* u.AA,flux=template[1].data.flux * u.adu)
            smoothfile=convolution_smooth(spec1, box1d_kernel) 
            template_finalwl=smoothfile.spectral_axis
            template_finalflux=smoothfile.flux.value
        
    
    if synth_template==0:
        print('You are assuming a template file using real data')
        multispec = SpectrumCollection.read(template)
        cut_size = int(trim * len(multispec[order].wavelength))
        #Array trimmed on each side 
        trimmed_array = multispec[order][cut_size:-cut_size]
        
        #fitting the blaze function
        model = models.Gaussian1D(max(trimmed_array.flux.value), np.nanmean(trimmed_array.wavelength).value, 30)
        g1_fit = fit_continuum(trimmed_array, model=model)
        g1_blaze = g1_fit(trimmed_array.wavelength)
        
        if synth_template==0:
            print('Next, we need to determine RV and barycentric correction')
            if instrument=='TRES':
                template_head=template[0].header
                template_bcv=float(template_head['BCV'])
                if printing ==1:
                    print('BCV=',template_bcv,' km/s')
                modified_wavelength=trimmed_array.wavelength*(1+((template_bcv)/vc)) #barycenter 
                customSimbad = Simbad()
                customSimbad.add_votable_fields('mk', 'rot', 'rv_value')
                result_table = customSimbad.query_object(template_file[0].header['OBJECT'])

                template_rv=result_table['RV_VALUE'][0]
                print('template RV = ',template_rv,' km/s')
                template_finalwl=modified_wavelength*(1-((template_rv)/vc)) 
                template_finalwl
               
                template_finalflux=trimmed_array.flux/g1_blaze
                #final arrays: template_finalwl,template_finalflux
    if synth_template==1:
        if printing ==1:
            print('Great, since you are using a synthetic template, there is need to perform RV or barycentric corrections!')
    if (instrument!='TRES'):
        abs_rv=np.nan
        print('WARNING --- the instrument ' +str(instrument)+' is not accounted for.')
    if printing ==1:
        fig=plt.figure(figsize=(15,8))
        plt.plot(template_finalwl,template_finalflux,label='template')
        plt.plot(target_finalwl,target_finaladu,label='target')
        plt.legend()
        plt.show()
    
    # Carry out the cross-correlation.
    # The RV-range is -30 - +100 km/s in steps of 0.6 km/s.
    # The first and last 50 points of the data are skipped.
    rv, cc = pyasl.crosscorrRV(target_finalwl.value, target_finaladu, 
                               template_finalwl.value, template_finalflux, 
                               minrv, maxrv, stepsize, skipedge=200)    
    if printing ==1:
        plt.plot(rv, cc/max(cc), 'b-', label="normalized unweighted CCF")
        #plt.xlim(-10,50)
        plt.grid(alpha=0.25)
        plt.title('Maximum CC value:'+str(np.argmax(cc))+', order='+str(order))
        plt.show()
    
    # Find the index of maximum cross-correlation function
    maxind = np.argmax(cc)
    if printing ==1:
        print('Maximum CC value:', np.argmax(cc))
        
    #print("Cross-correlation function is maximized at dRV = ", rv[maxind], " km/s")
    if printing ==1:
        if rv[maxind] > 0.0:
            print(" A red-shift with respect to the template")
        else:
            print(" A blue-shift with respect to the template")
        fig=plt.figure(figsize=(16,8))
        plt.plot(template_finalwl,template_finalflux,label='template',zorder=1)
        tmp_wavelength=target_finalwl*(1-((rv[maxind])/vc)) #barycenter 
        plt.plot(tmp_wavelength,target_finaladu,label='target',zorder=0)
        plt.legend()
        plt.grid(alpha=0.25)
        plt.show()
    
    abs_rv=rv[maxind]
    print('RV =',abs_rv,'km/s')
    template.close()
    target.close()
    return abs_rv


def abs_rv_complete(target,template,order_list=np.linspace(10,30,21),trim=0.15,synth_template=0,instrument='TRES',
           minrv=-30.,maxrv=100,stepsize=0.05,smoothing=0,kw=3,fit='spline', printing=1):
    RV_list = np.array([])
    for ORD in order_list:
        ORD =int(ORD)
        rv = abs_rv(target,template,ORD,trim,synth_template,instrument,
           minrv,maxrv,stepsize,smoothing,kw,fit, printing)
        RV_list = np.append(RV_list, rv)
    plt.plot(order_list, RV_list, color='black')
    plt.title('Radial Velocity vs. Order')
    plt.xlabel('Order')
    plt.ylabel('Radial Velocity (km/s)')
    median_rv = np.median(RV_list)
    plt.show()
    print('The median RV for this target is:', median_rv, 'km/s')
