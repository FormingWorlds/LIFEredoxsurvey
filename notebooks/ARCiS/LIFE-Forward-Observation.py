#!/usr/bin/env python
# coding: utf-8

# In[4]:


import warnings
warnings.simplefilter('ignore')

import numpy as np
import os
import matplotlib.pyplot as plt
from astropy import units as u
from matplotlib.pyplot import figure, show
import pandas
from spectres import spectres
import subprocess

import lifesim

#makes sure the working directory is unchanged by packages importation, if it is it fixes it
print(os.getcwd())
os.chdir('/Users/lorenzo/LIFE-ARCiS-Pipeline')
print(os.getcwd())

filename = input("Flux directory: ")
integration_time = float(input("Observation time (hours): "))
distance = float(input("Distance star-planet system (pc): "))
typ1 = input("Star Type (S, L or T): ")

subprocess.run(f"scp cesario@kapteyn.astro.rug.nl:/Users/users/cesario/BULK/Retrievals-on-forward/Forward-Files/{filename}/emis /Users/lorenzo/LIFE-ARCiS-Pipeline/MyFlux/{filename}.txt", shell=True)

emis = np.loadtxt(f"/Users/lorenzo/LIFE-ARCiS-Pipeline/MyFlux/{filename}.txt",comments='#')
lambdaa = emis[:,0] * u.micron
wavelength = lambdaa.to(u.m)
jansky_flux = emis[:,1] * u.Jansky

#wavelength in angstrom for the unit conversion
lambdaa = lambdaa.to(u.angstrom)

photon_flux = (jansky_flux.value * 1.51e3)/lambdaa.value * 10e13 * u.ph / u.m**3 / u.s

#integration time
intime = integration_time #hours
intime *= 3600

#star observing parameters
dist = distance #parsec

#STAR parameters

#TRAPPIST-1=0, AD LEONIS=1, SUN=2
if typ1 == 'T':
    star = 0
elif typ1 == 'L':
    star = 1
elif typ1 == 'S':
    star = 2

stemp = [2566,3390,5780][star] #Kelvin
sradius = [0.12,0.4,1][star] #Solar radii
slat = 0.78 #ecliptic latitude in radians 45deg = 0.78; 90deg = 1.57
orbit = [0.02,0.13,1][star] #AU
zodi = 3 #zodis, examples: 1,10,100,1000

#Instrument set-up, choose observing scenario
bus = lifesim.Bus()
#observing scenarios: optimistic, baseline, pessimistic
bus.data.options.set_scenario('baseline')
#spectral resolution
specres = 100
bus.data.options.set_manual(spec_res = specres)
#wavelength coverage
bus.data.options.set_manual(wl_min = 6)
bus.data.options.set_manual(wl_max = 16)
baseline_planet = True #optimizes baseline for planet

instrument = lifesim.Instrument(name='inst')
bus.add_module(instrument)
transm = lifesim.TransmissionMap(name='transm')
bus.add_module(transm)
#include the noise sources
exozodi = lifesim.PhotonNoiseExozodi(name='exo')
bus.add_module(exozodi)
localzodi = lifesim.PhotonNoiseLocalzodi(name='local')
bus.add_module(localzodi)
star_leak = lifesim.PhotonNoiseStar(name='star')
bus.add_module(star_leak)
#connect the modules
bus.connect(('inst', 'transm'))
bus.connect(('inst', 'exo'))
bus.connect(('inst', 'local'))
bus.connect(('inst', 'star'))
bus.connect(('star', 'transm'))

planet_flux = [wavelength, photon_flux]

SNR = 0
n = orbit * u.AU; n = n.to(u.pc) #radius of orbit
ang = (n.value/dist) * (360/(2*np.pi)) * 3600

while SNR<10:
    bins_snr, pflux, noyz = instrument.get_spectrum(temp_s=stemp, radius_s=sradius, lat_s=slat, z=zodi, distance_s=dist, angsep=ang, flux_planet_spectrum=planet_flux, integration_time=intime)
    SNR = np.sqrt(np.sum(bins_snr[1]**2))
    print(f"S/N: {SNR} \n Integration time: {intime/3600:.1f} hours \n")
    intime += 3600

print(f"S/N of observation: ",SNR,f" Integration time: {intime/3600:.1f} hours")

planet_flux_r = spectres(new_wavs=instrument.data.inst['wl_bin_edges'],spec_wavs=planet_flux[0].value,spec_fluxes=planet_flux[1].value,edge_mode=True)

#noise distribution
#calculate simulated noise
ratio = planet_flux_r / bins_snr[1]
noise = np.random.normal(0,ratio,size=pflux.shape)

bins = bins_snr[0] * u.m
bins = bins.to(u.micron)

width = (bins.value[-1] - bins.value[0])/ specres* u.micron
resolution = bins.value/width.value

#first we have to convert the emission spectra back to Jansky
preconv = bins_snr[0] * u.meter
conversion_lambda = preconv.to(u.angstrom)

mask = (preconv.to(u.micron).value > 6) & (preconv.to(u.micron).value < 16)

dat_emission = ((planet_flux_r * conversion_lambda.value)/(1.5e3 * 10e13))
dat_emission = dat_emission[mask]
dat_wavelength = preconv.to(u.micron)[mask]
dat_error = (noise * conversion_lambda.value)/(1.5e3 * 10e13)
dat_error = dat_error[mask]
dat_specres = resolution[mask]

data = np.column_stack((dat_wavelength.value, dat_emission, dat_error, dat_specres))
np.savetxt(f"ARCiS_fluxes/{filename}", data, fmt='%.6e', delimiter=' ')

local_file = f"ARCiS_fluxes/{filename}"
remote_path = "cesario@kapteyn.astro.rug.nl:/Users/users/cesario/BULK/Retrievals-on-forward/Fluxes"

subprocess.run(f"scp {local_file} {remote_path}", shell=True)

# In[ ]:




