# lifesim bug that changes working directory
import os
currentpath = os.getcwd()
import lifesim
os.chdir(currentpath)

# other imports
import configparser
import warnings
import numpy as np
from spectres import spectres
from astropy import units as u
import tkinter as tk
import matplotlib as mpl
#back-end for matplotlib to be able to output figure and a text output at the same time  
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show
import scienceplots
import scipy
from tkinter import scrolledtext
import os
import csv

from path import PATH
from functions import planet_blackbody, spectrum_photon, spectrum_Jansky, wavelength_converter, display_text
from lifesim_setup import signal_to_noise, wavel, widths, scenario

warnings.simplefilter('ignore')

config = configparser.ConfigParser(inline_comment_prefixes="#")
config.read(f"{PATH}/init_simulation.cfg")

#read-out key values from configuration file

# read PLANET and SIMULATION parameters

#distance to target system
distance_star = float(config['simulation']['simulation_distance-to-star'])
distance_planet = float(config['simulation']['simulation_distance-to-planet'])

#radius of planet's orbit
orbital_radius = float(config['planet']['planet_orbit'])
orbit = orbital_radius * u.AU; orbit = orbit.to(u.pc)

# truth value of baseline optimization based on target
if config['simulation']['simulation_baseline-optimization'] == 'True':
    baseline_planet = True
else:
    baseline_planet = False
    
# integration time
intime = float(config['simulation']['simulation_integration-time'])

# dust level
zodi = float(config['simulation']['simulation_zodi'])

#read STAR parameters
#radius
sradius = float(config['star']['star_radius']) #solar radii
#effective temperature
stemp = float(config['star']['star_temperature']) #kelvin
#ecliptic latitude
slat = float(config['star']['star_ecliptic-latitude'])


# check if user input a spectrum path in config file and load spectrum
if config['planet']['planet_insert-spectrum'] == 'True':
    manual = 'yes'
    #check unit of spectrum
    spec = config['planet']['planet_unit-emission']; unit_wavelength = config['planet']['planet_unit-wavelength']
    filepath = config['planet']['planet_pathtospectrum']; emission = np.loadtxt(filepath,comments='#')
    
    lambdaa = emission[:,0] #wavelength array
    photon_flux = emission[:,1] #irradiance array
    
    #check that given units for wavelength are correct and convert to meters
    lambdaa = wavelength_converter(lambdaa, unit_wavelength)
        
    if spec == 'ph/m2/s':
        planet_flux = spectrum_photon(lambdaa, photon_flux, widths)
    elif spec == 'Jy':
        planet_flux = spectrum_Jansky(lambdaa, photon_flux)
    else:
        raise ValueError("Invalid spectral irradiance unit for LIFEsim, check 'init_simulation.cfg' for accepted units [in brackets]")

# Planet spectra by blackbody function
elif config['planet']['planet_insert-spectrum'] != 'True':
    manual = 'no'
    pradius = float(config['planet']['planet_radius']); ptemp = float(config['planet']['planet_temperature'])
    #LIFEsim Blackbody function
    planet_flux = planet_blackbody(wavel,widths,ptemp,pradius,distance_planet)

#main function to calculate the signal-to-noise ratio of the simulation
SN, bins_snr, pflux, noise = signal_to_noise(planet_flux,stemp,sradius,slat,zodi,distance_star,orbit,intime)
    
if manual == 'yes':
    text = f"Observation of manually inserted spectrum. \n \n Filepath: {filepath} \n \n Wavelength Units: {unit_wavelength}. Instrument architecture: {scenario} \n \n Planet was located at a distance of {distance_planet} parsecs from the interferometer, at an orbit of {orbital_radius} AU from its host star. \n \n The host star has the following characteristics: radius = {sradius} Solar radii; temperature = {stemp} K \n \n The integration time of the observation is {intime} seconds, {intime/3600} hours. \n \n Signal-to-noise ratio of the observation: {SN}"
else:
    text = f"Observation of blackbody spectrum of a planet. \n \n Instrument architecture: {scenario} \n \n Planet's characteristics: radius = {pradius} Earth radii; temperature = {ptemp} K. \n \n Planet was located at a distance of {distance_planet} parsecs from the interferometer, at an orbit of {orbital_radius} AU from its host star. \n \n The host star has the following characteristics: radius = {sradius} Solar radii; temperature = {stemp} K \n \n The integration time of the observation is {intime} seconds, {intime/3600} hours.\n \n Signal-to-noise ratio of the observation: {SN}"


if config['plotting']['plotting_snr'] == 'True':
    fig, ax1 = plt.subplots(figsize=(8,6))
    ax2 = ax1.twinx()

    ax1.plot(bins_snr[0]*1e6,pflux,color='red')
    ax2.step(bins_snr[0]*1e6,bins_snr[1], color='black')

    ax1.set_ylabel(r"Photon flux [ph $m^{-3} s^{-1}$]",color='red')
    ax1.set_xlabel(r"$\lambda$ [$\mu$m]")
    ax2.set_ylabel(r"S/N")

    ax1.tick_params(axis='y', colors='red')

    ax1.set_title(r"Planet Spectrum for LIFEsim and Signal-to-noise ratio per $\lambda$ bin")

    ax1.grid()

    plt.savefig(f"{PATH}/outputdir/output-figure.pdf", format="pdf")
    show(block=False)

#display_text(text)

output_directory = f"{PATH}/outputdir"
filename = 'output.txt'
file_path = os.path.join(output_directory, filename)
with open(file_path, 'w') as file:
        file.write(text)

# save data to txt file
output_directory = f"{PATH}/outputdir"
filename = 'data.txt'
file_path = os.path.join(output_directory, filename)
with open(file_path, 'w') as file:
        writer = csv.writer(file)
        writer.writerow(['bins_snr_0','bins_snr_1','pflux'])
        for i in range(len(bins_snr[0])):
            string = [str(bins_snr[0][i]*1e6),str(bins_snr[1][i]),str(pflux[i])]
            writer.writerow(string)