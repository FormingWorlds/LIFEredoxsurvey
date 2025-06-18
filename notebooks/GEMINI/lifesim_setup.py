#!/usr/bin/env python
# coding: utf-8


import configparser
import lifesim
import warnings
import numpy as np
from spectres import spectres
from astropy import units as u

warnings.simplefilter('ignore')

from path import PATH

config = configparser.ConfigParser(inline_comment_prefixes="#")
config.read(f"{PATH}/init_simulation.cfg")

wavel = np.linspace(3.000e-6,20e-6,3000) * u.m
widths = 0.02e-6 * u.m

#setting up the LIFE Instrument
bus = lifesim.Bus()
scenario = config['instrument']['instrument_scenario']
bus.data.options.set_scenario(scenario)

if config['instrument']['instrument_scenario-override'] == 'True':
    bus.data.options.set_manual(diameter = float(config['instrument']['instrument_diameter']))
    bus.data.options.set_manual(wl_min = float(config['instrument']['instrument_wl-min']))
    bus.data.options.set_manual(wl_max = float(config['instrument']['instrument_wl-max']))

bus.data.options.set_manual(throughput = float(config['instrument']['instrument_photon-throughput']))
bus.data.options.set_manual(quantum_eff = float(config['instrument']['instrument_quantum-efficiency']))
bus.data.options.set_manual(spec_res = float(config['instrument']['instrument_spectral-resolution']))
bus.data.options.set_manual(bl_min = float(config['instrument']['instrument_minimum-baseline']))
bus.data.options.set_manual(bl_max = float(config['instrument']['instrument_maximum-baseline']))
bus.data.options.set_manual(ratio = float(config['instrument']['instrument_ratio-baseline']))
bus.data.options.set_manual(t_slew = float(config['instrument']['instrument_slewing-time']))
bus.data.options.set_manual(t_efficiency = float(config['instrument']['instrument_t-efficiency']))
bus.data.options.set_manual(image_size = int(float(config['instrument']['instrument_image-size'])))
bus.data.options.set_manual(wl_optimal = float(config['instrument']['instrument_optimal-wl']))

#set-up the modules
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

def signal_to_noise(planet_flux, star_temperature, star_radius, star_latitude, exozodi, distance, orbit, integration_time):
    """
    Function to calculate the signal-to-noise ratio of a target planet after simulating 
    an observation with LIFEsim. Orbit = radius of planet's orbit in pc.
    Returns the integrated signal to noise ratio (S/N) of the simulation, [wavelength, S/n per bin],
    planet's flux and the noise.
    """
    ang = (orbit.value/distance) * (360/(2*np.pi)) * 3600 #converts orbit to angular separation in arcsec
    bins_snr, pflux, noise = instrument.get_spectrum(temp_s=star_temperature, radius_s=star_radius, lat_s=star_latitude, z=exozodi, distance_s=distance, angsep=ang, flux_planet_spectrum=planet_flux, integration_time=integration_time)
    SN = np.sqrt(np.sum(bins_snr[1]**2))
    return SN, bins_snr, pflux, noise

