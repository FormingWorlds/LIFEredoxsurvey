#!/usr/bin/env python
# coding: utf-8

import configparser
import lifesim
import warnings
import os as os
import numpy as np
from spectres import spectres
from astropy import units as u
import tkinter as tk
from tkinter import scrolledtext

def planet_blackbody(wavelength, widths, temperature, radius, distance):
    """
    Calculates the blackbody emission spectrum of a planet in the right units for LIFEsim. 
    Wavelength and widths must be given with units of meters.
    """
    spec = lifesim.util.radiation.black_body(mode = 'planet', bins=wavelength.value, width=widths.value, temp = temperature, radius = radius, distance = distance)
    spec = spec / u.m**2 / u.s ; spec /= widths
    planet_flux = [wavelength, spec]
    return planet_flux

def spectrum_photon(wavelength, flux, width):
    """
    Takes the wavelength range (in meters) and photon flux from the planet (in ph per meter^2 second)
    and outputs the planet flux in the format usably by LIFEsim
    """
    #get the flux to the right units for LIFEsim
    spec = flux / u.m**2 / u.s ; spec /= width
    planet_flux = [wavelength, spec] #[meters, ph3 m-2 s-1]
    return planet_flux

def spectrum_Jansky(wavelength,flux):
    """
    Takes the wavelength range and flux (in Jansky) and converts it to units usable by LIFEsim
    """
    angstrom = wavelength.to(u.angstrom)
    pflux = (flux * 1.51e3)/angstrom.value * 10e13 * u.ph / u.m**3 / u.s
    planet_flux = [wavelength,pflux] #meter, ph3 m-2 s-1
    return planet_flux

def wavelength_converter(wavelength, unit_wavelength):
    """
    Takes given wavelength data, checks that the given unit is valid and converts it to meters for LIFEsim 
    if necessary
    """
    if unit_wavelength == 'm':
        lambdaa = wavelength * u.m
    elif unit_wavelength == 'micron':
        lambdaa = wavelength * u.micron
        lambdaa = lambdaa.to(u.m)
    else:
        raise ValueError("Invalid wavelength unit for LIFEsim, check 'init_simulation.cfg' for accepted units [in brackets]")
    return lambdaa

def display_text(text):
    root = tk.Tk()
    root.title("Text Output")

    text_area = scrolledtext.ScrolledText(root, wrap=tk.WORD, width=80, height=30, font=("Ariel", 15))
    text_area.pack(padx=15, pady=15)

    text_area.insert(tk.INSERT, text)

    root.mainloop()