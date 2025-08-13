#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""

import numpy as np
from scipy import constants


def H_to_eV(energy):
    # Conversion from Hartree to eV
    return energy/constants.physical_constants['electron volt-hartree relationship'][0]


def eV_to_H(energy):
    # Conversion from eV to Hartree
    return energy*constants.physical_constants['electron volt-hartree relationship'][0]


def H_to_kjmol(energy):
    # Conversion from Hartree to kJ / mol
    return energy*(constants.physical_constants['Hartree energy'][0] * 1e-3 * constants.Avogadro)


def kjmol_to_H(energy):
    # Conversion from kJ / mol to Hartree
    return energy/(constants.physical_constants['Hartree energy'][0] * 1e-3 * constants.Avogadro)


def au_to_angstrom(length):
    return length*(constants.physical_constants['atomic unit of length'][0] * 1e10)


def angstrom_to_au(length):
    return length/(constants.physical_constants['atomic unit of length'][0] * 1e10)


def cm_to_thz(freq):
    # Conversion from cm^-1 to THz
    return freq*(constants.physical_constants['speed of light in vacuum'][0] * 1e-10)


def thz_to_cm(freq):
    # Conversion from THz to cm^-1
    return freq/(constants.physical_constants['speed of light in vacuum'][0] * 1e-10)


def hartree_to_thz(freq):
    # Conversion from frequency in AU(angular) to THz(linear)
    return freq / (1e12 * constants.physical_constants['atomic mass unit-hartree relationship'][0] / constants.physical_constants['atomic mass unit-hertz relationship'][0])


def thz_to_hartree(freq):
    # Conversion from frequency in THz(linear) to AU(angular)
    return freq * (1e12 * constants.physical_constants['atomic mass unit-hartree relationship'][0] / constants.physical_constants['atomic mass unit-hertz relationship'][0])


def amu_to_me(mass):
    # Conversion from unified atomic mass unit to electron mass (mass unit in AU)
    return mass * (constants.physical_constants['unified atomic mass unit'][0] / constants.m_e)


def me_to_amu(mass):
    # Conversion from electron mass (mass unit in AU) to unified atomic mass unit
    return mass / (constants.physical_constants['unified atomic mass unit'][0] / constants.m_e)
