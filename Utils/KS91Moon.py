#!/usr/bin/env python

"""
RS 2015/07/15:  Krisciunas & Schaefer Moon Model

This module extends pyephem's Moon class to include calculations of sky
brightness from scattered moonlight.  The relevant reference is

    Krisciunas, K. & Schaefer, B. (1991), PASP 103, 1033.
"""

import ephem
import numpy as np


class KS91Moon(ephem.Moon):

    def sky_brightness(self, az, alt, k=0.2):
        """
        Implements the model of Krisciunas & Schaefer 1991, PASP 103, 1033.
        Might consider including Hayes & Latham 1975, ApJ 197, 593 for a model
        of the extinction coefficient, though maybe best to use historical
        values for your observatory instead.

        Inputs
            az, alt:  sky azimuth and elevation in *radians* (or ephem.Angles)
            k:  extinction coefficient in desired band in magnitudes
        Outputs
            Bmoon:  moon brightness (flux) in nano-Lamberts
        """

        # Check to see whether we've computed the moon's position
        if not hasattr(self, 'alt'):
            raise RuntimeError('sky brightness undefined until first compute')

        # If the moon is below the horizon, sky brightness is zero
        if self.alt < 0:
            return 0

        # Magnitude of the Moon as a function of its phase as given by pyephem.
        # This then gives the intensity of the moon on a suitable scale.
        I_moon = 10**(-0.4*(self.mag + 16.57))

        # Find the angle between the Moon's location and where we're looking
        rho_rad = ephem.separation((self.az, self.alt), (az, alt))
        rho_deg = np.degrees(rho_rad)

        # Scattering function:  Rayleigh + Mie components
        frho_ray = 2.3e+5*(1.06 + np.cos(rho_rad)**2)
        if rho_deg < 10:
            frho_mie = 6.2e+7/rho_deg**2
        else:
            frho_mie = 10**(6.15-rho_deg/40)
        frho = frho_ray + frho_mie

        # Airmass for the sky location at the given time
        XZ = 1.0/np.sqrt(1-0.96*np.sin(ephem.pi/2-alt)**2)
        XZm = 1.0/np.sqrt(1-0.96*np.sin(ephem.pi/2-self.alt)**2)

        # Finally, sky brightness including all the above ingredients
        Bmoon = frho*I_moon*(10**(-0.4*k*XZm))*(1 - 10**(-0.4*k*XZ))
        return Bmoon


def nL2Vmag(B):
    """Converts nanoLamberts to V magnitude per square arcsec"""
    return (20.7233 - np.log(B/34.08))/0.92104
