# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module computes various properties that can be extracted from lammps
calculations.

TODO: finish implementation
"""

import numpy as np
import scipy.integrate as sp_integrate


__author__ = "Kiran Mathew"
__email__ = "kmathew@lbl.gov"


class GreenKuboProperties(object):
    """
    Green-Kubo formulas for the transport coefficient relates a transport
    coefficient to a time integral of a time correlation function.
    """
    def __init__(self, lammpsrun):
        """
        Args:
             lammpsrun (LammpsRun)
        """
        self.lammpsrun = lammpsrun

    def get_integrated_correlation(self, array):
        """
        Compute the normalized autocorrelation and integrate it wrt time.

        Args:
            array (numpy.ndarray): input numpy array

        Returns:
            integrated autocorrelation
        """
        auto_corr_full = np.correlate(array, array, mode="full")
        auto_corr = auto_corr_full[auto_corr_full.size / 2:] / np.max(auto_corr_full)
        time = self.lammpsrun.traj_timesteps
        return sp_integrate.simps(auto_corr, time)

    @property
    def diffusivity(self):
        """
        Green-Kubo diffusivity from the velocity auto-correlation.

        Returns:
            [Dxx_mol1, Dyy_mol1, Dzz_mol1, ...]

        TODO: fix the units
        """
        return [self.get_integrated_correlation(self.lammpsrun.mol_velocity[:, i, dim])
                for i in range(self.lammpsrun.nmols)
                for dim in range(3)]

    @property
    def charge_flux(self):
        """
        net molecular charge_flux for each timestep
        J = sum(v_mol * charge_mol)

        Returns:
            numpy.ndarray(n_timesteps x 3)
        """
        mol_velocity = self.lammpsrun.mol_velocity
        mol_charges = self.lammpsrun.mol_charges
        jx = np.dot(mol_velocity[:, :, 0], mol_charges[:])
        jy = np.dot(mol_velocity[:, :, 1], mol_charges[:])
        jz = np.dot(mol_velocity[:, :, 2], mol_charges[:])
        jx = jx.reshape(jx.shape + (-1,))
        jy = jy.reshape(jy.shape + (-1,))
        jz = jz.reshape(jz.shape + (-1,))
        return np.concatenate((jx, jy, jz), axis=1)

    @property
    def electrical_conductivity(self):
        """
        Electrical conductivity from the auto-correlation of the charge flux.

        TODO: fix the units
        """
        mol_current = self.charge_flux
        kappa = [self.get_integrated_correlation(mol_current[:, dim])
                 for dim in range(3)]
        return kappa

    def viscosity(self, skip):
        """
        Computes Green-Kubo viscosity from the transverse component of the
        stress-stress correlation

        TODO: Fix units
        """
        if not self.lammpsrun.lammpslog.get('pxy'):
            print("no pressure data")
            raise KeyError
        else:
            nu = [
                self.get_integrated_correlation(np.array(self.lammpsrun.lammpslog[comp][skip:]))
                for comp in ['pxy', 'pxz', 'pyz', 'pxx', 'pyy', 'pzz']]
            return nu

    @property
    def nernst_einstein_conductivity(self):
        pass
