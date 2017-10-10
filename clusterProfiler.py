import tangos as db
import numpy as np
from scipy.interpolate import interp1d

class ClusterProfiler(object):

	def __init__(self, step):
		"""
		Assume that halo 1 in this step is the main cluster.  Then, create an interpolation scheme.

		Not doing a 2d interpolation scheme because the bins are not the same from time step to time step.
		"""

		#Obtain data
		clusterHalo = step.halos[0]
		redshifts, profiles, Rvir = clusterHalo.reverse_property_cascade('z()', 'gas_density_profile', 'Rvir')
		logInterpolationFunctions = []

		#At each time step, compute the profile.
		for z_index in range(len(redshifts)):
			xvalues = np.arange(len(profiles[z_index]))*0.1 / Rvir[z_index]
			logInterpolationFunctions.append(interp1d(np.log10(xvalues), np.log10(profiles[z_index]), bounds_error=False, fill_value=-np.inf))

		self.redshifts = redshifts
		self.logInterpolationFunctions = logInterpolationFunctions

	def _selectNearestFunction(self, redshift):
		"""
		Given a redshift, find the closest one for which we have data and return the function.
		"""

		return self.logInterpolationFunctions[np.argmin(np.abs(self.redshifts-redshift))]

	def computeGasDensity(self, distanceInRvir, redshift):
		"""
		Use interpolation functions to solve for gas density as a function of distance and redshift.
		"""

		#Profile data do not go forever.
		if redshift > self.redshifts[-1]:
			return 0

		#First, find the nearest redshift and get its interpolation function.
		interpFunct = self._selectNearestFunction(redshift)

		#Next, interpolate in log space.
		return 10**interpFunct(np.log10(distanceInRvir))
