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
		times, profiles, Rvir = clusterHalo.reverse_property_cascade('t()', 'gas_density_profile', 'Rvir')
		logInterpolationFunctions = []

		#At each time step, compute the profile.
		for t_index in range(len(times)):
			xvalues = np.arange(len(profiles[t_index]))*0.1 / Rvir[t_index]

			#Masking zeroes for better interpolation behavior.
			isZero = profiles[t_index] == 0
	
			#Also, let's get rid of kooky behavior within 0.05 Rvir
			sufficientlyFar = xvalues > 0.05

			logInterpolationFunctions.append(interp1d(np.log10(xvalues[(~isZero) & (sufficientlyFar)]), \
			np.log10(profiles[t_index][(~isZero) & (sufficientlyFar)]), bounds_error=False, \
			fill_value=(np.log10(profiles[t_index][~isZero][0]),-np.inf)))

		self.times = times
		self.logInterpolationFunctions = logInterpolationFunctions

	def _selectNearestFunction(self, time):
		"""
		Given a time, find the closest one for which we have data and return the function.
		"""

		return self.logInterpolationFunctions[np.argmin(np.abs(self.times-time))]

	def computeGasDensity(self, distanceInRvir, timeArr):
		"""
		Use interpolation functions to solve for gas density as a function of distance and time.
		"""

		if not hasattr(distanceInRvir, '__len__'):
			distanceInRvir = np.array([distanceInRvir])
		if not hasattr(timeArr, '__len__'):
			timeArr = np.full(len(distanceInRvir), timeArr)
		output = np.zeros(len(distanceInRvir))

		#First, find the nearest time and get its interpolation function.
		for d_index in range(len(distanceInRvir)):
			if timeArr[d_index] < self.times[-1]:
				output[d_index] = 0
			else:
				interpFunct = self._selectNearestFunction(timeArr[d_index])
				output[d_index] = 10**interpFunct(np.log10(distanceInRvir[d_index]))

		return output
