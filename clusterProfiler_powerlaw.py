import tangos as db
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import linregress

class ClusterProfiler(object):

	def __init__(self, step):
		"""
		Assume that halo 1 in this step is the main cluster.  Then, create an interpolation scheme.

		Not doing a 2d interpolation scheme because the bins are not the same from time step to time step.
		"""

		#Obtain data
		clusterHalo = step.halos[0]
		times, profiles, Rvir = clusterHalo.reverse_property_cascade('t()', 'gas_density_profile', 'Rvir')
		lineSlopes = []
		lineIntercepts = []

		#At each time step, compute the profile.
		for t_index in range(len(times)):
			xvalues = np.arange(len(profiles[t_index]))*0.1 / Rvir[t_index]

			#Masking zeroes for better interpolation behavior.
			isZero = profiles[t_index] == 0
	
			#Also, let's get rid of kooky behavior within 0.05 Rvir
			sufficientlyFar = xvalues > 0.05

			#Fit a power law (a line in log space)
			lineParameters = linregress(np.log10(xvalues[(~isZero) & (sufficientlyFar)]), np.log10(profiles[t_index][(~isZero) & (sufficientlyFar)]))
			lineSlopes.append(lineParameters[0])
			lineIntercepts.append(lineParameters[1])

		self.times = times
		self.lineSlopes = np.array(lineSlopes)
		self.lineIntercepts = np.array(lineIntercepts)
		self._generateInterpFunctions()

	def _generateInterpFunctions(self):
		"""
		Called by init.  Create interpolation formulas for line slopes and intercepts as a function of time.
		"""

		self.slopeOfT = interp1d(self.times, self.lineSlopes, bounds_error=False, fill_value=0)
		self.interceptOfT = interp1d(self.times, self.lineIntercepts, bounds_error=False, fill_value=0)

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
		output = 10**(self.slopeOfT(timeArr)*np.log10(distanceInRvir) + self.interceptOfT(timeArr))

		return output
