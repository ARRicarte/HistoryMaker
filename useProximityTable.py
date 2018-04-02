import cPickle as pickle
import numpy as np

class ProximityCalculator(object):

	def __init__(self, proximityPickle, mode='threshold', massType='Mstar', ratioThreshold=0.25):

		#Unpack proximity data
		with open(proximityPickle, 'r') as myfile:
			table = pickle.load(myfile)

		self.haloNumber = table['haloNumber']
		self.redshift = table['redshift']
		self.time = table['time']
		self.mass = table[massType]
		self.Rvir = table['Rvir']
		self.distanceMatrix = table['distanceMatrix']
		
		#Save kwargs
		self.mode = mode
		self.ratioThreshold = ratioThreshold

	def retraceProximity(self, haloNumbers, times):

		assert len(haloNumbers) == len(times)

		output = np.zeros(len(haloNumbers))
		for i in range(len(haloNumbers)):
			t_index = np.argmin(np.abs(self.time - times[i]))
			haloMatch = self.haloNumber[t_index] == haloNumbers[i]
			if any(haloMatch):
				if self.mode == 'threshold':
					relevanceMask = (self.haloNumber[t_index] != 1) & (self.haloNumber[t_index] != haloNumbers[i]) & \
					(self.mass[t_index]/self.mass[t_index][haloMatch] >= self.ratioThreshold)
					if any(relevanceMask):
						output[i] = np.min(self.distanceMatrix[t_index][haloMatch,relevanceMask])
					else:
						output[i] = np.inf
				elif self.mode == 'tidal':
					relevanceMask = (self.haloNumber[t_index] != 1) & (self.haloNumber[t_index] != haloNumbers[i])
					if any(relevanceMask):
						tidalScalar = self.Rvir[t_index][haloMatch] / self.distanceMatrix[t_index][haloMatch,relevanceMask] * \
						(self.mass[t_index][relevanceMask] / self.mass[t_index][haloMatch])**(1.0/3.0)
						output[i] = np.max(tidalScalar)
					else:
						output[i] = 0
			else:
				output[i] = np.nan
		return output

	def retraceProximityBetween(self, haloNumbers1, haloNumbers2, times1, times2):

		assert len(haloNumbers1) == len(times1)
		assert len(haloNumbers2) == len(times2)

		#Note:  Flipping upside-down because normal ordering is backwards.
		usedtimes = np.flipud(np.intersect1d(times1, times2))
		usedNumbers1 = haloNumbers1[np.in1d(times1, times2)]
		usedNumbers2 = haloNumbers2[np.in1d(times2, times1)]

		output = np.zeros(len(usedtimes))
		for i in range(len(usedtimes)):
			t_index = np.argmin(np.abs(self.time - usedtimes[i]))
			haloMatch1 = self.haloNumber[t_index] == usedNumbers1[i]
			haloMatch2 = self.haloNumber[t_index] == usedNumbers2[i]
			if (any(haloMatch1) & any(haloMatch2)):
				if self.mode == 'threshold':
					output[i] = self.distanceMatrix[t_index][haloMatch1,haloMatch2]
				elif self.mode == 'tidal':
					output[i] = self.Rvir[t_index][haloMatch1] / self.distanceMatrix[t_index][haloMatch1,haloMatch2] * \
					(self.mass[t_index][haloMatch2] / self.mass[t_index][haloMatch1])**(1.0/3.0)
			else:
				output[i] = np.nan
		return usedtimes, output

	def retraceClusterDistance(self, haloNumbers, times):

                assert len(haloNumbers) == len(times)

                output = np.zeros(len(haloNumbers))
                for i in range(len(haloNumbers)):
                        t_index = np.argmin(np.abs(self.time - times[i]))
                        haloMatch = self.haloNumber[t_index] == haloNumbers[i]
                        if any(haloMatch):
				output[i] = self.distanceMatrix[t_index][haloMatch,0]
                        else:
                                output[i] = np.nan
                return output

	def retraceVirialRadius(self, haloNumbers, times):

                assert len(haloNumbers) == len(times)

                output = np.zeros(len(haloNumbers))
                for i in range(len(haloNumbers)):
                        t_index = np.argmin(np.abs(self.time - times[i]))
                        haloMatch = self.haloNumber[t_index] == haloNumbers[i]
                        if any(haloMatch):
				output[i] = self.Rvir[t_index][haloMatch]
			else:
				output[i] = np.nan
		return output
