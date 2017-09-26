"""
ARR:  09.18.17

Plot the results of makeHistoryCollection.py
"""

import matplotlib.pyplot as plt
from tangos.examples import mergers
import numpy as np
from util import makeGaussianSmoothingKernel
import cPickle as pickle

class HistoryPlotter(object):

	def __init__(self, inputPickleName, showMergers=False, showDistance=True, smoothingWidth=3, \
		outputDirectory=None, majorMergerThreshold=0.25):
		"""
		Read in dictionary that includes all history data.
		"""

		#Read dictionary.
		with open(inputPickleName, 'r') as myfile:
			self.historyBook = pickle.load(myfile)

		#Save halo numbers for ease of access.  Only grabbing integers, since there may be other keys.
		self.haloNumbers = np.array([key for key in self.historyBook.keys() if isinstance(key, int)])

		#These are commonly used plotting options.
		self.showMergers = showMergers
		self.majorMergerThreshold = majorMergerThreshold
		self.showDistance = showDistance
		if outputDirectory is not None:
			if outputDirectory[-1] != '/':
				outputDirectory = outputDirectory + '/'
		self.outputDirectory = outputDirectory
		self._smoothingWidth = smoothingWidth
		self._smoothingKernel = makeGaussianSmoothingKernel(smoothingWidth)

	@property
	def smoothingWidth(self):
		return self._smoothingWidth

	@smoothingWidth.setter
	def smoothingWidth(self, value):
		self._smoothingWidth = value
		self._smoothingKernel = makeGaussianSmoothingKernel(value)

	def _addMergerMarkers(self, ax, haloNumber):
		"""
		Add merger bars to the plot.
		"""

		try:
			mergerTimes = self.historyBook[haloNumber]['mergerTimes']
			mergerRatios = self.historyBook[haloNumber]['mergerRatios']
		except KeyError:
			print "Warning: Merger information is not available for halo number {0}.".format(haloNumber)
			return

		for j in range(len(mergerTimes)):
			if mergerRatios[j] > self.majorMergerThreshold:
				color = 'r'
			else:
				color = 'k'
			ax.fill_between([mergerTimes[j][0],mergerTimes[j][1]], [-1e100,-1e100], [1e100,1e100], color=color, alpha=mergerRatios[j])

	def _addDistanceAxis(self, ax, haloNumber):
		"""
		Add a distance axis
		"""

		try:
			clusterDistance = self.historyBook[haloNumber]['clusterDistance']
		except KeyError:
			print "Warning: Cluster distance information not available for halo number {0}.".format(haloNumber)
			return

		ax_dist = ax.twinx()
		ax_dist.plot(self.historyBook[haloNumber]['time'], clusterDistance, lw=2, color='orange', ls='--')
		ax_dist.set_xlim(self.historyBook[haloNumber]['time'][0], self.historyBook[haloNumber]['time'][-1])
		ax_dist.set_ylim(0,3)
		ax_dist.set_ylabel(r"$D_\mathrm{Center}/R_\mathrm{vir,cluster}$", fontsize=16)

	def plotGrowth(self, haloNumberList):
		"""
		BHAR and SFR
		"""

		haloNumbers = np.atleast_1d(haloNumberList)
		for haloNumber in haloNumbers:
			fig, ax = plt.subplots()

			#Smooth data
			smooth_sfr = np.convolve(self.historyBook[haloNumber]['SFR'], self._smoothingKernel, mode='same')
			smooth_bhar = np.convolve(self.historyBook[haloNumber]['BHAR'], self._smoothingKernel, mode='same')
			time = self.historyBook[haloNumber]['time']

			#Let's also estimate the amount of gas mass depletion unexplained by SFR or BHAR.
			mgas = self.historyBook[haloNumber]['Mgas']
			dMdt_total = np.zeros(time.shape)
			dMdt_total[1:-1] = (mgas[2:] - mgas[:-2]) / (time[2:] - time[:-2])
			dMdt_total[0] = (mgas[1] - mgas[0]) / (time[1] - time[0])
			dMdt_total[-1] = (mgas[-1] - mgas[-2]) / (time[-1] - time[-2])
			dMdt_total /= 1e9

			#Plot
			ax.plot(time, smooth_sfr, lw=2, color='b', ls='-', label=r"$\dot{M}_*$")
			ax.plot(time, smooth_bhar/7e-4, lw=2, color='g', ls='-', label=r"$\dot{M}_\bullet$/7e-4")
			ax.plot(time, -1*dMdt_total, lw=2, color='r', ls='-', label=r"$-\dot{M}_g$")

			#Format
			ax.set_xlim(time[0], time[-1])
			ax.set_ylim(0, np.nanmax(smooth_sfr)*1.3)
			ax.set_xlabel(r'Age of the Universe [Gyr]', fontsize=16)
			ax.set_ylabel('Growth Rate [M$_\odot$ yr$^{-1}$]', fontsize=16)
			plt.figtext(0.15, 0.9, "#{0}:  $M_*$ = {1:1.1e} $M_\odot$".format(haloNumber, self.historyBook[haloNumber]['Mstar'][-1]), fontsize=12)
			ax.legend(loc='upper right', framealpha=0.5, fontsize=12)

			#Optional additions
			if self.showMergers:
				self._addMergerMarkers(ax, haloNumber)
			if self.showDistance:
				self._addDistanceAxis(ax, haloNumber)

			fig.tight_layout()

			#Save or show
			if self.outputDirectory is not None:
				fig.savefig(self.outputDirectory+'growth_halo{0}.png'.format(haloNumber))
			else:
				fig.show()
				raw_input("Please enter when finished.\n")
			plt.close()

	def plotSpecificGrowth(self, haloNumberList):
		"""
		Specific BHAR and SFR
		"""
		haloNumbers = np.atleast_1d(haloNumberList)
		for haloNumber in haloNumbers:
			fig, ax = plt.subplots()

			#Smooth data
			smooth_ssfr = np.convolve(self.historyBook[haloNumber]['SFR']/self.historyBook[haloNumber]['Mstar'], self._smoothingKernel, mode='same')
			smooth_sbhar = np.convolve(self.historyBook[haloNumber]['BHAR']/self.historyBook[haloNumber]['Mbh'], self._smoothingKernel, mode='same')
			time = self.historyBook[haloNumber]['time']
			
			#Plot
			ax.plot(time, smooth_ssfr, lw=2, color='b', ls='-', label=r"$\dot{M}_*/M_*$")
			ax.plot(time, smooth_sbhar, lw=2, color='g', ls='-', label=r"$\dot{M}_\bullet/M_\bullet$")

			#Format
			ax.set_yscale('log')
			ax.set_xlim(time[0], time[-1])
			ax.set_ylim(1e-13, 1e-6)
			ax.set_xlabel(r'Age of the Universe [Gyr]', fontsize=16)
			ax.set_ylabel('Specific Growth Rate [yr$^{-1}$]', fontsize=16)
                        plt.figtext(0.18, 0.9, "#{0}:  $M_*$ = {1:1.1e} $M_\odot$".format(haloNumber, self.historyBook[haloNumber]['Mstar'][-1]), fontsize=12)
                        ax.legend(loc='upper right', framealpha=0.5, fontsize=12)

			#Optional additions
			if self.showMergers:
				self._addMergerMarkers(ax, haloNumber)
			if self.showDistance:
				self._addDistanceAxis(ax, haloNumber)

                        fig.tight_layout()

			#Save or show
			if self.outputDirectory is not None:
				fig.savefig(self.outputDirectory+'specificGrowth_halo{0}.png'.format(haloNumber))
			else:
				fig.show()
				raw_input("Please enter when finished.\n")
			plt.close()

	def plotMass(self, haloNumberList):
		"""
		Mstar, Mbh, Mvir, Mgas
		"""
		haloNumbers = np.atleast_1d(haloNumberList)
                for haloNumber in haloNumbers:
                        fig, ax = plt.subplots()
        
                        #Get data.  It is assumed that these do not need smoothing.
			time = self.historyBook[haloNumber]['time']
			mvir = self.historyBook[haloNumber]['Mvir']
			mstar = self.historyBook[haloNumber]['Mstar']
			mbh = self.historyBook[haloNumber]['Mbh']
			mgas = self.historyBook[haloNumber]['Mgas']

			#Plot
			ax.plot(time, mvir, color='purple', lw=2, ls='-', label=r"$M_\mathrm{vir}$")
			ax.plot(time, mstar, color='b', lw=2, ls='-', label=r"$M_*$")
			ax.plot(time, mbh, color='g', lw=2, ls='-', label=r"$M_\bullet$")
			ax.plot(time, mgas, color='brown', lw=2, ls='-', label=r"$M_\mathrm{gas}$")

			#Format
			ax.set_yscale('log')
			ax.set_xlim(time[0], time[-1])
			ax.set_ylim(1e6, 1e15)
			ax.set_xlabel(r'Age of the Universe [Gyr]', fontsize=16)
			ax.set_ylabel(r'Mass [$M_\odot$]', fontsize=16)

			plt.figtext(0.15, 0.9, "#{0}:  $M_*$ = {1:1.1e} $M_\odot$".format(haloNumber, self.historyBook[haloNumber]['Mstar'][-1]), fontsize=12)
			ax.legend(loc='upper right', framealpha=0.5, fontsize=12, ncol=2)

			#Optional additions
			if self.showMergers:
				self._addMergerMarkers(ax, haloNumber)
			if self.showDistance:
				self._addDistanceAxis(ax, haloNumber)

                        fig.tight_layout()

			#Save or show
                        if self.outputDirectory is not None:
                                fig.savefig(self.outputDirectory+'mass_halo{0}.png'.format(haloNumber))
                        else:
                                fig.show()
				raw_input("Please enter when finished.\n")
                        plt.close()

	def makeAllPlots(self, haloNumberList=None):
		"""
		Make all of the plots for each of the halo numbers given.
		"""

		if haloNumberList is None:
			haloNumbers = self.haloNumbers
		else:
			haloNumbers = np.atleast_1d(haloNumberList)

		for haloNumber in haloNumbers:
			print "Plotting halo number {0}.".format(haloNumber)
			self.plotGrowth(haloNumber)
			self.plotSpecificGrowth(haloNumber)
			self.plotMass(haloNumber)
