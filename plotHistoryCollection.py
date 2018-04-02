"""
ARR:  09.18.17

Plot the results of makeHistoryCollection.py
"""

import matplotlib.pyplot as plt
from tangos.examples import mergers
import numpy as np
from util import makeGaussianSmoothingKernel, t2z
import cPickle as pickle

class HistoryPlotter(object):

	def __init__(self, inputPickleName, showMergers=True, showDistance=False, showPressure=True, smoothingWidth=3, \
		showStellarMass=True, showVirialMass=True, showBlackHoleMass=True, showGasMass=True, showColdMass=True, \
		showBHAR=True, showSFR=True, showProximity=False, proximityFile=None, showLegend=True, showLabel=True, \
		outputDirectory=None, majorMergerThreshold=0.25, minorMergerThreshold=0.1):
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
		self.minorMergerThreshold = minorMergerThreshold
		self.showDistance = showDistance
		self.showPressure = showPressure
		if proximityFile is not None:
			self.showProximity = showProximity
			with open(proximityFile, 'r') as myfile:
				self.proximityData = pickle.load(myfile)
		self.showStellarMass = showStellarMass
		self.showVirialMass = showVirialMass
		self.showBlackHoleMass = showBlackHoleMass
		self.showGasMass = showGasMass
		self.showColdMass = showColdMass
		self.showBHAR = showBHAR
		self.showSFR = showSFR
		self.showLegend = showLegend
		self.showLabel = showLabel
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
			if mergerRatios[j] < self.minorMergerThreshold:
				continue
			if mergerRatios[j] > self.majorMergerThreshold:
				color = 'r'
			else:
				color = 'k'
			ax.fill_between([mergerTimes[j][0],mergerTimes[j][1]], [-1e100,-1e100], [1e100,1e100], color=color, alpha=mergerRatios[j])

	def _addDistanceAxis(self, ax, haloNumber, xlim=None, ylim=None):
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
		if xlim is not None:
			ax_dist.set_xlim(xlim[0], xlim[1])
		else:
			ax_dist.set_xlim(self.historyBook[haloNumber]['time'][0], self.historyBook[haloNumber]['time'][-1])
		if ylim is not None:
			ax_dist.set_ylim(ylim[0], ylim[1])
		else:
			ax_dist.set_ylim(0,3)
		ax_dist.set_ylabel(r"$D_\mathrm{Center}/R_\mathrm{vir,cluster}$", fontsize=16)

	def _addPressureAxis(self, ax, haloNumber, xlim=None, ylim=None):
		"""
		Add a ram pressure axis
		"""

		try:
                        pressure = self.historyBook[haloNumber]['ramPressure']
                except KeyError:
                        print "Warning: Ram pressure information not available for halo number {0}.".format(haloNumber)
                        return

		ax_pres = ax.twinx()
		ax_pres.plot(self.historyBook[haloNumber]['time'], pressure, lw=2, color='orange', ls='--')
		if xlim is not None:
			ax_pres.set_xlim(xlim[0], xlim[1])
		else:
			ax_pres.set_xlim(self.historyBook[haloNumber]['time'][0], self.historyBook[haloNumber]['time'][-1])
		if ylim is not None:
			ax_pres.set_ylim(ylim[0], ylim[1])
		else:
			ax_pres.set_ylim(3e-15,1e-10)
		ax_pres.set_yscale('log')
		ax_pres.set_ylabel(r"$P_\mathrm{ram}$ [Pa]", fontsize=16)

	def _addProximityAxis(self, ax, haloNumber, xlim=None, ylim=None):
                """
                Add a proximity axis

		#THIS DOESN'T WORK.

		#Need to match these with the indices of the big matrices.
		plotHaloNumbers = self.historyBook[haloNumber]['haloNumber']
		plotTimes = self.historyBook[haloNumber]['t_slice']
		plotRedshifts = t2z(plotTimes)

		t_list = []
		proximity_list = []
		for plot_index in range(len(plotRedshifts)):
			z_index = np.argmin(np.abs(plotRedshifts[plot_index]-self.proximityData['redshift']))
			h_index = np.where(self.proximityData['haloNumber'][z_index] == plotHaloNumbers[plot_index])[0]
			matrixFilter = ((self.proximityData['haloNumber'][z_index] != 1) & (self.proximityData['haloNumber'][z_index] != plotHaloNumbers[plot_index]) & ())
			relevantDistances = self.proximityData['distanceMatrix'][z_index][h_index,matrixFilter]
			relevantMasses = self.proximityData['Mstar'][matrixFilter]
                """
		pass

	def plotGrowth(self, haloNumberList, showGas=False, xlim=None, ylim=None, ylim2=None):
		"""
		BHAR and SFR
		"""

		haloNumbers = np.atleast_1d(haloNumberList)
		for haloNumber in haloNumbers:
			fig, ax = plt.subplots()
			hasBH = 'Mbh' in self.historyBook[haloNumber].keys()

			#Smooth data
			smooth_sfr = np.convolve(self.historyBook[haloNumber]['SFR'], self._smoothingKernel, mode='same')
			if hasBH:
				smooth_bhar = np.convolve(self.historyBook[haloNumber]['BHAR'], self._smoothingKernel, mode='same')
			time = self.historyBook[haloNumber]['time']

			#Let's also estimate the amount of gas mass depletion.
			if showGas:
				mgas = self.historyBook[haloNumber]['Mgas']
				dMdt_total = np.zeros(time.shape)
				dMdt_total[1:-1] = (mgas[2:] - mgas[:-2]) / (time[2:] - time[:-2])
				dMdt_total[0] = (mgas[1] - mgas[0]) / (time[1] - time[0])
				dMdt_total[-1] = (mgas[-1] - mgas[-2]) / (time[-1] - time[-2])
				dMdt_total /= 1e9

			#Plot
			if self.showSFR:
				ax.plot(time, smooth_sfr, lw=2, color='b', ls='-', label=r"$\dot{M}_*$")
			if (hasBH) & (self.showBHAR):
				ax.plot(time, smooth_bhar/7e-4, lw=2, color='g', ls='-', label=r"$\dot{M}_\bullet$/7e-4")
			if showGas:
				ax.plot(time, -1*dMdt_total, lw=2, color='r', ls='-', label=r"$-\dot{M}_g$")

			#Format
			if xlim is not None:
				ax.set_xlim(xlim[0], xlim[1])
			else:
				ax.set_xlim(time[0], time[-1])
			if ylim is not None:
				ax.set_ylim(ylim[0], ylim[1])
			else:
				ax.set_ylim(0, np.nanmax(smooth_sfr)*1.3)
			ax.set_xlabel(r'Age of the Universe [Gyr]', fontsize=16)
			ax.set_ylabel('Growth Rate [M$_\odot$ yr$^{-1}$]', fontsize=16)
			if self.showLabel:
				plt.figtext(0.15, 0.9, "#{0}:  $M_*$ = {1:1.1e} $M_\odot$".format(haloNumber, self.historyBook[haloNumber]['Mstar'][-1]), fontsize=12)
			if self.showLegend:
				ax.legend(loc='upper right', framealpha=0.5, fontsize=12)

			#Optional additions
			if self.showMergers:
				self._addMergerMarkers(ax, haloNumber)
			if self.showDistance:
				self._addDistanceAxis(ax, haloNumber, xlim=xlim, ylim=ylim2)
			if self.showPressure:
                                self._addPressureAxis(ax, haloNumber, xlim=xlim, ylim=ylim2)

			fig.tight_layout()

			#Save or show
			if self.outputDirectory is not None:
				fig.savefig(self.outputDirectory+'growth_halo{0}.png'.format(haloNumber))
			else:
				fig.show()
				raw_input("Please enter when finished.\n")
			plt.close()

	def plotSpecificGrowth(self, haloNumberList, xlim=None, ylim=None, ylim2=None):
		"""
		Specific BHAR and SFR
		"""
		haloNumbers = np.atleast_1d(haloNumberList)
		for haloNumber in haloNumbers:
			fig, ax = plt.subplots()

			hasBH = 'Mbh' in self.historyBook[haloNumber].keys()
			#Smooth data
			if self.showSFR:
				smooth_ssfr = np.convolve(self.historyBook[haloNumber]['SFR']/self.historyBook[haloNumber]['Mstar'], self._smoothingKernel, mode='same')
			if (hasBH) & (self.showBHAR):
				smooth_sbhar = np.convolve(self.historyBook[haloNumber]['BHAR']/self.historyBook[haloNumber]['Mbh'], self._smoothingKernel, mode='same')
			time = self.historyBook[haloNumber]['time']
			
			#Plot
			ax.plot(time, smooth_ssfr, lw=2, color='b', ls='-', label=r"$\dot{M}_*/M_*$")
			if hasBH:
				ax.plot(time, smooth_sbhar, lw=2, color='g', ls='-', label=r"$\dot{M}_\bullet/M_\bullet$")

			#Format
			ax.set_yscale('log')
			if xlim is not None:
				ax.set_xlim(xlim[0], xlim[1])
			else:
				ax.set_xlim(time[0], time[-1])
			if ylim is not None:
				ax.set_ylim(ylim[0], ylim[1])
			else:
				ax.set_ylim(1e-13, 1e-6)
			ax.set_xlabel(r'Age of the Universe [Gyr]', fontsize=16)
			ax.set_ylabel('Specific Growth Rate [yr$^{-1}$]', fontsize=16)
			if self.showLabel:
				plt.figtext(0.18, 0.9, "#{0}:  $M_*$ = {1:1.1e} $M_\odot$".format(haloNumber, self.historyBook[haloNumber]['Mstar'][-1]), fontsize=12)
			if self.showlegend:
				ax.legend(loc='upper right', framealpha=0.5, fontsize=12)

			#Optional additions
			if self.showMergers:
				self._addMergerMarkers(ax, haloNumber)
			if self.showDistance:
				self._addDistanceAxis(ax, haloNumber, xlim=xlim, ylim=ylim2)
			if self.showPressure:
                                self._addPressureAxis(ax, haloNumber, xlim=xlim, ylim=ylim2)

                        fig.tight_layout()

			#Save or show
			if self.outputDirectory is not None:
				fig.savefig(self.outputDirectory+'specificGrowth_halo{0}.png'.format(haloNumber))
			else:
				fig.show()
				raw_input("Please enter when finished.\n")
			plt.close()

	def plotMass(self, haloNumberList, xlim=None, ylim=None, ylim2=None):
		"""
		Mstar, Mbh, Mvir, Mgas, Mcold
		"""
		haloNumbers = np.atleast_1d(haloNumberList)
                for haloNumber in haloNumbers:
                        fig, ax = plt.subplots()
			hasBH = 'Mbh' in self.historyBook[haloNumber].keys()

                        #Get data.  It is assumed that these do not need smoothing.
			time = self.historyBook[haloNumber]['time']
			mvir = self.historyBook[haloNumber]['Mvir']
			mstar = self.historyBook[haloNumber]['Mstar']
			if hasBH:
				mbh = self.historyBook[haloNumber]['Mbh']
			mgas = self.historyBook[haloNumber]['Mgas']
			mcold = self.historyBook[haloNumber]['Mcold']

			#Plot
			if self.showVirialMass:
				ax.plot(time, mvir, color='k', lw=2, ls='-', label=r"$M_\mathrm{vir}$")
			if self.showStellarMass:
				ax.plot(time, mstar, color='b', lw=2, ls='-', label=r"$M_*$")
			if (hasBH) & (self.showBlackHoleMass):
				ax.plot(time, mbh, color='g', lw=2, ls='-', label=r"$M_\bullet$")
			if self.showGasMass:
				ax.plot(time, mgas, color='brown', lw=2, ls='-', label=r"$M_\mathrm{gas}$")
			if self.showColdMass:
				ax.plot(time, mcold, color='indigo', lw=2, ls='-', label=r"$M_\mathrm{cold}$")

			#Format
			ax.set_yscale('log')
			if xlim is not None:
				ax.set_xlim(xlim[0], xlim[1])
			else:
				ax.set_xlim(time[0], time[-1])
			if ylim is not None:
				ax.set_ylim(ylim[0], ylim[1])
			else:
				ax.set_ylim(1e6, 1e15)
			ax.set_xlabel(r'Age of the Universe [Gyr]', fontsize=16)
			ax.set_ylabel(r'Mass [$M_\odot$]', fontsize=16)

			if self.showLabel:
				plt.figtext(0.15, 0.9, "#{0}:  $M_*$ = {1:1.1e} $M_\odot$".format(haloNumber, self.historyBook[haloNumber]['Mstar'][-1]), fontsize=12)
			if self.showLegend:
				ax.legend(loc='upper right', framealpha=0.5, fontsize=12, ncol=2)

			#Optional additions
			if self.showMergers:
				self._addMergerMarkers(ax, haloNumber)
			if self.showDistance:
				self._addDistanceAxis(ax, haloNumber, xlim=xlim, ylim=ylim2)
			if self.showPressure:
				self._addPressureAxis(ax, haloNumber, xlim=xlim, ylim=ylim2)

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
