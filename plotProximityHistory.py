import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
from useProximityTable import ProximityCalculator

class ProximityPlotter(object):

	def __init__(self, historyFile, proximityFile, mode='threshold', massType='Mstar', ratioThreshold=0.1):

		with open(historyFile, 'r') as myfile:
			self.historyBook = pickle.load(myfile)

		self.proximityCalculator = ProximityCalculator(proximityFile, mode=mode, massType=massType, ratioThreshold=ratioThreshold)

        def _addMergerMarkers(self, ax, haloNumber, minorMergerThreshold=0.1, majorMergerThreshold=0.25):
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
                        if mergerRatios[j] < minorMergerThreshold:
                                continue
                        if mergerRatios[j] > majorMergerThreshold:
                                color = 'r'
                        else:
                                color = 'k'
                        ax.fill_between([mergerTimes[j][0],mergerTimes[j][1]], [-1e100,-1e100], [1e100,1e100], color=color, alpha=mergerRatios[j])

	def plotProximity(self, haloNumber, savename=None, showLabel=True, showLegend=True):

		fig, ax = plt.subplots()
		taxis = self.historyBook[haloNumber]['t_slice']
		haloNumbers = self.historyBook[haloNumber]['haloNumber']
		proximity = self.proximityCalculator.retraceProximity(haloNumbers, taxis)	
		radii = self.proximityCalculator.retraceVirialRadius(haloNumbers, taxis)
		clusterDistance = self.proximityCalculator.retraceClusterDistance(haloNumbers, taxis)

		ax.plot(taxis, proximity, color='forestgreen', linestyle='-', linewidth=2, label='To Significant Neighbor')
		ax.plot(taxis, clusterDistance, color='orange', linestyle='-', linewidth=2, label='To Cluster Center')
		ax.plot(taxis, radii, color='darkturquoise', linestyle='-', linewidth=2, label=r'$R_{200}$')
		ax.set_yscale('log')
		ax.set_xlabel('Time [Gyr]', fontsize=16)
		ax.set_ylabel('Distance [kpc]', fontsize=16)
		ax.set_xlim(0,14)
		ax.set_ylim(1e1,5e3)
		if showLegend:
			ax.legend(frameon=False, fontsize=14, loc='lower right')
		if showLabel:
			plt.figtext(0.15, 0.9, "#{0}:  $M_*$ = {1:1.1e} $M_\odot$".format(haloNumber, self.historyBook[haloNumber]['Mstar'][-1]), fontsize=12)

		self._addMergerMarkers(ax, haloNumber)

		fig.tight_layout()

		if savename is None:
			fig.show()
		else:
			fig.savefig(savename)
			plt.close()

	def plotProximityBetween(self, haloNumber1, haloNumber2, savename=None, showLabel=True, showLegend=True):

                fig, ax = plt.subplots()
                taxis1 = self.historyBook[haloNumber1]['t_slice']
                taxis2 = self.historyBook[haloNumber2]['t_slice']
                haloNumbers1 = self.historyBook[haloNumber1]['haloNumber']
                haloNumbers2 = self.historyBook[haloNumber2]['haloNumber']
                overlaptime, proximity = self.proximityCalculator.retraceProximityBetween(haloNumbers1, haloNumbers2, taxis1, taxis2)
                radii = self.proximityCalculator.retraceVirialRadius(haloNumbers1, taxis1)
                clusterDistance = self.proximityCalculator.retraceClusterDistance(haloNumbers1, taxis1)

                ax.plot(overlaptime, proximity, color='forestgreen', linestyle='-', linewidth=2, label='To Halo {0}'.format(haloNumber2))
                ax.plot(taxis1, clusterDistance, color='orange', linestyle='-', linewidth=2, label='To Cluster Center')
                ax.plot(taxis1, radii, color='darkturquoise', linestyle='-', linewidth=2, label=r'$R_{200}$')
                ax.set_yscale('log')
                ax.set_xlabel('Time [Gyr]', fontsize=16)
                ax.set_ylabel('Distance [kpc]', fontsize=16)
                ax.set_xlim(0,14)
                ax.set_ylim(1e1,5e3)
                if showLegend:
                        ax.legend(frameon=False, fontsize=14, loc='lower right')
                if showLabel:
                        plt.figtext(0.15, 0.9, "#{0}:  $M_*$ = {1:1.1e} $M_\odot$".format(haloNumber1, self.historyBook[haloNumber1]['Mstar'][-1]), fontsize=12)

                self._addMergerMarkers(ax, haloNumber1)

                fig.tight_layout()

                if savename is None:
                        fig.show()
                else:
                        fig.savefig(savename)
                        plt.close()

	def makePlotDirectory(self, saveDirectory='./proximityPlots/'):

		fig, ax = plt.subplots()
		allHaloNumbers = [key for key in self.historyBook.keys() if isinstance(key, int)]
		for haloNumber in allHaloNumbers:
			print "Halo Number = {0}".format(haloNumber)
			self.plotProximity(haloNumber, savename=saveDirectory+'proximity_halo{0}.png'.format(haloNumber))
