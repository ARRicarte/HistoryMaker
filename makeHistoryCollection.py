"""
ARR: 09.15.17

Make a collection of histories and bind them all into a dictionary.
"""

import tangos as db
from tangos.live_calculation import NoResultsError
from getSuitableHalos import *
from stitched_merger_finder import *
from makeHistory import *
from clusterProfiler_powerlaw import *
import cPickle as pickle
import time
import smtplib
import constants
from email.mime.text import MIMEText

def createHistoryCollection(step, pickleName, maximumSkips=5, cutoffDistance=2, minStellarMass=1e8, contaminationTolerance=0.05, \
	minDarkParticles=1e4, requireBH=True, emailAddress=None, computeRamPressure=True, computeMergers=True, massForRatio='Mstar', \
	bhString="bh('BH_mass', 'max', 'BH_central')"):
	"""
	Create a dictionary of histories.

	:arg step - A timestep from which to start
	:arg pickleName - The name of the output file.  The suffix '.pkl' is not added by this function.

	:kwarg maximumSkips - The maximum number of skips stitched_reverse_property_cascade is allowed to take.
	:kwarg cutoffDistance - The maximum distance a SMBH is allowed to be from the center of its host.
	:kwarg minStellarMass - Minimum stellar mass of this collection.
	:kwarg contaminationTolerance - The maximum contamination allowed for the halo.
	:kwarg minDarkParticles - The minimum number of DM particles allowed in these galaxies.
	:kwarg requireBH - Whether or not we try to include galaxies without SMBHs.  (Not yet implemented).
	:kwarg emailAddress - The email address for an update when this function is finished.
	"""

	#Time the calculation
	t_start = time.time()

	#Obtain halos that meet the requirements.
	haloList = getSuitableHalos(step, requireBH=requireBH, minStellarMass=minStellarMass, contaminationTolerance=contaminationTolerance, \
	minDarkParticles=minDarkParticles)
	historyCollection = {}
	failedHaloNumbers = []

	#Loop through and find histories.
	for h_index in range(len(haloList)):
		haloNumber = haloList[h_index].halo_number
		print "Processing halo_number {0}, halo {1} of {2}.".format(haloNumber, h_index+1, len(haloList))
		try:
			historyBook = makeHistory(haloList[h_index], maximumSkips=maximumSkips, cutoffDistance=cutoffDistance, bhString=bhString)
			historyCollection[haloNumber] = historyBook
		except NoResultsError:
			#The galaxy lacks one of the items asked for, probably a BH.
			print "   FAILED"
			failedHaloNumbers.append(haloNumber)
			pass
		if computeMergers:
			mergerTimes, mergerRatios = stitched_merger_finder(haloList[h_index], maximumSkips=maximumSkips, \
                        cutoffDistance=cutoffDistance, massForRatio=massForRatio)
                        historyCollection[haloNumber]['mergerTimes'] = mergerTimes
                        historyCollection[haloNumber]['mergerRatios'] = mergerRatios

	if step.simulation.basename == 'h1.cosmo50':
		#Adding one new key:  The distance from the cluster center
		try:
			clusterCoordinates = historyCollection[1]['SSC']
			clusterRadius = historyCollection[1]['R200']
			print "Computing cluster distances."
			for h_index in range(len(haloList)):
				haloNumber = haloList[h_index].halo_number
				if (haloNumber == 1) | (haloNumber in failedHaloNumbers):
					continue
				else:
					displacement = historyCollection[haloNumber]['SSC'] - clusterCoordinates
					distance = np.sqrt(np.array([np.dot(displacement[:,i],displacement[:,i]) for i in range(displacement.shape[1])]))
					historyCollection[haloNumber]['clusterDistance'] = distance / clusterRadius
		except KeyError:
			#The cluster coordinates aren't in this set.
			pass

		#Adding another key:  ram pressure
		if computeRamPressure:
			cp = ClusterProfiler(step)
			clusterVelocity = historyCollection[1]['Vcom']
			for h_index in range(len(haloList)):
				haloNumber = haloList[h_index].halo_number
				if (haloNumber == 1) | (haloNumber in failedHaloNumbers):
					continue
				else:
					clusterDensity = cp.computeGasDensity(historyCollection[haloNumber]['clusterDistance'], \
					historyCollection[haloNumber]['time'])
					relativeVelocities = historyCollection[haloNumber]['Vcom'] - clusterVelocity
					relativeSpeedSquared = np.array([np.dot(relativeVelocities[:,i],relativeVelocities[:,i]) for i in range(relativeVelocities.shape[1])])
					historyCollection[haloNumber]['ramPressure'] = clusterDensity * relativeSpeedSquared * constants.M_sun / (constants.pc * 1e3)**3 * 1e6
	
	#Pickle the output
	historyCollection['failedHaloNumbers'] = failedHaloNumbers
	with open(pickleName, 'w') as myfile:
		pickle.dump(historyCollection, myfile)
	
	t_end = time.time()

	print "Process complete after {0:3.2f} hours.".format((t_end-t_start)/60/60)
	print "Saved to {0}.".format(pickleName)

	if emailAddress is not None:
		#Send an optional email alert.
		msg = MIMEText("Hey there!\n\nIt's me, your friend the HistoryMaker.  I'm just emailing to let you know that the collection you asked for is done.  It took me {0:3.2f} hours to complete.\n\nRegards,\nHM".format((t_end-t_start)/60/60))
		msg['From'] = "HistoryMaker"
		msg['To'] = emailAddress
		msg['Subject'] = "History Collection Complete"

		s = smtplib.SMTP('localhost')
		s.sendmail("HistoryMaker", [emailAddress], msg.as_string())
		s.quit()
