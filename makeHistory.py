import numpy as np
import tangos as db
from stitched_reverse_property_cascade import *

nbins = 2000
tmax_Gyr = 20.0

def bin_index(time):
        """
        Convert a time (Gyr) to a bin index in the histogram
        """

        index = int(nbins * time / tmax_Gyr)
        if index < 0:
                index = 0
        return index

def makeHistory(halo, bhString="bh('BH_central_distance', 'min', 'BH_central')", \
	maximumSkips=5, cutoffDistance=2):
	"""
	Track this halo as far back in time as possible.  Make arrays with the same resolution as
	mdot histograms.

	:arg halo - an input halo of type tangos.core.Halo

	:kwarg bhString - the selection of black hole to use for this reconstruction
	:kwarg maximumSkips - the maximum number of skips allowed when trying to reconstruct a history based on
        tracking the central black hole backwards in time
        :kwarg cutoffDistance - the maximum number of kpc that the central black hole is allowed to be from the 
        center of its host halo for tracking

	:returns historyBook - a dictionary of various pre-determined arrays
	"""

	#The existence of a black hole will add keys.
	hasBH = 'BH_central' in halo.keys()

	#These are the properties we will trace backwards in time.
	if halo.timestep.simulation.basename == 'cosmo25':
		rawProperties = ["t()", "halo_number()", "Mstar", "raw(SFR_histogram)", "Mvir", "radius(200)", "Mgas", "Mcold", "shrink_center"]
		if hasBH:
			rawBlackHoleProperties = ['BH_mass', 'raw(BH_mdot_histogram)', 'BH_central_distance']
			allRawProperties = rawProperties + [bhString+'.'+rawBHProp for rawBHProp in rawBlackHoleProperties]
			dictionaryNames = ["Time", "haloNumber", "Mstar", "SFR", "Mvir", "R200", "Mgas", "Mcold", "SSC", \
			"Mbh", "BHAR", "Dbh"]
		else:
			allRawProperties = rawProperties
			dictionaryNames = ["Time", "haloNumber", "Mstar", "SFR", "Mvir", "R200", "Mgas", "Mcold", "SSC"]
	elif halo.timestep.simulation.basename == 'h1.cosmo50':
		rawProperties = ["t()", "halo_number()", "Mstar", "raw(SFR_histogram)", "Mvir", "radius(200)", "Mgas", "MColdGas", "shrink_center", "Vcom"]
		if hasBH:
			rawBlackHoleProperties = ['BH_mass', 'raw(BH_mdot_histogram)', 'BH_central_distance']
			allRawProperties = rawProperties + [bhString+'.'+rawBHProp for rawBHProp in rawBlackHoleProperties]
			dictionaryNames = ["Time", "haloNumber", "Mstar", "SFR", "Mvir", "R200", "Mgas", "Mcold", "SSC", "Vcom", \
			"Mbh", "BHAR", "Dbh"]
		else:
			allRawProperties = rawProperties
			dictionaryNames = ["Time", "haloNumber", "Mstar", "SFR", "Mvir", "R200", "Mgas", "Mcold", "SSC", "Vcom"]

	#Get all the properties
	print "Querying database with a stitched_reverse_property_cascade."
	if hasBH:
		time, haloNumber, mstar, sfr, mvir, rvir, mgas, mcold, ssc, vel, mbh, bhar, dbh = stitched_reverse_property_cascade(halo, allRawProperties, \
		maximumSkips=maximumSkips, cutoffDistance=cutoffDistance)
	else:
		time, haloNumber, mstar, sfr, mvir, rvir, mgas, mcold, ssc, vel = stitched_reverse_property_cascade(halo, allRawProperties, \
                maximumSkips=maximumSkips, cutoffDistance=cutoffDistance)

	#For the histograms, some assembly is required.
	print "Processing histograms."	
	combinedSFR = np.zeros(bin_index(time[0]))

	for t_i, sfr_i in zip(time, sfr):
		#The start and end indices have overlap; don't worry.  Histograms go back a fixed time.
                end = bin_index(t_i)
                start = np.max((end - len(sfr_i), 0))

		#Raw SFR info is in solar masses per Gyr, for some reason.
		sfr_i = [val/1e9 for val in sfr_i]
		combinedSFR[start:end] = sfr_i

	if hasBH:
		combinedBHAR = np.zeros(bin_index(time[0]))
		tracedMbh = np.zeros(bin_index(time[0]))
		for t_i, m_i, bhar_i in zip(time, mbh, bhar):
			#The start and end indices have overlap; don't worry.  Histograms go back a fixed time.
			end = bin_index(t_i)
			start = np.max((end - len(bhar_i), 0))

			#Contingency in case a BH is detected in one step, but not a nearby one.
			combinedBHAR[start:end] = np.maximum(bhar_i[-(end-start):], combinedBHAR[start:end])

			#Retrace black hole mass with the resolution of the histogram.  Cannot account for BH mergers.
			cumulativeBHAR = np.cumsum(bhar_i) * tmax_Gyr * 1e9 / nbins
			cumulativeBHAR -= cumulativeBHAR[-1]
			tracedMbh[start:end] = cumulativeBHAR + m_i

	#This is a denser time axis than time, corresponding to the values in combinedBHAR and combinedSFR
	tracedTime = time[0] - np.arange(len(combinedSFR)-1,-1,-1)*tmax_Gyr/nbins

	#Trace stars.  Not doing the same thing as I do with Mbh due to stripping, accretion, and uncertainties with halo finding.
	tracedMstar = np.interp(tracedTime, np.flipud(time), np.flipud(mstar), left=0)

	#Get everything else to the same time resolution by interpolation.
	tracedMvir = np.interp(tracedTime, np.flipud(time), np.flipud(mvir), left=0)
	tracedR200 = np.interp(tracedTime, np.flipud(time), np.flipud(rvir), left=0)
	tracedMgas = np.interp(tracedTime, np.flipud(time), np.flipud(mgas), left=0)
	tracedMcold = np.interp(tracedTime, np.flipud(time), np.flipud(mcold), left=0)
	if hasBH:
		tracedDbh = np.interp(tracedTime, np.flipud(time), np.flipud(dbh), left=np.inf)

	#Interpolating each dimension of space separately, then combining.
	ssc = np.array(ssc)
	tracedx = np.interp(tracedTime, np.flipud(time), np.flipud(ssc[:,0]), left=np.inf)
	tracedy = np.interp(tracedTime, np.flipud(time), np.flipud(ssc[:,1]), left=np.inf)
	tracedz = np.interp(tracedTime, np.flipud(time), np.flipud(ssc[:,2]), left=np.inf)
	tracedCoordinates = np.vstack((tracedx, tracedy, tracedz))

	#The same must be done for velocities
	vel = np.array(vel)
        tracedvx = np.interp(tracedTime, np.flipud(time), np.flipud(vel[:,0]), left=0)
        tracedvy = np.interp(tracedTime, np.flipud(time), np.flipud(vel[:,1]), left=0)
        tracedvz = np.interp(tracedTime, np.flipud(time), np.flipud(vel[:,2]), left=0)
        tracedVelocities = np.vstack((tracedvx, tracedvy, tracedvz))

	#Combine everything into a dictionary
	if hasBH:
		historyBook = {"time": tracedTime, "haloNumber": np.array(haloNumber), "Mstar": tracedMstar, "SFR": combinedSFR, "Mvir": tracedMvir, \
		"R200": tracedR200, "Mgas": tracedMgas, "Mcold": tracedMcold, "SSC": tracedCoordinates, "Vcom": tracedVelocities, "t_slice": time, \
		"Mbh": tracedMbh, "BHAR": combinedBHAR, "Dbh": tracedDbh}
	else:
		historyBook = {"time": tracedTime, "haloNumber": np.array(haloNumber), "Mstar": tracedMstar, "SFR": combinedSFR, "Mvir": tracedMvir, \
                "R200": tracedR200, "Mgas": tracedMgas, "Mcold": tracedMcold, "SSC": tracedCoordinates, "Vcom": tracedVelocities, "t_slice": time}

	return historyBook

if __name__ == '__main__':
	simulationName = 'h1.cosmo50'
        sim = db.get_simulation(simulationName)
        #This is z=0.14
        startingStep = 59
        step = sim.timesteps[startingStep]

	book = makeHistory(step.halos[2])
