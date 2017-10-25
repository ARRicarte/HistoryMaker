"""
ARR:  10.25.17

Use the scheme of stitched_reverse_property_cascade to find halo mergers.
"""

import tangos as db
from tangos.live_calculation import NoResultsError
import numpy as np
from stitched_reverse_property_cascade import *

def stitched_merger_finder(halo, maximumSkips=5, cutoffDistance=2, massForRatio='Mstar'):
        """
        Given a halo and a list of properties, do a reverse property cascade and try to correct for missing halos
        by following central black holes.
        
        :arg halo - a halo of type tangos.core.Halo
        :arg propertyList - list of strings corresponding to halo keys

        :kwarg maximumSkips - the maximum number of skips allowed when trying to reconstruct a history based on
        tracking the central black hole backwards in time
        :kwarg cutoffDistance - the maximum number of kpc that the central black hole is allowed to be from the 
        center of its host halo for tracking
	:kwarg massForRatio - the key to use for mass ratios

        :returns mergerTimes - 2d array of merger times, since we only know the interval of merger times
	:returns mergerRatios - ratios taken with the mass specified
        """

	#First, get all progenitor halos in this roundabout way.
        times, halo_numbers = stitched_reverse_property_cascade(halo, ["t()", "halo_number()"], maximumSkips=maximumSkips, cutoffDistance=cutoffDistance)
	timesteps = halo.timestep.simulation.timesteps
	halos = []
	allTimes = np.array([step.time_gyr for step in timesteps])
	for t, h_number in zip(times, halo_numbers):
		matchedStep = timesteps[np.where(allTimes==t)[0][0]]
		all_halo_numbers, = matchedStep.gather_property('halo_number()')
		matchedHaloIndex = np.where(all_halo_numbers == h_number)[0][0]
		halos.append(matchedStep.halos[matchedHaloIndex])
	halos = np.array(halos)

	#Next, we're going to look at all of the halos along the main progenitor branch and see how many children there are.
	mergerTimes = []
	mergerRatios = []
	for halo in halos[:-1]:
		currentTime = halo.timestep.time_gyr
		previousTime = halo.timestep.previous.time_gyr
		relatives = halo['ptcls_in_common']
		if not hasattr(relatives, '__len__'):
			relatives = [relatives]
		relativeTimes = np.array([relative.timestep.time_gyr for relative in relatives])

		#This snippet breaks out if there are multiple parent nodes.  That means that there is probably a fake merger here.
		if halo.timestep.next is not None:
			nextTime = halo.timestep.next.time_gyr
			if np.sum(relativeTimes == nextTime) > 1:
				continue
		
		#Children are those halos which share particles and are in the previous step
		children = [halo for halo in relatives if (halo.timestep.time_gyr == previousTime)]

		#Only continue if there is more than one child.
		if len(children) == 1:
			continue

		#If you've made it this far, you can compute mass ratios and times.
		childMasses = np.array([child[massForRatio] for child in children])
		maximumMass = np.max(childMasses)
		mergerRatio = np.max(childMasses[childMasses!=maximumMass]) / maximumMass
		mergerTimes.append([previousTime,currentTime])
		mergerRatios.append(mergerRatio)
		
	return np.array(mergerTimes), np.array(mergerRatios)
