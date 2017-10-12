"""
ARR:  09.13.17

Attempts to correct for errors in reverse_property_cascade due to poor halo finding.
"""

import tangos as db
from tangos.live_calculation import NoResultsError
import numpy as np

def stitched_reverse_property_cascade(halo, propertyList, maximumSkips=5, cutoffDistance=2):
        """
        Given a halo and a list of properties, do a reverse property cascade and try to correct for missing halos
        by following central black holes.
        
        :arg halo - a halo of type tangos.core.Halo
        :arg propertyList - list of strings corresponding to halo keys

        :kwarg maximumSkips - the maximum number of skips allowed when trying to reconstruct a history based on
        tracking the central black hole backwards in time
        :kwarg cutoffDistance - the maximum number of kpc that the central black hole is allowed to be from the 
        center of its host halo for tracking

        :returns outputList - a list of properties going back in time, just like halo.reverse_property_cascade() 
        is supposed to return.  Note that this is indeed list instead of array format, due to shape inconsistencies
	 with raw histograms.
        """

        #Try to get as many values as the index of the halo's timestep.
        expectedLength = halo.timestep.simulation.timesteps.index(halo.timestep) + 1
        outputList = [[] for prop in propertyList]

        while True:
		try:
			#Do a reverse property cascade and append to output
			cascadedProperties = halo.reverse_property_cascade(*propertyList)
			for i in range(len(propertyList)):
				outputList[i].extend(cascadedProperties[i].tolist())
		except NoResultsError:
			#Missing properties that you wanted.
			break

                if len(outputList[0]) == expectedLength:
                        #You did it!
                        break
                else:
                        #Go to the last halo for which we have data.
                        if len(cascadedProperties[0]) == 1:
                                latestHalo = halo
                        else:
                                latestHalo = halo.calculate('earlier({0})'.format(len(cascadedProperties[0])-1))
                        problemHalo = latestHalo.previous
                        if problemHalo is None:
                                #That means the halo just didn't exist in the previous time step.  You should be done.
                                break
                        else:
                                try:
                                        #Let's find the most central black hole.  We'll track its halo history backwards.
                                        holeBeforeProblem = latestHalo.calculate("bh('BH_central_distance', 'min', 'BH_central')")
                                except NoResultsError:
                                        #Too bad, there are no central black holes to do this with.  Abort.
                                        break

				problemHole = holeBeforeProblem.previous
				if problemHole is None:
					#The BH just got seeded and there's nothing else to do.
					break
				if 'host_halo' in problemHole.keys():
					#Make sure there really is a kink in the tree going forward in time.
					relatedHalos = problemHalo['ptcls_in_common']
					if not hasattr(relatedHalos, '__len__'):
						relatedHalos = [relatedHalos]
					redshiftsOfChildHalos = [h.timestep.redshift for h in relatedHalos]
					if (redshiftsOfChildHalos.count(latestHalo.timestep.redshift) == 1) & (latestHalo in relatedHalos):
						#There was no problem with identifying the halo; something else went wrong.  Maybe a key you're after went missing.
						break

                                if holeBeforeProblem['BH_central_distance'] > cutoffDistance:
                                        #Alas, that's not much of a central black hole.  Abort.
                                        break

                                #Retrace its steps to before the problem
				try:
					holePastProblem = holeBeforeProblem.calculate('earlier(2)')
					hostPastProblem = holePastProblem['host_halo']
				except:
					#No previous time step, or the SMBH has no host.
					break
                                nSkips = 1
                                stitchFound = False
                                while (nSkips <= maximumSkips) & (holePastProblem.previous is not None):
                                        if holePastProblem['BH_central_distance'] <= cutoffDistance:
                                                #The gap has been breached!
                                                stitchFound = True
                                                break
                                        else:
                                                #This probably means that this SMBH is actually in a satellite that the halo finder did not detect.
                                                holePastProblem = holePastProblem.previous
						if holePastProblem is None:
							#The SMBH doesn't have previous time steps.  Tough.
							break
						else:
							hostPastProblem = holePastProblem['host_halo']
							nSkips += 1
                                if stitchFound:
                                        #Continue from the newly found halo.
                                        halo = hostPastProblem
                                else:
                                        #Stitching failed.  Just exit now.
                                        break

        return outputList
