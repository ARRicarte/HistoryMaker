"""
ARR: 09.15.17

Determines the halos which are suitable for tracking backwards in time.
"""

import tangos as db
import numpy as np
from util import crossmatch

def getSuitableHalos(step, minStellarMass=1e8, contaminationTolerance=0.05, minDarkParticles=1e4, requireBH=True):
        """
        Given a time step, return a list of halos.
        """

        #Making sure that we're in the zoom-in region, and not too contaminated.
	if contaminationTolerance is None:
		if requireBH:
			hid, mstar, ndm, mbh = step.gather_property('halo_number()', 'Mstar', 'NDM()', "bh().BH_mass")
			goodHalos = (mstar >= minStellarMass) & (ndm >= minDarkParticles) & (mbh > 0)
		else:
			hid, mstar, ndm = step.gather_property('halo_number()', 'Mstar', 'NDM()')
			goodHalos = (mstar >= minStellarMass) & (ndm >= minDarkParticles)
	else:
		if requireBH:
			hid, mstar, ndm, mbh, contamination = step.gather_property('halo_number()', 'Mstar', 'NDM()', "bh().BH_mass", 'contamination_fraction')
			goodHalos = (contamination < contaminationTolerance) & (mstar >= minStellarMass) & (ndm >= minDarkParticles) & (mbh > 0)
		else:
			hid, mstar, ndm, contamination = step.gather_property('halo_number()', 'Mstar', 'NDM()', 'contamination_fraction')
			goodHalos = (contamination < contaminationTolerance) & (mstar >= minStellarMass) & (ndm >= minDarkParticles)

        goodHaloNumbers = hid[goodHalos]
        goodHaloMasses = mstar[goodHalos]
        order = np.flipud(np.argsort(goodHaloMasses))
        finalHaloNumbers = goodHaloNumbers[order]

        #Because halo_number() does not always increase by 1, I must convert to indices.
        allHaloNumbers, = step.gather_property('halo_number()')
        finalHaloIndices = crossmatch(finalHaloNumbers, allHaloNumbers)[1]
	return [step.halos[i] for i in finalHaloIndices]
