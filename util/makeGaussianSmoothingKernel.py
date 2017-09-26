import numpy as np

'''
def makeGaussianSmoothingKernel(nbins):
	gaussian = lambda x: np.exp(-0.5*(float(x)/nbins)**2) / nbins / np.sqrt(2*np.pi)
	return [gaussian(x) for x in np.linspace(-3*nbins,3*nbins,num=6*nbins+1)]
'''

def makeGaussianSmoothingKernel(widthInBins, maxSigma=4):

        halfRangeInBins = np.floor(widthInBins * maxSigma)
        if halfRangeInBins < 1:
                return [1]
        else:
                gaussian = lambda x: np.exp(-0.5*(float(x)/widthInBins)**2) / widthInBins / np.sqrt(2*np.pi)
                binsSampled = np.linspace(-halfRangeInBins,halfRangeInBins,num=2*halfRangeInBins+1)
                return [gaussian(x) for x in binsSampled]
