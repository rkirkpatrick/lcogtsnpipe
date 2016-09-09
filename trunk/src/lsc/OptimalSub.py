#!/usr/bin/env python
import numpy as np
import scipy
from astropy.io import fits
from pyraf import iraf
import time
from os import system
import argparse
from lsc.lscastrodef import sextractor
import statsmodels.api as sm


import matplotlib.pyplot as plt

class OptimalSubtraction:
    '''Implementation of proper image subtraction from Zackey, Ofek, Gal-Yam 2016
       see https://arxiv.org/abs/1601.02655 for details'''


    def __init__(self, Nf, Rf, ArgsDict = {}):
        '''Take filenames and turn them into arrays and unpack arguments in ArgsDict'''

        self.ArgsDict = ArgsDict
        self.Nf, self.Rf = Nf, Rf
        self.MakeDefaults()

        

        #self.Nf, self.Rf = Nf, Rf

        self.CheckPSFsExist()

        self.CheckPSFfromIRAF()
        self.CheckAlign()

        self.N, self.R, self.Pn, self.Pr = CropImages(fits.getdata(self.Nf), fits.getdata(self.Rf), self.Pn, self.Pr)
        self.CheckBackground()
        self.CheckFlatten()

        self.sn = self.ArgsDict['NewNoise']
        self.sr = self.ArgsDict['RefNoise']

        self.CheckGainMatched()

    def MakeDefaults(self):
        '''Set up default optional arguments'''

        default = {
            'Align': False, 
            'Beta': 99999., 
            'Interactive': False, 
            'Flatten': False, 
            'Gamma': 99999., 
            'RefPSF': '', 
            'NewPSF': '', 
            'RefBackground': self.Rf, 
            'RefFWHM': 1.5, 
            'RefNoise': 5., 
            'RefReadNoise': 0., 
            'RefSaturation': 10e6, #GetSaturationCount(self.Rf), 
            'RefVariance': '', 
            'NewBackground': self.Nf, 
            'NewFWHM': 1.5, 
            'NewNoise': 5.,
            'NewReadNoise': 0., 
            'NewSaturation': 10e6, #GetSaturationCount(self.Nf), 
            'NewVariance': '', 
            'PSFfromIRAF': True, 
            'Show': True, 
            'SigmaX': 0.1, 
            'SigmaY': 0.1, 
            'Verbose': False
        }

        default.update(self.ArgsDict)
        self.ArgsDict = default


    def FindTransient(self, Threshold = 3.):
        try:
            self.Scorr_
        except AttributeError:
            self.Scorr()
        maxima = (self.Scorr_ == scipy.ndimage.maximum_filter(self.Scorr_, 4))
        PeaksIndex = np.where(maxima)
        Peaks = self.Scorr_[PeaksIndex]
        Transients = np.array([Peaks, PeaksIndex[0], PeaksIndex[1]]).transpose()
        TransientsSortIndex = np.argsort(Transients[:,0])
        Transients = Transients[TransientsSortIndex[::-1]]
        Transients = Transients[np.where(Transients[:,0] >= Threshold)]
        np.savetxt('transients.txt', Transients)


    def CheckAlign(self):
        '''Check if images need alignment'''

        if self.ArgsDict['Align']:
            self.RBackground = AlignImages(self.Nf, self.Rf)
            self.Rf = self.Nf.replace('.fits', '.ref.fits')


    def CheckBackground(self):
        try:
            self.NBackground = fits.getdata(self.ArgsDict['NewBackground'])
        except ValueError:
            self.NBackground = self.N
        try:
            self.RBackground = fits.getdata(self.ArgsDict['RefBackground'])
        except ValueError:
            self.RBackground = self.R


    def CheckFlatten(self):
        if self.ArgsDict['Flatten']:
            self.N = SubtractBackground(self.N)
            self.R = SubtractBackground(self.R)


    def CheckGainMatched(self):
        '''Check if the gain matching parameters have been found'''

        beta = self.ArgsDict['Beta']
        gamma = self.ArgsDict['Gamma']
        # if user hasn't supplied beta and gamma, fit them
        if beta == 99999 or gamma == 99999:
            print 'Gain matching not done; fitting manually'

            self.beta, self.gamma = self.MatchGain()

        else:
            self.beta = beta
            self.gamma = gamma


    def CheckPSFsExist(self):
        '''Check if user provided a psf and make one if necessary'''

        Prf = self.ArgsDict['RefPSF']
        Pnf = self.ArgsDict['NewPSF']
        if Pnf != '':
            self.Pnf = Pnf
        else:
            FWHM = self.ArgsDict['NewFWHM']
            Noise = self.ArgsDict['NewNoise']
            self.Pnf = FitPSF(self.Nf, FWHM, Noise, Verbose = self.ArgsDict['Verbose'], Show = self.ArgsDict['Show'])
            self.ArgsDict.update({'PSFfromIRAF': True})

        if Prf != '':
            self.Prf = Prf
        else:
            FWHM = self.ArgsDict['RefFWHM']
            Noise = self.ArgsDict['RefNoise']
            self.Prf = FitPSF(self.Rf, FWHM, Noise, Verbose = self.ArgsDict['Verbose'], Show = self.ArgsDict['Show'])
            self.ArgsDict.update({'PSFfromIRAF': True})
            

    def CheckPSFfromIRAF(self):
        '''Check if user specified IRAF psf or actual psf'''

        try:
            PSFfromIRAF = self.ArgsDict['PSFfromIRAF']
            if PSFfromIRAF:
                self.Pn = ExtractPSF(self.Pnf)
                self.Pr = ExtractPSF(self.Prf)
            else:
                self.Pn = fits.getdata(self.Pnf)
                self.Pr = fits.getdata(self.Prf)
        except KeyError:
            self.Pn = fits.getdata(self.Pnf)
            self.Pr = fits.getdata(self.Prf)


    def D(self, normalize = '', diagnostic = False):
        '''Calculate proper subtraction image and normalize to zero point of reference or target'''

        N, R, Pn, Pr, sn, sr = self.N, self.R, self.Pn, self.Pr, self.sn, self.sr

        #save convolved not gain matched images that will be subtracted
        if diagnostic:
            print 'Saving Deconvolved images'
            N_hat = np.fft.fft2(self.N)    
            R_hat = np.fft.fft2(self.R)
            Pr_hat = np.fft.fft2(self.Pr)
            Pn_hat = np.fft.fft2(self.Pn)
            Den = self.sr ** 2 * abs(Pn_hat) ** 2 + self.sn ** 2 * abs(Pr_hat) ** 2
            SqrtDen = np.sqrt(Den)
            DConvN = np.real(np.fft.ifft2(N_hat * Pr_hat / SqrtDen))
            hdu = fits.PrimaryHDU(DConvN)
            hdu.writeto('NDConvolve.fits', clobber = True)
            DConvR = np.real(np.fft.ifft2(R_hat * Pn_hat / SqrtDen))
            hdu = fits.PrimaryHDU(DConvR)
            hdu.writeto('RDConvolve.fits', clobber = True)


        # calculate D
        N = np.subtract(N, self.gamma)

        beta = self.beta

        N_hat = np.fft.fft2(N)
        R_hat = np.fft.fft2(R)
        Pr_hat = np.fft.fft2(Pr)
        Pn_hat = np.fft.fft2(Pn)

        Den = sr ** 2 * beta ** 2 * abs(Pn_hat) ** 2 + sn ** 2 * abs(Pr_hat) ** 2
        SqrtDen = np.sqrt(Den)

        D = np.fft.ifft2((Pr_hat * N_hat - beta * Pn_hat * R_hat) / SqrtDen)

        # apply user's normalization choice
        if normalize == 'i':
            DNormalize = D * beta / self.Fd()
        elif normalize == 't':
            DNormalize = D / self.Fd()
        else:
            DNormalize = D
        self.D_ = np.real(DNormalize)

        return np.real(DNormalize)


    def DNoise(self):
        a = self.N.shape
        deltaD = deltaD_hat * abs(np.divide(D_hat[a,a] * np.exp(2*np.pi*1j*np.dot(a,a)) - D_hat[0, 0], 2*np.pi*1j*D_hat))


    def Fd(self):
        '''Calculate the flux based zero point of D'''
        Fd = self.beta / np.sqrt(self.sn**2 + self.sr**2*self.beta**2)
        self.Fd_ = np.real(Fd)
        return np.real(Fd)


    def Flux(self, normalize = ''):
        '''Calculate transient Flux'''
        try:
            Flux = self.S_ / self.Fs_
        except AttributeError:
            Flux = self.S(normalize) / self.Fs()
        self.Flux_ = Flux
        return Flux


    def Fs(self):
        '''Calculate flux based zeropoint of S'''

        beta = self.beta
        Pr_hat = np.fft.fft2(self.Pr)
        Pn_hat = np.fft.fft2(self.Pn)
        sn, sr = self.sn, self.sr
        Den = sr ** 2 * beta ** 2 * abs(Pn_hat) ** 2 + sn ** 2 * abs(Pr_hat) ** 2
        Fs = np.sum(beta ** 2 * abs(Pn_hat) * abs(Pr_hat) / Den)
        self.Fs_ = Fs
        return Fs

    def Kernels(self):
        '''Calculate convolution kernels'''

        N, R, Pn, Pr, sn, sr = self.N, self.R, self.Pn, self.Pr, self.sn, self.sr
        N = np.subtract(N, self.gamma)

        N_hat = np.fft.fft2(N)    
        R_hat = np.fft.fft2(R)
        Pr_hat = np.fft.fft2(Pr)
        Pn_hat = np.fft.fft2(Pn)

        beta = self.beta
        Den = sr ** 2 * beta ** 2 * abs(Pn_hat) ** 2 + sn ** 2 * abs(Pr_hat) ** 2

        Kr_hat = (beta ** 2 * np.conj(Pr_hat) * abs(Pn_hat) ** 2) / Den
        Kr = np.real(np.fft.ifft2(Kr_hat))

        Kn_hat = (beta * np.conj(Pn_hat) * abs(Pr_hat) ** 2) / Den
        Kn = np.real(np.fft.ifft2(Kn_hat))
        return Kn, Kr


    def MakeCatalog(self, SortBy = 'magnitude'):
        '''Check for source catalog and make one if necessary'''
        try:
            cat = self.Catalog
            return cat

        except AttributeError:
            if SortBy == 'x':
                SortIndex = 0
            elif SortBy == 'y':
                SortIndex = 1
            elif SortBy == 'fwhm':
                SortIndex = 2
            elif SortBy == 'flux':
                SortIndex = 3
            elif SortBy == 'magnitude':
                SortIndex = 4
            elif SortBy == 'flag':
                SortIndex = 7
            else:
                SortIndex = 4

            sources = np.array(sextractor(self.Nf)).transpose()
            source_index = np.argsort(sources[:, SourceIndex])
            source_sorted = sources[source_index]

            self.Catalog = source_sorted
            return source_sorted

    def MatchGain(self):
        '''Call gain matching function'''

        NSat = self.ArgsDict['RefSaturation']
        RSat = self.ArgsDict['NewSaturation']
        beta, gamma, betaError, gammaError = FitB(self.N, self.R, self.Pn, self.Pr, self.sn, self.sr, 
                                                  NSaturation = NSat, RSaturation = RSat, Interactive = self.ArgsDict['Interactive'])

        return beta, gamma


    def Pd(self):
        '''Calculate PSF of D'''

        Pn_hat = np.fft.fft2(self.Pn)
        Pr_hat = np.fft.fft2(self.Pr)
        Den = self.sr ** 2 * self.beta ** 2 * abs(Pn_hat) ** 2 + self.sn ** 2 * abs(Pr_hat) ** 2
        SqrtDen = np.sqrt(Den)
        Pd_hat = self.beta * Pr_hat * Pn_hat / (self.Fd() * SqrtDen)
        Pd = np.fft.ifft2(Pd_hat)
        self.Pd_ = np.real(Pd)

        return np.real(Pd)


    def S(self):
        '''Calculate S array'''

        try:
            S_hat = self.Fd_ * np.fft.fft2(self.D_) * np.conj(np.fft.fft2(self.Pd_))
            S = np.fft.ifft2(S_hat)
        except AttributeError:
            S_hat = self.Fd() * np.fft.fft2(self.D()) * np.conj(np.fft.fft2(self.Pd()))
            S = np.fft.ifft2(S_hat)
        self.S_ = np.real(S)
        return np.real(S)

    def SaveD(self, Df, normalize = ''):
        '''Calculate and save proper subtraction image to database'''

        self.Df = Df
        self.D(normalize)
        hdu = fits.PrimaryHDU(np.real(self.D_))
        hdu.header = fits.getheader(self.Nf)
        hdu.header['PHOTNORM'] = normalize
        hdu.header['BETA'] = self.beta
        hdu.header['GAMMA'] = self.gamma
        hdu.writeto(self.Df, clobber = True)

    def SaveS(self, Sf, normalize = ''):
        '''Calculate and save S to database'''

        self.Sf = Sf
        self.S()
        hdu = fits.PrimaryHDU(np.real(self.S_))
        #hdu.header = fits.getheader(self.Nf)
        hdu.header['PHOTNORM'] = normalize
        hdu.header['BETA'] = self.beta
        hdu.header['GAMMA'] = self.gamma
        hdu.writeto(self.Sf, clobber = True)

    def SaveScorr(self, Scorrf, normalize = ''):
        '''Calculate and save S to database'''

        self.Scorrf = Scorrf
        self.Scorr()
        hdu = fits.PrimaryHDU(np.real(self.Scorr_))
        #hdu.header = fits.getheader(self.Nf)
        hdu.header['PHOTNORM'] = normalize
        hdu.header['BETA'] = self.beta
        hdu.header['GAMMA'] = self.gamma
        hdu.writeto(self.Scorrf, clobber = True)

    def SaveScorrThreshold(self, ScorrThreshf, Thresh = 3., normalize = ''):
        '''Save a Scorr such that pixel = 1 if > Thresh else pix = 0'''

        self.ScorrThreshf = ScorrThreshf
        self.Scorr()
        BrightIndex = np.where(self.Scorr_ >= Thresh)
        DarkIndex = np.where(self.Scorr_ < Thresh)
        ScorrThresh = np.copy(self.Scorr_)
        ScorrThresh[BrightIndex] = 1
        ScorrThresh[DarkIndex] = 0
        self.ScorrThresh = ScorrThresh        

        hdu = fits.PrimaryHDU(np.real(self.ScorrThresh))
        #hdu.header = fits.getheader(self.Nf)
        hdu.header['PHOTNORM'] = normalize
        hdu.header['BETA'] = self.beta
        hdu.header['GAMMA'] = self.gamma
        hdu.writeto(self.ScorrThreshf, clobber = True)


    def SaveSnoise(self, Snoisef, normalize = ''):
        '''Calculate and save S to database'''

        self.Snoisef = Snoisef
        self.Snoise()
        hdu = fits.PrimaryHDU(np.real(self.Snoise_))
        #hdu.header = fits.getheader(self.Nf)
        hdu.header['PHOTNORM'] = normalize
        hdu.header['BETA'] = self.beta
        hdu.header['GAMMA'] = self.gamma
        hdu.writeto(self.Snoisef, clobber = True)

    def SaveImageToWD(self):
        '''Save various images to working directory (testing only)'''
        Images = {'S.fits': self.S, 'Snoise.fits': self.Snoise, 'Scorr.fits': self.Scorr, 'D.fits': self.D, 'Flux.fits': self.Flux, 'Pd.fits': self.Pd}

        for element in Images:
            hdu = fits.PrimaryHDU(np.real(Images[element]))
            hdu.header=fits.getheader(self.Nf)
            hdu.writeto(element, clobber = True)


    def Scorr(self):
        '''Calculate Scorr'''

        try:
            Scorr = self.S_ / self.Snoise_
        except AttributeError:
            Scorr = self.S() / np.sqrt(self.Snoise())
        self.Scorr_ = np.real(Scorr)
        return np.real(Scorr)        


    def Snoise(self):
        '''Calculate the noise image for Scorr'''
        # this whole function needs optimization

        N, R, Pn, Pr, sn, sr = self.N, self.R, self.Pn, self.Pr, self.sn, self.sr
        N = np.subtract(N, self.gamma)

        N_hat = np.fft.fft2(N)    
        R_hat = np.fft.fft2(R)
        Pr_hat = np.fft.fft2(Pr)
        Pn_hat = np.fft.fft2(Pn)

        beta = self.beta
        Den = sr ** 2 * beta ** 2 * abs(Pn_hat) ** 2 + sn ** 2 * abs(Pr_hat) ** 2

        Kn, Kr = self.Kernels()

        Kn2_hat = np.fft.fft2(Kn ** 2)
        Kr2_hat = np.fft.fft2(Kr ** 2)

        if self.ArgsDict['NewVariance'] != '':
            NVarFilename = self.ArgsDict['NewVariance']
            EpsN = fits.getdata(NVarFilename)
            s = N.shape
            EpsN = EpsN[:s[0], :s[1]]
            V_Sn = np.fft.ifft2(Kn2_hat * np.fft.fft2(EpsN))

        else:
            NReadNoise = self.ArgsDict['NewReadNoise']
            EpsN = np.add(abs(self.NBackground), NReadNoise ** 2)
            V_Sn = np.fft.ifft2(Kn2_hat * np.fft.fft2(EpsN))

        if self.ArgsDict['RefVariance'] != '':
            RVarFilename = self.ArgsDict['RefVariance']
            EpsR = fits.getdata(RVarFilename)
            s = N.shape
            EpsR = EpsR[:s[0], :s[1]]
            V_Sr = np.fft.ifft2(Kr2_hat * np.fft.fft2(EpsR))

        else:
            RReadNoise = self.ArgsDict['RefReadNoise']
            EpsR = np.add(abs(self.RBackground), RReadNoise ** 2)
            V_Sr = np.fft.ifft2(Kr2_hat * np.fft.fft2(EpsR))

        xrms = self.ArgsDict['SigmaX']
        yrms = self.ArgsDict['SigmaY']

        GradNy, GradNx = np.gradient(np.fft.ifft2(np.fft.fft2(Kn) * N_hat))
        GradRy, GradRx = np.gradient(np.fft.ifft2(np.fft.fft2(Kr) * R_hat))

        Vr_ast = xrms ** 2 * GradRx ** 2 + yrms ** 2 * GradRy ** 2
        Vn_ast = xrms ** 2 * GradNx ** 2 + yrms ** 2 * GradNy ** 2

        Snoise = V_Sn + V_Sr + Vr_ast + Vn_ast
        self.Snoise_ = np.real(Snoise)
        return np.real(Snoise)


def AlignImages(NewFile, RefFile):
    '''Align reference image to new image when both images have WCS'''

    TempFile = NewFile.replace('.fits', '.ref.fits')
    iraf.wregister(RefFile, NewFile, TempFile)
    #system('cp {0} {1}'.format(TempFile, NewFile.replace('.fits', '.ref.fits')))
    AlignedRef = fits.getdata(TempFile)
    return AlignedRef


def CleanDataIndex(Data, SatCount = None, RmBackPix = True, Significance = 1.):
    '''Clean saturated and background pixels from dataset'''

    if SatCount is not None:
        UnSatPixIndex = np.where(Data < SatCount)
    else:
        UnSatPixIndex = np.arange(Data.size)

    if RmBackPix:
        Thresh = np.median(Data) +  Significance * np.std(Data)
        NotBackPixIndex = np.where(Data > Thresh)
    else:
        NotBackPixIndex = np.arange(Data.size)
    
    GoodPixInCommon = np.intersect1d(UnSatPixIndex, NotBackPixIndex)
    return GoodPixInCommon


def ConvertWeightToVariance(image):
    '''Convert swarp weight image to variance; taken from externaldata.py'''

    # the output of swarp give the normalized weight 
    # we want to store the variance image
    # we invert the normalized weight and we "de-normalized" this image
    hd = fits.getheader(image)
    ar = fits.getdata(image)
    hd2 = fits.getheader(image.replace('.fits', '.weight.fits'))
    ar2 = fits.getdata(image.replace('.fits', '.weight.fits'))
    variance = 1 / ar2
    #  this is to take in account that the weight is normalized
    variance *= (np.median(np.abs(ar - np.median(ar)))*1.48)**2/np.median(variance)
    varimg = image.replace('.fits', '.var.fits')
    _saturate = GetSaturationCount(image)
    ar = np.where(ar2 == 0, _saturate, ar)

    fits.writeto(varimg, variance, hd2, clobber=True)

    # put the saturation all values where the weight is zero 

def CropImages(N, R, Pn, Pr):
    '''Open and crop images and PSFs to the same size'''

    sN = N.shape
    #sR = R.shape
    #s = [np.min([sN[0], sR[0]]), np.min([sN[1], sR[1]])]
    #N = N[:s[0], :s[1]]
    #R = R[:s[0], :s[1]]

    Pn_extended = PadPSF(Pn, sN)
    Pr_extended = PadPSF(Pr, sN)

    return N, R, Normalize(Pn_extended), Normalize(Pr_extended)


def GetSaturationCount(filename):
    '''Get pixel saturation count for a given image'''

    sat = []
    Header = fits.getheader(filename)
    try:
        MaxLin = Header['MAXLIN']
        sat.append(MaxLin)
    except KeyError:
        pass
    try:
        Saturate = Header['SATURATE']
        sat.append(Saturate)
    except KeyError:
        pass
    if len(sat) < 1:
        print 'Saturation count not found in fits header'
        return None
    else:
        Saturation = min(sat)
        return Saturation


def ExtractPSF(Pf):
    '''Extract PSF array from iraf PSF file'''

    iraf.noao()
    iraf.digiphot()
    iraf.daophot(_doprint = 0)
    iraf.seepsf(Pf, 'temp.psf.fits')
    PSF = fits.open('temp.psf.fits')[0].data
    system('rm temp.psf.fits')

    return PSF


def FitB(NBackground, RBackground, Pn, Pr, sn, sr, NSaturation = None, RSaturation = None, Interactive = False):
    '''Fit gain matching (beta) and background (gamma) parameters'''

    center, d = NBackground.shape[0] / 2, NBackground.shape[0]
    a, b = center - d, center + d
    coords = [[a, b, a, b]]

    for star in coords:
        N = NBackground[star[0]:star[1], star[2]:star[3]]
        R = RBackground[star[0]:star[1], star[2]:star[3]]

        # trim PSFs to image size
        Pn = np.roll(Pn, Pn.shape[0] / 2, 0)
        Pn = np.roll(Pn, Pn.shape[1] / 2, 1)
        Pn = Pn[star[0]:star[1], star[2]:star[3]]
        Pn = np.roll(Pn, Pn.shape[0] / 2, 0)
        Pn = np.roll(Pn, Pn.shape[1] / 2, 1)

        Pr = np.roll(Pr, Pr.shape[0] / 2, 0)
        Pr = np.roll(Pr, Pr.shape[1] / 2, 1)
        Pr = Pr[star[0]:star[1], star[2]:star[3]]
        Pr = np.roll(Pr, Pr.shape[0] / 2, 0)
        Pr = np.roll(Pr, Pr.shape[1] / 2, 1)

        print 'x: {0}, y: {1}'.format(star[0] + d, star[2] + d)
        # iteratively solve for linear fit of beta and gamma
        if Interactive:
            beta, gamma, betaError, gammaError = InteractiveFit(N, R, Pn, Pr, sn, sr, NSaturation = NSaturation, RSaturation = RSaturation)
        else:
            beta, gamma, betaError, gammaError = IterativeSolve(N, R, Pn, Pr, sn, sr, NSaturation = NSaturation, RSaturation = RSaturation)

    return beta, gamma, betaError, gammaError


def SubtractBackground(data, NumStamps = 30, Border = 200):

    d = (data.shape[0] - 2 * Border) * (data.shape[1] - 2 * Border) / NumStamps ** 2
    LeftEdge = np.linspace(Border, data.shape[0] - Border, NumStamps)
    LeftEdge = [int(i) for i in LeftEdge]
    TopEdge = np.linspace(Border, data.shape[1] - Border, NumStamps)
    TopEdge = [int(i) for i in TopEdge]
    SmallImg = np.zeros(2 * [NumStamps])
    for iIndex, i in enumerate(TopEdge):
        for jIndex, j in enumerate(LeftEdge):
            #print i, i+d, j, j+d
            Stamp = data[i: i + d, j: j + d]
            Hist, BinEdge = np.histogram(Stamp, bins = 50)
            x = [(BinEdge[k] + BinEdge[k+1]) / 2 for k in range(len(BinEdge) - 1)]
            popt, cov = scipy.optimize.curve_fit(Gauss, x, Hist, p0 = [np.max(data), np.median(data), np.std(data)])
            #print popt
            SmallImg[jIndex, iIndex] = popt[1]
            #plt.imshow(Stamp)
            #plt.show()
            #plt.plot(x, Hist, x, Gauss(x, *popt))
            #plt.show()
    scale = [(data.shape[0] - 2 * d) / SmallImg.shape[0], (data.shape[1] - 2 * d) / SmallImg.shape[1]]
    print scale
    BackImg = scipy.ndimage.zoom(SmallImg, scale)
    s = [data.shape[0]/2 - BackImg.shape[0]/2, data.shape[1]/2 - BackImg.shape[1]/2]
    BackImgFull = np.rot90(np.fliplr(np.pad(BackImg, ((s[0], s[0]),(s[1], s[1])), 'constant')), 1)
#    BackImgFull = np.zeros(data.shape)
#    BackImgFull[Border:data.shape[0]-Border, Border:data.shape[1]-Border] = BackImg
    dataCorr = data - BackImgFull
    fits.writeto('back.fits', BackImgFull, clobber = True)
    fits.writeto('dataCorr.fits', dataCorr, clobber = True)

    plt.imshow(BackImg, interpolation = 'none')
    plt.show()
    plt.imshow(dataCorr, interpolation = 'none')
    plt.show()
    return dataCorr

def Gauss(x, a, b, c):
    '''Return a gaussian function'''

    return a * np.exp(-(x-b)**2/(2*c**2))


def FitPSF(ImageFile, FWHM = 5., Noise = 30., Verbose = True, Show = True, MaxCount = 15000.):
    '''Fit the PSF given an image file name'''

    if Verbose:
        verb = 'yes'
    else:
        verb = 'no'
    PSFisGood = False

    CoordsFile = ImageFile + '.coo'
    MagFile = ImageFile + '.mag'
    PSTFile = ImageFile + '.pst'
    PSFFile = ImageFile + '.psf'
    OpstFile = ImageFile + '.opst'
    GroupFile = ImageFile + '.group'
    SeeFile = ImageFile + '.see'


    while not PSFisGood:
        system('rm {} {} {} {} {} {} {}'.format(PSFFile + '.fits', CoordsFile, MagFile, PSTFile, OpstFile, GroupFile, SeeFile))

        try:
            # generate star catalog using daofind
            iraf.noao()
            iraf.digiphot()
            iraf.daophot(_doprint = 0)
            iraf.datapars.datamax = MaxCount
            iraf.datapars.fwhmpsf = FWHM
            iraf.datapars.sigma = Noise
            iraf.findpars.threshold = 5.
            iraf.datapars.datamin = 0
            iraf.datapars.datamax = MaxCount
            iraf.daofind(ImageFile, output = CoordsFile, verify = 'no', display = 'no', verbose = verb)

            # uncomment the raw_input line if daofind adds stars that do not exist in catalog
            # this gives you time to manually remove nonexistent stars that cause a bad psf fit
            # this is temporary until daofind works better with images coadded with swarp
            raw_input('Manually edit .coo file now if necessary; Press enter to continue ')

            # do aperture photometry
            a1, a2, a3, a4, = int(FWHM + 0.5), int(FWHM * 2 + 0.5), int(FWHM * 3 + 0.5), int(FWHM * 4 + 0.5)
            iraf.photpars.apertures = '{0},{1},{2}'.format(a2, a3, a4)
            iraf.centerpars.calgori = 'centroid'
            iraf.fitskypars.salgori = 'mode'
            iraf.fitskypars.annulus = 10
            iraf.fitskypars.dannulu = 10
            iraf.phot(ImageFile, CoordsFile, MagFile, verify = 'no', verbose = verb)

            # select PSF stars
            iraf.daopars.fitrad = a1
            iraf.daopars.nclean = 4
            iraf.daopars.varorder = 0
            iraf.daopars.recenter = 'yes'
            iraf.pstselect(ImageFile, MagFile, PSTFile, maxnpsf = 50, verify = 'no', verbose = verb)

            # make PSF
            iraf.psf(ImageFile, MagFile, PSTFile, PSFFile, OpstFile, GroupFile, verify = 'no', verbose = verb, interactive = 'no')

            # show psf to user for approval
            if Show:
                system ('rm {}'.format(SeeFile + '.fits'))
                iraf.seepsf(PSFFile, SeeFile)
                iraf.surface(SeeFile)
                PSFisGoodyn = raw_input('GoodPSF? y/n: ')
                if PSFisGoodyn == 'y':
                    PSFisGood = True
                else:
                    FWHMguess = raw_input('New FWHM: [{}] '.format(FWHM))
                    Noiseguess = raw_input('New Noise: [{}] '.format(Noise))
                    if FWHMguess != '':
                        FWHM = float(FWHMguess)
                    if Noiseguess != '':
                        Noise = float(Noiseguess)

            else:
                break

        except:
            print 'PSF fitting failed; try again with different parameters'
            FWHM = float(raw_input('New FWHM: '))
            Noise = float(raw_input('New Noise: '))

    return PSFFile


def gamma(beta, gammaPrime, sn, sr):
    '''Convert params in fourier space to image space'''

    return gammaPrime * np.sqrt(sn ** 2 + beta ** 2 * sr ** 2)

def InteractiveFit(N, R, Pn, Pr, sn, sr, NSaturation = None, RSaturation = None):
    N_hat = np.fft.fft2(N)    
    R_hat = np.fft.fft2(R)
    Pr_hat = np.fft.fft2(Pr)
    Pn_hat = np.fft.fft2(Pn)

    GoodFit = False
    beta = 1.
    gamma = 0.

    while not GoodFit:

        SqrtDen = np.sqrt(sn ** 2 * abs(Pr_hat) ** 2 + beta ** 2 * sr ** 2 * abs(Pn_hat) ** 2)
        Dn_hat = Pr_hat * N_hat / SqrtDen
        Dr_hat = Pn_hat * R_hat / SqrtDen
        Dn = np.real(np.fft.ifft2(Dn_hat))
        Dr = np.real(np.fft.ifft2(Dr_hat))
        DnFlatten = Dn.flatten()
        DrFlatten = Dr.flatten()

        # remove saturated pixels
        DnGoodPix = CleanDataIndex(DnFlatten, SatCount = NSaturation)
        DrGoodPix = CleanDataIndex(DrFlatten, SatCount = RSaturation)
        GoodPixInCommon = np.intersect1d(DnGoodPix, DrGoodPix)
        DnFlatten = DnFlatten[GoodPixInCommon]
        DrFlatten = DrFlatten[GoodPixInCommon]

        plt.plot(DrFlatten, DnFlatten, 'bo', DrFlatten, np.polyval([beta, gamma], DrFlatten), 'r-')
        plt.show(block = False)

        GoodFityn = raw_input('Good Fit? y/[n]/b: ')
        if GoodFityn == 'y':
            GoodFit = True
        elif GoodFityn == 'b':
            print 'Bad dataset: check your images for cosmic rays, saturated stars, registration errors and bad PSFs'
            quit()
        else:
            beta = float(input('New beta: [{}] '.format(beta)))
            gamma = float(input('New gamma: [{}] '.format(gamma)))
        plt.show()
    return beta, gamma, 0, 0 # functions expect error output


def IterativeSolve(N, R, Pn, Pr, sn, sr, NSaturation = None, RSaturation = None, interactive = False):
    '''Solve for linear fit iteratively'''

    BetaTolerance = 1e-8
    GammaPrimeTolerance = 1e-8
    beta = 1.
    gammaPrime = 0.
    beta0 = 10e5
    gammaPrime0 = 10e5
    i = 0
    MaxIteration = 10

    N_hat = np.fft.fft2(N)    
    R_hat = np.fft.fft2(R)
    Pr_hat = np.fft.fft2(Pr)
    Pn_hat = np.fft.fft2(Pn)

    while abs(beta - beta0) > BetaTolerance or abs(gammaPrime - gammaPrime0) > GammaPrimeTolerance:

        SqrtDen = np.sqrt(sn ** 2 * abs(Pr_hat) ** 2 + beta ** 2 * sr ** 2 * abs(Pn_hat) ** 2)
        Dn_hat = Pr_hat * N_hat / SqrtDen
        Dr_hat = Pn_hat * R_hat / SqrtDen
        Dn = np.real(np.fft.ifft2(Dn_hat))
        Dr = np.real(np.fft.ifft2(Dr_hat))
        DnFlatten = Dn.flatten()
        DrFlatten = Dr.flatten()

        # remove saturated pixels
        DnGoodPix = CleanDataIndex(DnFlatten, SatCount = NSaturation, Significance = 3.)
        DrGoodPix = CleanDataIndex(DrFlatten, SatCount = RSaturation, Significance = 3.)
        GoodPixInCommon = np.intersect1d(DnGoodPix, DrGoodPix)
        DnFlatten = DnFlatten[GoodPixInCommon]
        DrFlatten = DrFlatten[GoodPixInCommon]

        beta0, gammaPrime0 = beta, gammaPrime

        x = sm.add_constant(DrFlatten)
        y = DnFlatten

        
        RobustFit = sm.RLM(y, x).fit()
        Parameters = RobustFit.params
        Errors = RobustFit.bse
        beta = Parameters[1]
        gammaPrime = Parameters[0]
        betaError = Errors[1]
        gammaPrimeError = Errors[0]

        if i == MaxIteration: break
        i += 1
        print 'Iteration {}:'.format(i)
        print 'Beta = {0}, gamma = {1}'.format(beta, gamma(beta, gammaPrime, sn, sr))

    print 'Fit done in {} iterations'.format(i)

    plt.plot(DrFlatten, DnFlatten, 'bo', DrFlatten, RobustFit.fittedvalues, 'r-')
    plt.show()

    print 'Beta = ' + str(beta)
    print 'Gamma = ' + str(gamma(beta, gammaPrime, sn, sr))
    return beta, gamma(beta, gammaPrime, sn, sr), betaError, gammaPrimeError


def Normalize(Image):
    '''Normalize to sum = 1'''

    return Image / np.sum(Image)


def PadPSF(PSF, shape):
    '''Pad PSF and center at (0,0)'''

    p = PSF.shape
    s = shape
    PSF_extended = np.zeros(s)
    PSF_extended[s[0] / 2 - p[0] / 2 - 1:s[0] / 2 + p[0] / 2, s[1] / 2 - p[1] / 2 - 1:s[1] / 2 + p[1] / 2] = PSF
    PSF_extended = np.roll(PSF_extended, s[0] / 2, 0)
    PSF_extended = np.roll(PSF_extended, s[1] / 2, 1)

    return PSF_extended


def RemoveSaturatedFromCatalog(catalog):
    '''Remove stars with saturated pixels from catalog'''

    SafeIndex = np.where(catalog[:7].flatten() == 0)
    SafeSource = catalog[SafeIndex]
    return SafeSource


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform optimal image subtraction')
    parser.add_argument('-N', dest = 'input', help = 'New background subtracted image')
    parser.add_argument('-R', dest = 'template', help = 'Reference background subtracted image')
    parser.add_argument('-n', dest = 'input_PSF', default = '', help = 'New PSF')
    parser.add_argument('-r', dest = 'template_PSF', default = '', help = 'Reference PSF')
    parser.add_argument('-o', dest = 'output', help = 'Output image')
    parser.add_argument('--NewBackground', dest = 'NewBackground', help = 'New image with background')
    parser.add_argument('--RefBackground', dest = 'RefBackground', help = 'Reference image with background')
    parser.add_argument('--type', dest = 'type', default = 'D', help = 'Subtraction type')
    parser.add_argument('--threshold', dest = 'threshold', help = 'Sigma threshold for detection')
    
    parser.add_argument('--interactive', dest = 'interactive', action = 'store_true', default = False, help = 'Fit beta and gamma interactively')
    parser.add_argument('--verbose', dest = 'verbose', action = 'store_true', default = False, help = 'Show IRAF output in PSF fitting')
    parser.add_argument('-b', dest = 'beta', default = 99999, help = 'Gain matching parameter beta')
    parser.add_argument('-g', dest = 'gamma', default = 99999, help = 'Gain matching parameter gamma')
    parser.add_argument('--NewNoise', dest = 'NewNoise', default = 5, help = 'Standard deviation of new image background')
    parser.add_argument('--RefNoise', dest = 'RefNoise', default = 5, help = 'Standard deviation of ref image background')
    parser.add_argument('--flatten', dest = 'flatten', action = 'store_true', default = False, help = 'Fix Background locally')
    parser.add_argument('--normalize', default = '', dest = 'normalize', help = 'Normalize to which image')
    args = parser.parse_args()

    d = {
        'Flatten': args.flatten, 
        'NewPSF': args.input_PSF, 
        'RefPSF': args.input_PSF, 
        'Interactive': args.interactive, 
        'Beta': float(args.beta), 
        'Gamma': float(args.gamma), 
        'Verbose': args.verbose, 
        'NewBackground': args.NewBackground, 
        'RefBackground': args.RefBackground, 
        'NewNoise': float(args.NewNoise), 
        'RefNoise': float(args.RefNoise)
        }

    if args.type == 'D':
        OptimalSubtraction(args.input, args.template, d).SaveD(args.output, args.normalize)
    elif args.type == 'Scorr':
        OptimalSubtraction(args.input, args.template, d).SaveScorr(args.output, args.normalize)
    elif args.type == 'S':
        OptimalSubtraction(args.input, args.template, d).SaveS(args.output, args.normalize)
    elif args.type == 'Thresh':
        OptimalSubtraction(args.input, args.template, d).SaveScorrThreshold(args.output, normalize = args.normalize, Thresh = args.threshold)
    elif args.type == 'Find':
        OptimalSubtraction(args.input, args.template, d).FindTransient()
    else:
        pass

