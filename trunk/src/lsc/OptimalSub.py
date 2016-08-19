#!/usr/bin/env python
import numpy as np
import scipy
from astropy.io import fits
from pyraf import iraf
import time
from os import system
import argparse
from lsc.lscastrodef import sextractor
from lsc.lscastrodef import crossmatchxy
import statsmodels.api as sm

import matplotlib.pyplot as plt

class OptimalSubtraction:
    '''Implementation of proper image subtraction from Zackey, Ofek, Gal-Yam 2016
       see https://arxiv.org/abs/1601.02655 for details'''

    def __init__(self, Nf, Rf, Pnf, Prf, ArgsDict = {}):
        '''Take filenames and turn them into arrays for input into other functions'''

        self.ArgsDict = ArgsDict
        self.Nf, self.Rf = Nf, Rf
        self.Pnf, self.Prf = Pnf, Prf
        self.Pn = ExtractPSF(self.Pnf)
        self.Pr = ExtractPSF(self.Prf)
        self.NBackground, self.RBackground, self.Pn, self.Pr = CropImages(fits.getdata(self.Nf), fits.getdata(self.Rf), self.Pn, self.Pr)
        self.sn = np.std(self.NBackground)
        self.sr = np.std(self.RBackground)



    def D(self, normalize = '', diagnostic = False):
        '''Calculate proper subtraction image and normalize to zero point of reference or target'''

        NBackground, RBackground, Pn, Pr, sn, sr = self.NBackground, self.RBackground, self.Pn, self.Pr, self.sn, self.sr


        if diagnostic:
            print 'Saving Deconvolved images'
            N_hat = np.fft.fft2(self.NBackground)    
            R_hat = np.fft.fft2(self.RBackground)
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

        self.CheckGainMatched()

        # calculate D
        R = np.copy(RBackground)
        N = np.subtract(NBackground, self.gamma)

        N_hat = np.fft.fft2(N)
        R_hat = np.fft.fft2(R)
        Pr_hat = np.fft.fft2(Pr)
        Pn_hat = np.fft.fft2(Pn)

        beta = self.beta
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

    def S(self):
        '''Calculate S array'''

        self.CheckGainMatched()
        try:
            S_hat = self.Fd_ * np.fft.fft2(self.D_) * np.conj(np.fft.fft2(self.Pd_))
            S = np.fft.ifft2(S_hat)
        except AttributeError:
            S_hat = self.Fd() * np.fft.fft2(self.D()) * np.conj(np.fft.fft2(self.Pd()))
            S = np.fft.ifft2(S_hat)
        self.S_ = np.real(S)
        return np.real(S)

    def Flux(self):
        '''Calculate transient Flux'''
        self.CheckGainMatched()
        try:
            Flux = self.S_ / self.Fs_
        except AttributeError:
            Flux = self.S() / self.Fs()
        self.Flux_ = Flux
        return Flux


    def Snoise(self):
        '''Calculate the noise image for Scorr'''
        # this whole function needs optimization
        self.CheckGainMatched()

        NBackground, RBackground, Pn, Pr, sn, sr = self.NBackground, self.RBackground, self.Pn, self.Pr, self.sn, self.sr
        R = np.copy(RBackground)
        N = np.subtract(NBackground, self.gamma)

        N_hat = np.fft.fft2(N)    
        R_hat = np.fft.fft2(R)
        Pr_hat = np.fft.fft2(Pr)
        Pn_hat = np.fft.fft2(Pn)

        beta = self.beta
        Den = sr ** 2 * beta ** 2 * abs(Pn_hat) ** 2 + sn ** 2 * abs(Pr_hat) ** 2

        Kn_hat = (np.conj(Pn_hat) * abs(Pr_hat) ** 2) / Den
        Kr_hat = (np.conj(Pr_hat) * abs(Pn_hat) ** 2) / Den
        Kn2_hat = np.fft.fft2(np.fft.ifft2(beta * Kn_hat) ** 2)
        Kr2_hat = np.fft.fft2(np.fft.ifft2(beta ** 2 * Kr_hat) ** 2)

        try:
            VarFilename = self.ArgsDict['nVariance']
            EpsN = fits.getdata(VarFilename)
            V_Sn = np.fft.ifft2(Kr2_hat * np.fft.fft2(EpsN))
        except KeyError:
            V_Sn = np.fft.ifft2(Kn2_hat * N_hat)

        try:
            VarFilename = self.ArgsDict['rVariance']
            EpsR = fits.getdata(VarFilename)
            V_Sr = np.fft.ifft2(Kr2_hat * np.fft.fft2(EpsN))
        except KeyError:
            V_Sr = np.fft.ifft2(Kr2_hat * R_hat)

        f = open('transform')
        for line in f:
            if 'xrms' in line:
                xrms = float(line.split()[1])
            elif 'yrms' in line:
                yrms = float(line.split()[1])
        f.close()

        GradNy, GradNx = np.gradient(np.fft.ifft2(beta * Kn_hat * N_hat))
        GradRy, GradRx = np.gradient(np.fft.ifft2(beta ** 2 * Kr_hat * R_hat))

        Vr_ast = xrms ** 2 * GradRx ** 2 + yrms ** 2 * GradRy ** 2
        Vn_ast = xrms ** 2 * GradNx ** 2 + yrms ** 2 * GradNy ** 2

        Snoise = V_Sn + V_Sr #+ Vr_ast + Vn_ast
        self.Snoise_ = np.real(Snoise)
        return np.real(Snoise)


    def Scorr(self):
        '''Calculate Scorr'''
        try:
            Scorr = self.S_ / self.Snoise_
        except AttributeError:
            Scorr = self.S() / np.sqrt(self.Snoise())
        self.Scorr_ = np.real(Scorr)
        return np.real(Scorr)

    def Fs(self):
        '''Calculate flux based zeropoint of S'''
        self.CheckGainMatched()
        beta = self.beta
        Pr_hat = np.fft.fft2(self.Pr)
        Pn_hat = np.fft.fft2(self.Pn)
        sn, sr = self.sn, self.sr
        Den = sr ** 2 * beta ** 2 * abs(Pn_hat) ** 2 + sn ** 2 * abs(Pr_hat) ** 2
        Fs = np.sum(beta ** 2 * abs(Pn_hat) * abs(Pr_hat) / Den)
        self.Fs_ = Fs
        return Fs

    def Fd(self):
        '''Calculate the flux based zero point of D'''
        self.CheckGainMatched()
        Fd = self.beta / np.sqrt(self.sn**2 + self.sr**2*self.beta**2)
        self.Fd_ = np.real(Fd)
        return np.real(Fd)

    def Pd(self):
        self.CheckGainMatched()
        Pn_hat = np.fft.fft2(self.Pn)
        Pr_hat = np.fft.fft2(self.Pr)
        Den = self.sr ** 2 * self.beta ** 2 * abs(Pn_hat) ** 2 + self.sn ** 2 * abs(Pr_hat) ** 2
        SqrtDen = np.sqrt(Den)
        Pd_hat = self.beta * Pr_hat * Pn_hat / (self.Fd() * SqrtDen)
        Pd = np.fft.ifft2(Pd_hat)
        self.Pd_ = Pd

        return np.real(Pd)

    def CheckGainMatched(self):
        '''Check if the gain matching parameters have been found'''

        # if user hasn't supplied beta and gamma, fit them
        try:
            self.beta
            self.gamma
        except AttributeError:
            try:
                self.beta = self.ArgsDict['beta']
                self.gamma = self.ArgsDict['gamma']
            except KeyError:
                print 'Gain matching not done; fitting manually'
                #sources = np.array(sextractor(self.Nf)).transpose()
                #source_index = np.argsort(sources[:, 4])
                #source_sorted = sources[source_index]

                # add support for local gain matching once fitting improves
                self.beta, self.gamma = self.MatchGain()#catalog = source_sorted)



    def SaveDToDB(self, Df, normalize):
        '''Calculate and save proper subtraction image to database'''
        self.Df = Df
        self.DNormalize = self.D(normalize)
        hdu = fits.PrimaryHDU(np.real(self.DNormalize))
        hdu.header=fits.getheader(self.Nf)
        hdu.header['PHOTNORM'] = normalize
        hdu.writeto(self.Df, clobber = True)

    def MatchGain(self, catalog = 'center'):
        '''Call gain matching function'''
        beta, gamma, betaError, gammaError = FitB(self.NBackground, self.RBackground, self.Pn, self.Pr, self.sn, self.sr, catalog = catalog) # change to include catalog later
        return beta, gamma

    def SaveImageToWD(self):
        '''Save various images to working directory (testing only)'''
        Images = {'S.fits': self.S, 'Snoise.fits': self.Snoise, 'Scorr.fits': self.Scorr, 'D.fits': self.D, 'Flux.fits': self.Flux, 'Pd.fits': self.Pd}

        for element in Images:
            hdu = fits.PrimaryHDU(np.real(Images[element]))
            hdu.header=fits.getheader(self.Nf)
            hdu.writeto(element, clobber = True)


def AlignPSF(PSF, Coords):
    '''Roll PSF to correct for registration errors'''
    PSF = scipy.ndimage.interpolation.shift(PSF, Coords)
    return PSF
    

def RemoveSaturatedFromCatalog(catalog):
    '''Remove stars with saturated pixels from catalog'''

    SafeIndex = np.where(catalog[:7].flatten() == 0)
    SafeSource = catalog[SafeIndex]
    return SafeSource
        

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
    iraf.daophot()
    iraf.seepsf(Pf,'temp.psf.fits')
    PSF = fits.open('temp.psf.fits')[0].data
    system('rm temp.psf.fits')

    # normalize PSF
    PSF=PSF / np.sum(PSF)

    return PSF

def CropImages(N, R, Pn, Pr):
    '''Open and crop images and PSFs to the same size'''

    if N.shape[0] <= R.shape[0] & N.shape[1] <= R.shape[1]:
        s = N.shape
    elif N.shape[0] > R.shape[0] & N.shape[1] > R.shape[1]:
        s = R.shape
    else:
        s = [np.min([N.shape[0], R.shape[0]]), np.min([N.shape[1], R.shape[1]])]

    N = N[:s[0], :s[1]]
    R = R[:s[0], :s[1]]

    Pn_extended = PadPSF(Pn, N.shape)
    Pr_extended = PadPSF(Pr, N.shape)

    return N, R, Pn_extended, Pr_extended

def PadPSF(PSF, shape):
    '''Pad PSF and center at (0,0)'''

    p = PSF.shape
    s = shape
    PSF_extended = np.zeros(s)
    PSF_extended[s[0] / 2 - p[0] / 2 - 1:s[0] / 2 + p[0] / 2, s[1] / 2 - p[1] / 2 - 1:s[1] / 2 + p[1] / 2] = PSF
    PSF_extended = np.roll(PSF_extended, s[0] / 2, 0)
    PSF_extended = np.roll(PSF_extended, s[1] / 2, 1)

    return PSF_extended


def FitB(NBackground, RBackground, Pn, Pr, sn, sr, catalog = 'center'):
    '''Fit gain matching (beta) and background (gamma) parameters'''

    if catalog is 'center':
        center, d = NBackground.shape[0] / 2, 2000
        a, b = center - d, center + d
        coords = [[a, b, a, b]]

    elif catalog == 'grid':
        EdgeBuffer = 500
        resolution =  20
        d = 200
        xGrid = np.linspace(EdgeBuffer, NBackground.shape[0] - EdgeBuffer, num = resolution)
        yGrid = np.linspace(EdgeBuffer, NBackground.shape[1] - EdgeBuffer, num = resolution)
        x, y = np.meshgrid(xGrid, yGrid)
        centers = np.array([x.flatten(), y.flatten()]).transpose()
        coords = [[int(c[0] - d), int(c[0] + d), int(c[1] - d), int(c[1] + d)] for c in centers]

    else:
        # this doesn't work well; fix later or remove
        d = 100
        coords = [[int(star[0] - d), int(star[0] + d), int(star[1] - d), int(star[1] + d)] 
                  for star in catalog if star[0] > d and star[0] < NBackground.shape[0] - d and star[1] > d and star[1] < NBackground.shape[1] - d]
    g = []
    b = []
    be = []
    ge = []
    for star in coords:
        N = NBackground[star[0]:star[1], star[2]:star[3]]
        R = RBackground[star[0]:star[1], star[2]:star[3]]
        Pn = Pn[0:2 * d, 0:2 * d]
        Pr = Pr[0:2 * d, 0:2 * d]

        print 'x: {0}, y: {1}'.format(star[0] + d, star[2] + d)
        # iteratively solve for linear fit of beta and gamma
        beta, gamma, betaError, gammaError = IterativeSolve(N, R, Pn, Pr, sn, sr)
        if catalog == 'grid':
            g.append(gamma)
            b.append(beta)
            ge.append(gammaError)
            be.append(betaError)

    if catalog == 'center':
        return beta, gamma, betaError, gammaError

    elif catalog == 'grid':
        shape = NBackground.shape
        betaMatrix = ListToArray(b, shape)
        betaErrorMatrix = ListToArray(be, shape)
        gammaMatrix = ListToArray(g, shape)
        gammaErrorMatrix = ListToArray(ge, shape)
        return betaMatrix, gammaMatrix, betaErrorMatrix, gammaErrorMatrix

    else:
        return np.mean(beta), np.mean(gamma), np.mean(betaError), np.mean(gammaError)

def ListToArray(List, shape):
    '''Take a list and turn it into an interpolated array'''
    resolution = int(np.sqrt(len(List)))
    parMatrix = np.array(List).reshape(2 * [resolution]).transpose()
    parMatrixFull = scipy.misc.imresize(parMatrix, shape, interp = 'nearest')
    return parMatrixFull

def gamma(beta, gammaPrime, sn, sr):
    return gammaPrime * np.sqrt(sn ** 2 + beta ** 2 * sr ** 2)

def IterativeSolve(N, R, Pn, Pr, sn, sr):
    '''Solve for linear fit iteratively'''

    BetaTolerance = 1e-10
    GammaPrimeTolerance = 1e-10
    beta = 1.
    gammaPrime = 0.
    beta0 = 10e5
    gammaPrime0 = 10e5
    i = 0
    MaxIteration = 5

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

        beta0, gammaPrime0 = beta, gammaPrime

        x = sm.add_constant(DrFlatten)#[index])
        y = DnFlatten#[index]

        RobustFit = sm.OLS(y,x).fit()#RLM(y, x).fit()
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

    #plt.plot(DrFlatten, DnFlatten, 'bo', DrFlatten, RobustFit.fittedvalues, 'r-')
    #plt.show()
    print 'Beta = ' + str(beta)
    print 'Gamma = ' + str(gamma(beta, gammaPrime, sn, sr))
    return beta, gamma(beta, gammaPrime, sn, sr), betaError, gammaPrimeError


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform optimal image subtraction')
    parser.add_argument('-i', dest = 'input', help = 'Input image')
    parser.add_argument('-t', dest = 'template', help = 'Template image')
    parser.add_argument('-o', dest = 'output', help = 'Output image')
    args = parser.parse_args()

    OptimalSubtraction(args.input, args.template).SaveDtoDB(args.output)


