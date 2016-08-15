#!/usr/bin/env python
import numpy as np
import scipy
from astropy.io import fits
from pyraf import iraf
import time
from os import system
import argparse
from lsc.util import display_image
from lsc.lscastrodef import sextractor
import statsmodels.api as sm

import matplotlib.pyplot as plt

class OptimalSubtraction:
    '''Implementation of proper image subtraction from Zackey, Ofek, Gal-Yam 2016
       see https://arxiv.org/abs/1601.02655 for details'''

    def __init__(self, Nf, Rf, Pnf, Prf):

        self.Nf, self.Rf = Nf, Rf
        self.Pnf, self.Prf = Pnf, Prf
        self.Pn = ExtractPSF(self.Pnf)
        self.Pr = ExtractPSF(self.Prf)
        self.NBackground, self.RBackground, self.Pn, self.Pr = CropImages(fits.getdata(self.Nf), fits.getdata(self.Rf), self.Pn, self.Pr)
        self.sn = np.std(self.NBackground)
        self.sr = np.std(self.RBackground)



    def D(self, normalize = 't', diagnostic = False):
        '''Calculate proper subtraction image and normalize to zero point of reference or target'''

        NBackground, RBackground, Pn, Pr, sn, sr = self.NBackground, self.RBackground, self.Pn, self.Pr, self.sn, self.sr
        print '\nDoing subtractions...'


        if diagnostic:
            print 'Saving Deconvolved images'
            DConvN = np.real(np.fft.ifft2(N_hat * Pr_hat))
            hdu = fits.PrimaryHDU(DConvN)
            hdu.writeto('NDConvolve.fits', clobber = True)
            DConvR = np.real(np.fft.ifft2(R_hat * Pn_hat))
            hdu = fits.PrimaryHDU(DConvR)
            hdu.writeto('RDConvolve.fits', clobber = True)

        try:
            beta, gamma = self.beta, self.gamma
        except AttributeError:
            self.beta, self.gamma = FitB(self.NBackground, self.RBackground, self.Pn, self.Pr, self.sn, self.sr)
            beta, gamma = self.beta, self.gamma

        R = RBackground
        N = np.subtract(NBackground, gamma)

        N_hat = np.fft.fft2(N)    
        R_hat = np.fft.fft2(R)
        Pr_hat = np.fft.fft2(Pr)
        Pn_hat = np.fft.fft2(Pn)

        Den = sr ** 2 * beta ** 2 * abs(Pn_hat) ** 2 + sn ** 2 * abs(Pr_hat) ** 2
        SqrtDen = np.sqrt(Den)

        Fd = beta / np.sqrt(sn**2 + sr**2*beta**2)
        Pd_hat = beta * Pr_hat * Pn_hat / (Fd * SqrtDen)

        self.Pd = np.fft.ifft2(Pd_hat)
        self.Fd = Fd

        D = np.fft.ifft2((Pr_hat * N_hat - beta * Pn_hat * R_hat) / SqrtDen)

        if normalize == 'i':
            DNormalize = D * B / Fd
        elif normalize == 't':
            DNormalize = D / Fd
        else:
            DNormalize = D

        return DNormalize        

    def SaveDToDB(self, Df, normalize):
        '''Calculate and save proper subtraction image to database'''
        self.Df = Df
        self.DNormalize = self.D(normalize)
        hdu = fits.PrimaryHDU(np.real(self.DNormalize))
        hdu.header=fits.getheader(self.Nf)
        hdu.header['PHOTNORM'] = normalize
        hdu.writeto(self.Df, clobber = True)

    def SaveImageToWD(self):
        '''Save various images to working directory (testing only)'''
        Images = {'S.fits': self.S, 'Snoise.fits': self.Snoise, 'Scorr.fits': self.Scorr, 'D.fits': self.D, 'Flux.fits': self.Flux, 'Pd.fits': self.Pd}

        for element in Images:
            hdu = fits.PrimaryHDU(np.real(Images[element]))
            hdu.header=fits.getheader(self.Nf)
            hdu.writeto(element, clobber = True)

def ExtractPSF(Pf):
    '''Extract PSF array from iraf PSF file'''

    iraf.noao()
    iraf.digiphot()
    iraf.daophot()
    iraf.seepsf(Pf,'temp.psf.fits')
    PSF = fits.open('temp.psf.fits')[0].data
    system('rm temp.psf.fits')

    #normalize PSF
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


def FitB(NBackground, RBackground, Pn, Pr, sn, sr):
    '''Fit gain matching (beta) and background (gamma) parameters'''

    s, d = NBackground.shape[0] / 2, 500
    a, b = s - d, s + d
    N, R, Pn, Pr = NBackground[a:b, a:b], RBackground[a:b, a:b], Pn[0:2 * d, 0:2 * d], Pr[0:2 * d, 0:2 * d]

    BetaTolerance = 1e-5
    GammaTolerance = 1e-5
    beta = 1.
    gamma = 0.
    beta0 = 10e5
    gamma0 = 10e5
    i = 0
    b = []
    g = []
    MaxIteration = 4

    N_hat = np.fft.fft2(N)    
    R_hat = np.fft.fft2(R)
    Pr_hat = np.fft.fft2(Pr)
    Pn_hat = np.fft.fft2(Pn)

    #index = np.random.randint(0, N.size, size = int(1e5))

    while abs(beta - beta0) > BetaTolerance or abs(gamma - gamma0) > GammaTolerance and i < MaxIteration:

        SqrtDen = np.sqrt(sn ** 2 * abs(Pr_hat) ** 2 + beta ** 2 * sr ** 2 * abs(Pn_hat) ** 2)
        Dn_hat = Pr_hat * N_hat / SqrtDen
        Dr_hat = Pn_hat * R_hat / SqrtDen
        Dn = np.real(np.fft.ifft2(Dn_hat))
        Dr = np.real(np.fft.ifft2(Dr_hat))
        DnFlatten = Dn.flatten()
        DrFlatten = Dr.flatten()
        #intersect = np.intersect1d(np.where(abs(DnFlatten) < 10 * np.median(DnFlatten)), np.where(abs(DrFlatten) < 10 * np.median(DrFlatten)))
        #DnFlatten = DnFlatten[intersect]
        #DrFlatten = DrFlatten[intersect]

        beta0, gamma0 = beta, gamma

        x = sm.add_constant(DrFlatten)#[index])
        y = DnFlatten#[index]

        RobustFit = sm.RLM(y, x, M = sm.robust.norms.RamsayE()).fit()
        Parameters = RobustFit.params
        beta = Parameters[1]
        gamma = Parameters[0]

        i += 1
        print 'Iteration {}:'.format(i)
        print 'Beta = {0}, gamma = {1}'.format(beta, gamma * np.sqrt(sn**2 + beta ** 2 *sr**2))

        g.append(gamma)
        b.append(beta)

    print 'Fit done in {} iterations'.format(i)

    plt.plot(DrFlatten, DnFlatten, 'bo', DrFlatten, RobustFit.fittedvalues, 'r-')
    plt.show()
    print 'Beta = ' + str(beta)
    print 'Gamma = ' + str(gamma * np.sqrt(sn**2 + beta ** 2 *sr**2))

    return beta, gamma



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform optimal image subtraction')
    parser.add_argument('-i', dest = 'input', help = 'Input image')
    parser.add_argument('-t', dest = 'template', help = 'Template image')
    parser.add_argument('-o', dest = 'output', help = 'Output image')
    args = parser.parse_args()

    OptimalSubtraction(args.input, args.template).SaveDtoDB(args.output)








