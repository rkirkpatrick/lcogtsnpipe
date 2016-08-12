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
    def __init__(self, Nf, Rf, Pnf, Prf, Df, normalize):

        self.normalize = normalize
        self.Df = Df
        self.Nf, self.Rf = Nf, Rf
        self.Pnf, self.Prf = Pnf, Prf
        #self.RegisterImages()
        self.Pn = self.ExtractPSF(self.Pnf)
        self.Pr = self.ExtractPSF(self.Prf)
        self.CropImages()
        self.sn = np.std(self.Nnbs)
        self.sr = np.std(self.Rnbs)
        self.FitB() 
        self.N = np.subtract(self.Nnbs, self.gamma)
        self.R = np.copy(self.Rnbs)

        self.Subtract()

        self.SaveImageToDB()
        #self.SaveImageToWD()


    def ExtractPSF(self, Pf):

        iraf.noao()
        iraf.digiphot()
        iraf.daophot()
        iraf.seepsf(Pf,'temp.psf.fits')
        P=fits.open('temp.psf.fits')[0].data
        system('rm temp.psf.fits')

        #normalize PSF
        P=P/np.sum(P)

        return P

    def CropImages(self):
        N, R, Pn, Pr = fits.getdata(self.Nf), fits.getdata(self.Rf), self.Pn, self.Pr
        if N.shape[0] <= R.shape[0] & N.shape[1] <= R.shape[1]:
            s = N.shape
        elif N.shape[0] > R.shape[0] & N.shape[1] > R.shape[1]:
            s = R.shape
        else:
            s = [np.min([N.shape[0], R.shape[0]]), np.min([N.shape[1], R.shape[1]])]

        N = N[:s[0], :s[1]]
        R = R[:s[0], :s[1]]
        p = Pn.shape

        Pn_ext = np.zeros(s)
        Pn_ext[s[0]/2-p[0]/2-1:s[0]/2+p[0]/2, s[1]/2-p[1]/2-1:s[1]/2+p[1]/2] = Pn
        Pn_ext = np.roll(Pn_ext, s[0]/2, 0)
        Pn_ext = np.roll(Pn_ext, s[1]/2, 1)

        p = Pr.shape

        Pr_ext = np.zeros(s)
        Pr_ext[s[0]/2 - p[0] / 2 - 1:s[0] / 2 + p[0] / 2, s[1] / 2 - p[1] /2 -1 :s[1] / 2 + p[1] / 2] = Pr
        Pr_ext = np.roll(Pr_ext, s[0]/2, 0)
        Pr_ext = np.roll(Pr_ext, s[1]/2, 1)

        self.Nnbs, self.Rnbs, self.Pn, self.Pr = N, R, Pn_ext, Pr_ext

    def FitB(self):
        s, d = self.Nnbs.shape[0] / 2, 500
        a, b = s - d, s + d
        N, R, Pn, Pr = self.Nnbs[a:b, a:b], self.Rnbs[a:b, a:b], self.Pn[0:2*d, 0:2*d], self.Pr[0:2*d, 0:2*d]
        sn, sr = self.sn, self.sr

        beta_tol = 1e-5
        gamma_tol = 1e-5
        beta = 1.
        gamma = 0.
        beta0 = 10e5
        gamma0 = 10e5
        i = 0
        b = []
        g = []
        maxiter = 4

        N_hat=np.fft.fft2(N)    
        R_hat=np.fft.fft2(R)
        Pr_hat=np.fft.fft2(Pr)
        Pn_hat=np.fft.fft2(Pn)

        #index = np.random.randint(0, N.size, size = int(1e5))

        while abs(beta - beta0) > beta_tol or abs(gamma - gamma0) > gamma_tol:

            Dn_hat = Pr_hat*N_hat/np.sqrt(sn**2*abs(Pr_hat**2)+beta**2*sr**2*abs(Pn_hat)**2)
            Dr_hat = Pn_hat*R_hat/np.sqrt(sn**2*abs(Pr_hat**2)+beta**2*sr**2*abs(Pn_hat)**2)
            Dn = np.real(np.fft.ifft2(Dn_hat))
            Dr = np.real(np.fft.ifft2(Dr_hat))
            Dnx = Dn.flatten()
            Drx = Dr.flatten()
            #intersect = np.intersect1d(np.where(abs(Dnx) < 10 * np.median(Dnx)), np.where(abs(Drx) < 10 * np.median(Drx)))
            #Dnx = Dnx[intersect]
            #Drx = Drx[intersect]

            beta0, gamma0 = beta, gamma

            x = sm.add_constant(Drx)#[index])
            y = Dnx#[index]

            fi = sm.RLM(y, x, M = sm.robust.norms.RamsayE()).fit()
            par = fi.params
            beta = par[1]
            gamma = par[0]

            i += 1
            print 'Iteration {}:'.format(i)
            print beta, gamma
            #plt.plot(Drx, Dnx, 'bo')#, Drx, fi.fittedvalues, 'r-')
            #plt.show()
            g.append(gamma)
            b.append(beta)

            if i == maxiter: 
                beta = np.mean(b)
                gamma = np.mean(g)
                break
        #plt.plot(Drx, Dnx, 'bo', Drx, fi.fittedvalues, 'r-')
        #plt.show()

        print 'Fit done in {} iterations'.format(i)
        #plt.plot(range(len(b)), b, 'g', range(len(g)), g, 'r')
        #plt.show()

        gamma = gamma * np.sqrt(sn**2 + beta ** 2 *sr**2)
        print 'Beta = ' + str(beta)
        print 'Gamma = ' + str(gamma)
        self.beta = beta
        self.gamma = gamma




    def Subtract(self, diagnostic = False):
        N, R, Pn, Pr, sn, sr = self.N, self.R, self.Pn, self.Pr, self.sn, self.sr
        print '\nDoing subtractions...'
        N_hat = np.fft.fft2(N)    
        R_hat = np.fft.fft2(R)
        Pr_hat = np.fft.fft2(Pr)
        Pn_hat = np.fft.fft2(Pn)

        if diagnostic:
            print 'Saving Deconvolved images'
            DConvN = np.real(np.fft.ifft2(N_hat * Pr_hat))
            hdu = fits.PrimaryHDU(DConvN)
            hdu.writeto('NDConvolve.fits', clobber = True)
            DConvR = np.real(np.fft.ifft2(R_hat * Pn_hat))
            hdu = fits.PrimaryHDU(DConvR)
            hdu.writeto('RDConvolve.fits', clobber = True)

        B = self.beta

        # correct for astrometric noise

        #f = open('transform')
        #for line in f:
        #    if 'xrms' in line:
        #        xrms = float(line.split()[1])
        #    elif 'yrms' in line:
        #        yrms = float(line.split()[1])
        #f.close()

        Den = sr ** 2 * B ** 2 * abs(Pn_hat) ** 2 + sn ** 2 * abs(Pr_hat) ** 2
        SqrtDen = np.sqrt(Den)

        #Kn = (np.conj(Pn_hat) * abs(Pr_hat) ** 2) / Den
        #Kr = (np.conj(Pr_hat) * abs(Pn_hat) ** 2) / Den
        #print np.sum(np.fft.ifft2(Pr_hat * N_hat)), np.sum(N), np.sum(np.fft.ifft2(B * Pn_hat * R_hat)), np.sum(B * R)
        #Fdn = 1 / np.sqrt(sn ** 2 + sr ** 2 * B ** 2)
        #Fdr = B / np.sqrt(sn ** 2 + sr ** 2 * B ** 2)
        #print Fdn, Fdr
        #D_hat = (Pr_hat * N_hat - B * Pn_hat * R_hat) / SqrtDen
        #self.D = np.fft.ifft2(D_hat)

        #self.D = self.N - np.fft.ifft2(B * R_hat * Pn_hat /  Pr_hat)
        Fd = B / np.sqrt(sn**2 + sr**2*B**2)
        Pd_hat = B * Pr_hat * Pn_hat / (Fd * np.sqrt(sn**2*abs(Pr_hat**2)+sr**2*B**2*abs(Pn_hat**2)))
        D = np.fft.ifft2((Pr_hat * N_hat - B * Pn_hat * R_hat) / SqrtDen)
        self.Pd = np.fft.ifft2(Pd_hat)
        print B, 1, Fd

        normalize = self.normalize
        if normalize == 'i':
            self.D = D * B / Fd
        elif normalize == 't':
            self.D = D / Fd
        else:
            self.D = D

        #del Pd_hat, Fd#, D_hat

        #self.S = np.fft.ifft2(B * Kn * N_hat - B ** 2 * Kr * R_hat)

        #Kr2 = np.fft.fft2(np.fft.ifft2(B ** 2 * Kr) ** 2)
        #Kn2 = np.fft.fft2(np.fft.ifft2(B * Kn) ** 2)

        #GradNy, GradNx = np.gradient(np.fft.ifft2(B * Kn * N_hat))
        #GradRy, GradRx = np.gradient(np.fft.ifft2(B ** 2 * Kr * R_hat))

        #Vr_ast = xrms ** 2 * GradRx ** 2 + yrms ** 2 * GradRy ** 2
        #Vn_ast = xrms ** 2 * GradNx ** 2 + yrms ** 2 * GradNy ** 2

        #NBackground_hat = np.fft.fft2(self.Nnbs)
        #RBackground_hat = np.fft.fft2(self.Rnbs)
        #V_Sn = np.fft.ifft2(Kn2 * NBackground_hat)
        #V_Sr = np.fft.ifft2(Kr2 * RBackground_hat)
        #self.Snoise = np.sqrt(V_Sn + V_Sr + Vr_ast + Vn_ast)

        #self.Fs = abs(np.sum(B ** 2 * abs(Pn_hat) ** 2 * abs(Pr_hat) ** 2) / Den)
        #self.Flux = self.S / self.Fs
        #self.Scorr = self.S / self.Snoise
        
        print '\nDone!'

    def SaveImageToWD(self):
        Images = {'S.fits': self.S, 'Snoise.fits': self.Snoise, 'Scorr.fits': self.Scorr, 'D.fits': self.D, 'Flux.fits': self.Flux, 'Pd.fits': self.Pd}

        for element in Images:
            hdu = fits.PrimaryHDU(np.real(Images[element]))
            hdu.header=fits.getheader(self.Nf)
            hdu.writeto(element, clobber = True)

    def SaveImageToDB(self):
        hdu = fits.PrimaryHDU(np.real(self.D))
        hdu.header=fits.getheader(self.Nf)
        hdu.header['PHOTNORM'] = self.normalize
        hdu.writeto(self.Df, clobber = True)

    def RegisterImages(self):
        #NCoords = sextractor(self.Nf)[:2]
        #RCoords = sextractor(self.Rf)[:2]
        Coords = np.genfromtxt('tmpcoo') # not sure whether to trust this
        print Coords
        NCoords = Coords[:, 0:2]
        RCoords = Coords[:, 2:4]
        m = scipy.optimize.minimize(AffineResidual, [1, 0, 0, 1, 0, 0], args = (NCoords, RCoords))
        print m

def AffineResidual(AffineParams, NCoords, RCoords):
    a, b, c, d, tx, ty = AffineParams
    t1 = np.sum(NCoords[0] - a * RCoords[0] - b * RCoords[1] - tx * np.ones(RCoords.shape))
    t2 = np.sum(NCoords[1] - c * RCoords[0] - d * RCoords[1] - ty * np.ones(RCoords.shape))
    return t1 ** 2 + t2 ** 2


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform optimal image subtraction')
    parser.add_argument('-i', dest = 'input', help = 'Input image')
    parser.add_argument('-t', dest = 'template', help = 'Template image')
    parser.add_argument('-o', dest = 'output', help = 'Output image')
    #parser.add_argument('--diagnostics', default = False, dest = 'diagnostic', action = 'store_true', help = 'Generate diagnostic images')
    #parser.add_argument('--interactive', default = False, dest = 'interactive', action = 'store_true', help = 'Fit B and Gamma interactively')
    args = parser.parse_args()

    OptimalSubtraction(args.input, args.template, args.output)








