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
    def __init__(self, Nf, Rf, Df):

        self.Nf, self.Rf = Nf, Rf
        self.Pnf = self.Nf.replace('.fits', '.psf.fits')
        self.Prf = self.Rf.replace('.fits', '.psf.fits')
        self.Pn = self.ExtractPSF(self.Pnf)
        self.Pr = self.ExtractPSF(self.Prf)
        self.CropImages()
        self.sn = np.std(self.Nnbs)
        self.sr = np.std(self.Rnbs)
        FitB(self.Nnbs, self.Rnbs, self.Pn, self.Pr, self.sn, self.sr) 
        self.N = np.subtract(self.Nbns, self.gamma)
        self.R = np.copy(self.Rnbs)

        Subtract(self)
        self.Scorr = self.S / self.Snoise

        SaveImageToWD(self)
        SaveImageToDb(self)

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
        N, R, Pn, Pr = fits.getdata(self.Nf), fits.getdata(self.Nf), Pn, Pr
        if N.shape[0]<=R.shape[0] & N.shape[1]<=R.shape[1]:
            s=N.shape
        elif N.shape[0]>R.shape[0] & N.shape[1]>R.shape[1]:
            s=R.shape
        else:
            s=[np.min([N.shape[0],R.shape[0]]),np.min([N.shape[1],R.shape[1]])]

        N=N[:s[0],:s[1]]
        R=R[:s[0],:s[1]]
        p=Pn.shape

        Pn_ext=np.zeros(s)
        Pn_ext[s[0]/2-p[0]/2-1:s[0]/2+p[0]/2,s[1]/2-p[1]/2-1:s[1]/2+p[1]/2]=Pn
        Pn_ext=np.roll(Pn_ext,s[0]/2,0)
        Pn_ext=np.roll(Pn_ext,s[1]/2,1)

        p=Pr.shape

        Pr_ext=np.zeros(s)
        Pr_ext[s[0]/2-p[0]/2-1:s[0]/2+p[0]/2,s[1]/2-p[1]/2-1:s[1]/2+p[1]/2]=Pr
        Pr_ext=np.roll(Pr_ext,s[0]/2,0)
        Pr_ext=np.roll(Pr_ext,s[1]/2,1)

        self.N, self.R, self.Pn, self.Pr = N, R, Pn, Pr

    def FitB(self, N, R, Pn, Pr, sn, sr, interactive = False):
        N, R, Pn, Pr, sn, sr = self.N, self.R, self.Pn, self.Pr, self.sn, self.sr

        if interactive:
            beta, gamma = [1., 0.]

            N_hat=np.fft.fft2(N)    
            R_hat=np.fft.fft2(R)
            Pr_hat=np.fft.fft2(Pr)
            Pn_hat=np.fft.fft2(Pn)

            Dn_hat = Pr_hat*N_hat/np.sqrt(sn**2*abs(Pr_hat**2)+beta**2*sr**2*abs(Pn_hat)**2)
            Dr_hat = Pn_hat*R_hat/np.sqrt(sn**2*abs(Pr_hat**2)+beta**2*sr**2*abs(Pn_hat)**2)
            Dn = np.real(np.fft.ifft2(Dn_hat))
            Dr = np.real(np.fft.ifft2(Dr_hat))

            while True:
                eqn = Dn - (beta * Dr + gamma)

                print 'B = {0}, Background = {1}'.format(beta, gamma)
                import astropy.io.fits as f
                f.PrimaryHDU(eqn).writeto('TestFit.fits', clobber = True)
                print np.max(eqn), np.min(eqn)
                iraf.display('TestFit.fits', 1, zscale = 'no', zrange = 'yes', fill = True)
                system('rm TestFit.fits')

                CheckGamma = raw_input('Background ok? [y/n]: ')
                if CheckGamma is 'n':
                    gamma = float(raw_input('New background value: '))

                CheckB = raw_input('Scale ok? [y/n]: ')
                if CheckB is 'n':
                    beta = float(raw_input('New Scale value: '))
                print 'sn = {0}; sr = {1}'.format(sn,sr)
                if [CheckB, CheckGamma] == 2 * ['y']:
                    break
            gamma = gamma * np.sqrt(sn**2 + beta ** 2 *sr**2)

        else:

            beta_tol = 0.01
            gamma_tol = 0.001
            beta = 1.
            gamma = 0.
            beta0 = 10e5
            gamma0 = 10e5
            i = 0
            b = []
            g = []
            maxiter = 15

            N_hat=np.fft.fft2(N)    
            R_hat=np.fft.fft2(R)
            Pr_hat=np.fft.fft2(Pr)
            Pn_hat=np.fft.fft2(Pn)

            while abs(beta - beta0) > beta_tol or abs(gamma - gamma0) > gamma_tol:

                Dn_hat = Pr_hat*N_hat/np.sqrt(sn**2*abs(Pr_hat**2)+beta**2*sr**2*abs(Pn_hat)**2)
                Dr_hat = Pn_hat*R_hat/np.sqrt(sn**2*abs(Pr_hat**2)+beta**2*sr**2*abs(Pn_hat)**2)
                Dn = np.real(np.fft.ifft2(Dn_hat))
                Dr = np.real(np.fft.ifft2(Dr_hat))
                Dnx = Dn.flatten()
                Drx = Dr.flatten()
                #intersect = np.intersect1d(np.where(Dnx >  5 * np.std(Dnx)), np.where(Drx > 5 * np.std(Drx)))
                #Dnx = Dnx[intersect]
                #Drx = Drx[intersect]

                index = np.random.randint(0, Dnx.size, size = int(1e5))
                beta0, gamma0 = beta, gamma

                x = sm.add_constant(Drx[index])
                y = Dnx[index]

                print np.std(Dn - Dr)
                fi = sm.RLM(y, x, M = sm.robust.norms.TrimmedMean(c=np.std(Dn - Dr))).fit()
                par = fi.params

                beta = par[1]
                gamma = par[0]
                i += 1
                print 'Iteration {}:'.format(i)
                print beta, gamma
                #plt.plot(Dr.flatten()[index], Dn.flatten()[index], 'bo', Dr.flatten()[index], fi.fittedvalues, 'r-')
                #plt.show()
                g.append(gamma)
                b.append(beta)

                if i == maxiter: 
                    beta = np.mean(b)
                    gamma = np.mean(g)
                    break

            print 'Fit done in {} iterations'.format(i)
            plt.plot(range(len(b)), b, 'g', range(len(g)), g, 'r')
            plt.show()

            gamma = gamma * np.sqrt(sn**2 + beta ** 2 *sr**2)
        print 'Beta = ' + str(beta)
        print 'Gamma = ' + str(gamma)
        self.beta = beta
        self.gamma = gamma




    def Subtract(self, diagnostic = False):
        N, R, Pn, Pr, sn, sr = self.N, self.R, self.Pn, self.Pr, self.sn, self.sr
        print '1/4'
        N_hat = np.fft.fft2(N)    
        print '2/4'
        R_hat = np.fft.fft2(R)
        print '3/4'
        Pr_hat = np.fft.fft2(Pr)
        print '4/4'
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

        f = open('transform')
        for line in f:
            if 'xrms' in line:
                xrms = float(line.split()[1])
            elif 'yrms' in line:
                yrms = float(line.split()[1])
        f.close()

        Den = sr ** 2 * B ** 2 * abs(Pn_hat) ** 2 + sn ** 2 * abs(Pr_hat) ** 2 + 1e-8
        SqrtDen = np.sqrt(Den)
        Kn = (np.conj(Pn_hat) * abs(Pr_hat) ** 2) / Den
        Kr = (np.conj(Pr_hat) * abs(Pn_hat) ** 2) / Den

        D_hat = (Pr_hat * N_hat - B * Pn_hat * R_hat) / Den
        self.D = np.real(np.fft.ifft2(D_hat))
        del D_hat
        Fd = B / np.sqrt(sn**2 + sr**2*B**2)
        Pd_hat = B * Pr_hat * Pn_hat / (Fd * np.sqrt(sn**2*abs(Pr_hat**2)+sr**2*B**2*abs(Pn_hat**2)))
        self.Pd = np.fft.ifft2(Pd_hat)
        del Pd_hat

        self.S = np.fft.ifft2(B * Kn * N_hat - B ** 2 * Kr * R_hat)

        Kr2 = np.fft.fft2(np.fft.ifft2(B ** 2 * Kr) ** 2)
        Kn2 = np.fft.fft2(np.fft.ifft2(B * Kn) ** 2)

        GradNy, GradNx = np.gradient(np.fft.ifft2(B * Kn * N_hat))
        GradRy, GradRx = np.gradient(np.fft.ifft2(B ** 2 * Kr * R_hat))

        Vr_ast = xrms ** 2 * GradRx ** 2 + yrms ** 2 * GradRy ** 2
        Vn_ast = xrms ** 2 * GradNx ** 2 + yrms ** 2 * GradNy ** 2

        NBackground_hat = np.fft.fft2(Sub.Nnbs)
        RBackground_hat = np.fft.fft2(Sub.Rnbs)
        V_Sn = np.fft.ifft2(Kn2 * NBackground_hat)
        V_Sr = np.fft.ifft2(Kr2 * RBackground_hat)
        self.Snoise = np.sqrt(V_Sn + V_Sr + Vr_ast + Vn_ast)

        self.Fs = np.real(np.sum(B ** 2 * abs(Pn_hat) * abs(Pr_hat)) / Den)
        print 'Done!'

    def SaveImageToWD(self):
        Images = {'S.fits': self.S, 'Snoise.fits': self.Snoise, 'Scorr.fits': self.Scorr, 'D.fits': self.D, 'Fs.fits': self.Fs}

        for element in Images:
            hdu = fits.PrimaryHDU(Images[element])
            hdu.header=fits.getheader(self.Nf)
            hdu.writeto(element, clobber = True)

    def SaveImageToDB(self):
        hdu = fits.PrimaryHDU(np.real(Images[element]))
        hdu.header=fits.getheader(self.Nf)
        hdu.writeto(self.Df, clobber = True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform optimal image subtraction')
    parser.add_argument('-i', dest = 'input', help = 'Input image')
    parser.add_argument('-t', dest = 'template', help = 'Template image')
    parser.add_argument('-o', dest = 'output', help = 'Output image')
    #parser.add_argument('--diagnostics', default = False, dest = 'diagnostic', action = 'store_true', help = 'Generate diagnostic images')
    #parser.add_argument('--interactive', default = False, dest = 'interactive', action = 'store_true', help = 'Fit B and Gamma interactively')
    args = parser.parse_args()

    OptimalSubtraction(args.input, args.template, args.output)








