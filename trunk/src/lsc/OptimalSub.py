#!/usr/bin/env python
import numpy as np
from scipy.optimize import minimize
from scipy.stats import mode
from astropy.io import fits
from pyraf import iraf
import time
from os import system
import argparse
from lsc.util import display_image
from lsc.lscastrodef import sextractor
import statsmodels.api as sm

import matplotlib.pyplot as plt


def Subtract(Nf, Rf, Pnf, Prf, Df, diagnostic = False, interactive = False):
    #open fits files
    N_hdulist=fits.open(Nf)
    R_hdulist=fits.open(Rf)

    N=N_hdulist[0].data
    R=R_hdulist[0].data

    #fix and roll PSFs, crop the images to the same shape
    Pn=ExtractPSF(Pnf)
    Pr=ExtractPSF(Prf)

    N,R,Pn,Pr=CropImages(N,R,Pn,Pr)


    print 'Doing Fourier transforms...'

    sn=np.std(N)
    sr=np.std(R)

    print '1/4'
    N_hat=np.fft.fft2(N)    
    print '2/4'
    R_hat=np.fft.fft2(R)
    print '3/4'
    Pr_hat=np.fft.fft2(Pr)
    print '4/4'
    Pn_hat=np.fft.fft2(Pn)



    if diagnostic:
        system('rm {0} {1}'.format('NDConvolve.fits', 'RDConvolve.fits'))
        print 'Saving Deconvolved images'
        DConvN = np.real(np.fft.ifft2(N_hat * Pr_hat))
        hdu = fits.PrimaryHDU(DConvN)
        hdu.writeto('NDConvolve.fits', clobber = True)
        DConvR = np.real(np.fft.ifft2(R_hat * Pn_hat))
        hdu = fits.PrimaryHDU(DConvR)
        hdu.writeto('RDConvolve.fits', clobber = True)

    print 'Fitting B, gamma...'
    #fit B using robust fitting
    t0=time.time()

    B, gamma = FitB(N, R, Pn, Pr, sn, sr, interactive)

    print 'B = ' + str(B)
    print 'Gamma = ' + str(gamma)
    print 'Fit time: {} seconds'.format(str(round(time.time()-t0, 3)))
    print 'Subtracting images...'

    N=np.subtract(N,gamma)
    N_hat=np.fft.fft2(N)

    #subtract images
    D_hat=(Pr_hat*N_hat-B*Pn_hat*R_hat)/np.sqrt(sn**2*abs(Pr_hat**2)+sr**2*B**2*abs(Pn_hat**2))
    #D=np.real(np.fft.ifft2(D_hat))



        #import matplotlib.pyplot as plt
        #plt.scatter(DConvN.reshape(DConvN.size,1), DConvR.reshape(DConvR.size,1))
        #print np.polyfit(DConvN.reshape(1, DConvN.size), DConvR.reshape(1, DConvR.size), 1)
        #plt.show()

    Fd = B / np.sqrt(sn**2 + sr**2*B**2)
    Pd_hat = B * Pr_hat * Pn_hat / (Fd * np.sqrt(sn**2*abs(Pr_hat**2)+sr**2*B**2*abs(Pn_hat**2)))
    S_hat = Fd * D_hat * np.conj(Pd_hat)

    S = np.real(np.fft.ifft2(S_hat))


    hdu = fits.PrimaryHDU(S)
    hdu.header=N_hdulist[0].header
    hdu.writeto(Df)

    print 'Done!'
    return S



def FitB(N, R, Pn, Pr, sn, sr, interactive):
    simple = True

    if interactive:

        beta, gamma = [1., 0.]

        #a, b = np.shape(N)[0]/2 + 50, np.shape(N)[0]/2 - 50
        #N, R, Pn, Pr = N[a:b, a:b], R[a:b, a:b], Pn[a:b, a:b], Pr[a:b, a:b]

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
            #plt.imshow(eqn, cmap = 'gray', interpolation = 'none', vmin = np.percentile(eqn, 5), vmax = np.percentile(eqn, 95))
            #plt.show()

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

    elif simple:

        beta_tol = 0.01
        gamma_tol = 0.01
        beta = 1.
        gamma = 0.
        beta0 = 10e5
        gamma0 = 10e5
        i = 0
        b = []
        g = []
        maxiter = 7

        N_hat=np.fft.fft2(N)    
        R_hat=np.fft.fft2(R)
        Pr_hat=np.fft.fft2(Pr)
        Pn_hat=np.fft.fft2(Pn)

        while abs(beta - beta0) > beta_tol or abs(gamma - gamma0) > gamma_tol:
            i += 1
            print 'Iteration {}:'.format(i)


            Dn_hat = Pr_hat*N_hat/np.sqrt(sn**2*abs(Pr_hat**2)+beta**2*sr**2*abs(Pn_hat)**2)
            Dr_hat = Pn_hat*R_hat/np.sqrt(sn**2*abs(Pr_hat**2)+beta**2*sr**2*abs(Pn_hat)**2)
            Dn = np.real(np.fft.ifft2(Dn_hat)).flatten()
            Dr = np.real(np.fft.ifft2(Dr_hat)).flatten()
            intersect = np.intersect1d(np.where(Dn >  5 * np.std(Dn)), np.where(Dr > 5 * np.std(Dr)))
            Dn = Dn[intersect]
            Dr = Dr[intersect]

            index = np.random.randint(0, Dn.size, size = (Dn.size/4, 1))
                    
            beta0, gamma0 = beta, gamma

            #s = [int(c / 2) for c in Dn.shape]
            #d = 25
            #y = Dn[s[0] - d:s[0] + d, s[0] - d:s[0] + d].flatten()
            #x = sm.add_constant(Dr[s[0] - d:s[0] + d, s[0] - d:s[0] + d].flatten())
            x = sm.add_constant(Dr[index])
            y = Dn[index]

            
            fi = sm.RLM(y, x, M = sm.robust.norms.TrimmedMean(c=np.median(Dn - Dr) / 2)).fit()
            par = fi.params

            beta = par[1]
            gamma = par[0]
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


    else:
        out_fits = fits.PrimaryHDU(N)
        out_fits.writeto('TestFit.fits', clobber = True)
        xpix, ypix, fw, cl, cm, ell, bkg, fl = sextractor('TestFit.fits')
        system('rm TestFit.fits')
        coords = np.array([[int(xpix[i]), int(ypix[i]), cm[i]] for i in range(len(xpix))])
        coords = coords[np.argsort(coords[:,2])]
        bg = []
        d = 20
        Pr = Pr[:2*d,:2*d]
        Pn = Pn[:2*d,:2*d]
        Pr_hat = np.fft.fft2(Pr)
        Pn_hat = np.fft.fft2(Pn)

        for pair in coords:
            x, y = int(pair[0]), int(pair[1])
            border = range(0, 100)
            borderx = range(int(N.shape[0]) - 100, int(N.shape[0]))
            bordery = range(int(N.shape[1]) - 100, int(N.shape[1]))
            centerx = range(int(N.shape[0]) / 2 - 100, int(N.shape[0]) / 2 + 100)
            centery = range(int(N.shape[1]) / 2 - 100, int(N.shape[1]) / 2 + 100)
            forbiddenx = border + borderx + centerx
            forbiddeny = border + bordery + centery
            if x in forbiddenx or y in forbiddeny: continue

            Nsmall = N[x-d:x+d, y-d:y+d]
            Rsmall = R[x-d:x+d, y-d:y+d]

            N_hat = np.fft.fft2(Nsmall)    
            R_hat = np.fft.fft2(Rsmall)

            tol = .05
            beta = 1.
            gamma = 0.
            beta0 = 10e5
            gamma0 = 10e5
            i = 0
            b = []
            g = []
            maxiter = 20

            while abs(beta - beta0) > tol:
                g.append(gamma)
                b.append(beta)

                Dn_hat = Pr_hat*N_hat/np.sqrt(sn**2*abs(Pr_hat**2)+beta**2*sr**2*abs(Pn_hat)**2)
                Dr_hat = Pn_hat*R_hat/np.sqrt(sn**2*abs(Pr_hat**2)+beta**2*sr**2*abs(Pn_hat)**2)
                Dn = np.real(np.fft.ifft2(Dn_hat))
                Dr = np.real(np.fft.ifft2(Dr_hat))
                    

                beta0, gamma0 = beta, gamma
                par, cov = np.polyfit(Dr.flatten(), Dn.flatten(), 1, cov = True)
                beta, gamma = par
                i += 1
                if i == maxiter: break

            #print par, np.sqrt(cov[0,0]), np.sqrt(cov[1,1])
            #plt.scatter(Dr, Dn)
            #plt.show()
            sd = np.std(Dn - (beta * Dr + gamma))
            B_std = np.sqrt(cov[0,0])

            bg.append({'B': beta, 'gamma': gamma, 'std': sd, 'cov': B_std})

        sd_tol = np.percentile([i['std'] for i in bg], 10)
        B_std_tol = np.percentile([i['cov'] for i in bg], 10)
        #B_std_tol = 0.5
        print B_std_tol, sd_tol


        bg = [i for i in bg if i['std'] < sd_tol]
        bg = [i for i in bg if i['cov'] < B_std_tol]
        print len(bg)
        bg = sorted(bg, key = lambda k: k['cov'])
        #for i in bg: print i
        beta = [i['B'] for i in bg]
        gamma = [i['gamma'] for i in bg]
        sd = [i['std'] for i in bg]
        B_cov = [i['cov'] for i in bg]
        plt.scatter(sd, beta, color = 'red')
        plt.scatter(sd, gamma, color = 'blue')
        plt.scatter(beta, gamma, color = 'green')
        plt.scatter(B_cov, beta, color = 'black')
        plt.show()

        #beta_binned_data, beta_bins = np.histogram(beta, len(beta))
        #beta_index = np.argmax(beta_binned_data)
        #gamma_binned_data, gamma_bins = np.histogram(gamma, len(gamma))
        #gamma_index = np.argmax(gamma_binned_data)

        #plt.hist(beta, len(beta))
        #plt.show()
        #plt.hist(gamma, len(gamma))
        #plt.show()

        beta = np.median(beta)
        print np.sqrt(sn**2 + beta ** 2 *sr**2)
        gamma = np.median(gamma) * np.sqrt(sn**2 + beta ** 2 *sr**2)



    return beta, gamma



def CropImages(N,R,Pn,Pr):
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

    return N,R,Pn_ext,Pr_ext

def ExtractPSF(Pf):

    iraf.noao()
    iraf.digiphot()
    iraf.daophot()
    iraf.seepsf(Pf,'temp.psf.fits')
    P=fits.open('temp.psf.fits')[0].data
    system('rm temp.psf.fits')

    #normalize PSF
    P=P/np.sum(P)

    return P


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform optimal image subtraction')
    parser.add_argument('-i', dest = 'input', help = 'Input image')
    parser.add_argument('-t', dest = 'template', help = 'Template image')
    parser.add_argument('-o', dest = 'output', help = 'Output image')
    parser.add_argument('--diagnostics', default = False, dest = 'diagnostic', action = 'store_true', help = 'Generate diagnostic images')
    parser.add_argument('--interactive', default = False, dest = 'interactive', action = 'store_true', help = 'Fit B and Gamma interactively')
    args = parser.parse_args()

    Subtract(args.input, args.template, 
             args.input.replace('.fits', '.psf.fits'), 
             args.template.replace('.fits', '.psf.fits'), 
             args.output, 
             diagnostic = args.diagnostic,
             interactive = args.interactive)










