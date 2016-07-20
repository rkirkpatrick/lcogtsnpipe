#!/usr/bin/env python
"""
How to run it now:

sourceingest.py -e YYYYMMDD-YYYYMMDD -n OBJNAME 
runlsc.py -e YYYYMMDD-YYYYMMDD --obstype e92 e93

"""
import os
import sys

import MySQLdb
import numpy as np
from numpy import pi
import astropy.io.fits as pyfits
from astropy import wcs
from astropy.coordinates import Angle
import astropy.units as u
from pyraf import iraf

import lsc
from LCOGTingest import db_ingest

##############################################################################
#############################Ingest Functions#################################
##############################################################################
def choosemagnitude():
    if np.random.randint(2) == 0:
        return None
    else:
        i = np.random.uniform(15,20)
        return i

def comparefakesourcemagnitude(ofilep, magnitude):
    hdulist1 = pyfits.open(ofilep)
    header1 = hdulist1[0].header
    return(header1['FAKEMAG'] == magnitude)


##############################################################################
##############################Drop PSF########################################
##############################################################################
def gaussian2d(shape,x0,y0,sigma):
    """
    Creates an ndarray with a gaussian at x0, y0, with shape being the
    shape of the ndarray.
    
    The use of mgrid flips rows and columns for the creation of a 2d
    gaussian, so that's why in the first line x and y are switched.
    """
    y,x = np.mgrid[:shape[0],:shape[1]]
    G = ((1/(2*pi*sigma**2))*np.exp(-((x.astype(float)-x0)**2 +(y.astype(float)-y0)**2)/(2*float(sigma**2))))
    return G


def createpsf(shape,xpos,ypos,ifilep):
    from iraf import daophot

    ifileppsf = ifilep.replace(".fits",".psf.fits")
    iraf.noao.digiphot.daophot.seepsf(ifileppsf,"e93.psf.fits")
    blankarray = np.zeros(shape,dtype=float)
    HDU = pyfits.open("e93.psf.fits")

    residualdata = HDU[0].data
    residualshape = residualdata.shape
    ystart = int(ypos - residualshape[0]/2)
    xstart = int(xpos - residualshape[1]/2)
    blankarray[ystart:ystart+residualshape[0],xstart:xstart+residualshape[1]] = residualdata
    psf = blankarray

    HDU.close()
    del HDU
    lsc.util.delete("e93.psf.fits")

    return psf


def degs2coords(ra,dec,header):
    w = wcs.WCS(header)
    return w.all_world2pix(ra,dec,0)


#m =-2.5log(counts/exptime) + zp
#counts = amplitude * G.sum
def magnitude2amplitude(psf,magnitude,zeropoint,exptime):
    counts = float((10**((magnitude - zeropoint)/-2.5))*exptime)
    amplitude = float(counts / np.sum(psf,dtype=float))
    return amplitude


def FWHM2sigma(arcFWHM,header):
    """
    Converts arcsec FWHM to pixel sigma.
    """
    if 'PIXSCALE' in header:
        scale = float(header['PIXSCALE'])
    elif 'CCDSCALE' in header:
        if 'CCDXBIN' in header:
            scale = header['CCDSCALE'] * float(header['CCDXBIN'])
        elif 'CCDSUM' in header:
            scale = header['CCDSCALE'] * int(string.split(lsc.util.readkey3(hdr, 'CCDSUM'))[0])
	
    coordFWHM = arcFWHM/scale
    coordsig = coordFWHM/((2*(2*np.log(2))**.5))

    return coordsig


def source_drop(ifilep, ofilep, ra, dec, fwhm, zeropoint, magnitude=None,_psf=False):
    """
    Copies ifile into ofile. Drops a gaussian source into ofile, depending
    if magnitude exists. Then to all ofile 'FKSRC' is added to the header
    to describe the magnitude of the dropped source. If there is no
    magnitude, the 'FAKEMAG' header will be set as an empty string.
    """
    HDU = pyfits.open(ifilep)
    inimage = HDU[0].data
    hdr = HDU[0].header
    shape = inimage.shape
    if magnitude != None:
        xpos, ypos = degs2coords(ra,dec,hdr)	
        exptime = hdr["exptime"]
        coordsigma  = FWHM2sigma(fwhm,hdr)

        if _psf:
            psf = createpsf(shape,xpos,ypos,ifilep)
        else:
            psf = gaussian2d(shape, xpos, ypos, coordsigma)
        psf = psf * magnitude2amplitude(psf,magnitude,zeropoint,exptime)
        outimage = inimage + psf

    HDU[0].data = outimage
    HDU[0].header.set("FAKEMAG", magnitude,
                           "If #s, then magnitude. Else, no fake source")
    HDU.writeto(ofilep, clobber = True)
    HDU.close()
    del HDU

##############################################################################
description='Creates e92/e93 files with e92 files having a perfect gaussian '\
            'and e93 files having a psf based of the e91 file.'
usage='%(prog)s -e epoch [-f filename -m magnitude]'
version='%(prog)s 0.3'

if __name__ == "__main__":	
    from argparse import ArgumentParser
    parser = ArgumentParser(description=description,usage=usage,version=version)
    parser.add_argument("-f", "--filename", nargs='+',default=[],help="A "\
                        "list to search for specific filenames")
    parser.add_argument("-n", "--name", dest="name", default='', type=str,
                      help='-n object name   \t [%(default)s]')
    parser.add_argument("-m", "--magnitude",type=float, default=None, help=
                        "Use a specific magnitude when computing PSF")
    parser.add_argument("-e", "--epoch", dest="epoch", default='', type=str,
                      help='epoch to reduce  \t [%(default)s]')
    parser.add_argument("-F", "--force", dest="force", action="store_true")
    parser.add_argument("-p","--psf",dest="psf",action="store_true",
                        default=False, help="Use residuals when making PSF")
    args = parser.parse_args()
    _magnitude = args.magnitude
    _filename = args.filename
    _name = args.name
    _epoch = args.epoch
    _force = args.force
    _psf = args.psf

    if (_epoch == ''):
        sys.argv.append('--help')
        args = parser.parse_args()

    hostname,username, passwd, database = lsc.mysqldef.getconnection('lcogt2')
    conn = lsc.mysqldef.dbConnect(hostname, username, passwd, database)
    query_files = 'SELECT filename FROM photlco WHERE filename LIKE "%e91%"\
                   AND mag != 9999 AND filetype = 1 '

    epoch = "-e " + _epoch
    if '-' not in str(_epoch):
        query_files += ' AND dayobs = {0} '.format(_epoch)
    else:
        epoch1, epoch2 = _epoch.split('-')
        query_files += ' AND dayobs >= {0} AND dayobs <= {1} '.format(epoch1,epoch2)

    if _filename != []:
        query_files = lsc.mysqldef.queryfilenamelike(query_files,_filename)

    if _name != '':
        query_files += ' AND objname = "{0}" '.format(_name)
        name = ' -n ' + _name + ' '

    c = conn.cursor()
    c.execute(query_files)
    r = c.fetchall()
    rname = ()
    for row in r:
        rtuple = (row[0],)
        rname = rname + rtuple

    for ifile in rname:
        if _psf:
            ofile = ifile.replace('e91','e93')
        else:
            ofile = ifile.replace('e91','e92')
        if _magnitude == None:
            magnitude = choosemagnitude()
        else:
            magnitude = _magnitude

        query_join = 'SELECT * FROM photlco LEFT JOIN photlcoraw ON ' \
                     'photlco.filename = photlcoraw.filename WHERE ' \
                     'photlco.filename LIKE "%{0}%";'.format(ifile)

        conn.query(query_join)
        rcomb = conn.store_result()		
        rcomb = rcomb.fetch_row(how=2)
        rcomb = rcomb[0]
        filepath = rcomb["photlco.filepath"]
        ifilep = filepath + ifile
        ofilep = filepath + ofile

        preexistrow =  lsc.mysqldef.getlistfromraw(conn,'photlco','filename',ofile)
        if preexistrow:
            if _force == False and comparefakesourcemagnitude(ofilep,magnitude) == True:
                print ofile, "already created with magnitude of", magnitude
            else:
                lsc.mysqldef.deleteredufromarchive(ofile)
                print "Deleted old", ofile, "from archive"
                print "Dropping source into", ofile
                source_drop(ifilep, ofilep, rcomb["photlcoraw.cat_ra"],
                    rcomb["photlcoraw.cat_dec"],rcomb["photlco.fwhm"],
                    rcomb["photlco.z1"],magnitude,_psf)
                db_ingest(filepath,ofile,force=True)
        else:
                source_drop(ifilep, ofilep, rcomb["photlcoraw.cat_ra"],
                    rcomb["photlcoraw.cat_dec"],rcomb["photlco.fwhm"],
                    rcomb["photlco.z1"],magnitude,_psf)
                db_ingest(filepath,ofile,force=True)

    if _psf:
        os.system('lscingestredudata.py --obstype e93 -m ' + epoch + name)
    else:
        os.system('lscingestredudata.py --obstype e92 -m ' + epoch + name)

