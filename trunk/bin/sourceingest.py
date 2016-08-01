#!/usr/bin/env python
"""
How to run it now:

sourceingest.py -e YYYYMMDD-YYYYMMDD 
lscloop.py -e YYYYMMDD-YYYYMMDD --obstype e93 -s psf
lscloop.py -e YYYYMMDD-YYYYMMDD --obstype e93 -s psfmag -c
lscloop.py -e YYYYMMDD-YYYYMMDD --obstype e93 -n OBJNAME -s zcat --zcatnew
lscloop.py -e YYYYMMDD-YYYYMMDD --obstype e93 -s mag

If you use --move, you must use --RAS and --DEC in lscloop.py
-s psfmag. Otherwise, it won't chose the new source.

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


def findzeropoint(mag, psfmag, airmass):
    filter = 0.203
    zeropoint = mag - psfmag + (airmass * filter)
    return zeropoint


def sexa2deg(ra,dec):
    ra = Angle(ra,u.hour).degree
    dec = Angle(dec,u.degree).degree
    return ra, dec
##############################################################################
##############################Drop PSF########################################
##############################################################################

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


def source_drop(ifilep, ofilep, ra, dec, zeropoint, magnitude=None, _move=False):
    """
    Copies ifile into ofile. Drops a psf into ofile, depending
    if magnitude exists. Then to all ofile 'FAKEMAG' is added to the header
    to describe the magnitude of the dropped source. If there is no
    magnitude, the 'FAKEMAG' header will be set as an empty string.
    """
    HDU = pyfits.open(ifilep)
    inimage = HDU[0].data
    hdr = HDU[0].header
    shape = inimage.shape

    if _move == True:
        hdr['RA'] = hdr['CAT-RA'] = hdr['OFST-RA'] = ra
        hdr['DEC'] = hdr['CAT-DEC'] = hdr['OFST-DEC'] = dec
        HDU[0].header = hdr

        ra, dec = sexa2deg(ra,dec)
        print '#'*75
        print 'Use this RAS and DECS for lscloop.py -s psfmag:', ra, dec
        print '#'*75

    if magnitude != None:
        xpos, ypos = degs2coords(ra,dec,hdr)	
        exptime = hdr["exptime"]
        psf = createpsf(shape,xpos,ypos,ifilep)
        psf = psf * magnitude2amplitude(psf,magnitude,zeropoint,exptime)
        outimage = inimage + psf
    else:
        outimage = inimage

    HDU[0].data = outimage
    HDU[0].header.set("FAKEMAG", magnitude,
                      "If #s, then magnitude. Else, no fake source")
    HDU.writeto(ofilep, clobber = True)
    HDU.close()
    del HDU

##############################################################################
description='Creates e93 files where e93 files have an inserted psf based'\
            'off of the e91 file.'
usage='%(prog)s -e epoch [-f filename -m magnitude]'
version='%(prog)s 1.0'

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
    parser.add_argument("--move",dest="move",action="store_true",
                        default=False, help="Move the placed source")
    parser.add_argument("--ra", dest="ra", default='', type=str,
                        help="Manually input RA")
    parser.add_argument("--dec", dest="dec", default='', type=str,
                        help="Manually input DEC")
    parser.add_argument("--negdec", dest="negdec",action="store_true",
                        default=False, help="Add negative sign to DEC")

    args = parser.parse_args()
    _magnitude = args.magnitude
    _filename = args.filename
    _name = args.name
    _epoch = args.epoch
    _force = args.force
    _move = args.move
    _ra = args.ra
    _dec = args.dec
    _negdec = args.negdec

    if (_epoch == ''):
        sys.argv.append('--help')
        args = parser.parse_args()

    if(_move == True and (_ra == '' or _dec == '')):
        sys.argv.append('--help')
        args = parser.parse_args()

    if _negdec == True:
        _dec = "-" + _dec

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
    else:
        name = ''

    c = conn.cursor()
    c.execute(query_files)
    r = c.fetchall()
    rname = ()
    for row in r:
        rtuple = (row[0],)
        rname = rname + rtuple

    for ifile in rname:
        print '#'* 75, '\n', ifile
        ofile = ifile.replace('e91','e93')
        if _move == True:
            ofile = ofile.replace('.fits','.moved.fits')
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

        if _ra != '':
            ra = _ra
        else:
            ra = rcomb["photlcoraw.cat_ra"]

        if _dec != '':
            dec = _dec
        else:
            dec = rcomb["photlcoraw.cat_dec"]

        zeropoint = findzeropoint(rcomb["photlco.mag"],rcomb["photlco.psfmag"],rcomb["photlco.airmass"])
        print zeropoint
        preexistrow =  lsc.mysqldef.getlistfromraw(conn,'photlco','filename',ofile)
        if preexistrow:
            if _force == False and comparefakesourcemagnitude(ofilep,magnitude) == True:
                print ofile, "already created with magnitude of", magnitude
            else:
                lsc.mysqldef.deleteredufromarchive(ofile,"photlcoraw")
                lsc.mysqldef.deleteredufromarchive(ofile,"photlco")
                print "Deleted old", ofile, "from archive"
                print "Dropping source of magnitude", magnitude, " into", ofile, '\n'
                source_drop(ifilep, ofilep, ra, dec, zeropoint,
                            magnitude, _move)
                db_ingest(filepath,ofile,force=True)
        else:
                print "Dropping source into", ofile, '\n'
                source_drop(ifilep, ofilep, ra, dec, zeropoint,
                            magnitude, _move)
                db_ingest(filepath,ofile,force=True)


    print '\n',('#' * 75)
    command = 'lscingestredudata.py --obstype e93 -m ' + epoch + name
    print command, '\n'
    os.system(command)

