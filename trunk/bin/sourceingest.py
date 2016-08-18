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
from sourcesub import print_and_run_command

##############################################################################
################################Functions#####################################
##############################################################################
def query_for_files(conn, args):
    query = 'SELECT * FROM photlco WHERE filename LIKE "%e91%"\
                   AND z1 != 9999 AND filetype = 1 '

    if '-' not in str(args.epoch):
        query += ' AND dayobs = {0} '.format(args.epoch)
    else:
        epoch1, epoch2 = args.epoch.split('-')
        query += ' AND dayobs >= {0} AND dayobs <= {1} '.format(epoch1, epoch2)
    
    if args.filename != []:
        query = lsc.mysqldef.queryfilenamelike(query, args.filename)
    
    if args.name != '':
        query += ' AND objname = "{0}" '.format(args.name)

    resultSet = lsc.mysqldef.sqlquery(conn, query)

    return resultSet


def choosemagnitude(_magnitde=None):
    if _magnitude == None:
        if np.random.randint(2) == 0:
            mag = None
        else:
            mag = np.random.uniform(14,20)
    else:
        mag = _magnitude

    return mag


def comparefakemag(ofilep, magnitude):
    hdulist1 = pyfits.open(ofilep)
    header1 = hdulist1[0].header
    return(header1['FAKEMAG'] == magnitude)


def findzeropoint(row):
    zeropoint = 23
    if row['zcol1'] == "gr":
        zeropoint = row['z1']
    if row['zcol2'] == "gr":
        zeropoint = row['z2']

    print zeropoint

    return zeropoint


def get_ra_and_dec(row, args):
    if args.ras != None:
        ra = args.ras
    else:
        ra = row['ra0']

    if args.decs != None:
        dec = args.decs
    else:
        dec = row['dec0']

    ra = Angle(ra, unit=u.deg)
    sexa_ra = ra.to_string(unit=u.hour, sep=':')

    dec = Angle(dec, unit=u.deg)
    sexa_dec = dec.to_string(unit=u.degree, sep=':')

    return sexa_ra, sexa_dec


def get_parameters(row, args):
    sexa_ra, sexa_dec = get_ra_and_dec(row, args)
    zeropoint = findzeropoint(row)
    airmass = row['airmass']

    return sexa_ra, sexa_dec, zeropoint, airmass


def sexa2deg(ra,dec):
    ra = Angle(ra, u.hour).degree
    dec = Angle(dec, u.degree).degree
    return ra, dec


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
def magnitude2amplitude(psf, airmass, zeropoint, magnitude, exptime):
    filter = 0.203
    counts = float((10**((magnitude - zeropoint + (airmass*filter))/-2.5))*exptime)
    print counts
    amplitude = float(counts / np.sum(psf,dtype=float))
    return amplitude


def source_drop(ifilep, ofilep, magnitude, row, args):
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

    sexa_ra, sexa_dec, zeropoint, airmass = get_parameters(row, args)

    if args.ras != None or args.decs != None:
        hdr['RA'] = hdr['CAT-RA'] = hdr['OFST-RA'] = sexa_ra
        hdr['DEC'] = hdr['CAT-DEC'] = hdr['OFST-DEC'] = sexa_dec
        HDU[0].header = hdr

    ra, dec = sexa2deg(sexa_ra, sexa_dec)

    if magnitude != None:
        xpos, ypos = degs2coords(ra,dec,hdr)	
        exptime = hdr["exptime"]
        psf = createpsf(shape,xpos,ypos,ifilep)
        psf = psf * magnitude2amplitude(psf, airmass, zeropoint, magnitude, exptime)
        outimage = inimage + psf
    else:
        outimage = inimage

    HDU[0].data = outimage
    HDU[0].header.set("FAKEMAG", magnitude,
                      "If #s, then magnitude. Else, no fake source")
    HDU.writeto(ofilep, clobber = True)
    HDU.close()
    del HDU


def ingest_into_photlco(args):
    if args.name != '':
        name = ' -n ' + args.name + ''
    else:
        name = ''

    epoch = "-e " + args.epoch

    print '\n', ('#' * 75)
    command = 'lscingestredudata.py --obstype e93 -m ' + epoch + name
    print_and_run_command(command)


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
                        help='-n object name')
    parser.add_argument("-m", "--magnitude",type=float, default=None, help=
                        "Use a specific magnitude when computing PSF")
    parser.add_argument("-e", "--epoch", dest="epoch", default='', type=str,
                        help='epoch to reduce')
    parser.add_argument("-F", "--force", dest="force", action="store_true")
    parser.add_argument("--RAS", dest="ras", default=None, type=float,
                        help='Manually input right ascension in degrees')
    parser.add_argument("--DECS", dest="decs", default=None, type=float,
                        help='Manually input declination in degrees')


    args = parser.parse_args()
    _magnitude = args.magnitude
    _filename = args.filename
    _name = args.name
    _epoch = args.epoch
    _force = args.force
    _ras = args.ras
    _decs = args.decs

    if (_epoch == ''):
        sys.argv.append('--help')
        args = parser.parse_args()


    hostname,username, passwd, database = lsc.mysqldef.getconnection('lcogt2')
    conn = lsc.mysqldef.dbConnect(hostname, username, passwd, database)
    resultSet = query_for_files(conn, args)

    for row in resultSet:
        ifile = row['filename']
        ofile = ifile.replace('e91', 'e93')
        if _ras != None or _decs != None:
            ofile = ofile.replace('.fits','.moved.fits')
        filepath = row['filepath']
        ifilep = filepath + ifile
        ofilep = filepath + ofile
        print '#'* 75, '\n', ifile

        magnitude = choosemagnitude(args.magnitude)

        preexistrow =  lsc.mysqldef.getlistfromraw(conn,'photlco','filename',ofile)
        if preexistrow and _force == False and comparefakemag(ofilep, magnitude) == True:
            print ofile, "already created with magnitude of", magnitude
        else:
            lsc.mysqldef.deleteredufromarchive(ofile,"photlcoraw")
            lsc.mysqldef.deleteredufromarchive(ofile,"photlco")
            print "Dropping source of magnitude", magnitude, " into", ofile, '\n'
            source_drop(ifilep, ofilep, magnitude, row, args)
            db_ingest(filepath,ofile,force=True)


    ingest_into_photlco(args)


