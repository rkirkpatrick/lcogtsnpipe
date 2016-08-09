#!/usr/bin/env python

import os
import socket
import lsc
import numpy as np

description="Modified version of seatide_imagesub.py used to automatically "\
            "create and reduce fake source difference images."
usage = "%(prog)s [-e epoch -n objname -m magnitude]"

# Parse arguments
if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description=description, usage=usage)
    parser.add_argument("-n", "--name", dest="name", default='', type=str,
                        help='-n object name   \t [%(default)s]')
    parser.add_argument("-m", "--magnitude", nargs="+", default=[], help=
                        "Choose data from a specific magnitude", type=float)
    parser.add_argument("-e", "--epoch", dest="epoch", default='', type=str,
                        help='epoch to search for data \t [%(default)s]')
    parser.add_argument("-r", "--redo", dest="redo", default=False,
                        action="store_true", help="Redo psf on images and diff")
    parser.add_argument("-f","--fwhm", dest="fwhm", default=None,
                        help="Set the fwhm in the PSFs")
    parser.add_argument('-o','--optimal', dest="optimal", default=False,
                        action="store_true", help="Use optimal image differencing")

    args = parser.parse_args()
    _name = args.name
    _magnitude = args.magnitude
    _epoch = args.epoch
    _redo = args.redo
    _fwhm = args.fwhm
    _optimal = args.optimal

    # Change directory
    host = socket.gethostname()
    if host in ['intern1.lco.gtn']:
        os.chdir('/home/rkirkpatrick/seatide/tmp/')
    elif host in ['astrointern.lco.gtn']:
        os.chdir('/home/dguevel/snwd/tmp')
    # Try to connect to database
    try:
        hostname,username, passwd, database = lsc.mysqldef.getconnection('lcogt2')
        db = lsc.mysqldef.dbConnect(hostname, username, passwd, database)
    except:
        print 'Error: Could not connect to database'

    if _redo == False:
        # query photlco for available images
        query  = '''SELECT pr.objname, p.targetid, p.dayobs, p.filepath, p.wcs, p.filename, p.instrument '''
        query += '''FROM photlco AS p, photlcoraw AS pr '''
        query += '''WHERE pr.filename = p.filename '''
        query += '''AND pr.propid in ('CON2015B-002','CON2016A-004') '''
        query += '''AND p.filetype = 1 '''
        query += '''AND p.quality = 127 '''
        query += '''AND p.filename LIKE '%e91%' '''
        query += '''AND p.z1 != 9999 '''
        
    elif _redo == True:
        # query magcomparison for bad differences
        query = '''SELECT * FROM magcomparison AS p '''
        query += ''' WHERE (p.diffmag = 9999 '''
        query += ''' OR p.outlier = 1) '''

    if _epoch == '':
        pass
    elif '-' not in str(_epoch):
        query += ' AND p.dayobs = {0} '.format(_epoch)
    else:
        epoch1, epoch2 = _epoch.split('-')
        query += ' AND p.dayobs >= {0} AND p.dayobs <= {1} '.format(epoch1,epoch2)

    if _name != '':
        query += ' AND p.objname = "{0}" '.format(_name)

    query += '''ORDER by p.dayobs'''
    objects = lsc.mysqldef.sqlquery(db,query)

    if _magnitude == []:
        _magnitude = np.arange(14,20.5,0.5)

    if _redo == True and _fwhm == None:
        fwhm = " --fwhm=7 "
    elif _fwhm != None:
        fwhm = " --fwhm=" + str(_fwhm)
    else:
        fwhm = ''

    if _optimal == True:
        optimal = " --optimal "
    else:
        optimal = ''

    print objects

    for i,row in enumerate(objects):

        # pull new image from database
        objname = row['objname'] # name of object
        targetid = row['targetid']
        dayobs = row['dayobs']
        if _redo == True:
            instrument = "fl"
        else:
            instrument = row['instrument'][0:2]

        print ''
        print '''Checking for template images for %(dayobs)s target %(objname)s.''' %row

        # fix WCS if failed?
        if _redo == False:
            if row['wcs'] != 0:
                print '''No WCS, attempting to solve.'''
                command = 'lscloop.py -n ' + objname + ' -e ' + dayobs + ' -s wcs '
                print command
                os.system(command)

        # is there an LCOGT template?
        query  = '''SELECT filename, filepath, dayobs ,psf '''
        query += '''FROM photlco '''
        query += '''WHERE targetid=%s ''' %targetid
        query += '''AND instrument like '%fl%' '''
        query += '''AND filetype = 4 '''
        query += '''AND dayobs != '%s' ''' %dayobs
        query += '''AND filename LIKE '%e91%' '''
        query += '''ORDER BY dateobs, ut LIMIT 1 '''
        fileinfo = lsc.mysqldef.sqlquery(db,query)
        templatedate =''
        # if LCOGT template exists:
        if len(fileinfo)>0:
            templatesource = 'LCOGT'
            templatedate = fileinfo[0]['dayobs']
            print ''' Found LCOGT template from %s.''' %templatedate
            # does the template have a good PSF?
            if fileinfo[0]['psf'] == 'X':
                print ''' Making PSF for the template:'''
                command = 'lscloop.py --obstype e91 -n ' + objname + ' -e ' + templatedate + ' --filetype=4 -s psf' + fwhm
                print command
                os.system(command)
        else:
            # is there an LCOGT image that could be used as a template?
            query  = '''SELECT filename, filepath, dayobs, psf '''
            query += '''FROM photlco '''
            query += '''WHERE targetid=%s ''' %targetid
            query += '''AND instrument like '%fl%' '''
            query += '''AND filetype = 1 '''
            query += '''AND dayobs != %s ''' %dayobs
            query += '''AND filename LIKE '%e91%' '''
            query += '''ORDER BY dayobs, ut LIMIT 1 '''
            fileinfo = lsc.mysqldef.sqlquery(db,query)
            print fileinfo
            if len(fileinfo)>0:
                print ''' Found earlier LCOGT image to use as template. Will convert to template:'''
                templatesource = 'LCOGT'
                templatedate = fileinfo[0]['dayobs']
                if fileinfo[0]['psf'] == 'X' or _redo == True:
                    print ''' Making a psf for found image'''
                    command = 'lscloop.py --obstype e91 -n ' + objname + ' -e ' + templatedate + ' -s psf ' + fwhm
                    os.system(command)
                    print ''' Converting image to template:'''
                command = 'lscloop.py --obstype e91 -n ' + objname + ' -e ' + templatedate + ' -s template'
                print command
                os.system(command)
                print ''' Making PSF for the template:'''
                command = 'lscloop.py --obstype e91 -n ' + objname + ' -e ' + templatedate + ' --filetype=4 -s psf ' + fwhm
                print command
                os.system(command)

        if templatedate != '':
            print ''' Performing cosmic-ray rejection on template:'''
            command = 'lscloop.py --obstype e91 -n ' + objname + ' -e ' + templatedate + ' -s cosmic --filetype 4'
            print command
            os.system(command)

        # Create a redone psf for making fake sources
        if _redo == True:
            print "Redoing psf on e91"
            command = 'lscloop.py --obstype e91 -s psf ' + fwhm + ' -n ' + objname + ' -e ' + dayobs
            print command
            os.system(command)

        for i in _magnitude:
            inmag = i

            # create fake source images
            print '''Creating a fake source image with magnitude ''', inmag
            command = 'sourceingest.py --force -n ' + objname + ' -e ' + dayobs + ' -m ' + str(inmag)
            print command
            os.system(command)

            #create psf for fakesource image
            print ''' Creating psf for fake source image'''
            command = 'lscloop.py --obstype e93.fits -n ' + objname + ' -e ' + dayobs + ' -s psf ' + fwhm
            print command
            os.system(command)

            print ''' Performing cosmic-ray rejection on new image:'''
            command = 'lscloop.py --obstype e93.fits -n ' + objname + ' -e ' + dayobs + ' -s cosmic'
            print command
            os.system(command)

            # run image subtraction
            print ''' Running image subtraction:'''
            command = 'lscloop.py --obstype e93.fits e91 -n ' + objname + ' -e ' + dayobs + ' --tempdate ' + templatedate + ' -s diff --normalize t -T ' + instrument + optimal
            print command
            os.system(command)

            #run photometry on new diff image
            print '''Creating psf for diff e93 image'''
            command ='lscloop.py --obstype e93 --filetype=3 -s psf -n ' + objname + ' -e ' + dayobs
            print command
            os.system(command)

            print '''Calculate the instrumental magnitudes'''
            command = 'lscloop.py --obstype e93 --filetype=3 -s psfmag -x 1 -y 1 -c -n ' + objname + ' -e ' + dayobs
            print command
            os.system(command)

            print '''Calculate the zeropoint and color term'''
            command = 'lscloop.py --obstype e93 --filetype=3 -s zcat --zcatnew -n ' + objname + ' -e ' + dayobs
            print command
            os.system(command)

            print '''Calculate apparent magnitude'''
            command = 'lscloop.py --obstype e93 --filetype=3 -s mag -n ' + objname + ' -e ' + dayobs
            print command
            os.system(command)


            try:
                if _redo == False:
                    # read the apparent magnitude of the difference image
                    diffmagrow = lsc.mysqldef.getfromdataraw(db,'photlco','filename',diffname,'mag')
                    diffmag = diffmagrow[0]['mag']
                    print diffmag

                    query = 'SELECT * FROM magcomparison WHERE filename = "{0}" and inmag = {1}'.format(diffname, inmag)
                    preexistrow = lsc.mysqldef.sqlquery(db, query)
                    if preexistrow:
                        # update db with new data
                        row = preexistrow[0]
                        id = row['id']
                        lsc.mysqldef.updatevalue('magcomparison', 'outlier', False, id, filename0='id')
                        lsc.mysqldef.updatevalue('magcomparison', 'diffmag', diffmag, id, filename0='id')

                    else:
                        #put values into dictionary, then insert into magcomparison db
                        valuedict = {}
                        columns = ['targetid','filepath','filename','dayobs','objname','inmag','diffmag']
                        valuedict['targetid'] = targetid
                        valuedict['filepath'] = fpath
                        valuedict['filename'] = diffname
                        valuedict['dayobs'] = dayobs
                        valuedict['objname'] = objname
                        valuedict['inmag'] = inmag
                        valuedict['diffmag'] = diffmag
                        print valuedict
                        lsc.mysqldef.insert_values(db,'magcomparison',valuedict)

                if _redo == True:
                    # read the apparent magnitude of the difference image
                    diffname = row['filename'] # this row is from magcomparison
                    diffmagrow = lsc.mysqldef.getfromdataraw(db, 'photlco', 'filename', diffname, 'mag')
                    diffmag = diffmagrow[0]['mag']
                    print diffmag

                    # update db with new data
                    id = row['id']
                    lsc.mysqldef.updatevalue('magcomparison', 'outlier', False, id, filename0='id')
                    lsc.mysqldef.updatevalue('magcomparison', 'diffmag', diffmag, id, filename0='id')
                    
            except:
                print "Something went wrong with difference image photometry"


    # clean up
    os.system('rm -f *SDSS*.fits')
    os.system('rm -f *.list')
    os.system('rm -f *.txt')
    os.system('rm -f *.xml')
    os.system('rm -f substamplist')
    os.system('rm -f tmpcoo')


    #CREATE TABLE magcomparison( id bigint(20) UNSIGNED NOT NULL PRIMARY KEY AUTO_INCREMENT,
        #  targetid bigint(20) NOT NULL, filepath text DEFAULT NULL,
        #  filename varchar(100) DEFAULT NULL, dayobs text DEFAULT NULL,
        #  objname text DEFAULT NULL, inmag double DEFAULT 9999,
        #  diffmag double DEFAULT 9999, outlier boolean DEFAULT false);
