#!/usr/bin/env python

import os
import socket
import lsc
import numpy as np

description="Modified version of seatide_imagesub.py used to automatically "\
            "create and reduce fake source difference images."
usage = "%(prog)s [-e epoch -n objname -m magnitude]"
###############################################################################
def query_photlco_for_targets():
    # query photlco for available target images
    query  = '''SELECT pr.objname, p.targetid, p.dayobs, p.filepath, p.wcs, p.filename, p.instrument '''
    query += '''FROM photlco AS p, photlcoraw AS pr '''
    query += '''WHERE pr.filename = p.filename '''
    query += '''AND pr.propid in ('CON2015B-002','CON2016A-004') '''
    query += '''AND p.filetype = 1 '''
    query += '''AND p.quality = 127 '''
    query += '''AND p.filename LIKE '%e91%' '''
    query += '''AND p.z1 != 9999 '''

    return query


def query_magcomparison_for_redo_targets(_optimal):
    # query magcomparison for bad differences
    query = '''SELECT * FROM magcomparison AS p '''
    query += ''' WHERE (p.diffmag = 9999 '''
    query += ''' OR p.outlier = 1) '''

    if _optimal == True:
        query += ''' AND p.filename like '%optimal%' '''

    return query


def add_arguments_to_query(query, _epoch, _name):
    if _epoch == '':
        pass
    elif '-' not in str(_epoch):
        query += ' AND p.dayobs = {0} '.format(_epoch)
    else:
        epoch1, epoch2 = _epoch.split('-')
        query += ' AND p.dayobs >= {0} AND p.dayobs <= {1} '.format(epoch1,epoch2)

    if _name != '':
        query += ' AND p.objname = "{0}" '.format(_name)

    query += '''ORDER by p.dayobs '''
    return query


def choose_not_done_mags(db, diffname, _magnitude):
    magdonequery = "Select inmag from magcomparison where filename = '{0}' ".format(diffname)
    magdoneSet = lsc.mysqldef.sqlquery(db, magdonequery)
    magdone = []
    mags_to_be_done = []
    for magrow in magdoneSet:
        magdone.append(magrow['inmag'])

    for inmag in _magnitude:
        if inmag not in magdone:
            mags_to_be_done.append(inmag)

    return mags_to_be_done


def create_loop_command(dayobs, objname, obstype=[], filetype=1):
    if obstype != []:
        argument_obstype = ' --obstype '
        for obs in obstype:
            argument_obstype += obs + ' '
    else:
        argument_obstype = ''

    command = 'lscloop.py -e ' + dayobs + ' -n ' + objname + argument_obstype + ' --filetype ' + str(filetype)
    return command


def print_and_run_command(command):
    print(command)
    os.system(command)


def perform_lscloop(command, stage="psf", fwhm='', show='',instrument='fl', optimal='',tempdate='', _temptel='' ):
    #performs lscloop stages wcs, psf, psfmag, zcat, and mag
    if stage == 'wcs':
        command += ' -s wcs '
    elif stage == 'psf':
        command += ' -s psf ' + fwhm + show
    elif stage == 'cosmic':
        command += ' -s cosmic '
    elif stage == 'template':
        command += ' -s template '
    elif stage == 'psfmag':
        command += ' -s psfmag -c '
    elif stage == 'zcat':
        command += ' -s zcat --zcatnew '
    elif stage == 'mag':
        command += ' -s mag '
    elif stage == 'diff':
        if tempdate != '':
            argument_tempdate = ' --tempdate ' + tempdate
        else:
            argument_tempdate = ''
        if _temptel != '':
            argument_temptel = ' --temptel ' + _temptel[0:2]
        else:
            argument_temptel = ''

        command += ' -s diff --normalize t -T ' + instrument[0:2] + optimal + argument_tempdate + argument_temptel

    print_and_run_command(command)


def do_photometry(dayobs, objname, obstype=[], filetype=1, fwhm='', show=''):
    loop_command = create_loop_command(dayobs, objname, obstype, filetype)
    perform_lscloop(loop_command, stage="psf", fwhm=fwhm, show=show)
    perform_lscloop(loop_command, stage="psfmag")
    perform_lscloop(loop_command, stage="zcat")
    perform_lscloop(loop_command, stage="mag")


def query_for_previous_templates(db, row, _tempdate='',_temptel=''):
    query = '''SELECT * FROM photlco '''
    query += '''WHERE targetid=%s ''' % row['targetid']
    query += '''AND filetype = 4 '''
    query += '''AND filename NOT LIKE '%e93%' '''
    query += '''AND dayobs != '%s' ''' % row['dayobs']
    if _tempdate != '':
        query += '''AND dayobs = {} '''.format(_tempdate)
    if _temptel != '':
        query += '''AND instrument LIKE '%{}%' '''.format(_temptel)
    query += '''ORDER BY dateobs, ut LIMIT 1 '''
    print query
    fileinfo = lsc.mysqldef.sqlquery(db, query)
    return fileinfo


def query_for_possible_templates(db, row, _tempdate='',_temptel=''):
    # is there an LCOGT image that could be used as a template?
    query = '''SELECT * FROM photlco '''
    query += '''WHERE targetid=%s ''' % row['targetid']
    query += '''AND filetype = 1 '''
    query += '''AND filename NOT LIKE '%e93%' '''
    query += '''AND dayobs != %s ''' % row['dayobs']
    query += ''' AND wcs = 0 '''
    if _tempdate != '':
        query += ''' AND dayobs = {} '''.format(_tempdate)
    if _temptel != '':
        query += ''' AND instrument LIKE '%{}%' '''.format(_temptel)
    query += '''ORDER BY dayobs, ut LIMIT 2 '''
    fileinfo = lsc.mysqldef.sqlquery(db, query)
    return fileinfo


def items_from_temp(templates):
    tempdate = templates[0]['dayobs']
    temppsf = templates[0]['psf']
    objname = templates[0]['objname']
    tempname = templates[0]['filename']

    return tempdate, tempname, temppsf, objname


def make_template_psf(templates, tempfwhm, show):
    tempdate, tempname, temppsf, objname = items_from_temp(templates)
    if "temp" not in tempname:
        tempname = tempname.replace('.fits','.temp.fits')
    template_psf_command = create_loop_command(tempdate, objname, obstype=[tempname], filetype=4)
    print '''Making PSF for the template:'''
    perform_lscloop(template_psf_command, stage='psf', fwhm=tempfwhm, show=show)
    

def find_and_make_template(db, row, _redo=False, _tempdate='', _temptel='', tempfwhm='', show=''):

    previous_templates = query_for_previous_templates(db, row, _tempdate=_tempdate, _temptel=_temptel)
    # if LCOGT template exists and you don't want to redo the psfs:
    if len(previous_templates) > 0 and _redo == False:
        tempdate, tempname, temppsf, objname = items_from_temp(previous_templates)
        print '''Found LCOGT template from %s.''' % tempdate
        
        # does the template have a good PSF?
        if temppsf == 'X':
            make_template_psf(previous_templates, tempfwhm, show)

    # If old templates don't exist or you want to redo data
    else:
        new_templates = query_for_possible_templates(db,row,_tempdate=_tempdate, _temptel=_temptel)
        # If images exist for a new template
        if len(new_templates) > 0:
            print '''Found LCOGT image to use as template. Will convert to template:'''

            tempdate, filename, temppsf, objname = items_from_temp(new_templates)
            tempname = filename.replace('.fits', '.temp.fits')
            new_temp_loop_command = create_loop_command(tempdate, objname, obstype=[filename], filetype=1)

            if temppsf == 'X' or _redo == True:
                print '''Making a psf for found image'''
                perform_lscloop(new_temp_loop_command, stage='psf', fwhm=tempfwhm, show=show)

            print '''Converting image to template:'''
            perform_lscloop(new_temp_loop_command, stage='template')

            make_template_psf(new_templates, tempfwhm, show)
        else:
            print "Error! No images available with same object name!"
            tempdate, tempname = '', ''

    return tempdate, tempname


def define_different_obstypes(tempname, suffix):
    fake_img_obstype = ['e93']
    differencing_obstype = fake_img_obstype + [tempname]
    post_subtraction_obstype = [fake_img_obstype[0] + suffix]
    return fake_img_obstype, differencing_obstype, post_subtraction_obstype


def reduce_fake_source_image(dayobs, objname, inmag, obstype, fwhm, show):
    # Create fake source image
    source_ingest_command = 'sourceingest.py --force -n ' + objname + ' -e ' + dayobs + ' -m ' + str(inmag)
    print_and_run_command(source_ingest_command)

    # Create psf for fakesource image
    fake_source_command = create_loop_command(dayobs, objname, obstype=obstype)
    perform_lscloop(fake_source_command, stage='psf', fwhm=fwhm, show=show)

    # Perform cosmic ray rejection on fakesource image
    perform_lscloop(fake_source_command, 'cosmic')


def get_new_diffmag(db, diffname):
    # read the apparent magnitude of the difference image
    diffmagrow = lsc.mysqldef.getfromdataraw(db, 'photlco', 'filename', diffname, 'mag')
    diffmag = diffmagrow[0]['mag']

    return diffmag


def query_for_preexisting_row(db, diffname, inmag):
    query = 'SELECT * FROM magcomparison WHERE filename = "{0}" and inmag = {1}'.format(diffname, inmag)
    preexistrow = lsc.mysqldef.sqlquery(db, query)

    return preexistrow


def update_magcomparison(diffmag, id):
    lsc.mysqldef.updatevalue('magcomparison', 'outlier', False, id, filename0='id')
    lsc.mysqldef.updatevalue('magcomparison', 'diffmag', diffmag, id, filename0='id')


def insert_new_into_magcomparison(db, row, diffname, inmag, diffmag):
    valuedict = {}
    valuedict['targetid'] = row['targetid']
    valuedict['filepath'] = row['filepath']
    valuedict['filename'] = diffname
    valuedict['dayobs'] = row['dayobs']
    valuedict['objname'] = row['objname']
    valuedict['inmag'] = inmag
    valuedict['diffmag'] = diffmag
    valuedict['instrument'] = row['instrument']
    lsc.mysqldef.insert_values(db, 'magcomparison', valuedict)


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
                        action="store_true", help="Redo psf when diff is an outlier")
    parser.add_argument("-F","--fwhm", dest="fwhm", default=None,
                        help="Set the fwhm in the PSFs", type=float)
    parser.add_argument("--tempfwhm", dest="tempfwhm", default=None,
                        help="Set the fwhm in the template PSFs", type=float)
    parser.add_argument("--tempdate", dest="tempdate", default='', type=str,
                        help='--tempdate  template date \t [%(default)s]')
    parser.add_argument("--temptel", dest="temptel", default='', type=str,
                    help='--temptel  template instrument \t [%(default)s]')
    parser.add_argument('-o','--optimal', dest="optimal", default=False,
                        action="store_true", help="Use optimal image differencing")
    parser.add_argument('-f', '--force', dest='force', default=False,
                        action='store_true', help="Do subtraction on all specified images, even ones already done.")
    parser.add_argument('-s', '--show', dest='show', default=False,
                        action='store_true', help="Use show when creating PSFs")


    args = parser.parse_args()
    _name = args.name
    _magnitude = args.magnitude
    _epoch = args.epoch
    _redo = args.redo
    _fwhm = args.fwhm
    _optimal = args.optimal
    _force = args.force
    _tempdate = args.tempdate
    _temptel = args.temptel
    _tempfwhm = args.tempfwhm
    _show = args.show

    if _magnitude == []:
        _magnitude = np.arange(14,20.5,0.5)
    mags_to_be_done = _magnitude

    if _fwhm != None:
        fwhm = " --fwhm=" + str(_fwhm)
    else:
        fwhm = ''

    if _tempfwhm != None:
        tempfwhm = " --fwhm=" + str(_tempfwhm)
    else:
        tempfwhm = ''

    if _optimal == True:
        optimal = " --optimal "
    else:
        optimal = ''

    if _show == True:
        show = ' --show '
    else:
        show = ''

    # Query Database for images to be source injected and differenced
    try:
        hostname,username, passwd, database = lsc.mysqldef.getconnection('lcogt2')
        db = lsc.mysqldef.dbConnect(hostname, username, passwd, database)
    except:
        print 'Error: Could not connect to database'

    if _redo == True:
        query = query_magcomparison_for_redo_targets(_optimal)
    else:
        query = query_photlco_for_targets()

    query = add_arguments_to_query(query, _epoch, _name)
    print query
    objects = lsc.mysqldef.sqlquery(db, query)


    # Go through each image from the query and do injection and differencing
    for row in objects:
        # Store row info into variables
        objname = row['objname']
        targetid = row['targetid']
        dayobs = row['dayobs']
        instrument = row['instrument']
        filepath = row['filepath']
        filename = row['filename']

        suffix = lsc.myloopdef.difference_suffix(_temptel, _optimal)

        if _redo == False:
            diffname = filename.replace('e91.fits', 'e93' + suffix)

            # If false, remove from the mags to be done inmags that have already been done
            # So if true, the do all magnitudes specified
            if _force == False:
                mags_to_be_done = choose_not_done_mags(db, diffname, _magnitude)

        # If magnitudes have not been done yet for the image
        if len(mags_to_be_done) > 0:
            if _redo == True:
                # Redo photometry on original images
                do_photometry(dayobs, objname, obstype=['e91'], fwhm=fwhm, show=show)

            print '\n', '''Checking for template images for %(dayobs)s target %(objname)s.''' % row

            # fix WCS if failed?
            if _redo == False:
                if row['wcs'] != 0:
                    bad_wcs_loop_command = create_loop_command(dayobs, objname)
                    perform_lscloop(bad_wcs_loop_command, stage='wcs')

            tempdate, tempname = find_and_make_template(db, row, _redo, _tempdate, _temptel, tempfwhm, show)

            # If a template exists
            if tempdate != '':
                print '''Performing cosmic-ray rejection on template:'''
                template_cosmic_command = create_loop_command(tempdate, objname,filetype=4)
                perform_lscloop(template_cosmic_command, stage='cosmic')

                # Reduce the data for each magnitude
                for inmag in mags_to_be_done:
                    print '''Injecting fake source of magnitude''', inmag
                    fakeimgobs, diffobs, postobs = define_different_obstypes(tempname, suffix)
                    reduce_fake_source_image(dayobs, objname, inmag, fakeimgobs, fwhm, show)

                    # Run image subtraction
                    print '''Running difference imagint on fake image'''
                    diff_command = create_loop_command(dayobs, objname, obstype=diffobs)
                    perform_lscloop(diff_command, stage='diff', instrument=instrument, optimal=optimal, tempdate=tempdate, _temptel=_temptel)

                    #Do photometry on difference image
                    print '''Doing photometry on difference image'''
                    do_photometry(dayobs, objname, obstype=postobs, filetype=3)


                    try:
                        if _redo == False:
                            diffmag = get_new_diffmag(db, diffname)
                            print diffmag
                            preexistrow = query_for_preexisting_row(db, diffname, inmag)

                            if len(preexistrow) > 0:
                                # update db with new data
                                print ("Data exists for magnitude of " + str(inmag) + ", so magcomparsion will be updated")
                                update_magcomparison(diffmag, preexistrow[0]['id'])
                            else:
                                #put values into dictionary, then insert into magcomparison db
                                print("A row doesn't already exist. Inserting values into magcomparison.")
                                insert_new_into_magcomparison(db, row, diffname, inmag, diffmag)

                        if _redo == True:
                            diffname = row['filename'] # this row is from magcomparison
                            diffmag = get_new_diffmag(db, diffname)
                            print diffmag
                            print "Database is being updated with redone photometry"
                            update_magcomparison(diffmag, row['id'])

                    except:
                        print "Something went wrong with difference image photometry"
        else:
            print "No magnitudes left to be done for", filename


    # clean up
    os.system('rm -f *SDSS*.fits')
    os.system('rm -f *.list')
    os.system('rm -f *.txt')
    os.system('rm -f *.xml')
    os.system('rm -f substamplist')
    os.system('rm -f tmpcoo')


    #CREATE TABLE magcomparison( id bigint(20) UNSIGNED NOT NULL PRIMARY KEY AUTO_INCREMENT,
        #  targetid bigint(20) NOT NULL, instrument text DEFAULT NULL, filepath text DEFAULT NULL,
        #  filename varchar(100) DEFAULT NULL, dayobs text DEFAULT NULL,
        #  objname text DEFAULT NULL, inmag double DEFAULT 9999,
        #  diffmag double DEFAULT 9999, outlier boolean DEFAULT false);
