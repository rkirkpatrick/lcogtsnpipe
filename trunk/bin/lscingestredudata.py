#!/usr/bin/env python
description = "> Ingest reduced data "
usage = "%(prog)s  -e epoch -T telescope [--filestr filestr]"

import string
import re, sys
from argparse import ArgumentParser
import datetime
import lsc

if __name__ == "__main__":
    parser = ArgumentParser(usage=usage, description=description, version="%(prog)s 1.0")
    parser.add_argument("-T", "--telescope", dest="telescope", default='all', type=str,
                      help='-T telescope fts, ftn, all, ' + ', '.join(list(lsc.telescope0['all'])) + ', '.join(
                          list(lsc.site0)) + ' \t [%(default)s]')
    parser.add_argument("-e", "--epoch", dest="epoch", default='20121212', type=str,
                      help='-e epoch \t [%(default)s]')
    # parser.add_option("-f", "--force",dest="force",action="store_true")
    parser.add_argument("-f", "--force", dest="force", default='no', type=str,
                      help='force ingestion \t [no/yes/update] \n')
    parser.add_argument("-m", "--missing", dest="missing", action="store_true")
    parser.add_argument("-t", "--type", dest="type", default='raw', type=str,
                      help='type -t type \t [raw/reduced] \n')
    parser.add_argument("-n","--object", dest="object", default='', type=str,
                      help='type --object object \t [name] \n')
    parser.add_argument("--filestr",nargs="+",type=str,dest="filestr", default=[], help = '--filestr\
                         [e90,e91,e92]\t [%(default)s]\n')


    args = parser.parse_args()
    _filestr = args.filestr
    if args.type not in ['raw', 'redu']:  sys.argv.append('--help')
    if len(sys.argv) < 2:
        sys.argv.append('--help')
    if args.force not in ['no', 'yes', 'update']:
        sys.argv.append('--help')
    args = parser.parse_args()
    _telescope = args.telescope
    _object = args.object
    epoch = args.epoch
    _force = args.force
    _missing = args.missing
    if not _missing:
        _missing = False
    else:
        _force = True
    _type = args.type

    hostname, username, passwd, database = lsc.mysqldef.getconnection('lcogt2')
    conn = lsc.mysqldef.dbConnect(hostname, username, passwd, database)

    if not _missing:
        if '-' not in str(epoch):
            epoch0 = datetime.date(int(epoch[0:4]), int(epoch[4:6]), int(epoch[6:8]))
            pippo = lsc.mysqldef.getlistfromraw(conn, 'photlcoraw', 'dayobs', str(epoch0), '', '*', _telescope, _filestr)
        else:
            epoch1, epoch2 = string.split(epoch, '-')
            start = re.sub('-', '', str(datetime.date(int(epoch1[0:4]), int(epoch1[4:6]), int(epoch1[6:8]))))
            stop = re.sub('-', '', str(datetime.date(int(epoch2[0:4]), int(epoch2[4:6]), int(epoch2[6:8]))))
            pippo = lsc.mysqldef.getlistfromraw(conn, 'photlcoraw', 'dayobs', str(start), str(stop), '*',
                                                _telescope, _filestr)
        listingested = [i['filename'] for i in pippo]
    else:
        if '-' not in str(epoch):
            epoch0 = re.sub('-', '', str(datetime.date(int(epoch[0:4]), int(epoch[4:6]), int(epoch[6:8]))))
            pippo = lsc.mysqldef.getmissing(conn, epoch0, '', _telescope, 'photlco', _filestr)
        else:
            epoch1, epoch2 = string.split(epoch, '-')
            pippo = lsc.mysqldef.getmissing(conn, epoch1, epoch2, _telescope, 'photlco', _filestr)
            print pippo
            print 'here'
        listingested = [i['filename'] for i in pippo]
    
    listname = [i['objname'] for i in pippo]
#    print listingested
#    print listname

    if _object:
         for ii, img in enumerate(listingested):
             if listname[ii] == _object:
                 lsc.mysqldef.ingestredu([img], _force, 'photlco')
             else:
                 print listname[ii],ii
                 print 'not the right object'
    else:
        for img in listingested:
            lsc.mysqldef.ingestredu([img], _force, 'photlco')
