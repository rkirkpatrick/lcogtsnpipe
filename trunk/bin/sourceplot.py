#!/usr/bin/env python

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

import lsc
description = "Displays a plot of the data of fake sources within magcomparison"
usage = "%(prog)s -e epoch [-g graph -n objname -m magnitude]"
version = "0.1"

def updateoutliers(resultSet,_magnitude):
    # Group rows into lists based on inmag
    maggroup= []
    for inmag in _magnitude:
        maggroup.append([])

    for row in resultSet:
        ii = 0
        for inmag in _magnitude:
            if row['inmag'] == inmag:
                maggroup[ii].append(row)
            ii += 1
    mi = 0
    for inmagrow in maggroup:
        try:
            # Find outlier range in each inmag
            diffmagset = []
            for row in inmagrow:
                diffmagset.append(row["diffmag"])
            Q1 = np.percentile(diffmagset,25)
            Q3 = np.percentile(diffmagset,75)
            IQR = Q3 - Q1
            outlierrange = [Q1 - (1.5 * IQR), Q3 + (1.5 * IQR)]
            print "The outlier range for a magnitude of", _magnitude[mi], "is", outlierrange

            # Set the data outlier key accordingly & update db
            for row in inmagrow:
                if row["diffmag"] < outlierrange[0] or row["diffmag"] > outlierrange[1]:
                    row["outlier"] = True
                else:
                    row["outlier"] = False
                value = row['outlier']
                id = row['id']
                lsc.mysqldef.updatevalue('magcomparison', 'outlier', value, id, connection="lcogt2", filename0="id")
        except:
            print "No data for magnitude:", _magnitude[mi]
        mi += 1

    print "Outliers have been updated", '\n', '#' * 75


def len2shape(pltlength):
    #find different shapes that would work
    #choose shape with smallest difference between two parts
    pass



def makehistogram(oneset, inmag, mean, stddev):
    num_bins = 15
    # the histogram of the data
    n, bins, patches = plt.hist(oneset, num_bins, normed=1, facecolor='blue', alpha=0.5)
    # add a 'best fit' line
    y = mlab.normpdf(bins, mean, stddev)
    plt.plot(bins, y, 'r--')
    plt.xlabel("Difference Apparent Magnitudes")
    plt.ylabel('Probability')
    plt.title(r'Histogram of Input Magnitude:' + str(inmag))

    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)


def varstats(listofmagnitudes, datasets):
    means = []
    stddevs = []
    for i in range(len(listofmagnitudes)):
        means.append(np.mean(datasets[i]))
        stddevs.append(np.std(datasets[i]))

    return means, stddevs

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description=description, usage=usage, version=version)
    parser.add_argument("-n", "--name", dest="name", default='', type=str,
                        help='-n object name   \t [%(default)s]')
    parser.add_argument("-m", "--magnitude", nargs="+", default=[], help=
                        "Choose data from a specific magnitude", type=float)
    parser.add_argument("-e", "--epoch", dest="epoch", default='', type=str,
                        help='epoch to search for data \t [%(default)s]')
    parser.add_argument("-g", "--graph", dest="graph", type=str, default="scatter",
                        help="Choose type of graph [scatter, histogram, boxplot, varstats]")
    parser.add_argument("--bestfit",dest="bestfit", action="store_true",
                        default=False, help="Make a bestfit line in scatterplot")
    parser.add_argument("--yx", dest="yx", action="store_true",
                        default=False, help="Display a y = x graph with scatter")
    parser.add_argument("-b","--bad", dest="bad", action="store_true",
                        default=False, help="include images with bad diffmag from plots")
    parser.add_argument("-o","--keepout", dest="keepout", action="store_true",
                        default=False, help="Keep outliers in data")
    parser.add_argument("-u","--updateout", dest="updateout", action="store_true",
                        default=False, help="Update the outliers in the database")

    args = parser.parse_args()
    _name = args.name
    _magnitude = args.magnitude
    _epoch = args.epoch
    _graph = args.graph
    _bestfit= args.bestfit
    _yx = args.yx
    _bad = args.bad
    _keepout = args.keepout
    _updateout = args.updateout

    # Set query based on arguments
    query = "SELECT * from magcomparison "

    if (_epoch == ''):
        sys.argv.append('--help')
        args = parser.parse_args()

    epoch = "-e " + _epoch
    if '-' not in str(_epoch):
        query += ' WHERE dayobs = {0} '.format(_epoch)
    else:
        epoch1, epoch2 = _epoch.split('-')
        query += ' WHERE dayobs >= {0} AND dayobs <= {1} '.format(epoch1,epoch2)

    if _name != '':
        query += ' AND objname = "{0}" '.format(_name)

    if _bad is False:
        query += ' AND diffmag != 9999 '

    if _magnitude != []:
        if len(_magnitude) == 1:
            query += ' AND inmag = {0} '.format(_magnitude[0])
        else:
            for i in range(len(_magnitude)):
                if i == 0:
                    query += ' AND ( inmag = {0} '.format(_magnitude[i])
                else:
                    query += ' OR inmag = {0} '.format(_magnitude[i])
            query += ") "
    elif _magnitude == []:
        _magnitude = np.arange(14, 20.5, 0.5)
    _magnitude.sort()

    print "#" * 75

    # Connect to database, and query the db
    try:
        hostname,username, passwd, database = lsc.mysqldef.getconnection('lcogt2')
        db = lsc.mysqldef.dbConnect(hostname, username, passwd, database)
    except:
        print 'Error: Could not connect to database'

    resultSet = lsc.mysqldef.sqlquery(db,query)
    if _keepout is False:
        if _updateout is True:
            updateoutliers(resultSet,_magnitude)
        query += " AND outlier = FALSE"
        resultSet = lsc.mysqldef.sqlquery(db,query)


    # Set x(inmag) data and y(diffmag) data as ndarrays
    x = []
    y = []
    for row in resultSet:
        x.append(row['inmag'])
        y.append(row['diffmag'])
    x = np.array(x)
    y = np.array(y)

    # Create datasets for each inmag
    dataset = []
    for inmag in _magnitude:
        dataset.append([])
    for xi in range(len(x)):
        ii = 0
        for inmag in _magnitude:
            if x[xi] == inmag:
                dataset[ii].append(y[xi])
            ii += 1

    for i in range(len(_magnitude)):
        dataset[i] = np.array(dataset[i])

    # If user wants a scatterplot
    if _graph == "scatter":
        plt.scatter(x,y)
        # Create bestfit line in scatterplot
        if _bestfit == True:
            a = np.polyfit(x,y,1)
            plt.plot(x,a[0]*x+a[1],color="red",label="Best Fit")
            print "Best fit line data:"
            print "Slope =", str(a[0]) + ",", "y-int =", a[1]

        # Create y = x line in scatterplot
        if _yx == True:
            plt.plot(x,x,color="purple", label="y = x")
        # Create legend
        if _bestfit or _yx:
            plt.legend(loc="upper left")

    # If user wants a histogram
    elif _graph == "histogram":
        means, stddevs = varstats(_magnitude,dataset)
        for i in range(len(_magnitude)):
            makehistogram(dataset[i], _magnitude[i], means[i], stddevs[i])



    # If user wants a boxplot
    elif _graph == "boxplot":
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 7))

        # Rectangular box plot
        bplot = ax.boxplot(dataset)

        # Adding horizontal grid lines
        ax.yaxis.grid(True)
        ax.set_xticks([y + 1 for y in range(len(dataset))], )
        ax.set_ylabel("Difference Apparent Magnitudes")
        ax.set_xlabel("Input Apparent Magnitudes")

        # Add x-tick labels
        plt.setp(ax, xticks=[y + 1 for y in range(len(dataset))],
                 xticklabels=_magnitude)


    # If user wants variable statistics for each inmag
    elif _graph == "varstats":
        means, stddevs = varstats(_magnitude, dataset)
        for i in range(len(_magnitude)):
            print "For an input magnitude of", _magnitude[i]
            print "Mean =", means[i], "and Standard Deviation =", stddevs[i], '\n'


    else:
        sys.argv.append('--help')
        args = parser.parse_args()

    plt.show()