#!/usr/bin/env python
"""
Compute the expected alert times for a set of earthquakes based on observed
system delay times using the method of Behr et al. [2015].

Created on Jul 31, 2015

@author: behry
"""

import os
import sys
sys.path.append('../delays')
from delayeew import DelayEEW
import json
import numpy as np
from scipy.stats import scoreatpercentile


class NetworkInfo:
    """
    Load station information.
    """

    def __init__(self):
        self.lat = []
        self.lon = []
        self.chn = []
        self.nw = []
        self.nm = []
        self.lc = []
        self.excludes = ['RA.OGSI', 'RA.STBO', 'RA.STFL']
        self.size = 0
        self.networks = {}
        self.networks['xx'] = {'lat': [], 'lon': [], 'chn': [],
                               'nw': [], 'nm': [], 'lc': [],
                               'color': 'white', 'label': ''}

    def read(self, fn, sm=True):
        f = open(fn, 'r')
        for line in f.readlines():
            nt, st, chn, loc, lat, lon = line.split()
            if not sm:
                if chn[0:2] == 'HG' or chn[0:2] == 'HN':
                    continue
            ns = '%s.%s' % (nt, st)
            if ns in self.excludes:
                continue
            if ns not in self.networks['xx']['nm']:
                self.networks['xx']['nm'].append(ns)
                self.networks['xx']['nw'].append(nt)
                self.networks['xx']['chn'].append(chn)
                self.networks['xx']['lc'].append(loc)
                self.networks['xx']['lat'].append(float(lat))
                self.networks['xx']['lon'].append(float(lon))
                self.size += 1
        f.close()

    def get_networks(self):
        return self.networks


class DelayObs:
    """
    Load processing latencies.
    """
    def __init__(self, pkdelfn, envdelfn, octdelfn, magdelfn):
        self.pkdelfn = pkdelfn
        self.envdelfn = envdelfn
        self.octdelfn = octdelfn
        self.magdelfn = magdelfn

    def load(self, perstation=False):
        # envelopes
        fh = open(self.envdelfn)
        self.envdel = json.load(fh)
        fh.close()
        delays = []
        for st in self.envdel.keys():
            if perstation:
                delays.append(np.median(self.envdel[st]))
            else:
                delays += self.envdel[st]
        self.envdefault = {'xx':np.median(delays)}

        # picks
        fh = open(self.pkdelfn)
        self.pkdel = json.load(fh)
        fh.close()
        delays = []
        for st in self.pkdel.keys():
            if perstation:
                delays.append(np.median(self.pkdel[st]))
            else:
                delays += self.pkdel[st]
        self.pkdefault = {'xx':np.median(delays)}

        # associator
        a = np.load(self.octdelfn)
        self.ascdel = a['delays']

        # magnitude
        a = np.load(self.magdelfn)
        self.magdel = a['delays']

    def get_pick_delays(self):
        return self.pkdel

    def get_pick_default(self):
        return self.pkdefault

    def get_envelope_default(self):
        return self.envdefault

    def get_envelope_delays(self):
        return self.envdel

    def get_associator_delays(self):
        return self.ascdel

    def get_magnitude_delays(self):
        return self.magdel


def statistics(lat, lon, at, ttP, mapcoord):
    perc_min_1s = 16
    perc_max_1s = 84
    perc_min_2s = 2.5
    perc_max_2s = 97.5
    vals_1s = []
    vals_2s = []
    for idx in xrange(ttP.shape[1]):
        med = np.median(ttP[:, idx])
        lb1s = scoreatpercentile(ttP[:, idx], perc_min_1s)
        ub1s = scoreatpercentile(ttP[:, idx], perc_max_1s)
        lb2s = scoreatpercentile(ttP[:, idx], perc_min_2s)
        ub2s = scoreatpercentile(ttP[:, idx], perc_max_2s)
        val1s = (at[idx] - lb1s) / (ub1s - lb1s)
        vals_1s.append(val1s)
        val2s = (at[idx] - lb2s) / (ub2s - lb2s)
        vals_2s.append(val2s)
    # 1 sigma
    idx1s = np.where((np.array(vals_1s) <= 1.0) & (np.array(vals_1s) >= 0.0))
    idx1s_fast = np.where(np.array(vals_1s) < 0.0)
    idx1s_slow = np.where(np.array(vals_1s) > 1.0)

    # 2 sigma
    idx2s = np.where((np.array(vals_2s) <= 1.0) & (np.array(vals_2s) >= 0.0))
    idx2s_fast = np.where(np.array(vals_2s) < 0.0)
    idx2s_slow = np.where(np.array(vals_2s) > 1.0)
    print "%d events analysed." % lat.size
    print "%.1f%% lie within the %dth and %dth percentile (%d events)" % \
    ((float(idx1s[0].size) / len(vals_1s)) * 100,
     perc_min_1s, perc_max_1s, idx1s[0].size)
    print "--> %.1f%% are below the %dth percentile (%d events)" % \
    ((float(idx1s_fast[0].size) / len(vals_1s)) * 100, perc_min_1s,
     idx1s_fast[0].size)
    print "--> %.1f%% are above the %dth percentile (%d events)" % \
    ((float(idx1s_slow[0].size) / len(vals_1s)) * 100,
     perc_max_1s, idx1s_slow[0].size)
    print "%.1f%% lie within the %.1fth and %.1fth percentile (%d events)" % \
    ((float(idx2s[0].size) / len(vals_2s)) * 100, perc_min_2s, perc_max_2s,
     idx2s[0].size)
    print "--> %.1f%% are below the %.1fth percentile (%d events)" % \
    ((float(idx2s_fast[0].size) / len(vals_2s)) * 100, perc_min_2s,
     idx2s_fast[0].size)
    print "--> %.1f%% are above the %.1fth percentile (%d events)" % \
    ((float(idx2s_slow[0].size) / len(vals_2s)) * 100, perc_max_2s,
     idx2s_slow[0].size)

    if not any(x is None for x in mapcoord):
        try:
            import matplotlib.pyplot as plt
            from matplotlib.colors import Normalize
            from mpl_toolkits.basemap import Basemap
        except:
            print "To plot maps the Python packages matplotlib \
            and basemap have to be installed."
            raise()
        latmin, lonmin, latmax, lonmax = mapcoord
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111)
        m = Basemap(projection='merc', llcrnrlat=latmin,
                    urcrnrlat=latmax, llcrnrlon=lonmin,
                    urcrnrlon=lonmax, lat_ts=(latmin + latmax) / 2.,
                    resolution='h', ax=ax)
        cmap = plt.cm.get_cmap('RdBu_r')
        m.drawcoastlines()
        m.drawstates()
        m.drawcountries()
        x, y = m(lon, lat)
        sc = m.scatter(x, y, c=vals_2s, cmap=cmap,
                       norm=Normalize(vmin=-1.0, vmax=2.0), linewidths=(0,),
                       alpha=0.6)
        cb = plt.colorbar(sc, extend='both')
        cb.set_label('Alert time deviation')
        ax.set_title("Deviation from the 2.5th and 97.5th percentile")
        plt.show()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Compute the expected alert \
    times for a set of earthquakes based on observed system delay times using \
    the method of Behr et al. [2015]. Note that by default strong motion \
    channels are ignored as they tend to not contribute picks to most \
    (moderate) size earthquakes --unless they have been specifically configured \
    to do so -- and therefore introduce a bias towards \
    faster alert times. Also, by default events below 40 km are ignored \
    as the assumption of a constant P-wave velocity of 6.5 km/s is likely \
    not valid anymore for deep events.")
    parser.add_argument('netf', help="File path to ascii file containing the \
    station coordinates.")
    parser.add_argument('catf', help="File path to csv file containing the \
    hypocenter coordinates.")
    parser.add_argument('pkdelf', help="File path to json file containing \
    pick delays.")
    parser.add_argument('envdelf', help="File path to json file containing \
    envelope delays.")
    parser.add_argument('octdelf', help="File path to a npz file containing \
    origin creation time delays.")
    parser.add_argument('magdelf', help="File path to json file containing \
    magnitude computation delays.")
    parser.add_argument('-o', '--out', help="File path to store results.",
                        default='simulated_alert_times.npz')
    parser.add_argument('-a', '--archive', help="Use stored results instead of\
    running the simulation again.", action='store_true')
    parser.add_argument('--maxdepth', help='Define a cutoff depth below which \
    events are ignored.', default=40.0, type=float)
    parser.add_argument('--usesm',
                        help='Include strong motion channels [HN/HG].',
                        action='store_true')
    pm = parser.add_argument_group('Map options', 'Maps will only be plotted \
    if all of the following optional arguments are passed.')
    pm.add_argument('--latmin', help="Latitude of the map's lower left corner",
                    default=None, type=float)
    pm.add_argument('--lonmin', help="Longitude of the map's lower left corner",
                    default=None, type=float)
    pm.add_argument('--latmax', help="Latitude of the map's upper right corner",
                    default=None, type=float)
    pm.add_argument('--lonmax', help="Longitude of the map's upper right corner",
                    default=None, type=float)
    args = parser.parse_args()
    ni = NetworkInfo()
    ni.read(args.netf, sm=args.usesm)
    do = DelayObs(args.pkdelf, args.envdelf, args.octdelf, args.magdelf)
    do.load(perstation=True)
    elat, elon, edep, at = np.loadtxt(args.catf, unpack=True, delimiter=',',
                                      usecols=(2, 3, 4, 5))
    idx = np.where(edep <= args.maxdepth)
    elon = elon[idx]
    elat = elat[idx]
    edep = edep[idx]
    at = at[idx]
    de = DelayEEW()
    if not args.archive:
        lat, lon, dep, ttP, tstarget = \
                de.compute(ni, elon, elat, edep, vp=6.5, vs=3.5, nnst=6,
                           procdelay=True, nmaps=500, resultsfn=args.out,
                           latencies=do)
    else:
        a = np.load(args.out)
        ttP = a['ttP']
        lat = a['lat']
        lon = a['lon']
    statistics(lat, lon, at, ttP,
               (args.latmin, args.lonmin, args.latmax, args.lonmax))
