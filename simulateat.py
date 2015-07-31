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
import matplotlib.pyplot as plt


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

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Compute the expected alert \
    times for a set of earthquakes based on observed system delay times using \
    the method of Behr et al. [2015].")
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
    args = parser.parse_args()
    ni = NetworkInfo()
    ni.read(args.netf, sm=True)
    do = DelayObs(args.pkdelf, args.envdelf, args.octdelf, args.magdelf)
    do.load(perstation=True)
    elat, elon, edep, at = np.loadtxt(args.catf, unpack=True, delimiter=',',
                                      usecols=(2, 3, 4, 5))
    de = DelayEEW()
    if not args.archive:
        lat, lon, dep, ttP, tstarget = \
                de.compute(ni, elon, elat, edep, vp=6.5, vs=3.5, nnst=6,
                           procdelay=True, nmaps=500, resultsfn=args.out,
                           latencies=do)
