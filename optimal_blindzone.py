#!/usr/bin/env python
"""
Compute the optimal blindzone for the given seismicity and network geometry.
Created on Feb 23, 2015

@author: behry
"""
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import os
import pyproj
import sys
sys.path.append('../delays')
from delayeew import DelayEEW
from alerttimemap import AlertTimeMap
import json
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize


class NetworkInfo:
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
                               'color':'white', 'label':''}

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
    Dummy class to load processing latencies.
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

class DelayObsOpt:
    """
    Return optimal delays.
    """

    def __init__(self, networkinfo, nsamples):
        self.stations = []
        self.pkdel = {}
        self.envdel = {}
        self.ascdel = {}
        self.magdel = {}
        for _nt in networkinfo.networks.keys():
            self.stations += networkinfo.networks[_nt]['nm']
        for _st in self.stations:
            self.pkdel[_st] = 2.0 * np.random.random_sample(nsamples) + 0.1
            self.envdel[_st] = 2.0 * np.random.random_sample(nsamples) + 0.1
        self.ascdel = np.random.normal(0.7, 0.5, nsamples)
        self.magdel = np.random.normal(0.7, 0.5, nsamples)
        self.pkdefault = np.median(2.0 * np.random.random_sample(nsamples) + 0.1)
        self.envdefault = np.median(2.0 * np.random.random_sample(nsamples) + 0.1)

    def load(self):
        pass

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


def plot_bzreduction(fn, lat, lon, dep, at, bmap, cmap, vmin=0., vmax=15., extend='neither',
                     **kargs):
    cmap = plt.cm.ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax),
                             cmap=cmap)
    a = np.load(fn)
    ttP = a['ttP']
    olat = a['lat']
    olon = a['lon']
    odep = a['dep']
    vs = 3.5
    bz2 = at * at * vs * vs - dep * dep
    bz2 = np.where(bz2 > 0.0, bz2, 0.0)
    bz = np.sqrt(bz2)
    at = bz
    ttPmed = np.median(ttP, axis=0)
    obz2 = ttPmed * ttPmed * vs * vs - dep * dep
    obz2 = np.where(obz2 > 0.0, obz2, 0.0)
    obz = np.sqrt(obz2)
    redbz = bz - obz
    at = redbz
    x, y = bmap(lon, lat)
    for lt, ln, rbz in zip(lat, lon, redbz):
        cl = cmap.to_rgba(rbz)
        x, y = bmap(ln, lt)
        bmap.plot(x, y, ms=8, c=cl, marker='o', picker=5.)
    bmap.drawcoastlines(zorder=2)
    bmap.drawcountries(linewidth=1.0, zorder=2)
    bmap.drawstates(zorder=2)


def plot_bz(**kargs):
    elat, elon, edep, at = np.loadtxt(kargs['events'], unpack=True, delimiter=',',
                                      usecols=(2, 3, 4, 5))
    if 'geofilter' in kargs:
        msk = kargs['geofilter'].point_in_polygon(elon, elat)
        elon = elon[np.where(msk)]
        elat = elat[np.where(msk)]
        edep = edep[np.where(msk)]
        at = at[np.where(msk)]
    ni = NetworkInfo()
    ni.read(kargs['stations'])
    do = DelayObsOpt(ni, 1000)
    do.load()
    de = DelayEEW()
    if kargs['new']:
        lat, lon, dep, ttP, tstarget = \
        de.compute(ni, elon, elat, edep, vp=6.5, vs=3.5,
                   nnst=2, procdelay=True, nmaps=500, resultsfn=kargs['fout'],
                   latencies=do)
    print lat.shape
    if 'ax' not in kargs and 'fig' not in kargs:
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_axes([0.08, 0.35, 0.8, 0.6])
    else:
        ax = kargs['ax']
        fig = kargs['fig']
    if 'cax' not in kargs:
        cax = fig.add_axes([0.83, 0.35, 0.05, 0.6])
    else:
        cax = kargs['cax']

    if 'atmap' in kargs:
        am1 = AlertTimeMap(kargs['atmap'])
        am1.plot2D(fig, ax, procdelay=True, scale=True, blindzone=True)
    if kargs['at']:
        ttP = np.zeros((1, at.size))
        ttP[0] = at[:]
        am = AlertTimeMap(mapbox=kargs['mapbnd'], cmapname='RdBu_r', ax=ax,
                          ttP=ttP, lat=elat, lon=elon, dep=edep, tstarget=None)
        # am.plot_stations(fig, ax, ni)
    elif kargs['bzreduction']:
        am = AlertTimeMap(resultsfn=kargs['fout'], mapbox=kargs['mapbnd'],
                          cmapname='RdBu_r', ax=ax)
        plot_bzreduction(kargs['fout'], elat, elon, edep, at,
                         am.m, 'RdBu_r', vmin=kargs['vmin'],
                         vmax=kargs['vmax'])
    else:
        am = AlertTimeMap(resultsfn=kargs['fout'], mapbox=kargs['mapbnd'],
                          cmapname='RdBu_r', ax=ax)
        am.plot_stations(fig, ax, ni)
    if not kargs['bzreduction']:
        am.plot1D(vmin=kargs['vmin'], vmax=kargs['vmax'], blindzone=True)
    cb = ColorbarBase(cax, cmap='RdBu_r',
                      norm=Normalize(vmin=kargs['vmin'], vmax=kargs['vmax']))
    cb.set_label('Blind zone radius [km]')
    return am, cb

if __name__ == '__main__':
    pass
