#!/usr/bin/env python
"""
Created on Jun 15, 2015

@author: behry
"""
import matplotlib
matplotlib.use('WXAgg')
import matplotlib.pyplot as plt
import numpy as np
import pyproj
from collections import defaultdict
from scipy.stats import scoreatpercentile

txt_fontsize = 14


class MagErrVSTime:

    def __init__(self):
        self.fig = None
        self.ax = None
        self.debug = False
        self.fndict = {'Turkey':{'file':'./data/event_list_turkey.csv', 'correction':None},
                       'New Zealand':{'file':'./data/event_list_nz.csv', 'correction':None},
                       'Switzerland':{'file':'./data/event_list_ch.csv', 'correction':'%s+0.25'},
                       'Romania':{'file':'./data/event_list_romania.csv', 'correction':'0.83*((%s-0.8)/0.74)+0.17'},
                       'Patras':{'file':'./data/event_list_patras.csv', 'correction':None},
                       'Iceland':{'file':'./data/event_list_iceland.csv',
                                  'file_bb':'./data/event_list_iceland_bardarbunga.csv',
                                  'correction':None},
                       'California':{'file':'./data/event_list_ca.csv', 'correction':None}}
        self.g = pyproj.Geod(ellps='WGS84')

    def read_data(self):
        pass

    def stats(self, tdiff, merr):
        # bin the data
        t_bins = np.linspace(0, 60, 61)
        d = defaultdict(list)
        med = np.zeros(t_bins.size)
        ub = np.zeros(t_bins.size)
        lb = np.zeros(t_bins.size)
        for _dt, _merr in zip(tdiff, merr):
            idx = np.argmin(np.abs(t_bins - _dt))
            d[idx].append(_merr)
        for _i in xrange(t_bins.size):
            if not _i in d:
                med[_i] = np.nan
                ub[_i] = np.nan
                lb[_i] = np.nan
            # elif len(d[_i]) < 2:
            #    med[_i] = np.nan
            #    ub[_i] = np.nan
            #    lb[_i] = np.nan
            else:
                med[_i] = np.median(d[_i])
                ub[_i] = scoreatpercentile(d[_i], 84)
                lb[_i] = scoreatpercentile(d[_i], 16)
        return t_bins, med, ub, lb

    def mag_err(self, ax, depdisc=40., magdisc=0.0):
        minmag = 1e39
        maxmag = 0.
        nevents = 0
        allml = np.array([])
        allmvs = np.array([])
        alldist = np.array([])
        alldep = np.array([])
        alltdiff = np.array([])

        for _c in self.fndict.keys():
            fn = self.fndict[_c]['file']
            dep, tdiff, depvs, mvs, ml, lon, lat, lonvs, latvs = \
            np.loadtxt(fn, unpack=True, delimiter=',',
                       usecols=(4, 5, 8, 9, 1, 3, 2, 7, 6))
            if self.fndict[_c]['correction'] is not None:
                mvs = eval(self.fndict[_c]['correction'] % 'mvs')
            idx = np.where((dep <= depdisc) & (ml > magdisc))
            if idx[0].size < 2:
                continue
            nevents += idx[0].size
            minmag = min(ml[idx].min(), minmag)
            maxmag = max(ml[idx].max(), maxmag)
            az, baz, dist = self.g.inv(lon[idx], lat[idx], lonvs[idx], latvs[idx])
            dist /= 1000.
            ddep = np.abs(dep[idx] - depvs[idx])
            dist = np.sqrt(dist * dist + ddep * ddep)
            allml = np.hstack((allml, ml[idx]))
            allmvs = np.hstack((allmvs, mvs[idx]))
            alldist = np.hstack((alldist, dist))
            alldep = np.hstack((alldep, dep[idx]))
            alltdiff = np.hstack((alltdiff, tdiff[idx]))
            txt = r'shallow events ($<$ 40 km)'
        tbins, med, ub, lb = self.stats(alltdiff, allml - allmvs)
        ax.plot(tbins, med, ls='-', color='k', lw=2)
        ax.plot(tbins, ub, 'k--', lw=3)
        ax.plot(tbins, lb, 'k--', lw=3)
        ax.hlines(0, 0, 60)

if __name__ == '__main__':
    # fout = './plots/mag_comp.pdf'
    fig = plt.figure()
    ax = fig.add_subplot(111)
    me = MagErrVSTime()
    me.mag_err(ax)
    plt.show()

