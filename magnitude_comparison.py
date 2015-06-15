#!/usr/bin/env python
"""
Plot the magnitude comparisons in a composite figure.
Created on Nov 24, 2014

@author: behry
"""
import matplotlib
matplotlib.use('WXAgg')
import matplotlib.pyplot as plt
import numpy as np
import pyproj
from collections import defaultdict
from scipy.stats import scoreatpercentile
import ipdb

txt_fontsize = 14

class MagComp:

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

    def stats(self, ml, mvs):
        # bin the data
        mag_bins = np.linspace(2, 7, 11)
        d = defaultdict(list)
        med = np.zeros(mag_bins.size)
        ub = np.zeros(mag_bins.size)
        lb = np.zeros(mag_bins.size)
        for _ml, _mvs in zip(ml, mvs):
            idx = np.argmin(np.abs(mag_bins - _ml))
            d[idx].append(_mvs - _ml)
        for _i in xrange(mag_bins.size):
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
        return mag_bins, med, ub, lb

    def setup(self):
        """
        Set up the panels to plot in.
        """
        width = 18.0
        height = 18.0
        self.fig = plt.figure(1, (width, height))
        self.ax = []
        vmargin = 0.7 / width
        hmargin = 0.5 / width
        vpad = 0.7 / height
        shpad = 0.6 / width
        hpad = 0.5 / width
        spht = 3.2 / height
        spwd = 5.7 / width
        bpht = 3.5 / height
        bpwd = 6. / width
        # bottom row
        x = vmargin
        y = hmargin
        w = spwd
        h = spht
        self.ax.append(self.fig.add_axes((x, y, w, h)))
        x += shpad + spwd
        self.ax.append(self.fig.add_axes((x, y, w, h)))
        x += shpad + spwd
        self.ax.append(self.fig.add_axes((x, y, w, h)))
        plt.setp(self.ax[1].get_yticklabels(), visible=False)
        plt.setp(self.ax[2].get_yticklabels(), visible=False)

        # middle row
        x = vmargin
        w = bpwd
        h = bpht
        y += vpad + spht
        self.ax.append(self.fig.add_axes((x, y, w, h)))
        x += hpad + bpwd - 0.8 * hpad
        self.ax.append(self.fig.add_axes((x, y, w, h)))
        x += hpad + bpwd
        self.ax.append(self.fig.add_axes((x, y, spwd, spht)))

        # top row
        x = vmargin
        y += vpad + bpht
        w = spwd
        h = spht
        self.ax.append(self.fig.add_axes((x, y, w, h)))
        x += shpad + spwd
        self.ax.append(self.fig.add_axes((x, y, w, h)))
        x += shpad + spwd
        self.ax.append(self.fig.add_axes((x, y, w, h)))
        plt.setp(self.ax[7].get_yticklabels(), visible=False)
        plt.setp(self.ax[8].get_yticklabels(), visible=False)
        # plt.show()

    def plot_mag(self, ax, ml, mvs, dist, dep, xmin=2.0, xmax=7.0, fact=2000.0,
                 cbo='horizontal', marker='o', legend=True, minmag=2.0,
                 cb_aspect=20):
        dist = np.where(dist < 10.0, 10.0, dist)
        idx = np.where((dist < 100.0) & (ml >= minmag))
        if idx[0].size > 0:
            sc = ax.scatter(ml[idx], mvs[idx] - ml[idx], marker=marker, c=dep[idx],
                            s=fact / dist[idx], linewidths=0, alpha=0.5)
        idx = np.where((dist >= 100.0) & (ml >= minmag))
        if idx[0].size > 0:
            sc = ax.scatter(ml[idx], mvs[idx] - ml[idx], marker='+', c=dep[idx],
                            s=100.0, linewidths=1, alpha=0.5)
        if legend:
            small_error = np.where((dist < 100.) & (ml > minmag))
            good = np.where(np.abs(ml[small_error] - mvs[small_error]) <= 0.5)
            idx = np.where(ml >= minmag)
            mn = np.mean(mvs[idx] - ml[idx])
            std = np.std(mvs[idx] - ml[idx])
            if self.debug:
                print "Number of events: %d; Minimum magnitude: %.2f; Maximum magnitude: %.2f; Mean difference: %.2f; Std (difference): %.2f" % (ml[idx].size, ml[idx].min(), ml[idx].max(), mn, std)
                print "Number of events with location error <= 100: %d" % small_error[0].size
                print "---> out of these # of events within +-0.5 Ml: %d" % (good[0].size)
            cb = self.fig.colorbar(sc, orientation=cbo, ax=ax,
                                   aspect=cb_aspect)
            cb.set_label("Depth [km]")
            ax.text(5.5, -1.1, 'Location accuracy [km]',
                    horizontalalignment='center', fontsize=txt_fontsize)
            ax.scatter(4.8, -1.42, s=fact / 10, c='white')
            ax.text(4.8, -1.9, r'$\geq$ 10',
                    horizontalalignment='center')
            ax.scatter(5.4, -1.42, s=fact / 30.0, c='white')
            ax.text(5.4, -1.9, r'30', horizontalalignment='center')
            ax.scatter(5.95, -1.42, s=fact / 50, c='white')
            ax.text(5.95, -1.9, '50', horizontalalignment='center')
            ax.scatter(6.5, -1.42, s=100, c='black', marker='+')
            ax.text(6.5, -1.9, r'$\leq$ 100', horizontalalignment='center')
            mags, med, ub, lb = self.stats(ml[idx], mvs[idx])
            ax.plot(np.linspace(xmin, xmax + 0.1, 100),
                    np.zeros(100), 'k-')
            # ax.plot(np.linspace(xmin, xmax + 0.1, 100),
            #        np.zeros(100) + 0.5, 'k--')
            # ax.plot(np.linspace(xmin, xmax + 0.1, 100),
            #        np.zeros(100) - 0.5, 'k--')
            ax.plot(mags, med, 'k--')
            ax.plot(mags, ub, 'k--')
            ax.plot(mags, lb, 'k--')
            ax.set_xlabel('Ml')
            ax.set_ylabel('MVS - Ml')
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(-2, 2)

    def plot_mag_comp(self, ax, mtype='Ml', countryname='', panelnumber=''):
        fn = self.fndict[countryname]['file']
        dep, depvs, mvs, ml, lon, lat, lonvs, latvs = \
        np.loadtxt(fn, unpack=True, delimiter=',',
                   usecols=(4, 8, 9, 1, 3, 2, 7, 6))
        if self.fndict[countryname]['correction'] is not None:
            mvs = eval(self.fndict[countryname]['correction'] % 'mvs')
        az, baz, dist = self.g.inv(lon, lat, lonvs, latvs)
        dist /= 1000.
        ddep = np.abs(dep - depvs)
        dist = np.sqrt(dist * dist + ddep * ddep)
        if self.debug:
            print "Country: ", countryname
        ax.text(5.3, 1.5, countryname, horizontalalignment='left',
                verticalalignment='center', fontsize=txt_fontsize)
        ax.text(4.9, 1.5, panelnumber, horizontalalignment='right',
                verticalalignment='center', fontsize=txt_fontsize)
        self.plot_mag(ax, ml, mvs, dist, dep, cb_aspect=30)

    def shallow_deep_mag_comp(self, ax, depdisc=40., magdisc=0.0, shallow=True,
                              panelnumber=''):
        minmag = 1e39
        maxmag = 0.
        nevents = 0
        allml = np.array([])
        allmvs = np.array([])
        alldist = np.array([])
        alldep = np.array([])
        for _c in self.fndict.keys():
            fn = self.fndict[_c]['file']
            dep, depvs, mvs, ml, lon, lat, lonvs, latvs = \
            np.loadtxt(fn, unpack=True, delimiter=',',
                       usecols=(4, 8, 9, 1, 3, 2, 7, 6))
            if self.fndict[_c]['correction'] is not None:
                mvs = eval(self.fndict[_c]['correction'] % 'mvs')
            if shallow:
                idx = np.where((dep <= depdisc) & (ml > magdisc))
            else:
                idx = np.where((dep > depdisc) & (ml > magdisc))
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
        if shallow:
            txt = r'shallow events ($<$ 40 km)'
        else:
            txt = r'deep events ($\geq$ 40 km)'
        ax.text(4.4, 1.5, txt, horizontalalignment='left',
                verticalalignment='center', fontsize=txt_fontsize)
        ax.text(4.2, 1.5, panelnumber, horizontalalignment='right',
                verticalalignment='center', fontsize=txt_fontsize)
        self.plot_mag(ax, allml, allmvs, alldist, alldep, cbo='vertical')

    def add_bardarbunga(self, ax, depdisc=None, shallow=False, deep=False):
        """
        The events related to the volcanic activity at Bardarbunga volcano
        contain very low frequencies and MVS is therefore underestimated. To
        discriminate them from 'normal' crustal events we plot them with
        different symbols.
        """
        fn = self.fndict['Iceland']['file_bb']
        dep, depvs, mvs, ml, lon, lat, lonvs, latvs = \
        np.loadtxt(fn, unpack=True, delimiter=',',
                    usecols=(4, 8, -2, 1, 3, 2, 7, 6))
        if self.fndict['Iceland']['correction'] is not None:
            mvs = eval(self.fndict['Iceland']['correction'] % 'mvs')
        idx = [np.arange(dep.size)]
        if shallow and depdisc is not None:
            idx = np.where(dep <= depdisc)
            if idx[0].size < 2:
                return
        if deep and depdisc is not None:
            idx = np.where(dep > depdisc)
            if idx[0].size < 2:
                return
        az, baz, dist = self.g.inv(lon[idx], lat[idx], lonvs[idx], latvs[idx])
        dist /= 1000.
        ddep = np.abs(dep[idx] - depvs[idx])
        dist = np.sqrt(dist * dist + ddep * ddep)
        self.plot_mag(ax, ml[idx], mvs[idx], dist[idx], dep[idx], cbo='vertical',
                      marker='s', legend=False)

    def plot(self, fout):
        self.plot_mag_comp(self.ax[0], countryname='Turkey', panelnumber='g')
        self.plot_mag_comp(self.ax[1], countryname='Romania', panelnumber='h')
        self.plot_mag_comp(self.ax[2], countryname='Iceland', panelnumber='i')
        # self.add_bardarbunga(self.ax[2])
        for i in [1, 2]:
            self.ax[i].set_ylabel('')
            self.ax[i].set_yticklabels([])
        self.shallow_deep_mag_comp(self.ax[3], panelnumber='d')
        self.shallow_deep_mag_comp(self.ax[4], shallow=False, panelnumber='e')
        self.plot_mag_comp(self.ax[5], countryname='California', panelnumber='f')
        self.add_bardarbunga(self.ax[4], depdisc=40.0, deep=True)
        self.ax[4].set_ylabel('')
        self.ax[4].set_yticklabels([])
        self.plot_mag_comp(self.ax[6], countryname='Switzerland',
                           panelnumber='a')
        self.plot_mag_comp(self.ax[7], countryname='Patras',
                           panelnumber='b')
        self.plot_mag_comp(self.ax[8], countryname='New Zealand',
                           panelnumber='c')
        for i in [7, 8]:
            self.ax[i].set_ylabel('')
            self.ax[i].set_yticklabels([])
        self.fig.savefig(fout, dpi=300, bbox_inches='tight')
        plt.show()

if __name__ == '__main__':
    fout = './plots/mag_comp.pdf'
    mc = MagComp()
    mc.setup()
    mc.plot(fout)
