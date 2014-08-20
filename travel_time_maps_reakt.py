#!/usr/bin/env python
"""
A map showing the different P-wave travel times
for Romania, Switzerland, Greece, Turkey, Iceland
and New Zealand.

Created on Jun 13, 2014

@author: behry
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('WXAgg')
import wx
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from matplotlib.pyplot import cm
from mpl_toolkits.basemap import Basemap


def plot_ptt(datafn, bmap, cmap, dlevel, vmin=0., vmax=15., extend='max'):
    a = np.load(datafn)
    lat = a['lat']
    lon = a['lon']
    ttP = a['ttP']
    ttPmed = np.median(ttP, axis=2)
    xmin, ymin = bmap(lon.min(), lat.min())
    xmax, ymax = bmap(lon.max(), lat.max())
    nx = int(abs(xmax - xmin) / 5000.) + 1
    ny = int(abs(ymax - ymin) / 5000.) + 1
    dat, x, y = bmap.transform_scalar(ttPmed.T, lon, lat, nx, ny,
                                   returnxy=True, masked=True)
    dat_filt = np.where(dat > vmax, np.nan, dat)
    bmap.contourf(x, y, dat_filt, cmap=cmap,
                  levels=np.arange(vmin, vmax + dlevel, dlevel),
                  norm=Normalize(vmin=vmin, vmax=vmax), extend=extend,
                  zorder=1)


def plot_at(datafn, bmap, cmap, vmin=0., vmax=15., extend='neither'):
    lat, lon, at = np.loadtxt(datafn, unpack=True, delimiter=',',
                              usecols=(2, 3, 5))
    x, y = bmap(lon, lat)
    bmap.scatter(x, y, c=at, cmap=cmap,
                 norm=Normalize(vmin=vmin, vmax=vmax), zorder=1)

def main(traveltime=True, alerttime=False):
    lonmin, lonmax, latmin, latmax = 2.5, 37.5, 35.0, 48.0
    dlat = 2.0
    dlon = 2.0
    meridians = np.arange(2.0, 40, 4.0)
    parallels = np.arange(33, 48, 2.0)
    cmapname = 'RdBu_r'
    cmap = cm.get_cmap(cmapname)
    cmap.set_over('grey')
    extend = 'max'
    dlevel = 0.5
    datadir = './data/'

    if traveltime:
        vmin = 0.
        vmax = 15.
        fout = 'plots/travel_time_maps_reakt.png'
        cb_label = 'P-wave travel time to %d stations [s]' % 6

    if alerttime:
        vmin = 10.
        vmax = 60.
        fout = 'plots/alert_time_maps_reakt.png'
        cb_label = 'Initial alert time [s]'


    # fig = plt.figure(figsize=(12, 7))
    # without Iceland
    fig = plt.figure(figsize=(11, 7))
    ax = fig.add_axes([0.1, 0., .8, 1.0])

    # setup albers equal area conic basemap
    # lat_1 is first standard parallel.
    # lat_2 is second standard parallel.
    # lon_0,lat_0 is central point.
    if True:
        # m = Basemap(width=4000000, height=2000000,
        #            resolution='l', projection='aea', \
        #            lat_1=35., lat_2=48, lon_0=20, lat_0=42, ax=ax)
        # without Iceland
        m = Basemap(width=3500000, height=2000000,
                    resolution='l', projection='aea', \
                    lat_1=35., lat_2=48, lon_0=20, lat_0=42, ax=ax)
        m.drawmeridians(meridians, labels=[0, 0, 0, 1], color='lightgray',
                        linewidth=0.5, zorder=0)
        m.drawparallels(parallels, labels=[1, 0, 0, 0], color='lightgray',
                        linewidth=0.5, zorder=0)
        m.drawcoastlines(zorder=2)
        m.drawcountries(linewidth=1.0, zorder=2)
        m.fillcontinents('gray', zorder=0)

        if traveltime:
            # plot Romanian data
            resultsfn = os.path.join(datadir, 'ptt_ro_6stations.npz')
            plot_ptt(resultsfn, m, cmap, dlevel, vmax=vmax)

            # plot Greek data
            resultsfn = os.path.join(datadir, 'p_wave_tt_gr.npz')
            plot_ptt(resultsfn, m, cmap, dlevel, vmax=vmax)

            # plot Swiss data
            resultsfn = os.path.join(datadir, 'ptt_ch_6stations.npz')
            plot_ptt(resultsfn, m, cmap, dlevel, vmax=vmax)

            # plot Turkish data
            resultsfn = os.path.join(datadir, 'ptt_koeri_6stations.npz')
            plot_ptt(resultsfn, m, cmap, dlevel, vmax=vmax)

        if alerttime:
            # plot Romanian data
            resultsfn = os.path.join(datadir, 'event_list_romania.csv')
            plot_at(resultsfn, m, cmap, vmin=vmin, vmax=vmax)

            # plot Greek data
            resultsfn = os.path.join(datadir, 'event_list_patras.csv')
            plot_at(resultsfn, m, cmap, vmin=vmin, vmax=vmax)

            # plot Swiss data
            resultsfn = os.path.join(datadir, 'event_list_ch.csv')
            plot_at(resultsfn, m, cmap, vmin=vmin, vmax=vmax)

            # plot Turkish data
            resultsfn = os.path.join(datadir, 'event_list_turkey.csv')
            plot_at(resultsfn, m, cmap, vmin=vmin, vmax=vmax)

    # Create an inset for Iceland
    if False:
        ax_ice = fig.add_axes([0.64, 0.55, .25, .25])
        mi = Basemap(width=550000, height=500000,
                    resolution='l', projection='aea', \
                    lat_1=62.5, lat_2=67.5, lon_0=-19, lat_0=65,
                    ax=ax_ice)
        mi.drawcoastlines(zorder=2)
        mi.fillcontinents('gray', zorder=0)
        mi.drawmeridians(np.arange(-26, -12, 4), labels=[0, 0, 1, 0],
                        color='lightgray', linewidth=0.5, zorder=0)
        mi.drawparallels(np.arange(60, 70, 2), labels=[0, 1, 0, 0],
                        color='lightgray', linewidth=0.5, zorder=0)
        # plot Iceland data
        if traveltime:
            resultsfn = os.path.join(datadir, 'ptt_imo_6stations.npz')
            plot_ptt(resultsfn, mi, cmap, dlevel, vmin=vmin, vmax=vmax)
        # if alerttime:
        #    resultsfn = os.path.join(datadir, 'event_list_.csv')
        #    plot_at(resultsfn, mi, vmin=vmin, vmax=vmax)


    # Create an inset for New Zealand
    if True:
        ax_nz = fig.add_axes([0.05, 0.2, .4, .4])
        mnz = Basemap(width=1300000, height=1700000,
                    resolution='l', projection='aea', \
                    lat_1=-50., lat_2=-32, lon_0=172, lat_0=-41,
                    ax=ax_nz)
        mnz.drawcoastlines(zorder=2)
        mnz.fillcontinents('gray', zorder=0)
        mnz.drawmeridians(np.arange(164, 182, 4), labels=[0, 0, 1, 0],
                        color='lightgray', linewidth=0.5, zorder=0)
        mnz.drawparallels(np.arange(-51, -31, 2), labels=[0, 1, 0, 0],
                        color='lightgray', linewidth=0.5, zorder=0)
        # plot NZ data
        if traveltime:
            resultsfn = os.path.join(datadir, 'p_wave_tt_nz.npz')
            plot_ptt(resultsfn, mnz, cmap, dlevel, vmin=vmin, vmax=vmax)
        if alerttime:
            resultsfn = os.path.join(datadir, 'event_list_nz.csv')
            plot_at(resultsfn, mnz, cmap, vmin=vmin, vmax=vmax)

    cax = fig.add_axes([0.2, 0.07, 0.6, 0.02])
    cb = ColorbarBase(cax, cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax),
                      orientation='horizontal', extend=extend)
    cb.set_label(cb_label)
    fig.savefig(fout, dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    main()
