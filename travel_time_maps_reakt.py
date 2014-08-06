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
from scipy.stats import scoreatpercentile


def plot_ptt(datafn, bmap, vmin=0., vmax=15.):
    a = np.load(datafn)
    lat = a['lat']
    lon = a['lon']
    ttP = a['ttP']
    ttPmed = np.median(ttP, axis=2)
    xmin, ymin = m(lon.min(), lat.min())
    xmax, ymax = m(lon.max(), lat.max())
    nx = int(abs(xmax - xmin) / 5000.) + 1
    ny = int(abs(ymax - ymin) / 5000.) + 1
    dat, x, y = bmap.transform_scalar(ttPmed.T, lon, lat, nx, ny,
                                   returnxy=True, masked=True)
    dat_filt = np.where(dat > vmax, np.nan, dat)
    bmap.contourf(x, y, dat_filt, cmap=cmap,
                  levels=np.arange(vmin, vmax + dlevel, dlevel),
                  norm=Normalize(vmin=vmin, vmax=vmax), extend=extend,
                  zorder=1)

lonmin, lonmax, latmin, latmax = 2.5, 37.5, 35.0, 48.0
dlat = 2.0
dlon = 2.0
meridians = np.arange(2.0, 40, 4.0)
parallels = np.arange(33, 48, 2.0)
cmapname = 'RdBu_r'
cmap = cm.get_cmap(cmapname)
cmap.set_over('grey')
extend = 'max'
vmin = 0.
vmax = 15.
cb_label = 'P-wave travel time to %d stations [s]' % 6
dlevel = 0.5

fig = plt.figure(figsize=(12, 7))
ax = fig.add_axes([0.1, 0., .8, 1.0])

# setup albers equal area conic basemap
# lat_1 is first standard parallel.
# lat_2 is second standard parallel.
# lon_0,lat_0 is central point.
if True:
    m = Basemap(width=4000000, height=2000000,
                resolution='l', projection='aea', \
                lat_1=35., lat_2=48, lon_0=20, lat_0=42, ax=ax)
    m.drawmeridians(meridians, labels=[0, 0, 0, 1], color='lightgray',
                    linewidth=0.5, zorder=0)
    m.drawparallels(parallels, labels=[1, 0, 0, 0], color='lightgray',
                    linewidth=0.5, zorder=0)
    m.drawcoastlines(zorder=2)
    m.drawcountries(linewidth=1.0, zorder=2)
    m.fillcontinents('gray', zorder=0)

    # plot Romanian data
    datadir = '/home/behry/workspace/eew/delays/data/reakt_publication'
    resultsfn = os.path.join(datadir, 'ptt_ro_6stations.npz')
    plot_ptt(resultsfn, m, vmax=vmax)

    # plot Greek data
    datadir = '/home/behry/workspace/eew/delays/data/reakt_publication'
    resultsfn = os.path.join(datadir, 'p_wave_tt_gr.npz')
    plot_ptt(resultsfn, m, vmax=vmax)

    # plot Swiss data
    datadir = '/home/behry/workspace/eew/delays/data/reakt_publication'
    resultsfn = os.path.join(datadir, 'ptt_ch_6stations.npz')
    plot_ptt(resultsfn, m, vmax=vmax)

    # plot Turkish data
    datadir = '/home/behry/workspace/eew/delays/data/reakt_publication'
    resultsfn = os.path.join(datadir, 'ptt_koeri_6stations.npz')
    plot_ptt(resultsfn, m, vmax=vmax)

# Create an inset for Iceland
if True:
    ax_ice = fig.add_axes([0.64, 0.55, .25, .25])
    mi = Basemap(width=550000, height=500000,
                resolution='l', projection='aea', \
                lat_1=62.5, lat_2=67.5, lon_0=-19, lat_0=65,
                ax=ax_ice)
    mi.drawcoastlines(zorder=2)
    mi.fillcontinents('gray', zorder=0)
    # plot Iceland data
    datadir = '/home/behry/workspace/eew/delays/data/reakt_publication/'
    resultsfn = os.path.join(datadir, 'ptt_imo_6stations.npz')
    plot_ptt(resultsfn, mi, vmin=vmin, vmax=vmax)


# Create an inset for New Zealand
if True:
    ax_nz = fig.add_axes([0.05, 0.2, .4, .4])
    mnz = Basemap(width=1300000, height=1700000,
                resolution='l', projection='aea', \
                lat_1=-50., lat_2=-32, lon_0=172, lat_0=-41,
                ax=ax_nz)
    mnz.drawcoastlines(zorder=2)
    mnz.fillcontinents('gray', zorder=0)
    # plot NZ data
    datadir = '/home/behry/workspace/eew/delays/data/reakt_publication'
    resultsfn = os.path.join(datadir, 'p_wave_tt_nz.npz')
    plot_ptt(resultsfn, mnz, vmin=vmin, vmax=vmax)

cax = fig.add_axes([0.2, 0.07, 0.6, 0.02])
cb = ColorbarBase(cax, cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax),
                  orientation='horizontal', extend=extend)
cb.set_label(cb_label)
fout = 'data/reakt_publication/travel_time_maps_reakt.png'
fig.savefig(fout, dpi=300, bbox_inches='tight')
plt.show()
