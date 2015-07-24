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
from matplotlib.colors import Normalize, LogNorm
from matplotlib.pyplot import cm
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap
import ipdb

rcParams['axes.labelsize'] = 14
rcParams['xtick.labelsize'] = 14
rcParams['ytick.labelsize'] = 14

def plot_ptt(datafn, bmap, cmap, dlevel, vmin=0., vmax=15., extend='max'):
    a = np.load(datafn)
    lat = a['lat']
    lon = a['lon']
    ttP = a['ttP']
    if len(ttP.shape) > 3:
        ttPmed = np.median(ttP[0, :, :, :], axis=-1)
        lat = lat[0, :, 0]
        lon = lon[:, 0, 0]
    else:
        ttPmed = np.median(ttP, axis=-1)

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


def plot_at(datafn, bmap, cmap, vmin=0., vmax=15., extend='neither',
            blindzone=False, bzreduction=False, optimalbz=False,
            ssize=50, **kargs):
    lat, lon, dep, at0 = np.loadtxt(datafn, unpack=True, delimiter=',',
                              usecols=(2, 3, 4, 5))
    if blindzone:
        vs = 3.5
        bz2 = at0 * at0 * vs * vs - dep * dep
        bz2 = np.where(bz2 > 0.0, bz2, 0.0)
        bz = np.sqrt(bz2)
        at = bz
        a = np.load(kargs['optbz'])
        ttP = a['ttP']
        olat = a['lat']
        olon = a['lon']
        odep = a['dep']
        ttPmed = np.median(ttP, axis=0)
        try:
            obz2 = ttPmed * ttPmed * vs * vs - dep * dep
        except:
            ipdb.set_trace()
        obz2 = np.where(obz2 > 0.0, obz2, 0.0)
        obz = np.sqrt(obz2)
        if bzreduction:
            redbz = bz - obz
            at = redbz
        if optimalbz:
            at = obz

    x, y = bmap(lon, lat)
    bmap.scatter(x, y, c=at, cmap=cmap, s=ssize,
                 norm=LogNorm(vmin=vmin, vmax=vmax), zorder=1,
                 linewidths=0)
#                 norm=Normalize(vmin=vmin, vmax=vmax), zorder=1,


def main(traveltime=False, alerttime=False, blindzone=False,
         bzreduction=False, optimalbz=False):
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
    lbl_fontsize = 12

    cb_label = []
    if traveltime:
        vmin = 0.
        vmax = 15.
        # fout = 'plots/travel_time_maps_reakt.png'
        fout = 'plots/travel_time_maps_reakt.pdf'
        cb_label.append('P-wave travel time to %d stations [s]' % 6)
        cb_label.append('Approximate inter-station distance [km]')

    if alerttime:
        vmin = 10.
        vmax = 60.
        # fout = 'plots/alert_time_maps_reakt.png'
        fout = 'plots/alert_time_maps_reakt.pdf'
        cb_label.append('Initial alert time [s]')

    if blindzone:
        vmin = 10.
        vmax = 200.
        cb_int = 15
        # fout = 'plots/blind_zone_maps_reakt.png'
        fout = 'plots/blind_zone_maps_reakt.pdf'
        cb_label.append('Blind zone [km]')
        cb_label.append('Alert delay ($\Delta t_{alert}$) [s]')
        cb_label.append('Magnitude with positive EEW zone')

    if bzreduction:
        cmap = cm.get_cmap('RdBu')
        cmap.set_over('grey')
        vmin = 10
        vmax = 200
        cb_int = 15
        # fout = 'plots/blind_zone_reduction_reakt.png'
        fout = 'plots/blind_zone_reduction_reakt.pdf'
        cb_label.append('Blind zone reduction [km]')
        cb_label.append('Alert delay reduction [s]')

    if optimalbz:
        vmin = 10
        vmax = 200
        cb_int = 15
        # fout = 'plots/optimal_blind_zone_reakt.png'
        fout = 'plots/optimal_blind_zone_reakt.pdf'
        cb_label.append('Optimal blind zone [km]')
        cb_label.append('Alert delay ($\Delta t_{alert}$) [s]')
        cb_label.append('Magnitude with positive EEW zone')

    # get the damage zone, that is the zone for which intensity is >= 7.0
    # for magnitudes >= 5.0
    mags, dz = damage_zone()

    fig = plt.figure(figsize=(16, 7))
    # without Iceland
    # fig = plt.figure(figsize=(10, 7))
    ax = fig.add_axes([0.05, 0., .8, 1.0])

    # setup albers equal area conic basemap
    # lat_1 is first standard parallel.
    # lat_2 is second standard parallel.
    # lon_0,lat_0 is central point.
    if True:
        m = Basemap(width=5000000, height=2300000,
                    resolution='l', projection='aea', \
                    lat_1=35., lat_2=48, lon_0=20, lat_0=42, ax=ax)

        # without Iceland
        # m = Basemap(width=3000000, height=1800000,
        #            resolution='l', projection='aea', \
        #            lat_1=35., lat_2=48, lon_0=20, lat_0=42, ax=ax)
        m.drawmeridians(meridians, labels=[0, 0, 0, 1], color='lightgray',
                        linewidth=0.5, zorder=0, fontsize=lbl_fontsize)
        m.drawparallels(parallels, labels=[1, 0, 0, 0], color='lightgray',
                        linewidth=0.5, zorder=0, fontsize=lbl_fontsize)
        m.drawcoastlines(zorder=2)
        m.drawcountries(linewidth=1.0, zorder=2)
        m.fillcontinents('lightgray', zorder=0)

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

        if alerttime or blindzone or bzreduction or optimalbz:
            if bzreduction or optimalbz:
                blindzone = True
            # plot Romanian data
            resultsfn = os.path.join(datadir, 'event_list_romania.csv')
            optbz = os.path.join(datadir, 'optimal_blindzone_ro.npz')
            plot_at(resultsfn, m, cmap, vmin=vmin, vmax=vmax,
                    blindzone=blindzone, optbz=optbz, bzreduction=bzreduction,
                    optimalbz=optimalbz)

            # plot Greek data
            resultsfn = os.path.join(datadir, 'event_list_patras.csv')
            optbz = os.path.join(datadir, 'optimal_blindzone_gr.npz')
            plot_at(resultsfn, m, cmap, vmin=vmin, vmax=vmax,
                    blindzone=blindzone, optbz=optbz, bzreduction=bzreduction,
                    optimalbz=optimalbz)

            # plot Swiss data
            resultsfn = os.path.join(datadir, 'event_list_ch.csv')
            optbz = os.path.join(datadir, 'optimal_blindzone_ch.npz')
            plot_at(resultsfn, m, cmap, vmin=vmin, vmax=vmax,
                    blindzone=blindzone, optbz=optbz, bzreduction=bzreduction,
                    optimalbz=optimalbz)

            # plot Turkish data
            resultsfn = os.path.join(datadir, 'event_list_turkey.csv')
            optbz = os.path.join(datadir, 'optimal_blindzone_tr.npz')
            plot_at(resultsfn, m, cmap, vmin=vmin, vmax=vmax,
                    blindzone=blindzone, optbz=optbz, bzreduction=bzreduction,
                    optimalbz=optimalbz)

    darkgray = (107 / 255., 107 / 255., 107 / 255.)
    txt_fontsize = 14
    # add a panel for California
    if True:
        # ax_cal = fig.add_axes([0.62, 0., .31, 1.0])
        ax_cal = fig.add_axes([0.5, 0.54, .45, .35])
        mca = Basemap(width=700000, height=750000,
                    resolution='l', projection='aea', \
                    lat_1=31, lat_2=38., lon_0=-117.5, lat_0=34.5,
                    ax=ax_cal)
        # mca = Basemap(projection='merc', llcrnrlat=31, urcrnrlat=37.5,
        #              llcrnrlon=-121.5, urcrnrlon=-114, lat_ts=34.5,
        #              resolution='i', ax=ax_cal)
        mca.drawcoastlines(zorder=2)
        mca.drawcountries(zorder=2)
        mca.drawstates(zorder=2)
        mca.fillcontinents(color=darkgray, zorder=0)
        mca.drawmeridians([-119, -115], labels=[0, 0, 1, 0],
                          color='lightgray', linewidth=0.5, zorder=0,
                          fontsize=lbl_fontsize)
        mca.drawparallels([32, 34, 36], labels=[0, 1, 0, 0],
                          color='lightgray', linewidth=0.5, zorder=0,
                          fontsize=lbl_fontsize)
        xtxt, ytxt = mca(-120, 37.3)
        ax_cal.text(xtxt, ytxt, 'southern California', fontsize=txt_fontsize,
                    bbox=dict(facecolor='white', alpha=1.0))
        if traveltime:
            resultsfn = os.path.join(datadir, 'ptt_ca_6stations.npz')
            plot_ptt(resultsfn, mca, cmap, dlevel, vmin=vmin, vmax=vmax)
        if alerttime or blindzone or optimalbz:
            if bzreduction or optimalbz:
                blindzone = True
            resultsfn = os.path.join(datadir, 'event_list_ca.csv')
            optbz = os.path.join(datadir, 'optimal_blindzone_ca.npz')
            plot_at(resultsfn, mca, cmap, vmin=vmin, vmax=vmax,
                    blindzone=blindzone, optbz=optbz, bzreduction=bzreduction,
                    optimalbz=optimalbz, ssize=50)

    # Create an inset for Iceland
    if True:
        ax_ice = fig.add_axes([0.05, 0.63, .2, .25])
        mi = Basemap(width=550000, height=580000,
                    resolution='l', projection='aea', \
                    lat_1=62.5, lat_2=68.5, lon_0=-19, lat_0=65,
                    ax=ax_ice)
        mi.drawcoastlines(zorder=2)
        mi.fillcontinents(color=darkgray, zorder=0)
        mi.drawmeridians(np.arange(-26, -12, 4), labels=[0, 0, 1, 0],
                        color='lightgray', linewidth=0.5, zorder=0,
                        fontsize=lbl_fontsize)
        mi.drawparallels(np.arange(60, 70, 2), labels=[0, 1, 0, 0],
                        color='lightgray', linewidth=0.5, zorder=0,
                        fontsize=lbl_fontsize)
        xtxt, ytxt = mi(-23.5, 67)
        ax_ice.text(xtxt, ytxt, 'Iceland', fontsize=txt_fontsize)
        # plot Iceland data
        if traveltime:
            resultsfn = os.path.join(datadir, 'ptt_imo_6stations.npz')
            plot_ptt(resultsfn, mi, cmap, dlevel, vmin=vmin, vmax=vmax)
        if alerttime or blindzone or optimalbz:
            if bzreduction or optimalbz:
                blindzone = True
            resultsfn = os.path.join(datadir, 'event_list_iceland_all.csv')
            optbz = os.path.join(datadir, 'optimal_blindzone_is.npz')
            plot_at(resultsfn, mi, cmap, vmin=vmin, vmax=vmax,
                    blindzone=blindzone, optbz=optbz, bzreduction=bzreduction,
                    optimalbz=optimalbz)


    # Create an inset for New Zealand
    if True:
        ax_nz = fig.add_axes([-0.04, 0.12, .5, .45])
        mnz = Basemap(width=1300000, height=1700000,
                    resolution='l', projection='aea', \
                    lat_1=-50., lat_2=-32, lon_0=172, lat_0=-41,
                    ax=ax_nz)
        mnz.drawcoastlines(zorder=2)
        mnz.fillcontinents(color=darkgray, zorder=0)
        mnz.drawmeridians(np.arange(164, 182, 4), labels=[0, 0, 0, 1],
                        color='lightgray', linewidth=0.5, zorder=0,
                        fontsize=lbl_fontsize)
        mnz.drawparallels(np.arange(-51, -31, 2), labels=[1, 0, 0, 0],
                        color='lightgray', linewidth=0.5, zorder=0,
                        fontsize=lbl_fontsize)
        xtxt, ytxt = mnz(165.5, -35)
        ax_nz.text(xtxt, ytxt, 'New Zealand', fontsize=txt_fontsize)
        # plot NZ data
        if traveltime:
            resultsfn = os.path.join(datadir, 'p_wave_tt_nz.npz')
            plot_ptt(resultsfn, mnz, cmap, dlevel, vmin=vmin, vmax=vmax)
        if alerttime or blindzone or optimalbz:
            if bzreduction or optimalbz:
                blindzone = True
            resultsfn = os.path.join(datadir, 'event_list_nz.csv')
            optbz = os.path.join(datadir, 'optimal_blindzone_nz.npz')
            plot_at(resultsfn, mnz, cmap, vmin=vmin, vmax=vmax,
                    blindzone=blindzone, optbz=optbz, bzreduction=bzreduction,
                    optimalbz=optimalbz)

    cb_fontsize = 14
    lbl_fontsize = 14
    lbl_pad = 10
    cax = fig.add_axes([0.87, 0.1, 0.01, 0.8])
    if traveltime:
        cb = ColorbarBase(cax, cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax),
                          orientation='vertical', extend=extend)
        cb.set_label(cb_label[0], fontsize=cb_fontsize, labelpad=lbl_pad)
        cb.ax.tick_params(labelsize=lbl_fontsize)
        cax2 = fig.add_axes([.95, 0.1, 0.01, 0.8])
        cb2 = ColorbarBase(cax2, cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax),
                           orientation='vertical', extend=extend)
        cb2.set_label(cb_label[1], fontsize=cb_fontsize, labelpad=lbl_pad)
        cb2.set_ticks(np.arange(1.5, 16.5, 1.5))
        cb2.set_ticklabels(station_distance(np.arange(1.5, 16.5, 1.5), 8.0, 6.5, 6))
        cb2.ax.tick_params(labelsize=lbl_fontsize)
    if (blindzone or optimalbz) and not bzreduction:
        cb = ColorbarBase(cax, cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax),
                          orientation='vertical', extend=extend)
        ticks = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200])
        cb.set_ticks(ticks)
        cb.set_ticklabels(ticks)
        cb.set_label(cb_label[0], fontsize=cb_fontsize, labelpad=lbl_pad)
        cb.ax.tick_params(labelsize=lbl_fontsize)
        cax2 = fig.add_axes([0.95, 0.1, 0.01, 0.8])
        cax3 = fig.add_axes([1.03, 0.1, 0.01, 0.8])
        cb2 = ColorbarBase(cax2, cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax),
                           orientation='vertical', extend=extend)
        cb2.set_label(cb_label[1], fontsize=cb_fontsize, labelpad=lbl_pad)
        cb2.set_ticks(ticks)
        cb2.set_ticklabels([int(x / 6.5 + 0.5) for x in ticks])
        cb2.ax.tick_params(labelsize=lbl_fontsize)
        cb3 = ColorbarBase(cax3, cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax),
                           orientation='vertical', extend=extend)
        tklbl = []
        for _bz in ticks:
            idx = np.argmin(np.abs(dz - _bz))
            tklbl.append(mags[idx])
        cb3.set_label(cb_label[2], fontsize=cb_fontsize, labelpad=lbl_pad)
        cb3.set_ticks(ticks)
        cb3.set_ticklabels(tklbl)
        cb3.ax.tick_params(labelsize=lbl_fontsize)
    if bzreduction:
        cb = ColorbarBase(cax, cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax),
                          orientation='vertical', extend=extend)
        ticks = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200])
        cb.set_ticks(ticks)
        cb.set_ticklabels(ticks)
        cb.set_label(cb_label[0], fontsize=cb_fontsize, labelpad=lbl_pad)
        cb.ax.tick_params(labelsize=lbl_fontsize)
        cax2 = fig.add_axes([.95, 0.1, 0.01, 0.8])
        cb2 = ColorbarBase(cax2, cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax),
                           orientation='vertical', extend=extend)
        cb2.set_label(cb_label[1], fontsize=cb_fontsize, labelpad=lbl_pad)
        cb2.set_ticks(ticks)
        cb2.set_ticklabels([int(x / 6.5 + 0.5) for x in ticks])
        cb2.ax.tick_params(labelsize=lbl_fontsize)

    fig.savefig(fout, dpi=300, bbox_inches='tight')
    plt.show()

def station_distance(traveltime, depth, vp, nst):
    dist = traveltime * vp
    return [int(sdist + 0.5) for sdist in np.sqrt(dist * dist - depth * depth)]

def alert_time(blindzone, depth, vp):
    return [int(at + 0.5) for at in np.sqrt(blindzone * blindzone + depth * depth) / vp]

def damage_zone():
    dz = []
    mags = np.arange(5.0, 9.0, 0.1)
    c0 = 2.085
    c1 = 1.428
    c2 = -1.402
    c4 = 0.078
    m1 = -0.209
    m2 = 2.042
    Rhyp = np.arange(1, 300, 0.5)
    dist_correct = c4 * np.log(Rhyp / 50.)
    S = 0.0
    for mag in mags:
        Rm = m1 + m2 * np.exp(mag - 5)
        I = c0 + c1 * mag + c2 * np.log(np.sqrt(Rhyp ** 2 + Rm ** 2)) + S
        I = np.where(Rhyp > 50., I + dist_correct, I)
        idx = np.where(I >= 7)
        dz.append(Rhyp[idx][-1])
    return (mags, dz)
if __name__ == '__main__':
    main(blindzone=True)
