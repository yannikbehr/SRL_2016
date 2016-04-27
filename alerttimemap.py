#!/usr/bin/env python
"""
Make a plot of P-wave travel time to 6 stations in the network.

Created on Jun 24, 2013

@author: behry
"""
import sys
import numpy as np
import matplotlib
try:
    matplotlib.use('WXAgg')
    import wx
except:
    print "WX package for Python not installed"
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from matplotlib.pyplot import cm
from mpl_toolkits.basemap import Basemap
from scipy.stats import scoreatpercentile
from copy import copy
import warnings
import pyproj
import ipdb


class AlertTimeMap:

    def __init__(self, resultsfn=None, mapbox=(45, 48.5, 5, 12, 47),
                 cmapname='RdBu_r', ax=None, ** kargs):
        llcrnrlat, urcrnrlat, llcrnrlon, urcrnrlon, lat_ts = mapbox
        self.m = Basemap(projection='merc', llcrnrlat=llcrnrlat,
                         urcrnrlat=urcrnrlat, llcrnrlon=llcrnrlon,
                         urcrnrlon=urcrnrlon, lat_ts=lat_ts,
                         resolution='h', ax=ax)
        self.cmapname = cmapname
        self.extend = 'max'
        self.load_data(resultsfn, **kargs)
        self.vs = 3.5

    def load_data(self, resultsfn, **kargs):
        if resultsfn is not None:
            a = np.load(resultsfn)
            self.tstarget = a['tstarget']
            self.ttP = a['ttP']
            self.dims = self.ttP.shape
            if len(self.dims) > 2:
                self.lat = a['lat'][:, :, 0]
                self.lon = a['lon'][:, :, 0]
                self.dep = a['dep'][0, 0, :]
                # Keep also original matrices to compute
                # earthquake to site distance more easily
                self.elat = a['lat']
                self.elon = a['lon']
                self.edep = a['dep']
            else:
                self.lat = a['lat']
                self.lon = a['lon']
                self.dep = a['dep']
        else:
            if 'ttP' in kargs:
                self.ttP = kargs['ttP']
                self.lat = kargs['lat']
                self.lon = kargs['lon']
                self.dep = kargs['dep']
                self.tstarget = kargs['tstarget']
        self.x, self.y = self.m(self.lon, self.lat)

    def smooth(self, ttPmed):
        xmin, ymin = self.m(self.lon.min(), self.lat.min())
        xmax, ymax = self.m(self.lon.max(), self.lat.max())
        nx = int((xmax - xmin) / 5000.) + 1
        ny = int((ymax - ymin) / 5000.) + 1
        lon_lin = np.linspace(self.lon.min(), self.lon.max(),
                              self.lon.shape[1])
        lat_lin = np.linspace(self.lat.min(), self.lat.max(),
                              self.lat.shape[0])
        dat, x, y = self.m.transform_scalar(ttPmed.T, lon_lin, lat_lin,
                                            nx, ny, returnxy=True, masked=True)
        return x, y, dat

    def plot1D(self, scale=True, meridians=np.arange(5, 12, 2),
               parallels=np.arange(44, 49, 2), vmin=0., vmax=15,
               blindzone=False):
        """
        Show alert times for selected events.
        """
        ttPmed = np.median(self.ttP, axis=0)
        if blindzone:
            bz2 = ttPmed * ttPmed * self.vs * self.vs - self.dep * self.dep
            bz2 = np.where(bz2 > 0.0, bz2, 0.0)
            bz = np.sqrt(bz2)
            ttPmed = bz

        cmap = cm.ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax),
                                 cmap=self.cmapname)
        for lt, ln, delay in zip(self.lat, self.lon, ttPmed):
            cl = cmap.to_rgba(delay)
            x, y = self.m(ln, lt)
            self.m.plot(x, y, ms=8, c=cl, marker='o', picker=5.)
        self.m.drawmeridians(meridians, labels=[0, 0, 0, 1], color='lightgray',
                             linewidth=0.5, zorder=0)
        self.m.drawparallels(parallels, labels=[1, 0, 0, 0], color='lightgray',
                             linewidth=0.5, zorder=0)
        self.m.drawcoastlines(zorder=2)
        self.m.drawcountries(linewidth=1.0, zorder=2)
        self.m.drawstates(zorder=2)

    def plot2D(self, fig, ax, nnst=6, vs=3.5, procdelay=False,
               scale=False, boxin=(45.4, 48.3, 5.6, 11.1),
               geofilter=None, boxout=(45, 48.5, 5, 12, 47), target=None,
               meridians=np.arange(5, 12, 2), parallels=np.arange(44, 49, 2),
               blindzone=False, clevels=None, clevels_only=False,
               clevel_colours=['lightgray', 'gray'], cbshrink=0.8,
               smooth=True):

            self.m.ax = ax
            if target is not None:
                # If target is given compute S-wave traveltime from every
                # potential epicenter to the target site
                tln, tlt = target
                lat = self.elat.ravel()
                lon = self.elon.ravel()
                dep = self.edep.ravel()
                tstarget = np.zeros(self.dims)
                tstmp = np.zeros(lon.size)
                g = pyproj.Geod(ellps='WGS84')
                idx = 0
                for _lat, _lon, _dep in zip(lat, lon, dep):
                    azt, bazt, distt = g.inv(_lon, _lat, tln, tlt)
                    distt /= 1000.
                    tstmp[idx] = np.sqrt(distt * distt + _dep * _dep) / vs
                    idx += 1
                tstarget[:] = tstmp.reshape(self.elon.shape)

            latmin, latmax, lonmin, lonmax = boxin
            cmap = cm.get_cmap(self.cmapname)
            cmap.set_over('grey')
            if procdelay:
                self.vmin = 4.
                self.vmax = 25.
                self.cb_label = 'Time since origin time [s]'
                unit = 's'
                dlevel = 0.5
            if target is not None:
                self.extend = 'both'
                self.vmin = -10.
                self.vmax = 60.
                cmap = cm.get_cmap(self.cmapname)
                cmap.set_over('grey')
                self.cb_label = 'Lead time [s]'
                unit = 's'
                dlevel = 0.5
            if not procdelay and target is None:
                self.vmin = 0.
                self.vmax = 15.
                self.cb_label = 'P-wave travel time to %d stations [s]' % nnst
                unit = 's'
                dlevel = 0.5
            if blindzone:
                self.vmin = 22.
                self.vmax = 55.
                self.cb_label = 'Blind zone [km]'
                unit = 'km'
                dlevel = 0.5

            # Mask points outside of polygon
            if geofilter is not None:
                if self.dims > 2:
                    rowidx = 0
                    colidx = 0
                    idx = 0
                    ydim = self.lat.shape[1]
                    for _lat, _lon in zip(self.lat.ravel(), self.lon.ravel()):
                        rowidx = idx / ydim
                        colidx = idx - rowidx * ydim
                        idx += 1
                        if not geofilter.point_in_polygon([_lon], [_lat])[0]:
                            self.ttP[:, rowidx, colidx, :] = np.nan
                else:
                    idx = 0
                    for _lat, _lon in zip(self.lat, self.lon):
                        if not geofilter.point_in_polygon([_lon], [_lat])[0]:
                            self.ttP[:, idx] = np.nan
                        idx += 1

            ttPmed = np.median(self.ttP, axis=0)
            ttPmed = np.median(ttPmed, axis=-1)
            print "The minimum alert time is: ", np.ma.masked_invalid(ttPmed).min()
            if target is not None:
                print tstarget.min(), tstarget.max()
                ttPmed = np.median((tstarget - self.ttP), axis=0)
                ttPmed = np.median(ttPmed, axis=-1)

            if blindzone:
                depth = np.median(self.dep)
                bz2 = ttPmed * ttPmed * self.vs * self.vs - depth * depth
                bz2 = np.where(bz2 > 0.0, bz2, 0.0)
                bz = np.sqrt(bz2)
                ttPmed = bz
                print "The minimum blindzone is: ", ttPmed.min()
            if smooth:
                x, y, ttPmed = self.smooth(ttPmed)
            else:
                x, y = self.x, self.y
            if not clevels_only:
                cf = self.m.contourf(x, y, ttPmed, cmap=cmap,
                                     levels=np.arange(self.vmin,
                                                      self.vmax + dlevel, dlevel),
                                     norm=Normalize(vmin=self.vmin,
                                                    vmax=self.vmax),
                                     extend=self.extend)

            if target is not None:
                xt, yt = self.m(tln, tlt)
                self.m.plot(xt, yt, marker='o', ms=14, color='black')

            # Add contour lines
            if clevels is not None:
                for _lev, _col in zip(clevels, clevel_colours):
                    cs = self.m.contour(x, y, ttPmed,
                                        colors=_col, levels=[_lev],
                                        linestyles='solid', linewidths=3)
                    with warnings.catch_warnings(record=True):
                        plt.clabel(cs, fmt="%d " + unit, fontsize=12,
                                   colors='black')

            if scale:
                cb = fig.colorbar(cf, ax=ax, extend=self.extend,
                                  orientation='vertical',
                                  spacing='uniform', shrink=cbshrink)
                cb.set_label(self.cb_label)
            self.m.drawmeridians(meridians, labels=[0, 0, 0, 1], color='lightgray',
                            linewidth=0.5, zorder=0)
            self.m.drawparallels(parallels, labels=[1, 0, 0, 0], color='lightgray',
                            linewidth=0.5, zorder=0)
            self.m.drawcoastlines(zorder=2)
            self.m.drawcountries(linewidth=1.0, zorder=2)
            self.m.drawstates(zorder=2)

            return self.m

    def plot_stations(self, fig, ax, networkinfo, interactive=False,
                      networklabels=False):
        # Plot station locations
        dataX = []
        dataY = []
        values = []
        networks = networkinfo.get_networks()
        for net in networks.keys():
            stlat = networks[net]['lat']
            stlon = networks[net]['lon']
            names = networks[net]['nm']
            color = networks[net]['color']
            label = networks[net]['label']
            nwcode = np.array(networks[net]['nw'])
            for _ln, _lt, nm, nt in zip(stlon, stlat, names, nwcode):
                x, y = self.m(_ln, _lt)
                values.append('%s' % (nm))
                dataX.append(x)
                dataY.append(y)
                self.m.plot(x, y, color=color, marker='^', mec='black', ms=6,
                            picker=5.)
            if networklabels:
                self.m.plot(x, y, color=color, marker='^', mec='black', ms=5,
                       label=label, ls='None')
        if networklabels:
                ax.legend(numpoints=1, borderpad=0.2, fontsize='small',
                          markerscale=1.5)

        # add an interactive picker
        if interactive:
            try:
                # add a pop-up window showing the station and its value
                tooltip = wx.ToolTip(tip='')
                tooltip.Enable(False)
                tooltip.SetDelay(0)
                fig.canvas.SetToolTip(tooltip)

                def onMotion(event):
                    line2d = event.artist
                    x = line2d.get_xdata()[0]
                    y = line2d.get_ydata()[0]
                    found = False
                    for i in xrange(len(dataX)):
                        radius = 5
                        if abs(x - dataX[i]) < radius and abs(y - dataY[i]) < radius:
                            tip = '%s' % values[i]
                            tooltip.SetTip(tip)
                            tooltip.Enable(True)
                            found = True
                            break
                    if not found:
                        tooltip.Enable(False)
                fig.canvas.mpl_connect('pick_event', onMotion)
            except Exception, e:
                print "Cannot add wx.ToolTip: ", e

    def error_plot(self, ax_lb, ax_ub, cax, cborientation='vertical',
                   smooth=True):
        # plot the error map
        ttP_lb = np.zeros((self.dims[1::]))
        ttP_ub = ttP_lb.copy()
        for _i1 in xrange(self.dims[1]):
            for _i2 in xrange(self.dims[2]):
                for _i3 in xrange(self.dims[3]):
                    ttP_lb[_i1, _i2, _i3] = scoreatpercentile(self.ttP[:, _i1, _i2, _i3], 16)
                    ttP_ub[_i1, _i2, _i3] = scoreatpercentile(self.ttP[:, _i1, _i2, _i3], 84)

        mlb = copy(self.m)
        mlb.ax = ax_lb
        mub = copy(self.m)
        mub.ax = ax_ub

        cmap = cm.get_cmap(self.cmapname)
        cmap.set_over('grey')
        if smooth:
            x, y, ttP_lb = self.smooth(ttP_lb[:, :, 0])
        else:
            x, y, ttP_lb = self.x, self.y, ttP_lb[:, :, 0]

        mlb.contourf(x, y, ttP_lb, cmap=cmap,
                     levels=np.arange(self.vmin, self.vmax + 0.5, 0.5),
                     norm=Normalize(vmin=self.vmin, vmax=self.vmax),
                     extend=self.extend)

        if smooth:
            x, y, ttP_ub = self.smooth(ttP_ub[:, :, 0])
        else:
            x, y, ttP_ub = self.x, self.y, ttP_ub[:, :, 0]
        mub.contourf(x, y, ttP_ub, cmap=cmap,
                     levels=np.arange(self.vmin, self.vmax + 0.5, 0.5),
                     norm=Normalize(vmin=self.vmin, vmax=self.vmax),
                     extend=self.extend)
        mlb.drawcoastlines(zorder=2)
        mlb.drawcountries(linewidth=1.0, zorder=2)
        mub.drawcoastlines(zorder=2)
        mub.drawcountries(linewidth=1.0, zorder=2)
        cb = ColorbarBase(cax, cmap=cmap, norm=Normalize(vmin=self.vmin,
                                                         vmax=self.vmax),
                         orientation=cborientation, extend=self.extend)
        cb.set_label(self.cb_label)
        return mlb, mub


if __name__ == '__main__':
    pass
