#!/usr/bin/env python
"""
Plot the cumulative distribution functions of delay measurements.
Created on Aug 31, 2015

@author: behry
"""
import matplotlib
# matplotlib.use('Agg')
import json
import numpy as np
from pandas import Series, DataFrame
import pandas as pd
import matplotlib.pyplot as plt
import string
pd.options.display.mpl_style = 'default'


def plot_cdf(ax, data, cl, lbl, binmin=0, binmax=30, dbin=1.0):
    data = np.array(data)
    dmin = data.min()
    dmax = data.max()
    if dmin >= binmin and dmax <= binmax:
        bins = np.arange(binmin, binmax + dbin, dbin)
    elif dmin < binmin and dmax <= binmax:
        bins = np.r_[dmin, np.arange(binmin, binmax + dbin, dbin)]
    elif dmin >= binmin and dmax > binmax:
        bins = np.r_[np.arange(binmin, binmax + dbin, dbin), dmax]
    else:
        bins = np.r_[dmin, np.arange(binmin, binmax + dbin, dbin), dmax]
    hist, edges = np.histogram(data, bins=bins, density=False)
    x = edges[0:-1] + np.diff(edges) / 2.
    y = np.cumsum(hist.astype(float)) / np.sum(hist)
    ax.plot(x, y, color=cl)
    # ax.hist(data, cumulative=True, bins=bins, normed=True,
    #        histtype='step', color=_cl)
    lines.append(plt.Line2D([], [], lw=2.0, linestyle='-', color=cl))
    labels.append('%s (%d)' % (lbl, len(data)))

pkd_fn = './data/pick_delays_ch.txt'
dl = {'centaur':('./data/ch_dataloggers/centaur.txt', '#A60628'),
      'fixnet':('./data/ch_dataloggers/fixnet.txt', '#A60628'),
      'mobile-lte':('./data/ch_dataloggers/mobile-lte.txt', '#467821'),
      'taurus':('./data/ch_dataloggers/taurus.txt', '#467821')}

df = None
all = {}
fig = plt.figure(figsize=(16, 8))

fh = open(pkd_fn)
pkdel = json.load(fh)
fh.close()

# Different dataloggers
ax1 = fig.add_subplot(1, 3, 1)
lines = []
labels = []
allstats = []
stats_no_info = []
for _k in ['centaur', 'taurus']:
    _fn, _cl = dl[_k]
    data = []
    stats = map(string.strip, open(_fn).readlines())
    for _s in stats:
        _sk = 'CH.%s' % _s
        if _sk not in pkdel:
            if _sk not in stats_no_info:
                stats_no_info.append(_sk)
        else:
            data.extend(pkdel[_sk])
            if _sk not in allstats:
                allstats.append(_sk)
            else:
                print "Station %s in more than one list" % _sk
    plot_cdf(ax1, data, _cl, _k)
    all[_k] = Series(data)

data = []
for _s in pkdel.keys():
    if _s not in allstats:
        net, stat = _s.split('.')
        if net in ['8D', 'CH']:
            data.extend(pkdel[_s])

plot_cdf(ax1, data, '#7A68A6', 'other')
all['other'] = Series(data)
ax1.set_xlim(0, 30)
ax1.set_ylim(0, 1.0)
ax1.set_xlabel('Pick delay [s]')
ax1.set_ylabel('Percentile')
ax1.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
ax1.set_yticklabels([r'$0^{th}$', r'$25^{th}$', r'$50^{th}$',
                    r'$75^{th}$', r'$100^{th}$'])
ax1.legend(lines, labels, loc='lower right', ncol=1)

# Different connections
ax2 = fig.add_subplot(1, 3, 2)
lines = []
labels = []
allstats = []
for _k in ['fixnet', 'mobile-lte']:
    _fn, _cl = dl[_k]
    data = []
    stats = map(string.strip, open(_fn).readlines())
    for _s in stats:
        _sk = 'CH.%s' % _s
        if _sk not in pkdel:
            if _sk not in stats_no_info:
                stats_no_info.append(_sk)
        else:
            data.extend(pkdel[_sk])
            if _sk not in allstats:
                allstats.append(_sk)
            else:
                print "Station %s in more than one list" % _sk
    plot_cdf(ax2, data, _cl, _k)
    all[_k] = Series(data)

data = []
for _s in pkdel.keys():
    if _s not in allstats:
        net, stat = _s.split('.')
        if net in ['8D', 'CH']:
            data.extend(pkdel[_s])

plot_cdf(ax2, data, '#7A68A6', 'other')
all['other'] = Series(data)

ax2.set_xlim(0, 30)
ax2.set_ylim(0, 1.0)
ax2.set_xlabel('Pick delay [s]')
ax2.set_ylabel('Percentile')
ax2.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
ax2.set_yticklabels([r'$0^{th}$', r'$25^{th}$', r'$50^{th}$',
                    r'$75^{th}$', r'$100^{th}$'])
ax2.legend(lines, labels, loc='lower right', ncol=1)

# CH vs foreign stations
ax3 = fig.add_subplot(1, 3, 3)
lines = []
labels = []
allstats = []
ch_data = []
f_data = []
for _k, _d in pkdel.iteritems():
    net, stat = _k.split('.')
    if net in ['CH', '8D']:
        ch_data.extend(_d)
    else:
        f_data.extend(_d)
plot_cdf(ax3, ch_data, '#A60628', 'Swiss stations')
all['swiss'] = Series(ch_data)
plot_cdf(ax3, f_data, '#467821', 'Foreign stations')
all['foreign'] = Series(f_data)

for _s in stats_no_info:
    print "No pick delay information for %s" % _s

ax3.set_xlim(0, 30)
ax3.set_ylim(0, 1.0)
ax3.set_xlabel('Pick delay [s]')
ax3.set_ylabel('Percentile')
ax3.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
ax3.set_yticklabels([r'$0^{th}$', r'$25^{th}$', r'$50^{th}$',
                    r'$75^{th}$', r'$100^{th}$'])
ax3.legend(lines, labels, loc='lower right', ncol=1)
fig.savefig('./plots/pick_delay_cdf_ch.png', bbox_inches='tight')
df = DataFrame(all)
plt.show()
