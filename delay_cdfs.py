#!/usr/bin/env python
"""
Plot the cumulative distribution functions of delay measurements.
Created on Aug 31, 2015

@author: behry
"""
import matplotlib
matplotlib.use('WxAgg')
import json
import numpy as np
from pandas import Series, DataFrame
import pandas as pd
import matplotlib.pyplot as plt
pd.options.display.mpl_style = 'default'

def main(fout=None):
    fns = {'Turkey': ('./data/pick_delays_tr.txt', 'green', 'o'),
           'Switzerland': ('./data/pick_delays_ch.txt', 'blue', 's'),
           'Southern California': ('./data/pick_delays_ca.txt', 'red', '-'),
           'Greece': ('./data/pick_delays_gr.txt', 'cyan', '^'),
           'New Zealand': ('./data/pick_delays_nz.txt', 'black', '-.'),
           'Iceland': ('./data/pick_delays_is.txt', 'brown', '--'),
           'Romania': ('./data/pick_delays_ro.txt', 'orange', '*')}

    bw = True
    df = None
    all = {}
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    lines = []
    labels = []
    binmin = 0
    binmax = 30
    dbin = 1.0
    for _k in ['Southern California', 'Switzerland', 'New Zealand', 'Romania', 'Greece', 'Iceland', 'Turkey']:
        _fn, _cl, _m = fns[_k]
        fh = open(_fn)
        pkdel = json.load(fh)
        fh.close()
        data = []
        for _nst, _del in pkdel.iteritems():
            data.extend(_del)
        data = np.array(data)
        if _k == 'New Zealand':
            data -= 8.0
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

        if bw:
            hist, edges = np.histogram(data, bins=bins, density=False)
            x = edges[0:-1] + np.diff(edges) / 2.
            y = np.cumsum(hist.astype(float)) / np.sum(hist)
            if _m in ['-', '--', '-.']:
                ax.plot(x, y, color='k', linestyle=_m, lw=1.5)
                lines.append(plt.Line2D([], [], lw=1.5, linestyle=_m, color='k'))
            else:
                ax.plot(x, y, color='k', marker=_m, linestyle='-', lw=0.3)
                lines.append(plt.Line2D([], [], lw=0.3, marker=_m, color='k',
                                        linestyle='-'))
        else:
            ax.hist(data, cumulative=True, bins=bins, normed=True,
                    histtype='step', color=_cl)
            lines.append(plt.Line2D([], [], lw=2.0, linestyle='-', color=_cl))
        labels.append(_k)
        all[_k] = Series(data)
    ax.set_xlim(0, 30)
    ax.set_xlabel('Pick delay [s]')
    ax.set_ylabel('Percentile')
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels([r'$0^{th}$', r'$25^{th}$', r'$50^{th}$',
                        r'$75^{th}$', r'$100^{th}$'])
    ax.legend(lines, labels, loc='lower right', ncol=2)
    if fout is not None:
        fig.savefig(fout, bbox_inches='tight')
    df = DataFrame(all)
    plt.show()

if __name__ == '__main__':
    fout = './plots/pick_delay_cdf.pdf'
    main(fout)