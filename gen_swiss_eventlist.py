#!/usr/bin/env python

import os
import sys
sys.path.append('../reports')
import json
import glob
from obspy import UTCDateTime
from scvsmaglog_report_parser import parser as scvs_parser

datdir = '/home/sysop/VS_tests/data/iaspei2013_with_hg'
datdir2 = '/home/sysop/VS_tests/data/iaspei_no_data_latency/'
print '#origin time,true mag,true lat,true lon,true depth,alert time,1st estimate: lat, lon, dep, mag,last estimate:lat, lon, dep, mag'
for dirn in glob.glob(os.path.join(datdir, '*'))[:]:
    evid = os.path.basename(dirn)
    reps = glob.glob(os.path.join(dirn, '*_report.txt'))
    dtls = os.path.join(dirn, evid + '_details.txt')
    if not os.path.isfile(dtls) or len(reps) < 1:
        print "problems with event %s" % evid
        continue
    else:
        f = open(dtls, 'r')
        f.readline()
        line = f.readline()
        a = line.split()
        true_ot = UTCDateTime(a[4])
        lat,lon = map(float,a[0:2])
        true_dep = float(a[2])
        true_loc = (lat, lon)
        true_mag = float(a[3])
        f.close()
        if len(reps) > 1:
            # print reps
            min_tdiff = 1e39
            for _rr in reps:
                _t = scvs_parser(_rr, trueot=true_ot, trueloc=true_loc, format='new')
                mag, lt, ln, dep, ct, ot, td, ddiff, lkh, nstorig, nstmag = _t[0]
                _tdiff = abs(UTCDateTime(ct) - true_ot)
                if _tdiff < min_tdiff:
                    min_tdiff = _tdiff
                    cts = _t
        else:
            cts = scvs_parser(reps[0], trueot=true_ot, trueloc=true_loc, format='new')
        mags, lats, lons, deps, ct_s, ots, tds, ddiffs, lkh, nos, nms = cts[0]
        mage, late, lone, depe, ct_e, ote, tde, ddiffe, lkh, noe, nme = cts[-1]
        del1 = UTCDateTime(ct_s)
        # The playback runs with station delays have to be made
        # in realtime mode. Therefore origin times are not the original
        # ones anymore. To get an idea of the origin time error I compare
        # the origin time error for the same events from historic playbacks.
        reps2 = glob.glob(os.path.join(datdir2, evid, '*_report.txt'))
        if len(reps2) > 1:
            # print reps
            min_tdiff = 1e39
            for _rr in reps2:
                _t = scvs_parser(_rr, trueot=true_ot, trueloc=true_loc, format='new')
                mag, lt, ln, dep, ct, ot, td, ddiff, lkh, no, nm = _t[0]
                _tdiff = abs(UTCDateTime(ct) - true_ot)
                if _tdiff < min_tdiff:
                    min_tdiff = _tdiff
                    cts2 = _t
        else:
            cts2 = scvs_parser(reps2[0], trueot=true_ot, trueloc=true_loc, format='new')
        ot_del = UTCDateTime(cts2[-1, 5]) - UTCDateTime(true_ot)
        fmt = '{:s},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f}, {:.2f}, {:.2f}'
        print fmt.format(true_ot, true_mag, lat, lon, true_dep, del1 - UTCDateTime(ote) + ot_del, lats,lons,deps,mags,late,lone,depe,mage)

fh = open('/home/behry/workspace/eew/reports/data/swiss_nuclear_events.txt')
repdir = '/home/sysop_sc3rt/.seiscomp3/log/VS_reports/'
for _l in fh.readlines():
    if _l.startswith('#'):
        continue
    a = _l.split()
    true_lat = float(a[1])
    true_lon = float(a[2])
    true_loc = (true_lat, true_lon)
    true_dep = float(a[3])
    true_mag = float(a[4])
    true_ot = UTCDateTime(a[5])
    bn = a[-2]
    evid = a[-1].rstrip()
    if bn != 'None':
        repfn = os.path.join(repdir, bn)
        # only read new format reports that contain the number of
        # stations used for the origins
        try:
            ct = scvs_parser(repfn, trueot=true_ot, trueloc=true_loc, format='new')
        except:
            ct = scvs_parser(repfn, trueot=true_ot, trueloc=true_loc, format='old')
        mags, lts, lns, deps, cts, ots, tds, ddiffs, lkhs, nos, nms = ct[0]
        mage, lte, lne, depe, cte, ote, tde, ddiffe, lkhe, noe, nme = ct[-1]
        tdiff= UTCDateTime(cts) - true_ot
        print fmt.format(true_ot, true_mag,true_lat,true_lon, true_dep,tdiff,lts,lns,deps,mags,lte,lne,depe,mage)
