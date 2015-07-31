#!/usr/bin/env python

import sys
sys.path.append('../reports')
from pyrep import Reporter

fin = 'data/vssc3_reports_cit.txt'
r = Reporter()
evts = r.read_report(fin)

fh = open('data/event_list_ca.csv', 'w')
print >> fh, '#origin time,true mag,true lat,true lon,true depth,alert time,1st estimate: lat, lon, dep, mag,last estimate:lat, lon, dep, mag'
for _e in evts['match']:
    id, true_ot, true_lon, true_lat, true_mag, dep, tdiff, mag1, mag2, lat1, \
    lon1, lat2, lon2, lh1, lh2, dep1, dep2 = _e
    print >> fh, ','.join(map(str, (true_ot, true_mag, true_lat, true_lon, dep,
                           tdiff, lat1, lon1, dep1, mag1, lat2, lon2,
                           dep2, mag2)))
fh.close()

