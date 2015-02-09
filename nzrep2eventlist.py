#!/usr/bin/env python

matched = []
lines = open('data/vs_reports_nz_2014.txt').readlines()
for _l in lines:
    if _l.startswith('Match:'):
        matched.append(_l)

fh = open('data/event_list_nz.csv','w')
print >>fh, '#origin time,true mag,true lat,true lon,true depth,alert time,1st estimate: lat, lon, dep, mag,last estimate:lat, lon, dep, mag'
for _l in matched:
    tmp = _l.replace('Match: ','')
    idx = tmp.rfind(';')
    a = tmp.replace(';',',')
    print >>fh, a[0:idx]
fh.close()
                    
