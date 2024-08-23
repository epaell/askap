#!/usr/bin/env python3
import numpy as np
from astroquery.utils.tap.core import Tap
import datetime
from astropy.time import Time
from astropy.coordinates import EarthLocation, Angle, SkyCoord
from astropy import units as u
import sys

if len(sys.argv) != 2:
    sys.exit("Usage: %s sbids.txt" %(sys.argv[0]))

casda_tap = Tap(url="https://casda.csiro.au/casda_vo_tools/tap")

# Get the list of sbids to find
sbids = []
for line in open(sys.argv[1]):
    if len(line)<2:
        continue
    if line[0]=="#":
        continue
    sbids.append(int(line))

job = casda_tap.launch_job_async("SELECT sbid,obs_start,obs_end,event_date,event_type,obs_program,project_code FROM casda.observation_event")
results_table = job.get_results()
results_table.sort('event_date')

# Find the latest event for each SBID
sbid_status = {}
sbid_time = {}
for item in results_table:
    sbid = int(item['sbid'])
    event = item['event_type']
    sbid_status[sbid] = event
    sbid_time[sbid] = item["event_date"]

sbid_list = list(sbid_status.keys())
sbid_list.sort()

# Go through the list and out the results
for sbid in sbids:
    if sbid in sbid_list:
        print("%d %s" %(sbid, sbid_status[sbid]))
