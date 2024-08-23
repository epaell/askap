#!/usr/bin/env python
import numpy as np
from astroquery.utils.tap.core import Tap
import datetime
from astropy.time import Time
from astropy.coordinates import EarthLocation, Angle, SkyCoord
from astropy import units as u
import sys

casda_tap = Tap("https://casda.csiro.au/casda_vo_tools/tap")

sbids = []
for line in open("sbids.txt"):
    if len(line)<2:
        continue
    if line[0]=="#":
        continue
    sbids.append(int(line))

for sbid in sbids:
    job = casda_tap.launch_job_async("SELECT sbid,obs_start,obs_end,event_date,event_type,obs_program,project_code FROM casda.observation_event WHERE sbid='%d'" %(sbid))
    results_table = job.get_results()
    results_table.sort('event_date')
    if len(results_table):
        print(results_table[0])
