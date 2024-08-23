#!/usr/bin/env python
import numpy as np
from astroquery.utils.tap.core import Tap
import datetime
from astropy.time import Time
from astropy.coordinates import EarthLocation, Angle, SkyCoord
from astropy import units as u

proj_code = 'AS201'

casda_tap = Tap("https://casda.csiro.au/casda_vo_tools/tap")

job = casda_tap.launch_job_async("SELECT sbid,obs_start,obs_end,event_date,event_type,obs_program,project_code FROM casda.observation_event WHERE project_code='%s'" %(proj_code))
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

sbids = list(sbid_status.keys())
sbids.sort()

released = []
unreleased = []
for sbid in sbids:
    if sbid_status[sbid] in ["UPDATED", "RELEASED"]:
        released.append(sbid)
    elif sbid_status[sbid] in ["REJECTED", "DELETED"]:
        continue
    else:
        unreleased.append(sbid)

print("Summary for project code: %s" %(proj_code))
for sbid in np.unique(unreleased):
    if sbid_status[sbid] in ["REJECTED"]:
        print("%d rejected" %(sbid))

# Where are we observing from?
latitude = Angle("-26:41:46.0", unit=u.deg)
longitude = Angle("116:38:13.0", unit=u.deg)
observing_location = EarthLocation(lat=latitude, lon=longitude)

t_now = Time(datetime.datetime.utcnow(), scale='utc', location=observing_location)
for sbid in np.unique(unreleased):
    if sbid_status[sbid] in ["REJECTED", "DELETED"]:
        continue
    event_time = Time(sbid_time[sbid].replace("T", " ").replace("Z", ""))
    elapsed_time_days = (t_now - event_time).sec / 60.0 / 60.0 / 24.0
    print("%d not yet released in CASDA: status = %s %.1f days ago" %(sbid, sbid_status[sbid], elapsed_time_days))