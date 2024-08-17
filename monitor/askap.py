#!/usr/bin/env python
from urllib.request import urlopen
import json 
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
import numpy as np
import sys
import os
from time import sleep
from astropy.table import Table, vstack
from astropy.io import ascii

def read_known_footprints():
    footprints = {}
    footprints["VAST"] = ["square_6x6", 1.05, 45.0]
    footprints["EMU"] = ["closepack36", 0.90, 45.0]
    footprints["RACS"] = ["closepack36", 0.90, 45.0]
    
    fields = []
    for fname in ["vast_fields.csv", "emu_fields.csv", "racs_fields.csv"]:
        sst_data = Table.read(fname, format='csv')
#        sst_data = sst_data[np.where(sst_data["SBID"]>0)]
        if len(fields) == 0:
            fields = sst_data
        else:
            fields = vstack([fields, sst_data])
    return fields, footprints

def get_footprint(fields, footprints, field_name):
    field = fields[np.where(fields["FIELD_NAME"]==field_name)]
    if len(field) == 0:
        return 0.0, 0.0, "", 0.0, 0.0
    # FIELD_NAME,RA_HMS,DEC_DMS,RA_DEG,DEC_DEG,GAL_LONG,GAL_LAT,POL_AXIS
    field = field[0]
    sst = field_name.split("_")[0]
    footprint = footprints[sst][0]
    pitch = footprints[sst][1]
    rotation = footprints[sst][2]
    return field["RA_DEG"], field["DEC_DEG"], footprint, pitch, field["POL_AXIS"] + rotation
    
def get_askap_data():
    url = "https://prod-api.vlbi.atnf.tools/get_status/askap"
    response = urlopen(url)
    data = response.read().decode("utf-8")
    return data, json.loads(data)

def dump_raw(data):
    keys = ['infoTime', 'antennasList', 'template', 'configuration', 'azimuth', 'state', 'alias', 'antennaName', 'declinationICRF', 'rightAscensionICRF', 'progress', 'errors', 'antennasIn', 'elevation', 'id', 'duration', 'stateError', 'weather', 'status']
    for key in keys:
        print(key, " : ", data[key])
        print()

def median_ra_dec(data):
    ra_str = np.array(data['rightAscensionICRF'])
    dec_str = np.array(data['declinationICRF'])
    ra_str = ra_str[np.where(ra_str != None)]   # Remove dodgy values
    dec_str = dec_str[np.where(ra_str != None)]
    ra_str = ra_str[np.where(dec_str != None)]
    dec_str = dec_str[np.where(dec_str != None)]
    coords = SkyCoord(Angle(ra_str, unit=u.deg), Angle(dec_str, unit=u.deg), frame='fk5')
    ra_deg = np.median(coords.ra.deg)
    dec_deg = np.median(coords.dec.deg)
    coord = SkyCoord(Angle(ra_deg, unit=u.deg), Angle(dec_deg, unit=u.deg), frame='fk5')
    return coord

def find_missing(data):
    ants = data['antennasList']
    nant = data['antennasIn']
    missing_ants = []
    if nant != 36:
        for ant in range(1,37):
            if ant in ants:
                continue
            missing_ants.append("ak%02d" %(ant))
    return missing_ants

def get_states(data):
    state = data['state']
    ant_states = {}
    ant_states["Idle"] = state.count("Idle")
    ant_states["Stowed"] = state.count("Stowed")
    ant_states["Slewing"] = state.count("Slewing")
    ant_states["Tracking"] = state.count("Tracking")
    obs_status = "undefined"
    if ant_states["Tracking"] > 18:
        obs_status = "Tracking"
    elif ant_states["Stowed"] > 18:
        obs_status = "Stowed"
    elif ant_states["Idle"] > 18:
        obs_status = "Idle"
    elif ant_states["Slewing"] > 18:
        obs_status = "Slewing"
    return obs_status, ant_states

def dump_summary(fields, footprints, start_time, current_time, field_name, sbid, coord, duration, progress, missing_ants, obs_status, ant_states, status, weather, footprint):
    print("Current time: %s" %(current_time))
    print("Currently observing: %s" %(field_name))
    print("Scheduling block: %d" %(sbid))
    print("Coordinate: %s" %(coord.to_string(style='hmsdms')))
    print("Total duration %.1f s" %(duration))
    print("Progress: %.1f %%" %(progress))
    print("Observing status: %s" %(obs_status))
#    print("%d tracking; %d slewing; %d idle; %d stowed" %(ant_states["Tracking"], ant_states["Slewing"], ant_states["Idle"], ant_states["Stowed"]))
    if len(missing_ants) == 0:
        print("All antennas available")
    else:
        print("Missing antennas: %s" %(", ".join(missing_ants)))

    fra_deg, fdec_deg, ffootprint, fpitch, frotation = get_footprint(fields, footprints, field_name)
    print("ffootprint=",ffootprint)
    sc = SkyCoord(Angle(fra_deg, unit=u.deg), Angle(fdec_deg, unit=u.deg), frame='fk5')
    
    print("Footprint:")
    print("  Footprint   : %s" %(footprint["name"]))
    print("  Pitch       : %.2f" %(footprint["pitch"]))
    print("  Rotation    : %.1f" %(footprint["rotation"]))
    if ffootprint == "":
        print("Unknown DB footprint")
    else:
        print("DB Footprint:")
        print("  Coordinate  : %s" %(sc.to_string(style='hmsdms')))
        print("  Footprint   : %s" %(ffootprint))
        print("  Pitch       : %.2f" %(fpitch))
        print("  Rotation    : %.1f" %(frotation))
#    print("status: %s" %(status))
#    print("stateError: %s" %(stateError))
    print("Weather conditions")
    print("  Temperature: %.1f C" %(weather["temperature"]))
    print("  Humidty    : %.1f %%" %(weather["humidity"]))
    print("  Pressure   : %.1f millibars" %(weather["pressure"]))
    print("  WindSpeed. : %.1f kph" %(weather["windSpeed"]))
    print("  StormFlag. : %s" %(weather["stormFlag"]))

fields, footprints = read_known_footprints()
summary = True
log = True
while True:
    response, data = get_askap_data()
    if log:
        log_file = open("log.txt", "at")
        log_file.write(response)
        log_file.write("\n")
        log_file.close()

    os.system('clear')
    if summary == False:
        dump_raw(data)
        sys.exit(0)

    # Extract what appears to be more useful information from the data provided
    sbid = int(data['id'])
    field_name = data['alias']
    current_time = data['infoTime']
    coord = median_ra_dec(data)

    progress = data['progress']
    duration = data['duration']
    stateError = data['stateError']
    status = data['status']
    weather = data['weather']
    footprint = data["footprint"]
    start_time = data["startTime:"]

    missing_ants = find_missing(data)
    obs_status, ant_states = get_states(data)
    
    dump_summary(fields, footprints, start_time, current_time, field_name, sbid, coord, duration, progress, missing_ants, obs_status, ant_states, status, weather, footprint)
    sleep(60)