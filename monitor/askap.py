#!/usr/bin/env python3
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

def get_askap_data():
    url = "https://prod-api.vlbi.atnf.tools/get_status/askap"
    response = urlopen(url)
    data = response.read().decode("utf-8")
    return data, json.loads(data)

def dump_raw(data):
    for key1 in data.keys():
        dataval = data[key1]
        if type(dataval) is dict:
            for key2 in data[key1].keys():
                print("  %s : " %(key2), dataval[key2])
        else:
            print("%s : " %(key1), dataval)

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
    missing = []
    if nant != 36:
        for ant in range(1,37):
            if ant in ants:
                continue
            missing_ants.append("ak%02d" %(ant))
            missing.append(ant-1)
    return missing_ants, missing

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

def flatten(data):
    fdata = {}
    for key1 in data.keys():
        dataval = data[key1]
        if type(dataval) is dict:
            for key2 in data[key1].keys():
                fdata[key2] = dataval[key2]
        else:
            fdata[key1] = dataval
    missingAnts,missing = find_missing(fdata)
    fdata["missingAnts"] = missingAnts
    for i in range(len(fdata["azimuth"])):
        if fdata["azimuth"][i] is None:
            data["azimuth"][i] = np.nan
    # For antennas that aren't listed assume they are bad and so mark the antenna-dependent data
    # as bad.
    for i in missing:
        fdata["azimuth"][i] = np.nan
        fdata["elevation"][i] = np.nan
        fdata["polarisation"][i] = np.nan

    fdata["medianAzimuth"] = np.nanmedian(fdata["azimuth"])
    obs_status, ant_states = get_states(fdata)
    fdata["obsStatus"] = obs_status
    fdata["antStates"] = ant_states
    fdata["coord"] = median_ra_dec(fdata)
    fdata["sbid"] = int(fdata["id"])
    
    # Convert these key value to float
    floatdefs = {}
    floatdefs["duration"] = 0.0
    floatdefs["progress"] = 0.0
    floatdefs["pol_axis_angle"] = 0.0
    floatdefs["temperature"] = 0.0
    floatdefs["humidity"] = 0.0
    floatdefs["pressure"] = 0.0
    floatdefs["windSpeed"] = 0.0
    floatdefs["rotation"] = 0.0
    floatdefs["pitch"] = 0.0
    for key in floatdefs.keys():
        if key in fdata:
            if fdata[key] == 'null':
                data[key] = floatdefs[key]
            else:
                fdata[key] = float(fdata[key])
        else:
            fdata[key] = floatdefs[key]
    keyValues = ["infoTime", "obsStatus", "field_name", "coord", "startTime", "missingAnts", "name"]
    for key in keyValues:
        if key in fdata:
            continue
        else:
            fdata[key] = "Undefined"
    return fdata

def dump_summary(fdata):
    print("Current time (UTC): %s" %(fdata["infoTime"]))
    print("Observing status: %s" %(fdata["obsStatus"]))
    # Observing parameters only make sense when not in Stowed state
    if (fdata["obsStatus"] in ["Stowed"]) == False:
        print("Currently observing: %s" %(fdata["field_name"]))
        print("Scheduling block: %d" %(fdata["sbid"]))
        print("Coordinate: %s" %(fdata["coord"].to_string(style='hmsdms')))
        print("Start time: %s" %(fdata["startTime"]))
        print("Total duration %.1f s" %(fdata["duration"]))
        print("Progress: %.1f %%" %(fdata["progress"]))

        if len(fdata["missingAnts"]) == 0:
            print("All antennas available")
        else:
            print("Missing antennas: %s" %(", ".join(fdata["missingAnts"])))

        print("Footprint:")
        print("  Footprint    : %s" %(fdata["name"]))
        print("  Pitch        : %.2f" %(fdata["pitch"]))
        print("  Rotation     : %.1f" %(fdata["rotation"]+fdata["pol_axis_angle"]))

    print("Weather conditions")
    print("  Temperature: %.1f C" %(fdata["temperature"]))
    print("  Humidty    : %.1f %%" %(fdata["humidity"]))
    print("  Pressure   : %.1f millibars" %(fdata["pressure"]))
    print("  WindSpeed. : %.1f kph" %(fdata["windSpeed"]))
    print("  StormFlag. : %s" %(fdata["stormFlag"]))

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
    fdata = flatten(data)    
    dump_summary(fdata)
    sleep(60)