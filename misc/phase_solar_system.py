#!/usr/bin/env python

from casacore.tables import *
import numpy as np
import astropy.constants as const
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, EarthLocation, Angle, SkyCoord
from astropy.coordinates import get_body
import glob
import os
import sys

def phase_rotate_to_target(observing_location, ms, ms_new, target, reverse, data_column):
    """ 
    Perform phase rotation
      observing_location - EarthLocation where observatory is located
      ms          - measurement set to work on
      ms_new      - results are written to this measurementset
      target      - target solar system object
      data_column - the data column to work on e.g. 'DATA', 'CORRECTED_DATA'
    """
    if ms == ms_new:
        sys.exit("Error: destination = source")
    # Remove old version of the destination if it exists
    if os.path.exists(ms_new):
        print("Destination already exists, removing previous destination")
        os.system("rm -fr %s" %(ms_new))
    # Copy over the the original data
    print("Copying data to %s" %(ms_new))
    os.system("cp -R %s %s" %(ms, ms_new))

    t = table(ms_new, readonly=True, ack=False)
    t1 = taql("select from $t where (ANTENNA1==0) && (ANTENNA2==0)")
    vis_time = t1.getcol('TIME')
    t1.close()
    t.close()

    t = Time(vis_time/60.0/60.0/24.0, format='mjd', scale='utc')
    print("Determining location of %s for all integrations ..." %(target))
    target_sc = get_body(target, t, observing_location)
    sc = SkyCoord(target_sc.to_string(style='hmsdms'), frame='fk5')

    mset = table(ms_new, readonly=False, ack=False)
    field_id = np.unique(mset.getcol("FIELD_ID"))
    if len(field_id) != 1:
        sys.exit("Can't work with multi-field measurement sets!")
    field_id = field_id[0]
    tp = table("%s/FIELD" %(ms_new), readonly=False, ack=False)
    ms_phase = tp.getcol('PHASE_DIR')
    ra0, dec0 = ms_phase[field_id][0]
    freqs = mset.SPECTRAL_WINDOW.getcell('CHAN_FREQ', 0)
    nchan = len(freqs)
    lambdas = const.c.value / freqs
    filtered = taql("select * from $mset where not FLAG_ROW and ANTENNA1 <> ANTENNA2")
    nvis = len(filtered)
    nant = len(np.unique([filtered.col("ANTENNA1"), filtered.col("ANTENNA2")]))
    nbl = int(nant*(nant-1)/2)
    nint = len(sc)
    print("Total antennas: %d" %(nant))
    print("Total baselines: %d" %(nbl))
    print("Total integrations: %d" %(nint))
    assert(nint*nbl==nvis)
    tindex = 0
    pfrac = -1
    for row in range(0, nvis, nbl):
        # Show simple progress
        frac = int(10.0*(tindex/nint+0.005))*10
        if pfrac != frac:
            pfrac = frac
            print("Processed %d%%" %(pfrac))
        flags = filtered.getcol('FLAG', startrow=row, nrow=nbl)
        uvw = filtered.getcol('UVW', startrow=row, nrow=nbl)
        data = np.complex128(filtered.getcol(data_column, startrow=row, nrow=nbl))
        data[flags] = np.nan
    
        # Calculate rotated uvw values
        if reverse == False:    # Phase to Sun
            new_uvw = rotateuvw(uvw, sc[tindex].ra.rad, sc[tindex].dec.rad, ra0, dec0)
        else:                   # Phase back to field
            new_uvw = rotateuvw(uvw, ra0, dec0, sc[tindex].ra.rad, sc[tindex].dec.rad)

        # Calculate phase offset
        new_data = woffset(data, uvw.T[2], new_uvw.T[2], lambdas)
        filtered.putcol('UVW', new_uvw, startrow=row, nrow = nbl)
        filtered.putcol(data_column, new_data, startrow=row, nrow = nbl)
        tindex += 1
    filtered.close()
    # Don't update the new direction because it is already the original position.
#    ms_phase[0][0][0] = sc[0].ra.rad
#    ms_phase[0][0][1] = sc[0].dec.rad
#    tp.putcol("DELAY_DIR", ms_phase)
#    tp.putcol("PHASE_DIR", ms_phase)
#    tp.putcol("REFERENCE_DIR", ms_phase)
    tp.close()
    print("Done.")


def rotateuvw(uvw, ra, dec, ra0, dec0):
    """
    We calculate new uvw values based on existing uvw values. Whilst this has the effect
    of propagating any existing uvw errors, it has the benefit of being mathematically
    self-consistent.

    Adopted from matrix equation 4.1, in Thompson, Moran, Swenson (3rd edition).
    Let (uvw) = r(ra, dec) * (xyz), then this formula is: r(ra, dec) * r^-1(ra0, dec0)
    """
    u, v, w = uvw.T
    uvwprime = np.empty_like(uvw)

    uvwprime[:, 0] = (
        u * np.cos(ra - ra0)
        + v * np.sin(dec0) * np.sin(ra - ra0)
        - w * np.cos(dec0) * np.sin(ra - ra0)
    )
    uvwprime[:, 1] = (
        -u * np.sin(dec) * np.sin(ra - ra0)
        + v * (np.sin(dec0) * np.sin(dec) * np.cos(ra - ra0) + np.cos(dec0) * np.cos(dec))
        + w * (np.sin(dec0) * np.cos(dec) - np.cos(dec0) * np.sin(dec) * np.cos(ra - ra0))
    )
    uvwprime[:, 2] = (
        u * np.cos(dec) * np.sin(ra - ra0)
        + v * (np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0))
        + w * (np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0))
    )
    return uvwprime

def woffset(data, oldw, neww, lambdas):
    offset = -2j * np.pi * (neww - oldw)
    phase = np.empty_like(data)
    ofs = np.tile(offset, (phase.shape[1],1)).T / lambdas
    for pol in range(data.shape[2]):
        phase[:, :, pol] = ofs #tmp

    return data * np.exp(phase)

# Phase shift a measurement set to a solar system object
# Note that the phase shifting retains the same measurement set direction (since a moving target doesn't have a single direction).

# List of valid targets that can be shifted to
valid_targets = ["sun", "moon", "mercury", "venus", "mars", "jupiter", "saturn", "uranus", "neptune"]
valid_destinations = ["field", "target"]
valid_columns = ["DATA", "CORRECTED_DATA", "MODEL_DATA"]

# Observatory location (change for your favourite observatory - currently ASKAP)
latitude = Angle("-26:41:46.0", unit=u.deg)
longitude = Angle("116:38:13.0", unit=u.deg)
observing_location = EarthLocation(lat=latitude, lon=longitude)

# Basic command-line parameter checks
if len(sys.argv) != 5:
    sys.exit("Usage: %s file.ms target datacolumn destination\nWhere:\n\ttarget=[%s]\n\tdatacolumn=[%s]\n\tdestination=[%s]" %(sys.argv[0], " | ".join(valid_targets), " | ".join(valid_columns), " | ".join(valid_destinations)))
    
ms = sys.argv[1]
if ms[-3:] != ".ms":
    sys.exit("Measurement set requires '.ms' extension")

target = sys.argv[2]
if (target in valid_targets) == False:
    sys.exit("Target must be one of:  %s" %(",".join(valid_targets)))

# Which data column to work on e.g. DATA, CORRECTED_DATA
data_column = sys.argv[3]
if (data_column in valid_columns) == False:
    sys.exit("datacolumn must be one of:  %s" %(",".join(valid_columns)))

# Shift to "target" or shift back to "field"
destination = sys.argv[4]
if (destination in valid_destinations) == False:
    sys.exit("destination must be one of: %s" %(",".join(valid_destinations)))

# The name of the target ms (just append the target name to the ms)
ms_new = ms.replace(".ms", "_%s.ms" %(target))

if destination in ["target"]:
    # Phase rotate from the field to the target solar system object
    phase_rotate_to_target(observing_location, ms, ms_new, target, False, data_column)
else:
    # Phase rotate from the target solar system object back to the field
    ms = ms.replace(".ms", "_%s.ms" %(target))
    ms_new = ms_new = ms.replace(".ms", "_field.ms")

    phase_rotate_to_target(observing_location, ms, ms_new, target, True, data_column)
