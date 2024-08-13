#!/usr/bin/env python

import sys
import os
import warnings
import numpy as np
import glob
from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u

db_base_path = "/Users/len067/Desktop/aces/calibration/askap_surveys"
survey = "RACS"
epoch = 9 # RACS-low3

warnings.filterwarnings("ignore")

#if len(sys.argv) != 3:
#    sys.exit("Usage:\n\t%s [NVSS|ICRF|VLASS] sbid" %(sys.argv[0]))
#ref = sys.argv[1]
db_path = '%s/%s/db/epoch_%d/' %(db_base_path, survey, epoch)
sbid = int(sys.argv[1])
if len(sys.argv)==3:
    maxd_from_beam = float(sys.argv[2])

# Read field data
# INDEX,SRC,FIELD_NAME,SBID,SCAN,CAL_SBID,STATE,RA_HMS,DEC_DMS,RA_DEG,DEC_DEG,GAL_LONG,GAL_LAT,
# POL_AXIS,SCAN_START,SCAN_LEN,SCAN_TINT,NBEAMS_I,MOSAIC_TIME,NPIXELS_V,
# PSF_MAJOR,PSF_MINOR,PSF_ANGLE,MINIMUM,MAXIMUM,MED_RMS_uJy,MODE_RMS_uJy,STD_RMS_uJy,MIN_RMS_uJy,MAX_RMS_uJy,
# SELAVY_TIME,NUM_SELAVY,SELECT,DEFECT,ANOMALY,COMMENT,MinUV

db_data = Table.read('%s/field_data.csv' %(db_path))
field = db_data[np.where(db_data["SBID"]==sbid)][0]
field_name = field["FIELD_NAME"]

# Read beam positions for current field
beam_inf = Table.read('%s/beam_inf_%d-%s.csv' %(db_path,sbid, field_name))
beam_sc = SkyCoord(Angle(beam_inf["RA_DEG"], unit=u.deg),Angle(beam_inf["DEC_DEG"], unit=u.deg), frame='fk5')

print("Field name: %s" %(field_name))
print("PSF: %.1f x %.1f at %.1f deg" %(field["PSF_MAJOR"], field["PSF_MINOR"], field["PSF_ANGLE"]))
print("RMS: %0.1f uJy/beam" %(field["MODE_RMS_uJy"]))