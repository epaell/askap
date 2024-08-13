#!/usr/bin/env python
import glob
import os
import sys
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
from multiprocessing import Process

db_base_path = '/Users/len067/Desktop/aces/calibration/askap_surveys'
max_concurrent = 4

def my_process(process):
    db_data = Table.read('%s/RACS/db/epoch_9/field_data.csv' %(db_base_path))

    flist = glob.glob("cat/selavy*beam00*components.xml")
    flist.sort()
    sbids = []
    for fname in flist:
        sbid = int(fname.split('.')[3][2:])
        sbids.append(sbid)
    for index in range(len(sbids)):
        if (index % max_concurrent) == process:
            sbid = sbids[index]
            field = db_data[np.where(db_data["SBID"]==sbid)][0]
            os.system("./fit_astrometry.py %d .xml" %(sbid))
            os.system("./plot_fit.py %d .xml" %(sbid))
            os.system("./apply_fit.py %d .xml" %(sbid))

# Read field data
# INDEX,SRC,FIELD_NAME,SBID,SCAN,CAL_SBID,STATE,RA_HMS,DEC_DMS,RA_DEG,DEC_DEG,GAL_LONG,GAL_LAT,
# POL_AXIS,SCAN_START,SCAN_LEN,SCAN_TINT,NBEAMS_I,MOSAIC_TIME,NPIXELS_V,
# PSF_MAJOR,PSF_MINOR,PSF_ANGLE,MINIMUM,MAXIMUM,MED_RMS_uJy,MODE_RMS_uJy,STD_RMS_uJy,MIN_RMS_uJy,MAX_RMS_uJy,
# SELAVY_TIME,NUM_SELAVY,SELECT,DEFECT,ANOMALY,COMMENT,MinUV

db_data = Table.read('%s/field_data.csv' %(db_path))

flist = glob.glob("cat/selavy*beam00*components.xml")
flist.sort()
sbids = []
for fname in flist:
    sbid = int(fname.split('.')[3][2:])
    sbids.append(sbid)

processes = []
for m in range(max_concurrent):
    p = Process(target=my_process, args=(m,))
    p.start()
    processes.append(p)

for p in processes:
   p.join()
