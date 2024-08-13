#!/usr/bin/env python
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import os

db_base_path = "/Users/len067/Desktop/aces/calibration/askap_surveys"
survey = "RACS"
epoch = 9 # RACS-low3

# Read field data
# INDEX,SRC,FIELD_NAME,SBID,SCAN,CAL_SBID,STATE,RA_HMS,DEC_DMS,RA_DEG,DEC_DEG,GAL_LONG,GAL_LAT,
# POL_AXIS,SCAN_START,SCAN_LEN,SCAN_TINT,NBEAMS_I,MOSAIC_TIME,NPIXELS_V,
# PSF_MAJOR,PSF_MINOR,PSF_ANGLE,MINIMUM,MAXIMUM,MED_RMS_uJy,MODE_RMS_uJy,STD_RMS_uJy,MIN_RMS_uJy,MAX_RMS_uJy,
# SELAVY_TIME,NUM_SELAVY,SELECT,DEFECT,ANOMALY,COMMENT,MinUV
db_path = "%s/%s/db/epoch_%d/field_data.csv" %(db_base_path, survey, epoch)
if os.path.exists(db_path) == False:
    sys.exit("Unable to load database")
field_data = Table.read(db_path, format='csv')

fig, axes = plt.subplots(2,2, figsize=(12,8))

astro = Table.read("check_all.csv")
decs = []
for index in range(len(astro)):
    sbid = astro[index]["sbid"]
    field = field_data[np.where(field_data["SBID"]==sbid)][0]
    dec = field["DEC_DEG"]
    decs.append(dec)
    if astro[index]["xfstd"]>0.25:
        print("%d,%s,high sigma in RA:%.3f" %(sbid, field["FIELD_NAME"], astro[index]["xfstd"]))
    if astro[index]["yfstd"]>0.59:
        print("%d,%s,high sigma in Dec:%.3f" %(sbid, field["FIELD_NAME"], astro[index]["yfstd"]))
decs = np.array(decs)

dd = []
psf_maj = []
psf_min = []
for d in range(-87,50,5):
    ddata = field_data[np.where(field_data["DEC_DEG"]>d-2.5)]
    ddata = ddata[np.where(ddata["DEC_DEG"]<=d+2.5)]
    dd.append(d)
    psf_maj.append(np.median(ddata["PSF_MAJOR"]))
    psf_min.append(np.median(ddata["PSF_MINOR"]))

dd = np.array(dd)
psf_maj = np.array(psf_maj)/40.0
psf_min = np.array(psf_min)/40.0

ox_std = astro["xostd"]
oy_std = astro["yostd"]
fx_std = astro["xfstd"]
fy_std = astro["yfstd"]

ox_mean = astro["xomean"]
oy_mean = astro["yomean"]
fx_mean = astro["xfmean"]
fy_mean = astro["yfmean"]

axes[0,0].plot(decs, ox_std, label="original", c="grey",marker="+", ls="None")
axes[0,0].plot(decs, fx_std, label="corrected", c="black", marker="+", ls="None")
axes[0,0].plot(dd, psf_min, label="PSF/2SNR", c="red")

axes[0,1].plot(decs, oy_std, label="original", c="grey", marker="+", ls="None")
axes[0,1].plot(decs, fy_std, label="corrected", c="black", marker="+", ls="None")
axes[0,1].plot(dd, psf_maj, label="PSF/2SNR", c="red")

axes[1,0].plot(decs, ox_mean, label="original", c="grey",marker="+", ls="None")
axes[1,0].plot(decs, fx_mean, label="corrected", c="black", marker="+", ls="None")

axes[1,1].plot(decs, oy_mean, label="original", c="grey", marker="+", ls="None")
axes[1,1].plot(decs, fy_mean, label="corrected", c="black", marker="+", ls="None")

axes[0,0].set_xlabel('Declination (degrees)')
axes[0,0].set_ylabel('$\sigma$dRA (arcsec)')
axes[0,1].set_xlabel('Declination (degrees)')
axes[0,1].set_ylabel('$\sigma$dDec (arcsec)')

axes[1,0].set_xlabel('Declination (degrees)')
axes[1,0].set_ylabel('mean dRA (arcsec)')
axes[1,1].set_xlabel('Declination (degrees)')
axes[1,1].set_ylabel('mean dDec (arcsec)')

if True:
    axes[0,0].set_xlim(-90,+70)
    axes[0,0].set_ylim(0.0,2.0)
    axes[0,1].set_xlim(-90,+70)
    axes[0,1].set_ylim(0.0,2.0)

    axes[1,0].set_xlim(-90,+70)
    axes[1,0].set_ylim(-2.0,0.5)
    axes[1,1].set_xlim(-90,+70)
    axes[1,1].set_ylim(-2.0,0.5)

axes[0,0].legend()
axes[0,1].legend()
axes[1,0].legend()
axes[1,1].legend()
plt.show()
plt.close()