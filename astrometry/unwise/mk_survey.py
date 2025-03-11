#!/usr/bin/env python

from astropy.coordinates import Angle, SkyCoord
from astropy import units as u
from astroquery import vizier
from astropy.table import Table, unique, vstack
import numpy as np
import time
import sys
import os

# Extracts per-beam catalogues and merges them into a single field catalogue containing all unwise sources.
def field_survey(survey, ref_cat, ref_sc, field_name, beam_sc, radius_deg: float=5.):
    unwise = "II/363/unwise"
    field_cat = None
    for beam in range(36):
        print(f"Extracting {survey} table for field {field_name} and beam:{beam}")
        if os.path.exists(f"{survey}/{survey}_{field_name}_beam{beam:02d}.fits") == False:
            seps = beam_sc[beam].separation(ref_sc).deg
            beam_cat = ref_cat[np.where(seps<radius_deg)]
            if len(beam_cat) == 0:
                break
            beam_cat.write(f"{survey}/{survey}_{field_name}_beam{beam:02d}.fits", format="fits", overwrite=False)
        else:
            print(f"Found cached file, reading {survey}_{field_name}_beam{beam:02d}.fits")
            beam_cat = Table.read(f"{survey}/{survey}_{field_name}_beam{beam:02d}.fits", format="fits")
        print("%d rows downloaded" %(len(beam_cat)))
        print("Merging into field catalogue")
        if beam == 0:
            field_cat = beam_cat
        else:
            field_cat = unique(vstack([field_cat, beam_cat]))

    if field_cat != None:
        field_cat.write(f"{survey}/{survey}_{field_name}.fits", format="fits", overwrite=True)
        print("Cleaning cached beam catalogues")
        os.system(f"rm -fr {survey}/{survey}_{field_name}_beam*fits")
    return

survey = "FIRST"

os.system("mkdir {survey}")
if survey == "NVSS":
    print("Reading NVSS catalogue")
    ref_cat = Table.read('ref/NVSS_vizier.fits')
    ref_cat = ref_cat[np.where(ref_cat["MajAxis"] < 50.0)]
    print("Found %d sources read" %(len(ref_cat)))
    ref_sc = SkyCoord(Angle(ref_cat["RAJ2000"], unit=u.deg), Angle(ref_cat["DEJ2000"], unit=u.deg), frame='fk5')
elif survey == "VLASS":
    print("Reading VLASS catalogue")
    ref_cat = Table.read('ref/vlass_isolated.fits')
    print("Found %d sources read" %(len(ref_cat)))
    ref_sc = SkyCoord(Angle(ref_cat["RA_1"], unit=u.deg), Angle(ref_cat["DEC_1"], unit=u.deg), frame='fk5')
    ref_cat.rename_column('RA_1', 'RAJ2000')
    ref_cat.rename_column('DEC_1', 'DEJ2000')
elif survey == "SUMSS":
    print("Reading SUMSS catalogue")
    ref_cat = Table.read('ref/SUMSS_MGPS2.fits')
    print("Found %d sources read" %(len(ref_cat)))
    ref_sc = SkyCoord(Angle(ref_cat["RAJ2000"], unit=u.deg), Angle(ref_cat["DEJ2000"], unit=u.deg), frame='fk5')
elif survey == "ICRF":
    print("Reading ICRF catalogue")
    ref_cat = Table.read('ref/icrf.csv')
    print("Found %d sources read" %(len(ref_cat)))
    ref_sc = SkyCoord(Angle(ref_cat["RA"], unit=u.deg), Angle(ref_cat["Dec"], unit=u.deg), frame='fk5')
    ref_cat.rename_column('RA', 'RAJ2000')
    ref_cat.rename_column('Dec', 'DEJ2000')
elif survey == "VLBI":
    print("Reading VLBI catalogue")
    ref_cat = Table.read('ref/rfc_cat_compact.csv')
    print("Found %d sources read" %(len(ref_cat)))
    ref_sc = SkyCoord(Angle(ref_cat["RAJ2000"], unit=u.deg), Angle(ref_cat["DEJ2000"], unit=u.deg), frame='fk5')
elif survey == "FIRST":
    print("Reading FIRST catalogue")
    ref_cat = Table.read('ref/first.fits')
    ref_cat = ref_cat[np.where(ref_cat["DEJ2000"] < 50.0)]
    ref_cat = ref_cat[np.where(ref_cat["Maj"] < 15.0)]
    ref_cat = ref_cat[np.where(ref_cat["Fpeak"] > 1.0)]
    print("Found %d sources read" %(len(ref_cat)))
    ref_sc = SkyCoord(Angle(ref_cat["RAJ2000"], unit=u.deg), Angle(ref_cat["DEJ2000"], unit=u.deg), frame='fk5')

for col in ref_cat.colnames:
    if col in ['RAJ2000', 'DEJ2000']:
        continue
    ref_cat.remove_column(col)

field_beams = Table.read("closepack36_beams.fits")
print("Read %d beam positions" %(len(field_beams)))

unique_fields = np.unique(field_beams["FIELD_NAME"])

if False: # Set to True if you only want to create catalogues for "bad" fields.
    unique_fields = [ "RACS_1146-37",  "RACS_1232+09",  "RACS_1650-41",  "RACS_1701-37",  "RACS_1710-32",  "RACS_1714-46",  
        "RACS_1724-28",  "RACS_1727-37",  "RACS_1731-23",  "RACS_1735-32",  "RACS_1735-37",  "RACS_1748-18",  "RACS_1748-28",  "RACS_1748-41",  "RACS_1754-23",  
        "RACS_1800-14",  "RACS_1800-32",  "RACS_1811-18",  "RACS_1812-28",  "RACS_1816-09",  "RACS_1817-23",  "RACS_1821-14",  "RACS_1824-32",  "RACS_1833-18",  
        "RACS_1836-28",  "RACS_1837+00",  "RACS_1837+04",  "RACS_1837-04",  "RACS_1837-09",  "RACS_1839-23",  "RACS_1843-14",  "RACS_1859+00",  "RACS_1859+04",  
        "RACS_1859+09",  "RACS_1859-04",  "RACS_1905+14",  "RACS_1918+18",  "RACS_1920+04",  "RACS_1920+09",  "RACS_1925+23",  "RACS_1927+14",  "RACS_1936+41",  
        "RACS_1941+18",  "RACS_1948+23",  "RACS_1948+28",  "RACS_2004+32",  "RACS_2004+37",  "RACS_2004+41",  "RACS_2032+41",  "RACS_2044+46",  "RACS_2114+46"]

for field_name in unique_fields:
    print("Processing %s" %(field_name))
    if os.path.exists(f"{survey}/{survey}_{field_name}.fits") == True:
        continue
    beam_inf = field_beams[np.where(field_beams["FIELD_NAME"]==field_name)]
    if len(beam_inf) != 36:
        print("Missing beams")
        continue
    field_survey(survey, ref_cat, ref_sc, field_name, SkyCoord(Angle(beam_inf["RA_DEG"], unit=u.deg),Angle(beam_inf["DEC_DEG"], unit=u.deg), frame='fk5'))
