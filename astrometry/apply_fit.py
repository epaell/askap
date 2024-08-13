#!/usr/bin/env python
import os
import numpy as np
import sys
import pickle
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle, match_coordinates_sky
import glob
import astropy.coordinates as ac
import astropy.table as at
import warnings

def apply_fit(flist, sbid, field_name, bref, fit):
    flist.sort()
#    print(f'field {field_name}')

    # Depending on the format work out what columns to use
    fit_deg = fit / 3600.0  # Convert fitted offset to degrees
    # get the combined list of coordinates, fluxes, peaks, rms and major axis
    for beam in range(36):             # for each beam
        incat = flist[beam]
#        print(f"Reading {incat}: ", end="")
        cat = at.Table.read(incat)  # read the catalogue
        dra = Angle(np.zeros(len(cat), dtype=np.float128) + fit_deg[0][beam], u.deg)
        ddec = Angle(np.zeros(len(cat), dtype=np.float128) + fit_deg[1][beam], u.deg)
        sc = SkyCoord(Angle(cat['col_ra_deg_cont'], u.deg), Angle(cat['col_dec_deg_cont'], u.deg), frame='fk5')
        sc_ofs = sc.spherical_offsets_by(dra, ddec)
        cat['col_ra_deg_cont'] = sc_ofs.ra.deg
        cat['col_dec_deg_cont'] = sc_ofs.dec.deg
#        print("%d sources %.3f %.3f" %(len(cat), fit_deg[0][beam]*3600.0, fit_deg[1][beam]*3600.0))
        cname = incat.split("/")[-1]
        fout = "fitted/%s" %(cname.replace('.xml', '_fit.xml'))
        cat.write(fout, format='votable', overwrite=True)

warnings.filterwarnings("ignore")

sbid = int(sys.argv[1])
ext = '.xml'
if len(sys.argv)==3:
    ext = sys.argv[2]   # either ".xml" (if starting with original catalogue) or "_fit.xml" (if applying to already corrected catalogue)

bref = 20
pfile = f'pickle/beamwise_fitted_rel_SB{sbid}_{bref}.pkl'

[bref, fit] = np.load(pfile, allow_pickle=True, encoding='bytes')
flist = glob.glob("cat/*SB%d*.components%s" %(sbid, ext))
fname = flist[0]
field_name = fname.split(".")[2]
apply_fit(flist, sbid, field_name, bref, fit)
