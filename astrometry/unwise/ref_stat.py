#!/usr/bin/env python

import sys
import os
import warnings
import numpy as np
import glob
from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u


# --- Define Astrometric matched ---

def ast(c1, c2, radius, nneigh):
    """
    Matches two catalogues within a certain angular separation
    - Uses input of ra/dec in degrees and radius in arcsec
    - Output match has best matches for catalogue 1 and may have sources that
        are matched to the same source in catalogue 2
    - Returns matched (within radius) skycoord array and distances
    :param c1: Catalogue 1
    :param c2: Catalogue 2
    :param radius: Angular separation (degrees) required for match
    :param nneigh: Which closest neighbor to search for: 1 (1st) or 2 (2nd).
    :return: c1_out Cat 1 coords of matched sources
             c2_out Cat 2 coords of matched sources
             d2_out Source angular separation (degrees)
             id_match_1 Indices in Cat 1
             id_match_2 Indices in Cat 2
    """

    id_m, d2_m, d3_m = c1.match_to_catalog_sky(c2, nthneighbor=nneigh)

    mk = d2_m.value <= radius / 3600.
    c1_out = c1[mk]
    c2_out = c2[id_m[mk]]
    d2_out = d2_m.value[mk]

    id_match_1 = np.arange(len(c1))[mk]
    id_match_2 = id_m[mk]

    return c1_out, c2_out, d2_out, id_match_1, id_match_2


def match_self(c1, radius, cmp = 'gt'):
    """
    Matches the same catalogues within a certain angular separation
    - Uses input of ra/dec in degrees and radius in arcsec
    - if type==gt then only keep those sources with matches father away,
    - if lt then chooses nearby matches
    """

    id_m, d2_m, d3_m = c1.match_to_catalog_sky(c1, nthneighbor=2)
    if cmp == 'lt':
        mk = d2_m.value <= radius / 3600.
    elif cmp == 'gt':
        mk = d2_m.value >= radius / 3600.

    c1_out = c1[mk]
    d2_out = d2_m[mk]
    id_match = np.arange(len(c1))[mk]

    return c1_out, d2_out, id_match

def weighted_mean_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    mean = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average(np.power(values-mean, 2.0), weights=weights)
    return (mean, np.sqrt(variance))

def compare(ref_sc, cat_list, dxs, dys, max_sep):
    ares_min = 0.2   # Minimum ratio of int / peak
    ares_max = 1.2   # Maximum ratio of int / peak
    snr_min = 20.0   # Minimum SNR
    rad_match = 15.0 # Maximum radius for matching (arcsec)
    rad_isol = 20.0  # Search radius to ensure isolated source (arcsec)

    # Work out if this is a selavy catalogue or something else (and so which columns to use)
    cat_file = cat_list[0]
    if cat_file.find("selavy") != -1:
        cat_type = "selavy"
    else:
        cat_type = "default"

    # Get the field name
    field_name_pos = cat_file.find("RACS_")
    assert(field_name_pos != -1)
    field_name = cat_file[field_name_pos:field_name_pos+12]

    # Get the SBID
    sbid_pos = cat_file.find("SB")
    assert(sbid_pos != -1)
    cdata = cat_file[sbid_pos+2:].split(".")
    sbid = int(cdata[0])

    all_dx = []
    all_dy = []
    for beam in range(len(cat_list)):
        cat = cat_list[beam]
        #    print("Reading catalogue: %s" %(comparison_cat))
        comparison_cat = cat
        comp_cat = Table.read(comparison_cat)
        if cat_type in ["selavy"]:
            comp_cat.rename_column('col_ra_deg_cont', "ra")
            comp_cat.rename_column('col_dec_deg_cont', "dec")
            comp_cat.rename_column('col_flux_int', "int_flux")
            comp_cat.rename_column('col_flux_peak', "peak_flux")
            comp_cat.rename_column('col_rms_image', "local_rms")
            comp_cat.rename_column('col_rms_fit_gauss', "residual_std")
        
        comp_sc = SkyCoord(Angle(comp_cat["ra"], unit=u.deg), Angle(comp_cat["dec"], unit=u.deg), frame='fk5')
        comp_sc = comp_sc.spherical_offsets_by(dxs[beam], dys[beam])
       
        if max_sep != None:
            mean_sc = SkyCoord(Angle(np.median(comp_sc.ra.deg), unit=u.deg), Angle(np.median(comp_sc.dec.deg), unit=u.deg), frame='fk5')
            seps = mean_sc.separation(comp_sc).deg
            comp_sc = comp_sc[np.where(seps<max_sep)]
            comp_cat = comp_cat[np.where(seps<max_sep)]

        med_rms = np.median(comp_cat["local_rms"])

        # Find isolated sources
        c1, d2, id_match = match_self(comp_sc, rad_isol)
        comp_sc = comp_sc[id_match]
        comp_cat = comp_cat[id_match]

        # Filter out weak sources and sources that may be resolved
        filt_snr = (comp_cat["peak_flux"] / comp_cat["local_rms"] >= snr_min)
        filt_r_lo = (comp_cat["int_flux"] / comp_cat["peak_flux"] >= ares_min)
        filt_r_hi = (comp_cat["int_flux"] / comp_cat["peak_flux"] <= ares_max)
        comp_cat = comp_cat[filt_snr & filt_r_lo & filt_r_hi]
        comp_sc = comp_sc[filt_snr & filt_r_lo & filt_r_hi]
        #    print("Found %d suitable sources" %(len(comp_sc)))
        if len(comp_sc) == 0:
#            print("Skipping %s because no suitable sources" %(comparison_cat))
            continue

        #    print("Cross-matching against alternate catalogue ...")
        c_r_comp, c_srv_c, d2_c_r, id_r_comp, id_srv_c = ast(comp_sc, ref_sc, rad_match, 1)

        if len(comp_sc[id_r_comp]) == 0:
#            print("Skipping %s because no matching sources" %(comparison_cat))
            continue

        dra, ddec = ref_sc[id_srv_c].spherical_offsets_to(comp_sc[id_r_comp])
        # X and Y errors in arcsec
        dx = dra.arcsec
        dy = ddec.arcsec
        if len(all_dx) == 0:
            all_dx = dx
            all_dy = dy
        else:
            all_dx = np.append(all_dx, dx)
            all_dy = np.append(all_dy, dy)
#        print("Beam%02d,%s,%d %.3f+/-%.3f %.3f+/-%.3f" %(beam, comp_type, len(dx), np.mean(dx), np.std(dx), np.mean(dy), np.std(dy)))
    # selavy-image.i.RACS_2049+25.SB45261.cont.taylor.0.restored.conv.components.xml
    cdata = comparison_cat.split(".")

    return field_name, sbid, len(all_dx), np.mean(all_dx), np.median(all_dx), np.std(all_dx), np.mean(all_dy), np.median(all_dy), np.std(all_dy)

warnings.filterwarnings("ignore")

offset_file = ""
if len(sys.argv) != 39:
    sys.exit("Usage:\n\t%s [NVSS|ICRF|VLASS|VLBI|FIRST|SUMSS] catalogues offset_file" %(sys.argv[0]))
ref = sys.argv[1]
cat_list = sys.argv[2:-1]
offset_file = sys.argv[-1]
radius = None

cat_list.sort()

if  (ref in ["NVSS", "VLASS", "SUMSS", "ICRF", "VLBI", "FIRST"]) == False:
    sys.exit("Usage:\n\t%s [NVSS|ICRF|VLASS|VLBI|FIRST|SUMSS] catalogues" %(sys.argv[0]))
    

if ref == "NVSS":
    print("Reading NVSS catalogue")
    ref_cat = Table.read('ref/NVSS_vizier.fits')
    ref_cat = ref_cat[np.where(ref_cat["MajAxis"] < 50.0)]
    print("Found %d sources read" %(len(ref_cat)))
    ref_sc = SkyCoord(Angle(ref_cat["RAJ2000"], unit=u.deg), Angle(ref_cat["DEJ2000"], unit=u.deg), frame='fk5')
elif ref == "VLASS":
    print("Reading VLASS catalogue")
    ref_cat = Table.read('ref/vlass_isolated.fits')
    print("Found %d sources read" %(len(ref_cat)))
    ref_sc = SkyCoord(Angle(ref_cat["RA_1"], unit=u.deg), Angle(ref_cat["DEC_1"], unit=u.deg), frame='fk5')
elif ref == "SUMSS":
    print("Reading SUMSS catalogue")
    ref_cat = Table.read('ref/SUMSS_MGPS2.fits')
    print("Found %d sources read" %(len(ref_cat)))
    ref_sc = SkyCoord(Angle(ref_cat["RAJ2000"], unit=u.deg), Angle(ref_cat["DEJ2000"], unit=u.deg), frame='fk5')
elif ref == "ICRF":
    print("Reading ICRF catalogue")
    ref_cat = Table.read('ref/icrf.csv')
    print("Found %d sources read" %(len(ref_cat)))
    ref_sc = SkyCoord(Angle(ref_cat["RA"], unit=u.deg), Angle(ref_cat["Dec"], unit=u.deg), frame='fk5')
elif ref == "VLBI":
    print("Reading VLBI catalogue")
    ref_cat = Table.read('ref/rfc_cat_compact.csv')
    print("Found %d sources read" %(len(ref_cat)))
    ref_sc = SkyCoord(Angle(ref_cat["RAJ2000"], unit=u.deg), Angle(ref_cat["DEJ2000"], unit=u.deg), frame='fk5')
elif ref == "FIRST":
    print("Reading FIRST catalogue")
    ref_cat = Table.read('ref/first.fits')
    ref_cat = ref_cat[np.where(ref_cat["DEJ2000"] < 50.0)]
    ref_cat = ref_cat[np.where(ref_cat["Maj"] < 15.0)]
    ref_cat = ref_cat[np.where(ref_cat["Fpeak"] > 1.0)]
    print("Found %d sources read" %(len(ref_cat)))
    ref_sc = SkyCoord(Angle(ref_cat["RAJ2000"], unit=u.deg), Angle(ref_cat["DEJ2000"], unit=u.deg), frame='fk5')

shift = Table.read(offset_file, format="csv")
dxs = Angle(shift["DXS"].value - shift["DXS"].value, u.arcsec)
dys = Angle(shift["DYS"].value - shift["DYS"].value, u.arcsec)
field_name, sbid, onx, odx_mean, odx_med, odx_std, ody_mean, ody_med, ody_std = compare(ref_sc, cat_list, dxs, dys, radius)

dxs = -Angle(shift["DXS"].value, u.arcsec)
dys = -Angle(shift["DYS"].value, u.arcsec)
field_name, sbid, nx, dx_mean, dx_med, dx_std, dy_mean, dy_med, dy_std = compare(ref_sc, cat_list, dxs, dys, radius)
fout = open(f"stats/SB{sbid}.{field_name}.ref.{ref}.stat.csv", "wt")
fout.write("FIELD_NAME,SBID,ON,ODX_MEAN,ODX_MEDIAN,ODX_STD,ODY_MEAN,ODY_MEDIAN,ODY_STD,N,DX_MEAN,DX_MEDIAN,DX_STD,DY_MEAN,DY_MEDIAN,DY_STD\n")
fout.write(f"{field_name},{sbid},{onx},{odx_mean:.3f},{odx_med:.3f},{odx_std:.3f},{ody_mean:.3f},{ody_med:.3f},{ody_std:.3f},{nx},{dx_mean:.3f},{dx_med:.3f},{dx_std:.3f},{dy_mean:.3f},{dy_med:.3f},{dy_std:.3f}\n")
fout.close()

