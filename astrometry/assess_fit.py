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

def compare(cat_list, beam_sc, maxd_beam_to_beam = 1.0, maxd_from_beam=0.75):
    ares_min = 0.2   # Minimum ratio of int / peak
    ares_max = 1.2   # Maximum ratio of int / peak
    snr_min = 20.0   # Minimum SNR
    rad_match = 15.0 # Maximum radius for matching (arcsec)
    rad_isol = 30.0  # Search radius to ensure isolated source (arcsec)

    
    # Pre-read all of the catalogues
    cats = []
    for cat_file in cat_list:
        cats.append(Table.read(cat_file))

    all_dx = []
    all_dy = []
    all_weights = []
    for beam1 in range(len(cat_list)-1):
        # Set up the reference catalogue (beam 1)
        ref_cat_file = cat_list[beam1]
        ref_cat = cats[beam1]
        ref_sc = SkyCoord(Angle(ref_cat["col_ra_deg_cont"], unit=u.deg), Angle(ref_cat["col_dec_deg_cont"], unit=u.deg), frame='fk5')

        beam1_sc = beam_sc[beam1]
        seps1 = beam1_sc.separation(ref_sc).deg
        ref_sc = ref_sc[np.where(seps1<maxd_from_beam)]
        ref_cat = ref_cat[np.where(seps1<maxd_from_beam)]

        if len(ref_sc) < 2:
#            print("Skipping beam %d: ra=%.3f,dec=%.3f" %(beam1, mean_sc1.ra.deg, mean_sc1.dec.deg))
            continue

        # Find isolated sources
        c1, d2, id_match = match_self(ref_sc, rad_isol)
        ref_sc = ref_sc[id_match]
        ref_cat = ref_cat[id_match]
        ref_snr = ref_cat["col_flux_int"] / ref_cat["col_rms_image"]

        # Filter out weak sources and sources that may be resolved
        filt_snr = (ref_snr >= snr_min)
        filt_r_lo = (ref_cat["col_flux_int"] / ref_cat["col_flux_peak"] >= ares_min)
        filt_r_hi = (ref_cat["col_flux_int"] / ref_cat["col_flux_peak"] <= ares_max)
        ref_cat = ref_cat[filt_snr & filt_r_lo & filt_r_hi]
        ref_sc = ref_sc[filt_snr & filt_r_lo & filt_r_hi]
        ref_snr = ref_snr[filt_snr & filt_r_lo & filt_r_hi].value

        
        if len(ref_sc) < 2:
#            print("Skipping beam %d B" %(beam1))
            continue
        
        for beam2 in range(beam1+1, len(cat_list)):
            beam2_sc = beam_sc[beam2]
            if beam1_sc.separation(beam2_sc).deg > maxd_beam_to_beam: # Don't bother checking if the beams are highly separated
                continue
            comp_cat = cats[beam2]

            comp_sc = SkyCoord(Angle(comp_cat["col_ra_deg_cont"], unit=u.deg), Angle(comp_cat["col_dec_deg_cont"], unit=u.deg), frame='fk5')

            seps2 = beam2_sc.separation(comp_sc).deg
            comp_sc = comp_sc[np.where(seps2<maxd_from_beam)]
            comp_cat = comp_cat[np.where(seps2<maxd_from_beam)]

            if len(comp_sc) < 2:
                continue

            # Find isolated sources
            c1, d2, id_match = match_self(comp_sc, rad_isol)
            comp_sc = comp_sc[id_match]
            comp_cat = comp_cat[id_match]
            comp_snr = comp_cat["col_flux_int"] / comp_cat["col_rms_image"]

            # Filter out weak sources and sources that may be resolved
            filt_snr = (comp_snr >= snr_min)
            filt_r_lo = (comp_cat["col_flux_int"] / comp_cat["col_flux_peak"] >= ares_min)
            filt_r_hi = (comp_cat["col_flux_int"] / comp_cat["col_flux_peak"] <= ares_max)
            comp_cat = comp_cat[filt_snr & filt_r_lo & filt_r_hi]
            comp_sc = comp_sc[filt_snr & filt_r_lo & filt_r_hi]
            comp_snr = comp_snr[filt_snr & filt_r_lo & filt_r_hi].value
            #    print("Found %d suitable sources" %(len(comp_sc)))
            if len(comp_sc) < 2:
                continue

            #    print("Cross-matching against alternate catalogue ...")
            c_r_comp, c_srv_c, d2_c_r, id_r_comp, id_srv_c = ast(comp_sc, ref_sc, rad_match, 1)

            if len(comp_sc[id_r_comp]) < 2:
                continue
            dra, ddec = ref_sc[id_srv_c].spherical_offsets_to(comp_sc[id_r_comp])
            # X and Y errors in arcsec
            dx = dra.arcsec
            dy = ddec.arcsec
            if len(all_dx) == 0:
                all_dx = dx
                all_dy = dy
                all_weights = np.sqrt(np.power(ref_snr[id_srv_c], 2.0), np.power(comp_snr[id_r_comp], 2.0))
            else:
                all_dx = np.append(all_dx, dx)
                all_dy = np.append(all_dy, dy)
                all_weights = np.append(all_weights, np.sqrt(np.power(ref_snr[id_srv_c], 2.0), np.power(comp_snr[id_r_comp], 2.0)))
#        if beam == 1:
#            break
    # selavy-image.i.RACS_2049+25.SB45261.cont.taylor.0.restored.conv.components.xml
    cdata = cat_list[0].split(".")
    sbid = cdata[3][2:]
    field_name = cdata[2]

#    print("%s,%s,%s,%d %.3f+/-%.3f %.3f+/-%.3f" %(comp_type, field_name, sbid, len(all_dx), np.mean(all_dx), np.std(all_dx), np.mean(all_dy), np.std(all_dy)))
    if len(all_dx) == 0:
        return 0,0.0,0.0,0.0,0.0
    # Get a weighted mean and std to take into consideration the brightness of the source
    xmean, xstd = weighted_mean_and_std(all_dx, all_weights)
    ymean, ystd = weighted_mean_and_std(all_dy, all_weights)
    return len(all_dx), xmean, xstd, ymean, ystd

warnings.filterwarnings("ignore")

#if len(sys.argv) != 3:
#    sys.exit("Usage:\n\t%s [NVSS|ICRF|VLASS] sbid" %(sys.argv[0]))
#ref = sys.argv[1]
sbid = int(sys.argv[1])

if len(sys.argv)==3:
    maxd_from_beam = float(sys.argv[2])
    maxd_beam_to_beam = 1.0 # Maximum distance between beams to consider

if len(sys.argv)==4:
    maxd_beam_to_beam = float(sys.argv[2]) # Maximum distance between beams to consider
    maxd_from_beam = float(sys.argv[3])

db_base_path = "/Users/len067/Desktop/aces/calibration/askap_surveys"

cat_list = glob.glob("cat/*SB%d*beam*components.xml" %(sbid))
cat_list.sort()
cdata = cat_list[0].split(".")
field_name = cdata[2]

# Read beam positions for current field
beam_inf = Table.read('%s/RACS/db/epoch_9/beam_inf_%d-%s.csv' %(db_base_path, sbid, field_name))
beam_sc = SkyCoord(Angle(beam_inf["RA_DEG"], unit=u.deg),Angle(beam_inf["DEC_DEG"], unit=u.deg), frame='fk5')

nso, xo_mean,xo_std, yo_mean,yo_std = compare(cat_list, beam_sc, maxd_beam_to_beam, maxd_from_beam)

cat_list = glob.glob("fitted/*SB%d*beam*components_fit.xml" %(sbid))
cat_list.sort()
nsf, xf_mean,xf_std, yf_mean,yf_std = compare(cat_list, beam_sc, maxd_beam_to_beam, maxd_from_beam)
flag=" "
if (xf_std > xo_std) or (yf_std > yo_std):
    flag = "*"
print("%d,%s,%4d,%6.3f,%6.3f,%6.3f,%6.3f,%4d,%6.3f,%6.3f,%6.3f,%6.3f,%s" %(sbid, field_name,nso,xo_mean,xo_std, yo_mean,yo_std,nsf, xf_mean,xf_std, yf_mean,yf_std,flag))