#!/usr/bin/env python

import os
import numpy as np
import sys
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle, match_coordinates_sky
import glob
from scipy.optimize import minimize
import astropy.coordinates as ac
import astropy.table as at
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import warnings

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

# Read catalogue and extract the necessary quantities
def process_sbid(catalogue_paths, ref_beam):
    ares_min = 0.2   # Minimum ratio of int / peak
    ares_max = 1.8   # Maximum ratio of int / peak
    snr_min = 10.0   # Minimum SNR
    rad_match = 10.0 # Maximum radius for matching (arcsec)
    rad_isol = 30.0  # Search radius to ensure isolated source (arcsec)

    path_name = catalogue_paths[0]
    cat_file = path_name.split("/")[-1]
    # Work out if this is a selavy catalogue or something else (and so which columns to use)
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

    if cat_type in ["selavy"]:
        output_prefix = f"SB{sbid}.{field_name}.selavy."
    else:
        output_prefix = f"SB{sbid}.{field_name}."

    srcs = []
    fluxes = []
    pkfluxes = []
    im_rms = []
    fit_rms = []
    # get the combined list of coordinates, fluxes, peaks, rms and major axis
    for incat in catalogue_paths:             # for each beam
        bpos = incat.find("beam")
        beam = int(incat[bpos+4:bpos+6])
        cat = at.Table.read(incat)  # read the catalogue
        if (cat_type in ["selavy"]):
            cat.rename_column('col_ra_deg_cont', "ra")
            cat.rename_column('col_dec_deg_cont', "dec")
            cat.rename_column('col_flux_int', "int_flux")
            cat.rename_column('col_flux_peak', "peak_flux")
            cat.rename_column('col_rms_image', "local_rms")
            cat.rename_column('col_rms_fit_gauss', "residual_std")

        comp_sc = SkyCoord(Angle(cat["ra"], unit=u.deg), Angle(cat["dec"], unit=u.deg))

        # Find isolated sources
        c1, d2, id_match = match_self(comp_sc, rad_isol)
        comp_sc = comp_sc[id_match]
        cat = cat[id_match]

        # Filter out weak sources and sources that may be resolved
        filt_snr = (cat["int_flux"] / cat["local_rms"] >= snr_min)
        filt_r_lo = (cat["int_flux"] / cat["peak_flux"] >= ares_min)
        filt_r_hi = (cat["int_flux"] / cat["peak_flux"] <= ares_max)
        cat = cat[filt_snr & filt_r_lo & filt_r_hi]
        comp_sc = comp_sc[filt_snr & filt_r_lo & filt_r_hi]
        
        ra = cat['ra']         
        dec = cat['dec']
        srcs.append(SkyCoord(ra, dec, unit='deg'))
        fluxes.append(cat['int_flux'])
        pkfluxes.append(cat['peak_flux'])
        im_rms.append(cat['local_rms'])

    # srcs, fluxes, pkfluxes,im_rms and maj all have 36 (beams) x nsrcs
    # At this point have srcs, fluxes, pkfluxes and im_rms
    
    # get the SNR for each source in each beam
    snr = [np.array(f)/np.array(i) for f,i in zip(fluxes,im_rms)]

    # Construct dictionaries keyed by the beam pair
    # Each has lists of indices, separations, position angles, fluxes, SNRs, major axes, and deconvolved major axes
    # do both upper and lower triangles
    max_sep = rad_match * u.arcsec
    indices = {}
    seps = {}
    pas = {}
    fluxdict = {}
    snrdict = {}
    majdict = {}
    majdcdict = {}
    sources = {}
    for j in range(36):
        for i in range(36):
            if i != j:
                key = f'{j:02d}-{i:02d}'
                q = match_coordinates_sky(srcs[j], srcs[i], nthneighbor=1, storekdtree='kdtree_sky')
                idx, sep2d, dist3d = q
                sel = (sep2d < max_sep)
                indices[key] = idx[sel]
                seps[key] = sep2d[sel]
                sj = srcs[j][sel]
                si = srcs[i][idx[sel]]
                pas[key] = sj.position_angle(si)
                fluxdict[key] = [fluxes[j][sel],fluxes[i][idx[sel]]]
                snrdict[key] = [snr[j][sel],snr[i][idx[sel]]]
                sources[key] = [sj,si]

    # indices[j,i], seps[j,i], pas[j,i], fluxdict[j,i], snrdist[j,i] and sources[j,i] are dictionaries holding cross-matched value pairs between beams j and i
    
    # For each of the beam pairs, calculate stats relating to offsets (median, std, ncross-matches)
    dxs = np.zeros([36,36])
    dys = np.zeros([36,36])
    dxs_std = np.zeros([36,36])
    dys_std = np.zeros([36,36])
    ns = np.zeros([36,36], dtype=int)
    for b in range(36):
        for i in range(36):
            if i != b:
                key = f'{b:02d}-{i:02d}'
                dx = seps[key].arcsec * np.sin(pas[key])
                dy = seps[key].arcsec * np.cos(pas[key])
                dxs[b,i] = np.median(dx) #dx.mean()
                dys[b,i] = np.median(dy) #dy.mean()
                dxs_std[b,i] = np.std(dx) #dx.std()
                dys_std[b,i] = np.std(dy) #dy.std()
                ns[b,i] = len(dx)

    # Some dxs and dys will have nans in them if there were no fits.

    # Objective function to be minimized
    def ofunc(p, ds, wt):
        r, rm = pfunc(p, ds)
        diffsq = (r.T - rm)**2
        ws = (diffsq*wt).sum()/wt.sum()     # Weight by SNR
        return ws

    def pfunc(p, ds):
        n = ds.shape[0]                     # Shape of the offsets array
        r = np.zeros([n,n])                 # start with zero offsets
        for j in range(n):                  # For each beam "j"
            for i in range(n):              # check against beam "i"
                r[i,j] = ds[i,j] - p[j]     # apply the test offset p to all of the offsets
        rm = np.median(r, axis=1)           # Find the median offset
        return r, rm

    fitted_rel = []
    sl = slice(0,36)
    wt = ns[sl,sl]      # wt has the number of cross-matched sources in a 36 x 36 array
    
    # iterate through each axis i.e. dx and then dy
    for d in [dxs,dys]:
        ds = d[sl,sl]           # ds contains the distance for the current axis from each beam to each other beam [36 x 36]
#        p = ds[ref_beam - sl.start,:] # not sure why Dave had this
        p = ds[ref_beam,:]          # Initial guess for offsets (assume beam ref_beam has 0 offset)
        ds[np.isnan(ds)] = 0.   # For offsets that have nans set them to 0
        popt = minimize(ofunc, p, args=(ds, wt), method='Nelder-Mead')
        fitted_rel.append(pfunc(popt.x, ds)[1])
    dxs = np.array(fitted_rel)[0,:]
    dys = np.array(fitted_rel)[1,:]
    dxs -= dxs[ref_beam]
    dys -= dys[ref_beam]

    if cat_type in ["selavy"]:
        fout = open(f"SB{sbid}.{field_name}.selavy.shifts.csv", "wt")
    else:
        fout = open(f"SB{sbid}.{field_name}.shifts.csv", "wt")
    fout.write("BEAM,DXS,DYS\n")
    for beam in range(36):
        fout.write("%d,%f,%f\n" %(beam, -dxs[beam], -dys[beam]))
    fout.close()
    
    return

warnings.filterwarnings("ignore")

flist = sys.argv[1:]
flist.sort()
process_sbid(catalogue_paths=flist, ref_beam=20)
