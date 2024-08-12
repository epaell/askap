#!/usr/bin/env python

import os
import numpy as np
import sys
import pickle
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle, match_coordinates_sky
import glob
from scipy.optimize import minimize
import astropy.coordinates as ac
import astropy.table as at
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import warnings

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

def get_cols(fmt):
    col = {}
    if fmt == 'xml':
        col['ra'] = 'col_ra_deg_cont'
        col['dec'] = 'col_dec_deg_cont'
        col['flux'] = 'col_flux_int'
        col['pkflux'] = 'col_flux_peak'
        col['rms'] = 'col_rms_image'
        col['fitrms'] = 'col_rms_fit_gauss'
        col['chi2'] = 'col_chi_squared_fit'
        col['maj'] = 'col_maj_axis'
        col['maj_deconv'] = 'col_maj_axis_deconv'
    elif fmt == 'vot':
        col['ra'] = 'ra'
        col['dec'] = 'dec'
        col['flux'] = 'int_flux'
        col['pkflux'] = 'peak_flux'
        col['rms'] = 'local_rms'
        col['fitrms'] = 'fit_rms'
        col['chi2'] = 'chi_squared'
        col['maj'] = 'psf_a'
        col['maj_deconv'] = 'maj_axis_deconv'
    return col


# # Read catalogue 
# and extract the necessary quantities
def process_sbid(flist, sbid, field_name):
    flist.sort()
    print(f'field {field_name}')

    # Check the catalogue format
    fmt = [fname.split('.')[-1] for fname in flist]
    # if non-homogenous then bail out with an error
    if len(set(fmt)) > 1:
        sys.exit('Different formats in the files: ', fmt)
    # Depending on the format work out what columns to use
    col = get_cols(fmt[0])

    ares_min = 0.2   # Minimum ratio of int / peak
    ares_max = 1.2   # Maximum ratio of int / peak
    snr_min = 20.0   # Minimum SNR
    rad_match = 10.0 # Maximum radius for matching (arcsec)
    rad_isol = 30.0  # Search radius to ensure isolated source (arcsec)
    bref = 20        # Reference beams against which to compare offsets

    srcs = []
    fluxes = []
    pkfluxes = []
    im_rms = []
    fit_rms = []
    chi2 = []
    maj = []
    maj_deconv = []
    # get the combined list of coordinates, fluxes, peaks, rms and major axis
    for incat in flist:             # for each beam
        bpos = incat.find("beam")
        beam = int(incat[bpos+4:bpos+6])
        print(f"Beam {beam}: ", end="")
        cat = at.Table.read(incat)  # read the catalogue
        
        comp_sc = SkyCoord(Angle(cat["col_ra_deg_cont"], unit=u.deg), Angle(cat["col_dec_deg_cont"], unit=u.deg))

        # Find isolated sources
        c1, d2, id_match = match_self(comp_sc, rad_isol)
        comp_sc = comp_sc[id_match]
        cat = cat[id_match]

        # Filter out weak sources and sources that may be resolved
        filt_snr = (cat["col_flux_int"] / cat["col_rms_image"] >= snr_min)
        filt_r_lo = (cat["col_flux_int"] / cat["col_flux_peak"] >= ares_min)
        filt_r_hi = (cat["col_flux_int"] / cat["col_flux_peak"] <= ares_max)
        cat = cat[filt_snr & filt_r_lo & filt_r_hi]
#         comp_sc = comp_sc[filt_snr & filt_r_lo & filt_r_hi]
        
        print("%d sources" %(len(cat)))
        ra = cat[col['ra']]         
        dec = cat[col['dec']]
        srcs.append(SkyCoord(ra, dec, unit='deg'))
        fluxes.append(cat[col['flux']])
        pkfluxes.append(cat[col['pkflux']])
        im_rms.append(cat[col['rms']])
        # fit_rms.append(cat[col['fitrms']])
        # chi2.append(cat[col['chi2']])
        maj.append(cat[col['maj']])
        # maj_deconv.append(cat[col['maj_deconv']])

    # srcs, fluxes, pkfluxes,im_rms and maj all have 36 (beams) x nsrcs
#    print("len(srcs)=",len(srcs), "len(srcs[0])=",len(srcs[0]))
    # At this point have srcs, fluxes, pkfluxes and im_rms
    
    # get the SNR for each source in each beam
    snr = [np.array(f)/np.array(i) for f,i in zip(fluxes,im_rms)]

    # Construct disctionaries keyed by the beam pair
    # Each has lists of indices, separations, position angles, fluxes, SNRs, major axes, and deconvolved major axes
    # do both upper and lower triangles
    max_sep = 10.0 * u.arcsec
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
                # print(key)
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
                majdict[key] = [maj[j][sel],maj[i][idx[sel]]]
                # majdcdict[key] = [maj_deconv[j][sel],maj_deconv[i][idx[sel]]]
                sources[key] = [sj,si]
                # print(key, len(indices[key]))
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
#                print("%d-%d %5.1f %5.1f" %(i, b, dxs[i,b], dys[i,b]))
    pfile = f'medians/beamwise_medians_SB{sbid}.pkl'
    with open(pfile, 'wb') as f:
        pickle.dump([dxs,dys,dxs_std,dys_std,ns], f)

#    for beam in range(36):
#        if beam == bref:
#            continue
#        print("%d-20: %.3f,%.3f,%d" %(beam,dxs[beam,bref],dys[beam,bref],ns[beam,bref]))
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
#        p = ds[bref - sl.start,:] # not sure why Dave had this
        p = ds[bref,:]          # Initial guess for offsets (assume beam bref has 0 offset)
        ds[np.isnan(ds)] = 0.   # For offsets that have nans set them to 0
        popt = minimize(ofunc, p, args=(ds, wt), method='Nelder-Mead')
        # print(popt.x)
        fitted_rel.append(pfunc(popt.x, ds)[1])

    fitted_rel = np.array(fitted_rel)

    for beam in range(36):
        print("%2d: %6.3f,%6.3f" %(beam, fitted_rel[0][beam], fitted_rel[1][beam]))

    pfile1 = f'pickle/beamwise_fitted_rel_SB{sbid}_{bref}.pkl'
    with open(pfile1, 'wb') as f:
        pickle.dump([bref, fitted_rel], f)
        
    return len(flist), bref, fitted_rel

warnings.filterwarnings("ignore")
sbid = int(sys.argv[1])
ext = '.xml'
if len(sys.argv)==3:
    ext = sys.argv[2]   # either ".xml" (if starting with original catalogue) or "_fit.xml" (if applying to already corrected catalogue)

flist = glob.glob("cat/*SB%d*.components%s" %(sbid, ext))
if len(flist) != 36:
    sys.exit("SB%d missing beams (found %d)" %(sbid, len(flist)))
flist.sort()
fname = flist[0]
field_name = fname.split(".")[2]
#field_name = "RACS"
nf, bref, fitted_rel = process_sbid(flist, sbid, field_name)
