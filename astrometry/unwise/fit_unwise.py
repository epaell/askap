#!/usr/bin/env python

from time import time
from collections import defaultdict
from pathlib import Path
from dataclasses import dataclass
from itertools import combinations
from typing import Tuple
import os
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame, concatenate, match_coordinates_sky, search_around_sky, Angle
from astropy.coordinates import match_coordinates_sky
from astropy.table import Table, unique, vstack
import matplotlib.pyplot as plt
from astroquery import vizier
import sys

Paths = tuple[Path, ...]

# Extracts per-beam catalogues and merges them into a single field catalogue containing all unwise sources.
def load_unwise(field_name):
    field_cat = None
    if os.path.exists(f"../unwise/unwise_{field_name}.fits") == True:
        field_cat = Table.read(f"../unwise/unwise_{field_name}.fits", format="fits")
        print("%d rows loaded" %(len(field_cat)))
        return field_cat
    return None 

@dataclass
class Catalogue:
    sbid: int
    field_name: str
    beam: int
    cat_type: str
    table: Table
    path: Path
    centre: SkyCoord
    offset = (0,0)
    
    def __repr__(self):
        return f"Catalogue(beam={self.beam}, sbid={self.sbid}, field_name={self.field_name} cat_type={self.cat_type} table={len(self.table)} sources, path={self.path})"

    def filter(self):
        """Filter radio components out of an aegean radio catalogue
        based on their distance to neighbouring components and compactness. 

        Args:
            table (Table): Aegean radio component catalogue

        Returns:
            np.ndarray: Boolean array of components to keep. 
        """
        sky_coord = SkyCoord(self.table["ra"], self.table["dec"], unit=(u.deg, u.deg), frame='fk5')
    
        isolation_mask = sky_coord.match_to_catalog_sky(sky_coord, nthneighbor=2)[1] > (0.01 * u.deg)
    
        ratio = self.table["int_flux"] / self.table["peak_flux"]
        ratio_mask = (0.8 < ratio) & (ratio < 1.2)

        table_mask = isolation_mask & ratio_mask
        self.table = self.table[table_mask]

Catalogues = tuple[Catalogue, ...]

def estimate_skycoord_centre(sky_positions: SkyCoord, final_frame: str = "fk5") -> SkyCoord:
    
    xyz_positions = sky_positions.cartesian.xyz
    xyz_mean_position = np.mean(xyz_positions, axis=1)

    mean_position = SkyCoord(*xyz_mean_position, representation_type="cartesian").transform_to(final_frame)

    return mean_position

def filter_table(table: Table):
    sky_coord = SkyCoord(table["ra"], table["rec"], unit=(u.deg, u.deg))
    
    isolation_mask = sky_coord.match_to_catalog_sky(sky_coord, nthneighbor=2)[1] > (0.01 * u.deg)
    
    ratio = table["int_flux"] / table["peak_flux"]
    ratio_mask = (0.8 < ratio) & (ratio < 1.2)

    return isolation_mask & ratio_mask
                         
def load_catalogue(catalogue_path: Path) -> Catalogue:
    """Load a beam catalogue astropy table

    Args:
        catalogue_path (Path): Path to load catalogue from

    Returns:
        Catalogue: Loaded catalogue
    """
    path_name = str(catalogue_path.name) if isinstance(catalogue_path, Path) else catalogue_path
    cat_file = path_name.split("/")[-1]
    # Work out if this is a selavy catalogue or something else (and so which columns to use)
    if cat_file.find("selavy") != -1:
        cat_type = "selavy"
    else:
        cat_type = "default"

    # Get the beam number for the catalogue
    beam_pos = cat_file.find("beam")
    assert(beam_pos != -1)
    beam = int(cat_file[beam_pos+4:beam_pos+6])
    
    # Get the field name
    field_name_pos = cat_file.find("RACS_")
    assert(field_name_pos != -1)
    field_name = cat_file[field_name_pos:field_name_pos+12]
    
    # Get the SBID
    sbid_pos = cat_file.find("SB")
    assert(sbid_pos != -1)
    cdata = cat_file[sbid_pos+2:].split(".")
    sbid = int(cdata[0])

    bdata = field_beams[np.where(field_beams["FIELD_NAME"]==field_name)]
    centre = SkyCoord(Angle(bdata["RA_DEG"][beam], u.deg), Angle(bdata["DEC_DEG"][beam], u.deg), frame='fk5')
        
    print(f"Loading {catalogue_path}")
    table = Table.read(catalogue_path)
    if cat_type in ["selavy"]:
        table.rename_column('col_ra_deg_cont', "ra")
        table.rename_column('col_dec_deg_cont', "dec")
        table.rename_column('col_flux_int', "int_flux")
        table.rename_column('col_flux_peak', "peak_flux")

    cat = Catalogue(field_name=field_name, sbid=sbid, beam=beam, cat_type=cat_type, table=table, path=catalogue_path, centre=centre)
    cat.filter()
    return cat
    
    
def load_catalogues(catalogue_paths: Paths) -> Catalogues:
    """Load in all of the catalgues"""    
    return [load_catalogue(catalogue_path=catalogue_path) for catalogue_path in catalogue_paths]
  

def make_sky_coords(table):
    if isinstance(table, Catalogue):
        table = table.table
    sky_pos = SkyCoord(table["ra"], table["dec"], unit=(u.deg, u.deg))
    return sky_pos

def make_skycoord(table):
    return SkyCoord(table["ra"], table["dec"], unit=(u.deg, u.deg), frame='fk5')

def add_offset_to_coords(sky_coords, offset):
    new_coords = SkyCoord([sc.spherical_offsets_by(-offset[0]*u.arcsec, -offset[1]*u.arcsec) for sc in sky_coords])

    return new_coords

def add_offset_to_coords_skyframeoffset(sky_coords, offset):
    d_ra = (np.zeros_like(sky_coords) -offset[0])*u.arcsec
    d_dec = (np.zeros_like(sky_coords) -offset[1])*u.arcsec
    
    new_coords = sky_coords.spherical_offsets_by(d_ra, d_dec)
    return new_coords

def add_offset_to_cata(cata, offset):
    cata_table = cata.table.copy()
    
    sky_coords = make_skycoord(table=cata_table)
    new_coords = add_offset_to_coords_skyframeoffset(sky_coords, offset)
    
    cata_table["ra"] = new_coords.ra.deg
    cata_table["dec"] = new_coords.dec.deg

    from copy import deepcopy
    
    new_cata = deepcopy(cata)

    new_cata.table = cata_table
    new_cata.offset = offset
    return new_cata

def add_offset_to_coords_skyframeoffset2(sky_coords, offset):
    d_ra = -offset[0]*u.arcsec
    d_dec = -offset[1]*u.arcsec
    
    new_coords = sky_coords.spherical_offsets_by(d_ra, d_dec)
    return new_coords

@dataclass
class OffsetGridSpace:
    ra_offsets: float
    dec_offsets: float
    beam: int
    seps: np.ndarray
    n_sources: int

def find_minimum_offset_space(offset_space):
    minimum_sep = None
    minimum_ra = None
    minimum_dec = None
    for (dec, ra, sep) in zip(offset_space.dec_offsets, offset_space.ra_offsets, offset_space.seps.flatten()):
        if minimum_sep is None or minimum_sep > sep:
            minimum_sep = sep
            minimum_ra = ra
            minimum_dec = dec

    return minimum_ra, minimum_dec, minimum_sep

def get_offset_space(cata, unwise_table, window, plot=True, radius_deg: float=5.):
    unwise_sky = SkyCoord(unwise_table["RAJ2000"], unwise_table["DEJ2000"], unit=(u.deg, u.deg))
    
    shifted_table = cata.table.copy()
    cata_sky = SkyCoord(shifted_table["ra"], shifted_table["dec"], unit=(u.deg, u.deg))

    # The range of RA and Dec searched
    ra_extent = np.abs(window[1]-window[0])
    dec_extent = np.abs(window[3]-window[2])
    # The number of bins in RA and Dec
    ra_bins = int(ra_extent / window[4])
    dec_bins = int(dec_extent / window[4])

    ras = np.linspace(window[0], window[1], ra_bins)
    decs = np.linspace(window[2], window[3], dec_bins)

    coords = []
    
    n_sources = len(cata_sky)
    n_delta = n_sources * len(decs) * len(ras)
    
    broadcast_d_ra = np.zeros(n_delta)
    broadcast_d_dec = np.zeros(n_delta)
    
    for d_idx, dec in enumerate(decs):
        for r_idx, ra in enumerate(ras):
            i = len(coords)
            j = i + 1
            broadcast_d_ra[i * n_sources:j * n_sources] = ra
            broadcast_d_dec[i * n_sources:j * n_sources] = dec
            
            coords.append(cata_sky)
    collected_coords = concatenate(coords)
    
    shifted_sky = add_offset_to_coords_skyframeoffset2(collected_coords, (broadcast_d_ra, broadcast_d_dec))
    
    matches = match_coordinates_sky(shifted_sky, unwise_sky, nthneighbor=1, storekdtree=True)
    
    accumulated_seps = np.zeros((dec_bins, ra_bins))
    
    results = {}
    minimum_key = None
    minimum_seps = None
    pair_decra = []
    
    for d_idx, dec in enumerate(decs):
        for r_idx, ra in enumerate(ras):
            start_idx = (d_idx*len(ras) + r_idx)
            end_idx = start_idx + 1
            k = (d_idx, r_idx)
            v = np.sum(matches[1][start_idx*n_sources:end_idx*n_sources].value) / n_sources
    
            seps = (dec, ra, v)
            pair_decra.append((dec, ra))
            results[k] = seps
            accumulated_seps[k] = v

    array_decra = np.array(pair_decra)
    
    return OffsetGridSpace(dec_offsets=array_decra[:,0], ra_offsets=array_decra[:,1], beam=cata.beam, n_sources=n_sources, seps=accumulated_seps)
        
def plot_offset_grid_space(fname, offset_grid_space, window):
    min_ra, min_dec, min_sep = find_minimum_offset_space(offset_grid_space)

    fig, ax = plt.subplots(1,1)

    cim = ax.imshow(offset_grid_space.seps, extent=(window[0], window[1], window[2], window[3]), origin="lower")
    
    ax.grid()
    ax.axhline(min_dec, ls="--", color="white")
    ax.axvline(min_ra, ls="--", color="white")
    
    ax.set(xlabel="Delta RA (arcsec)", ylabel="Delta Dec (arcsec)", title=f"Beam {offset_grid_space.beam} ({min_ra:.4f} {min_dec:.4f}) arcsec {offset_grid_space.n_sources} beam srcs")
    cbar = fig.colorbar(cim, label="Summed offsets (Degrees)")

    fig.tight_layout()
    plt.savefig(fname, dpi=150)
    plt.close()

field_beams = Table.read("closepack36_beams.fits")
field_centres = Table.read("closepack36_fields.fits")

flist = sys.argv[1:]
flist.sort()
catas = load_catalogues(flist)
field_name = catas[0].field_name
sbid = catas[0].sbid
if (catas[0].cat_type in ["selavy"])==True:
    fout = open(f"SB{sbid}.{field_name}.unwise.selavy.shifts.csv", "wt")
else:
    fout = open(f"SB{sbid}.{field_name}.unwise.shifts.csv", "wt")
fout.write("BEAM,DXS,DYS\n")

unique_fields = np.unique(field_beams["FIELD_NAME"])
beam_inf = field_beams[np.where(field_beams["FIELD_NAME"]==field_name)]
beam_sc = SkyCoord(Angle(beam_inf["RA_DEG"], unit=u.deg), Angle(beam_inf["DEC_DEG"], unit=u.deg), frame='fk5')

unwise_field_cat = load_unwise(field_name)
niter = 6
zoom = 4.0
for beam in range(36):
    width = 25.0
    delta = 5.0
    min_ra = min_dec = 0
    for iteration in range(niter):
        window = (min_ra-width, min_ra+width, min_dec-width, min_dec+width, delta)
#        print(window)
        offset_results = get_offset_space(catas[beam], unwise_field_cat, window=window)
        min_ra, min_dec, min_sep = find_minimum_offset_space(offset_results)
        if iteration==0:
            plot_offset_grid_space(f"SB{sbid}.{field_name}.beam{beam:02d}.offset.grid.png", offset_results, window=window)
        width /= zoom
        delta /= zoom
    print(f"Beam:{beam:02d} : {min_ra:6.3f},{min_dec:6.3f}")
    fout.write("%d,%f,%f\n" %(beam, min_ra, min_dec))
os.system(f"montage SB{sbid}*beam*.png -geometry +6+6 SB{sbid}.{field_name}.offset_grid.png")
os.system(f"rm SB{sbid}*beam*.png")
fout.close()
