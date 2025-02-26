#!/usr/bin/env python

import sys
from argparse import ArgumentParser
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from typing import NewType
import glob
import os

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord, SkyOffsetFrame, concatenate, match_coordinates_sky, search_around_sky, Angle
from astropy.coordinates import match_coordinates_sky
from astropy.table import Table

Paths = tuple[Path, ...]
MatchMatrix = NewType("MatchMatrix", np.ndarray)

@dataclass
class Offset:
    """Contains offsets in the RA and Dec directions in arcsec"""
    ra: float = 0.0
    """Offset in RA direction"""
    dec: float = 0.0
    """Offset in Dec direction"""

@dataclass
class Catalogue:
    """Represent a per-beam ASKAP component catalogue"""
    sbid: int
    field_name: str
    beam: int
    cat_type: str
    """The type of catalogue"""
    table: Table
    """The table loaded"""
    path: Path
    """Original path to the loaded catalogue"""
    centre: SkyCoord
    """Beam centre derived from coordinates of componetns in catalogue"""
    fixed: bool = False
    """Indicates whether beam has been fixed into a place"""
    offset: Offset = field(default_factory=Offset)
    """Per beam offsets, if known, in arcsec"""
    
    def __repr__(self) -> str:
        return f"Catalogue(beam={self.beam}, sbid={self.sbid}, field_name={self.field_name} cat_type={self.cat_type} table={len(self.table)} sources, path={self.path}, fixed={self.fixed})"

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

# Used to determine true beam centres for RACS-low3, RACS-mid1, RACS-high1 and RACS-mid2 observations
field_beams = Table.read("closepack36_beams.fits")
field_centres = Table.read("closepack36_fields.fits")

@dataclass
class Match:
    """Components around matching Catalogue 1 to Catalogue 2"""
    sky_pos_1: SkyCoord
    """Sky positions from catalogue 1"""
    sky_pos_2: SkyCoord
    """Sky positions from catalogue 2"""
    matches: tuple
    """The indicies of the matches"""
    match_1: SkyCoord
    """The sky-coordinate of a match in catalogue 1"""
    match_2: SkyCoord
    """The sky-coordinate of a match in catalogue 2"""
    n: int
    """Number of matches"""
    offset_mean: tuple[float,float]
    """Mean of the offset in arcseconds in the RA and Declination directions"""
    offset_std: tuple[float,float]
    """Std of the offset in arcseconds in the RA and Declination directions"""
    err_ra: np.ndarray
    """Difference in RA coordinates between matches"""
    err_dec: np.ndarray
    """Different in Dec coordinates between matches"""

@dataclass
class BeamPair:
    """Represents a stage in the alignment process"""
    fixed_beam_idx: int
    """The idx of the catalogue that will not change"""
    shift_beam_idx: int
    """The idx of the catalogue that will be shifted"""
    matches: Match
    """The result of the cross match"""

@dataclass
class StepInfo:
    """Statistics around the step in the alignmnet process"""
    accumulated_seps: float
    """The total separation amoung matched sources"""
    number_of_matches: int
    """The total number of matches"""


def calculate_matches(catalogue_1: Catalogue, catalogue_2: Catalogue, sep_limit_arcsecond: float=9) -> Match:
    """Match a pair of catalogues to identify the sources in common. 

    Args:
        catalogue_1 (Catalogue): The first loaded catalogue
        catalogue_2 (Catalogue): The second loaded catalogue
        sep_limit_arcsecond (float, optional): The separation limit for a match, in arcseconds. Defaults to None.

    Returns:
        Match: The result of the matching
    """

    sky_pos_1 = make_sky_coords(catalogue_1)
    sky_pos_2 = make_sky_coords(catalogue_2)

    matches = search_around_sky(
        sky_pos_1, sky_pos_2, seplimit=sep_limit_arcsecond*u.arcsec
    )
    match_1 = sky_pos_1[matches[0]]
    match_2 = sky_pos_2[matches[1]]

    # Extract the offsets of positions as angular offsets of the sphere
    deltas = match_1.spherical_offsets_to(match_2)
    err_ra = deltas[0].to(u.arcsec).value
    err_dec = deltas[1].to(u.arcsec).value
    
    mean_ra, mean_dec = np.mean(err_ra), np.mean(err_dec)
    std_ra, std_dec = np.std(err_ra), np.std(err_dec)
    
    return Match(
        sky_pos_1=sky_pos_1, 
        sky_pos_2=sky_pos_2, 
        matches=matches, 
        match_1=match_1, 
        match_2=match_2, 
        n=len(match_2), 
        offset_mean = (mean_ra, mean_dec),
        offset_std=(std_ra, std_dec),
        err_ra=err_ra,
        err_dec=err_dec
    )

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
  

def make_sky_coords(table: Table | Catalogue) -> SkyCoord:
    """Create the sky-coordinates from a catalogue table

    Args:
        table (Table | Catalogue): Loaded table or catalogue

    Returns:
        SkyCoord: Sky-positions loaded
    """
    table = table.table if isinstance(table, Catalogue) else table
    sky_pos = SkyCoord(table["ra"], table["dec"], unit=(u.deg, u.deg))
    return sky_pos

def make_catalogue_matrix(catalogues: Catalogues) -> MatchMatrix:
    """Match each catalogue to each other

    Args:
        catalogues (Catalogues): Collection of beamwise component catalogues

    Returns:
        MatchMatrix: Matrix of matches
    """
    no_catas = len(catalogues)
    match_matrix = np.zeros((no_catas, no_catas))

    combos = list(combinations(list(range(len(catalogues))), 2))
    
    print("Generating sky-posiitons in matrix")
    sky_positions = [make_sky_coords(table=catalogue) for catalogue in catalogues]

    for (b1, b2) in combos:
#        print(f"Matching {b1} to {b2}")
        sky_pos_1, sky_pos_2 = sky_positions[b1], sky_positions[b2]
        
        match_results = search_around_sky(
            sky_pos_1, sky_pos_2, seplimit=9*u.arcsec
        )
        match_matrix[b1,b2] = len(match_results[0])

    print(f"Have matched {len(combos)}")
    return match_matrix

def plot_match_matrix(matrix: MatchMatrix, output_path: None |Path = None) -> Path:
    """Plot the match matrix from the beam-wise matching

    Args:
        matrix (MatchMatrix): The beame to beam number of matches
        output_path (None | Path, optional): Location to write image to. If None 'match_matrix.png' is used. Defaults to None.

    Returns:
        Path: Path of new plot 
    """
    print("Plotting match matrix")
    output_path = Path("match_matrix.png") if output_path is None else output_path
    fig, (ax1, ax2) = plt.subplots(1,2)
    
    cim = ax1.imshow(matrix)
    fig.colorbar(cim, label="N")
    ax1.set(
        xlabel="Beam no.",
        ylabel="Beam no.",
        title="Source matches"
    )
    
    ax2.hist(matrix[matrix > 0].flatten(), bins=20)
    ax2.set(
        xlabel="Number of matches",
        ylabel="Count"
    )
    
    fig.tight_layout()
    fig.savefig(fname=output_path)
    print(f"Have created {output_path=}")
    return output_path

def make_and_plot_match_matrix(catalogues: Catalogues, plot_path: None | Path = None) -> tuple[np.ndarray, Path]:
    """Run the making and plotting of the match matrix"""

    matrix = make_catalogue_matrix(catalogues=catalogues)
#    plot_path = plot_match_matrix(matrix=matrix, output_path=plot_path)
    
    return matrix, plot_path


def set_seed_catalogues(catalogues: Catalogues, match_matrix: MatchMatrix) -> Catalogues:
    """Select a beam to fix into place so others are matched to it. 
    This is done by identifying the beam with the most matches to 
    other beams. 

    Args:
        catalogues (Catalogues): The collection of beam catalogues to consider
        match_matrix (MatchMatrix): The beam-to-beam sky-match result set

    Returns:
        Catalogues: The same as the input catalouges, with the exception of a fixed beam
    """
    sum_matrix = np.sum(match_matrix, axis=0)
    idx = np.argmax(sum_matrix)
    
    catalogues[idx].fixed = True

    assert len([catalogue.fixed for catalogue in catalogues if catalogue.fixed]) == 1, "Too many seeds"
    
    return catalogues

def find_next_pair(catalogues: Catalogues) -> BeamPair | None:
    """Identify a pair of beams that will form a step in the deshifter

    Args:
        catalogues (Catalogues): Collection of beam cataloues to consider

    Returns:
        BeamPair | None: The pair of beam catalogues for this step. If there are no catalogues to shift None is returned.
    """

    assert any([catalogue.fixed for catalogue in catalogues]), "There are no fixed catalogues"

    # Split the catalogues into groups of deshifted and ones to shift
    fixed_beam_idxs = [idx for idx, cata in enumerate(catalogues) if cata.fixed]
    candidate_beam_idxs =  [idx for idx, cata in enumerate(catalogues) if not cata.fixed] 

    if len(candidate_beam_idxs) == 0:
        return None
    
    ideal_fixed_beam_idx = None
    ideal_shift_beam_idx = None
    current_best_match = None

    for fixed_beam_idx in fixed_beam_idxs:
        fixed_beam_cata = catalogues[fixed_beam_idx]

        for candidate_beam_idx in candidate_beam_idxs:
            candidate_beam_cata = catalogues[candidate_beam_idx]
            matches = calculate_matches(
                catalogue_1=fixed_beam_cata, catalogue_2=candidate_beam_cata
            )
            
            if current_best_match is None or matches.n > current_best_match.n:
                current_best_match = matches
                ideal_fixed_beam_idx = fixed_beam_idx
                ideal_shift_beam_idx = candidate_beam_idx
#                print(f"Update {ideal_fixed_beam_idx=} {ideal_shift_beam_idx=}")
                
    return BeamPair(
        fixed_beam_idx=ideal_fixed_beam_idx, 
        shift_beam_idx=ideal_shift_beam_idx, 
        matches=current_best_match
    )

def _select_random_index(max_index: int) -> int:
    from random import randint
    return randint(a=0, b=max_index-1)

def add_offset_to_coords_skyframeoffset(
    sky_coords: SkyCoord, offset: tuple[float, float]
) -> SkyCoord:
    """Add offsets to sky coordinate offsets. This attempts to 
    be consistent with the `spherical_offsets_to` astropy function
    and adds the angular offsets appropriately on the sphere.

    Args:
        sky_coords (SkyCoord): The base set of coordinates to shift
        offset (tuple[float, float]): The angular offsets from `spherical_offsets_to`

    Returns:
        SkyCoord: The shifted sky positions
    """
    # NOTE: Spheres are hard and confuse me. I am not particularly convinced
    # that simply subtracting/adding d(RA) and d(Dec) from sets is correct.
    # I trust the astropy more than I. This function attempts to do the 
    # reverse of the `spherical_offsets_to` method. 
    
    # The shift needs to be an array of same shape
    d_ra = (np.zeros_like(sky_coords) - offset[0])*u.arcsec
    d_dec = (np.zeros_like(sky_coords) - offset[1])*u.arcsec
    
    new_coords = sky_coords.spherical_offsets_by(
        d_ra, d_dec
    )
    return new_coords

def add_offset_to_catalogue(
    catalogue: Catalogue, offset: tuple[float, float]
) -> Catalogue:
    """Add offsets to a catalogue and its table. 
    Args:
        catalogue (Catalogue): The catalogue object to shift
        offset (tuple[float, float]): The angular units to shift by

    Returns:
        Catalogue: The shifted catalogue
    """
    cata_table = catalogue.table.copy()
    
    sky_coords = make_sky_coords(table=cata_table)
    new_coords = add_offset_to_coords_skyframeoffset(sky_coords, offset)
    
    cata_table["ra"] = new_coords.ra.deg
    cata_table["dec"] = new_coords.dec.deg

    # Trust no one
    from copy import deepcopy
    new_cata = deepcopy(catalogue)

    new_cata.table = cata_table
    
    # Update the offset in case many passes over
    summed_offsets = Offset(
        ra=offset[0] + new_cata.offset.ra,
        dec=offset[1] + new_cata.offset.dec,
    )
    new_cata.offset = summed_offsets

    return new_cata


def calculate_catalogue_jitter(catalogues: Catalogues, sep_limit_arcsecond: float=9) -> StepInfo:
    """Calculate global statistics of all matches across all catalogues

    Args:
        catalogues (Catalogues): The set of catalogues to consider
        sep_limit_arcsecond (float, optional): The separation limit to condider, in arcseconds. Defaults to 9.

    Returns:
        StepInfo: Information of separations across all matches
    """
    # TODO: The is similar to the make match matrix function, but
    # extra statistics are accumulated. Potentially they can be merged
    
    num_catalogues = len(catalogues)
    seps = 0
    no_matches = 0
    combos = list(
        combinations(
            list(range(num_catalogues)), 
            2
        )
    )

    for (b1, b2) in combos:
        cata_1, cata_2 = catalogues[b1], catalogues[b2]
        
        sky_pos_1, sky_pos_2 = make_sky_coords(cata_1), make_sky_coords(cata_2)
        
        match_results = search_around_sky(
            sky_pos_1, sky_pos_2, seplimit=sep_limit_arcsecond*u.arcsec
        )
        seps += np.sum(match_results[2])
        no_matches += len(match_results[2])
    
    return StepInfo(accumulated_seps=seps, number_of_matches=no_matches)

def round_header(step: int, stats: StepInfo) -> None:
    print(f"Round {step}, {stats}")


def plot_iterative_shift_stats(step_statistics: list[StepInfo], output_path: Path | None = None) -> Path:
    """Plot the progression of matching statistics over rounds. The list of 
    ``StepInfo`` is assumed ot be in orfer.

    Args:
        step_statistics (list[StepInfo]): Colleciton of statistics gathered from the iterative convergence
        output_path (Path | None, optional): The output path. If None `stats_step_info.png` will be used. Defaults to None.

    Returns:
        Path: The output image path
    """

    output_path = output_path if output_path else Path("stats_step_info.png")

    seps = np.array([s.accumulated_seps.value for s in step_statistics])
    srcs = np.array([s.number_of_matches for s in step_statistics])
    
    fig, (ax, ax2) = plt.subplots(1,2)
    
    ax.plot(seps)
    ax.axvline(36, ls="--", label="New Seed Beam")
    ax.set(
        xlabel="Step",
        ylabel="Accumulated Separations (deg)"
    )
    ax.legend()
    ax.grid()
    
    ax2.plot(seps/srcs*3600)
    ax2.set(
        xlabel="Step",
        ylabel="Average Separation (arcsec)"
    )
    ax2.grid()
    
    fig.tight_layout()
    fig.savefig(fname=output_path)
    
    return output_path

def reseed_initial_fixed_catalogue(catalogues: Catalogues) -> Catalogues:
    for cata in catalogues:
            cata.fixed = False
    catalogues[
        _select_random_index(max_index=len(catalogues))
        ].fixed = True

    return catalogues    

def perform_iterative_shifter(
    catalogues: Catalogues, 
    passes: int = 1, 
    gather_statistics: bool = True,
    output_prefix: str | None = None
) -> Catalogues:
    """Attempt to shift catalogues to a common reference frame. A seed catalogue
    is selected, then a catalogue at a time is selected and aligned. This may be
    repeated a number of times. 

    Args:
        catalogues (Catalogues): Catalogues that should be aligned
        passes (int, optional): How many passes over all catalogues should be performed. Defaults to 1.
        gather_statistics (bool, optional): Whether statistcs across the convergence should be collected. This can be time consuming as all catalogues are matched to one another. Defaults to True.
        output_prefix (str | None, optional): The prefix to attach to output products. Defaults to None.

    Returns:
        Catalogues: The shifted catalogues
    """
    
    print(f"Shifting {len(catalogues)}")
    step_statistics = []
    for step in range(len(catalogues) * passes):
        pair_match = find_next_pair(catalogues)

        # This is triggered if everything is already matched
        # should some iterative procedure be invoked 
        if pair_match is None:
            catalogues = reseed_initial_fixed_catalogue(catalogues=catalogues)
            continue
            
        new_catalogue = add_offset_to_catalogue(
                catalogue=catalogues[pair_match.shift_beam_idx], 
                offset=pair_match.matches.offset_mean
            )
        new_catalogue.fixed = True

        catalogues[pair_match.shift_beam_idx] = new_catalogue
        
        if gather_statistics:
            total_seps = calculate_catalogue_jitter(catalogues=catalogues)
            step_statistics.append(total_seps)
            round_header(step=step, stats=total_seps)
        else:
            print(f"Shifted in round {step}")
        
    return catalogues

def save_catalogue_shift_positions(ref_beam, catalogues: Catalogues, output_path: Path | None = None) -> Path:
    dxs = np.zeros((36), dtype=np.float64)
    dys = np.zeros((36), dtype=np.float64)
    
    output_path = output_path if output_path else Path("shifts.csv")
    for catalogue in catalogues:
        dxs[catalogue.beam] = catalogue.offset.ra
        dys[catalogue.beam] = catalogue.offset.dec
    dxs -= dxs[ref_beam]
    dys -= dys[ref_beam]
    fout = open(output_path, "wt")
    fout.write("BEAM,DXS,DYS\n")
    for beam in range(36):
        fout.write("%d,%f,%f\n" %(beam, dxs[beam], dys[beam]))
    fout.close()
    
    return output_path

def beam_wise_shifts(
    catalogue_paths: Paths,
    ref_beam: int
) -> Catalogues:
    """Load in a set of catalogues and attempt to align them
    onto an internally consisten positional reference frame

    Args:
        catalogue_paths (Paths): The set of fits component cataloges to load
        output_prefix (str | None, optional): The prefix to use for output products. If None the default names are used. Defaults to None.

    Returns:
        Catalogues: The catalogues that have been shifted
    """
    
    print(f"Will be processing {len(catalogue_paths)} catalogues")
    catalogues: Catalogues = load_catalogues(catalogue_paths=catalogue_paths)
    field_name = catalogues[0].field_name
    sbid = catalogues[0].sbid
    if catalogues[0].cat_type in ["selavy"]:
        output_prefix = f"SB{sbid}.{field_name}.selavy."
    else:
        output_prefix = f"SB{sbid}.{field_name}."
    
    match_matrix: MatchMatrix
#    match_matrix_plot = Path(output_prefix+'match_matrix.png') if output_prefix else Path("match_matrix.png")
#    match_matrix, _ = make_and_plot_match_matrix(catalogues=catalogues, plot_path=match_matrix_plot)
    match_matrix = make_catalogue_matrix(catalogues=catalogues)


    catalogues = set_seed_catalogues(catalogues=catalogues, match_matrix=match_matrix)    
    catalogues = perform_iterative_shifter(
        catalogues=catalogues, 
        passes=1, 
        gather_statistics=True,
        output_prefix=output_prefix
    )

    shift_path = Path(output_prefix + "shifts.csv") if output_prefix else Path("shifts.csv")
    save_catalogue_shift_positions(ref_beam=ref_beam, catalogues=catalogues, output_path=shift_path)

    return catalogues

flist = sys.argv[1:]
flist.sort()
cats = beam_wise_shifts(catalogue_paths=flist, ref_beam=20)
#os.system("./plot_fit.py SB57171.RACS_1110-51.shifts.csv")

#flist = glob.glob("catalogues/*SB57171*.xml")
#flist.sort()
#beam_wise_shifts(catalogue_paths=flist, ref_beam=20)
#os.system("./plot_fit.py SB57171.RACS_1110-51.selavy.shifts.csv")
