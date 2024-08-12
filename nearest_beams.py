#!/usr/bin/env python

import sys
import astropy.io.fits as fits
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u

if len(sys.argv) != 4:
    sys.exit("Usage:\n   %s ASKAP_fits_image RA Dec" %(sys.argv[0]))
    
im = sys.argv[1]
ra_str = sys.argv[2]
dec_str = sys.argv[3]

max_sep = 1.0 # Max separation from beam centre

# Check for HMSDMS format
if ra_str.find("h") != -1:
    src_dir = SkyCoord(Angle(ra_str, unit=u.hourangle), Angle(dec_str, unit=u.deg), frame='fk5')
# Check for ":" separators
elif ra_str.find(":") != -1:
    src_dir = SkyCoord(ra='{}h{}m{}s'.format(*ra_str.split(':')), dec='{}d{}m{}s'.format(*dec_str.split(':')), frame='fk5')
# Otherwise treat as degrees for RA and Dec
else:
    src_dir = SkyCoord(Angle(float(ra_str), unit=u.deg), Angle(float(dec_str), unit=u.deg), frame='fk5')

hdu = fits.open(im)
beams = hdu[1].data
bra = []
bdec = []
for beam in beams:
    bra.append(beam[2])
    bdec.append(beam[3])
beam_sc = SkyCoord(Angle(bra, u.deg), Angle(bdec, u.deg), frame='fk5')
sep = src_dir.separation(beam_sc).deg
print("Beams nearest to: %s" %(src_dir.to_string(style='hmsdms')))
for index in range(len(beams)):
    if sep[index] > max_sep:
        continue
    print("beam %02d : separation %.3f deg (beam %d position = %s)" %(index, sep[index], index, beam_sc[index].to_string(style='hmsdms')))
