#!/usr/bin/env python

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord, Angle
from astropy.table import Table
import sys

# SB57171.RACS_1110-51.shifts.csv
# BEAM,DXS,DYS
shift_file = sys.argv[1]
sdata = shift_file.split(".")
sbid = int(sdata[0][2:])
field_name = sdata[1]

output_path = shift_file.replace(".csv", ".png")

shifts = Table.read(shift_file, format='csv')

# Used to determine true beam centres for RACS-low3, RACS-mid1, RACS-high1 and RACS-mid2 observations
field_beams = Table.read("closepack36_beams.fits")
field_centres = Table.read("closepack36_fields.fits")

field = field_centres[np.where(field_centres["FIELD_NAME"]==field_name)][0]
field_sc = SkyCoord(Angle(field["RA_DEG"], u.deg), Angle(field["DEC_DEG"], u.deg), frame='fk5')

bdata = field_beams[np.where(field_beams["FIELD_NAME"]==field_name)]
beam_sc = SkyCoord(Angle(bdata["RA_DEG"], u.deg), Angle(bdata["DEC_DEG"], u.deg), frame='fk5')
boffsets = np.array(field_sc.spherical_offsets_to(beam_sc))

dxs = shifts["DXS"].value
dys = -shifts["DYS"].value

fig, axes = plt.subplots(1,2, figsize=(12,8))
kwargs1 = {'angles':'uv', 'scale_units':'xy', 'pivot':'mid', 'scale':2, 'headwidth':4, 'headlength':3, 'headaxislength':2, 'width':0.01}

axes[0].quiver(-boffsets[0,:], boffsets[1,:], dxs, dys, **kwargs1)
#axes[0].quiver(-boffsets[0,:], boffsets[1,:], dxs, dys, **kwargs1)

# plot a single red arrow for scale
scale_pos = [2.,-2.5]
scale_len = 1.0
dx, dy = 0.4, 0.2
box_x = [scale_pos[0] - dx, scale_pos[0] + dx, scale_pos[0] + dx, scale_pos[0] - dx, scale_pos[0] - dx]
box_y = [scale_pos[1] - dy, scale_pos[1] - dy, scale_pos[1] + dy, scale_pos[1] + dy, scale_pos[1] - dy]
axes[0].quiver(scale_pos[0],scale_pos[1],scale_len,0, color='r', **kwargs1)
axes[0].plot(box_x, box_y, 'r')
axes[0].text(scale_pos[0] + 0.5, scale_pos[1], f'{scale_len} arcsec', color='r', fontsize=8)
axes[0].plot()
axes[0].plot(-boffsets[0,:], boffsets[1,:], 'og', ms=15, alpha = 0.3)
axes[0].set_xlim(-4,4)
axes[0].set_ylim(-4,4)
axes[0].set_aspect('equal')
axes[0].set_xlabel('Offset (degrees)')
axes[0].set_ylabel('Offset (degrees)')

mdx = np.mean(dxs)
mdy = np.mean(dys)
ndxs = dxs - mdx
ndys = dys - mdy
axes[1].set_title(f'SBID: {sbid} ({field_name}) Mean subtracted ({mdx:.1f}",{mdy:.1f}")')
axes[1].plot(ndxs, ndys, marker="+", ls="None")
axes[1].set_xlim(-5,5)
axes[1].set_ylim(-5,5)
axes[1].set_aspect('equal')
axes[1].set_xlabel('Offset (arcsec)')
axes[1].set_ylabel('Offset (arcsec)')
plt.savefig(output_path)
