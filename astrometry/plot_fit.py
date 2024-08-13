#!/usr/bin/env python

import os
import numpy as np
import sys
import pickle
from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky
import glob
import astropy.coordinates as ac
import astropy.table as at
import matplotlib.pyplot as plt
from astropy.table import Table

db_base_path = "/Users/len067/Desktop/aces/calibration/askap_surveys"
survey = "RACS"
epoch = 9 # RACS-low3
db_path = "%s/%s/db/epoch_%d/field_data.csv" %(db_base_path, survey, epoch)
if os.path.exists(db_path) == False:
    sys.exit("Unable to load database")
field_data = Table.read(db_path, format='csv')

ext = '.xml'
sbid = int(sys.argv[1])
if len(sys.argv) == 3:
    ext = sys.argv[2]   # either ".xml" (if starting with original catalogue) or "_fit.xml" (if applying to already corrected catalogue)

field_data = field_data[np.where(field_data["SBID"]==sbid)]
if len(field_data) == 0:
    sys.exit("unable to find sbid in database")
field_data = field_data[0]
cal_sbid = field_data["CAL_SBID"]
field_name = field_data["FIELD_NAME"]

bref = 20

# Read the fitted offsets
[bref, fitted_rel] = np.load(f'./pickle/beamwise_fitted_rel_SB{sbid}_{bref}.pkl', allow_pickle=True, encoding='bytes')
# bref = reference beam
# fitted_rel[0] = RA offsets for each beam (arcsec)
# fitted_rel[1] = dec offsets for each beam (arcsec)

boffsets = [[-2.75, -2.16506], [-1.75, -2.16506], [-0.75, -2.16506], [0.25, -2.16506], [1.25, -2.16506], [2.25, -2.16506],
    [-2.25, -1.29904], [-1.25, -1.29904], [-0.25, -1.29904], [0.75, -1.29904], [1.75, -1.29904], [2.75, -1.29904], 
    [-2.75, -0.433013], [-1.75, -0.433013], [-0.75, -0.433013], [0.25, -0.433013], [1.25, -0.433013], [2.25, -0.433013], 
    [-2.25, 0.433013], [-1.25, 0.433013], [-0.25, 0.433013], [0.75, 0.433013], [1.75, 0.433013], [2.75, 0.433013], 
    [-2.75, 1.29904], [-1.75, 1.29904], [-0.75, 1.29904], [0.25, 1.29904], [1.25, 1.29904], [2.25, 1.29904], 
    [-2.25, 2.16506], [-1.25, 2.16506], [-0.25, 2.16506], [0.75, 2.16506], [1.75, 2.16506], [2.75, 2.16506]]
boffsets = np.array(boffsets)
fig, axes = plt.subplots(1,1, figsize=(8,8))
kwargs1 = {'angles':'uv', 'scale_units':'xy', 'pivot':'mid', 'scale':2, 'headwidth':4, 'headlength':3, 'headaxislength':2, 'width':0.01}

axes.quiver(-boffsets[:,0], boffsets[:,1], fitted_rel[0], fitted_rel[1], **kwargs1)
for beam in range(36):
    print("Beam %02d : dx=%6.3f, dy=%6.3f" %(beam, fitted_rel[0][beam], fitted_rel[1][beam]))
#    print(fitted_rel[0], fitted_rel[1])
axes.quiver(-boffsets[:,0], boffsets[:,1], fitted_rel[0], fitted_rel[1], **kwargs1)

# plot a single red arrow for scale
scale_pos = [2.,-2.5]
scale_len = 1.0
dx, dy = 0.4, 0.2
box_x = [scale_pos[0] - dx, scale_pos[0] + dx, scale_pos[0] + dx, scale_pos[0] - dx, scale_pos[0] - dx]
box_y = [scale_pos[1] - dy, scale_pos[1] - dy, scale_pos[1] + dy, scale_pos[1] + dy, scale_pos[1] - dy]
axes.quiver(scale_pos[0],scale_pos[1],scale_len,0, color='r', **kwargs1)
axes.plot(box_x, box_y, 'r')
axes.text(scale_pos[0] + 0.5, scale_pos[1], f'{scale_len} arcsec', color='r', fontsize=8)
axes.plot()
axes.plot(-boffsets[:,0], boffsets[:,1], 'og', ms=15, alpha = 0.3)
axes.set_xlim(-3.2,3.2)
axes.set_ylim(-2.8,2.8)
axes.set_aspect('equal')
axes.set_xlabel('Offset (degrees)')
axes.set_ylabel('Offset (degrees)')
axes.set_title(f'Cal SBID: {cal_sbid} SBID: {sbid} ({field_name})   offsets')
plt.savefig("./png/SB%d_%s.png" %(sbid, field_name))
