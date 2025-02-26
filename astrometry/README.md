# astrometry
Tools to check for astrometric variations between per-beam ASKAP field catalogues.

# Tools

fit_astrometry.py sbid - examines beam-wise cataloues for the given SBID and generates offsets to align beams

plot_fit.py sbid - plot the offset corrections for the specified SBID

apply_fit.py sbid - takes offset corrections and applies to the beam-wise catalogues for the specified SBID
    
assess_fit.py sbid maxd maxr - assess the fit by comparing all beams seperated by no more than maxd degrees (1.0 to use only nearest beams) and using all sources no further than maxr degrees from the beam centre (0.9 deg works well). The script will generate a CSV-format line with the following information:
    
sbid,field_name,nso,xo_mean,xo_std,yo_mean,yo_std,nsf,xf_mean,xf_std,yf_mean,yf_std,flag

where the first set of values pertain to the original catalogues ("o") i.e. nso is the number of sources used in the assessment, xo_mean is the mean offset in RA, xo_std is the standard deviation of the offsets in RA, yo_mean is the mean offset in Dec, yo_std is the standard deviation of the offsets in declination. The second set of values pertain to the fitted catalogues ("f").

db.py sbid - print out some parameters extracted from the survey database for the specified sbid

run.py sbid - fits astrometry, applies it, plots it and assesses the result for a single SBID.

process_all.py - fits astrometry, applies it and plots the result for all sbids using multiprocessing.

plot_all.py - Reads the file check_all.csv (which has all a header and results from assess_fit.py appended to it) and plots a summary of the stats.

unwise/. - Contains tools that attempt to reference beam catalogues to UNWISE (based on Tims code).

# Setup

db_base_path should be set to point to the survey database in assess_fit.py, db.py, plot_all.py, plot_fit.py and process_all.py

max_concurrent in process_all.py should be set for the number of cores available to maximise processing speed

Directory structure:

./cat - location where beamwise catalogues are stored

./fitted - location where corrected beamwise catalogues are stored

./medians - location where median fits are stored

./png - location where per-sbid correction plots are stored

./pickle - location where corrections are stored

