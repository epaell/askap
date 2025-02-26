# astrometry
Tools to check for astrometric variations between per-beam ASKAP field catalogues and with respect to UNWISE (specifically for RACS-low3).

# Tools

align_beams.py <catalogues> - a method of examining beam-wise catalogues for relative shifts (based on Daves code). Result is written to an offset file.

align_beams_alt.py <catalogues> - a fairly slow method of examining beam-wise catalogues for relative shifts (based on Tims code). Result is written to an offset file.

fit_unwise.py <catalogues> - an iterative method of comparing beam catalogues against UNWISE to correct for astrometric shifts. Result is written to an offset file. The appropriate UNWISE catalogue for a field will be downloaded if it is not already present.

plot_fit.py offset_file - plot the offset corrections for the specified file for each of the beams

plot_ref_fit.py offset_file - Same as plot_fit.py but references all offsets to beam 20 (to allow direct comparisons across methods)

self_stat.py <catalogues> offset_file - do a self check across overlapping beams to measure mean residual shifts and std (before and after applying the offsets)

ref_stat.py reference <catalogues> offset_file - do a check against external reference catalogues to measure mean residual shifts and std (before and after applying the offsets). 
   
ref/* - Reference catalogues used by ref_stat.

closepack36_beams.fits - a table containing the beam positions for all of the available RACS-low3 (EMU/RACS-mid1/RACS-high1 use the same footprints) observations.

closepack36_fields.fits - a table containing the field pointing centre for all of the available RACS-low3 (EMU/RACS-mid1/RACS-high1 use the same footprints) observations.