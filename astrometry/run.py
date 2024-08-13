#!/usr/bin/env python
import os
import sys

sbid = int(sys.argv[1])
os.system("./fit_astrometry.py %d" %(sbid))
os.system("./apply_fit.py %d" %(sbid))
os.system("./plot_fit.py %d" %(sbid))
os.system("./assess_fit.py %d 1.0 0.75" %(sbid))
