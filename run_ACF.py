# This is the wrapper that calls KeplerACF.py

import numpy as np
import glob
from Kepler_ACF import corr_run

id_list = ['3544595','10124866']
lc_file = np.array(glob.glob('/Users/angusr/angusr/Kepler/%s/kplr*_llc.fits' %id_list[0]))

# corr_run(id_list[0], lc_file[0], 0)
quarters = range(0, 18)

def load_data(lc_file):
    hdulist = pyfits.open(lc_file)
    tbdata = hdulist[1].data
    x = tbdata["TIME"]
    y = tbdata["PDCSAP_FLUX"]
    yerr = tbdata["PDCSAP_FLUX_ERR"]
    q = tbdata["SAP_QUALITY"]

    # remove nans
    n = np.isfinite(x)*np.isfinite(y)*np.isfinite(yerr)*(q==0)

    # median normalise
    y = y/np.median(y[n]) - 1.
    yerr = yerr/np.median(y[n])

    return x[n], y, yerr

for i, file in enumerate(lc_file):
    corr_run(id_list[0], file, quarters[i])
