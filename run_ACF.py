# This is the wrapper that calls KeplerACF.py
import numpy as np
import glob
from Kepler_ACF import corr_run
import pyfits

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
    y[n] = y[n]/np.median(y[n]) - 1.
    yerr[n] = yerr[n]/np.median(y[n])
    return x[n], y[n], yerr[n]

def runacf(id_list, lc_file, join_quarters = False):
    # run ACF on individual quarters
    if join_quarters == False:
        for i, lc in enumerate(lc_file):
            time, flux, flux_err = load_data(lc)
            corr_run(time, flux, flux_err, id_list, q = i)
    # join quarters together
    else:
        time, flux, flux_err = load_data(lc_file[0])
        for i in range(1, len(lc_file)):
            time = np.concatenate((time, load_data(lc_file[i])[0]))
            flux = np.concatenate((flux, load_data(lc_file[i])[1]))
            flux_err = np.concatenate((flux_err, load_data(lc_file[i])[2]))
        corr_run(time, flux, flux_err, id_list[0], q = 'all')

# check consistency of period measurements
def check(DIR):
    nq = 18
    p = np.empty(nq); p_err = np.empty(nq)
    for i in range(nq):
        data = np.genfromtxt("%s/%sresult.txt"%(DIR, i))
        p[i] = data[0]
        p_err[i] = data[1]
    print "mean = ", np.mean(p)
    print "median = ", np.median(p)
    print np.genfromtxt("%s/allresult.txt"%DIR)
    for i in range(len(p)):
        print p[i], p_err[i]

id_list = ["11904151"]
lc_file = np.array(glob.glob('/Users/angusr/angusr/Kepler/%s/kplr*_llc.fits' %id_list[0]))
runacf(id_list[0], lc_file, join_quarters = False)
check("/Users/angusr/angusr/Kepler")

