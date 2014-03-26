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

# define directory where results will be saved
savedir = "/Users/angusr/angusr/ACF2" # directory in which results are saved

# index fits files
lc_files = np.array(glob.glob('/Users/angusr/angusr/data2/all_Qs/kplr*_llc.fits'))

# quarter list
quarters = ["2009131105131", "2009166043257", "2009259160929", "2009350155506", \
        "2010078095331", "2010174085026", "2010265121752", "2010355172524", \
        "2011073133259", "2011177032512", "2011271113734", "2012004120508", \
        "2012088054726", "2012179063303", "2012277125453", "2013011073258", \
        "2013098041711"]

# record kid and quarter and run ACF
for i, lc_file in enumerate(lc_files):
    q = quarters.index(str(lc_file[48:61]))
    kid = lc_file[38:47]
    time, flux, flux_err = load_data(lc_file)
    corr_run(time, flux, flux_err, kid, q, savedir)

# # load list of targets
# id_list = np.genfromtxt("/Users/angusr/Python/Gyro/data/astero_targets.txt").T

# # Join quarters together
# # loop over stars
# for i, kid in enumerate(id_list):
#     # load each quarter
#     lc_files = np.array(glob.glob('/Users/angusr/angusr/data2/all_Qs/kplr*%s*_llc.fits'\
#             %int(kid)))
#     # initialise the time and flux arrays
#     time, flux, flux_err = load_data(lc_files[0])
#     # loop over quarters
#     for i in range(1, len(lc_files)):
#         time = np.concatenate((time, load_data(lc_files[i])[0]))
#         flux = np.concatenate((flux, load_data(lc_files[i])[1]))
#         flux_err = np.concatenate((flux_err, load_data(lc_files[i])[2]))
#     # run ACF, once per star
#     corr_run(time, flux, flux_err, kid, 'all', savedir)
