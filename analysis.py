# Read in and assess all period measurement data

import numpy as np
import matplotlib.pyplot as pl

# load 'all' data
data = np.genfromtxt("/Users/angusr/angusr/ACF2/results.txt").T
all_KIDs = data[0]
all_ps = data[1]
all_perrs = data[2]

# load target list
KIDs = np.genfromtxt("/Users/angusr/Python/Gyro/data/astero_targets.txt").T

# tuning params
m = .2 # periods must lie within what fraction of the median?
f = .7 # what fraction of total measurements required to be consistent?
l = 0. # what number of quarters must be available?
d = 100. # within what fraction of the median must the residual p measurement lie?

accept = []; accept_p = []; accept_perr = []
for i, k in enumerate(KIDs):
    qdata = np.genfromtxt("/Users/angusr/angusr/ACF2/%s_results.txt"%(int(k))).T
    ind_KIDs = qdata[0]
    ps = qdata[1]
    perrs = qdata[2]
    med = np.median(ps[ps>0])
    a = (ps[ps>0]<med+(m*med))*(med-(m*med)<ps[ps>0])
    fraction = float(len(ps[a]))/float(len(ps[ps>0]))
    if fraction > f and len(ps[ps>0]) > l and np.sqrt((med-all_ps[i])**2)<d*med:
        print "accept"
        accept.append(k)
        accept_p.append(all_ps[i])
        accept_perr.append(all_perrs[i])
        print k, all_ps[i], all_perrs[i]

print len(accept)
np.savetxt("/Users/angusr/angusr/ACF2/periods.txt", \
        np.transpose((accept, accept_p, accept_perr)))
