# Read in and save all period measurement data

import numpy as np
import matplotlib.pyplot as pl
import glob

# load target list
KIDs = np.genfromtxt("/Users/angusr/Python/Gyro/data/astero_targets.txt").T

# load 'all'
data = np.zeros((len(KIDs), 3))
data[:,0] = KIDs
for i, k in enumerate(KIDs):
    ps = np.genfromtxt("/Users/angusr/angusr/ACF2/%s_all_result.txt"%int(k)).T
    data[i,1:] = ps
    np.savetxt("/Users/angusr/angusr/ACF2/results.txt", data)

# Load period data for each target
nq = 14
qdata = np.zeros((len(KIDs), nq, 3))
qdata[:,:,0] = KIDs[:, None]
for i, k in enumerate(KIDs):
    qs = []
    K = len(str(k))
    p_data = np.array(glob.glob("/Users/angusr/angusr/ACF2/%s_*_result.txt"%int(k)))
    [qs.append(p[25+K:27+K]) for p in p_data]
    try:
        qs.remove('al')
    except:
        "list.remove(x)"
    for j in range(len(qs)):
        if qs[j][1] == '_':
            qs[j] = float(qs[j][0])
        else: qs[j] = float(qs[j])
    qs = np.array(qs)
    qs, p_data = zip(*sorted(zip(qs, p_data)))
    ps = np.zeros((nq, 2))
    for q in range(len(qs)):
        ps[q,:] = np.genfromtxt(p_data[q]).T
    qdata[i,:,1:] = ps
    np.savetxt("/Users/angusr/angusr/ACF2/%s_results.txt"%int(k), qdata[i])

for i in range(len(KIDs)):
    print data[i,:]
    print qdata[i,:,:]
    print qdata[i,:,1]
    periods = qdata[i,:,1])
    print np.median(qdata[i,:,1])
    raw_input('enter')
