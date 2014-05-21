# Read in and save all period measurement data

import numpy as np
import matplotlib.pyplot as pl
import glob

# load 'all'
def load_all(KIDs, DIR):
    data = np.zeros((len(KIDs), 3))
    data[:,0] = KIDs
    for i, k in enumerate(KIDs):
        ps = np.genfromtxt("%s/%s_all_result.txt"%(DIR, int(k))).T
        data[i,1:] = ps
    np.savetxt("%s/results.txt"%DIR, data)

# Load period data for each target
def load_each(KIDs, DIR, Kepler = True):
    nq = 16
    nresults = 800
    KIDs = KIDs[:nresults]
    qdata = np.zeros((len(KIDs), nq, 3))
    KIDs = np.array(KIDs)
    qdata[:,:,0] = KIDs[:, None]
    for i, k in enumerate(KIDs):
        p_data = np.array(glob.glob("%s/%s_*_result.txt"%(DIR, int(k))))
        p_data = np.concatenate((p_data[7:16], p_data[0:7]))
        qs = []
        K = len(str(k))
        if Kepler == True:
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
        else:
            qs = np.array(range(1, 17))
        ps = np.zeros((nq, 2))
        for q in range(len(qs)):
            try:
                ps[q,:] = np.genfromtxt(p_data[q]).T
            except:
                pass
        qdata[i,:,1:] = ps
        np.savetxt("%s/%s_results.txt"%(DIR, int(k)), qdata[i])

if __name__ == "__main__":
    # load target list
#     KIDs = np.genfromtxt("/Users/angusr/Python/Gyro/data/astero_targets.txt").T
    KIDs = range(1, 1001)

    # set directories and load results
#     DIR = "/Users/angusr/angusr/ACF2"
#     DIR = "/Users/angusr/angusr/injections"
    DIR = "/Users/angusr/angusr/Suz_simulations"
#     load_all(KIDs, DIR)
    load_each(KIDs, DIR, Kepler = False)

#     for i in range(len(KIDs)):
#         data = np.genfromtxt("%s/%s_all_result.txt"%(DIR, i)
#         qdata = np.genfromtxt(
#         print data[i,:]
#         print qdata[i,:,:]
#         print qdata[i,:,1]
#         periods = qdata[i,:,1]
#         print np.median(qdata[i,:,1])
#         raw_input('enter')
