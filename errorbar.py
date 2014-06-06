import numpy as np
from match import match

# load data
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/new_matched.txt").T
KID = data[0]
period = data[1]
period_err = data[2]

# load ACF results
data = np.genfromtxt("/Users/angusr/angusr/ACF2/results.txt").T

data = match(KID, data)
nKID = data[0]
nperiod = data[1]
nperiod_err = data[2]

n = 0
for i in range(len(KID)):
    if period_err[i] > nperiod_err[i]:
        print KID[i], nKID[i]
        print period[i], nperiod[i]
        print period_err[i], nperiod_err[i], '\n'
        n+=1
    else:
        data[2][i] = nperiod_err[i]

print n, 'out of', len(KID)

np.savetxt("/Users/angusr/Python/Gyro/data/p_errs.txt", data)
