import numpy as np
from match import match

# load data
data = np.genfromtxt("/Users/angusr/Python/Gyro/data/new_matched.txt").T
KID = data[0]
period = data[1]
period_err = data[2]

# load ACF results
data2 = np.genfromtxt("/Users/angusr/angusr/ACF2/results.txt").T
KID2 = data2[0]
period2 = data2[1]
period_err2 = data2[2]

data3 = match(KID, data2)
KID3 = data3[0]
period3 = data3[1]
period_err3 = data3[2]

n = 0
for i in range(len(KID)):
    if period_err[i] > period_err3[i]:
        print KID[i], KID3[i]
        print period[i], period3[i]
        print period_err[i], period_err3[i], '\n'
        n+=1
    else:
        data[2][i] = period_err3[i]

print n, 'out of', len(KID)

np.savetxt("/Users/angusr/Python/Gyro/data/p_errs.txt", data.T)
