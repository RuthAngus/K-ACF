# Attempt to rewrite period_plots

import numpy as np
import matplotlib.pyplot as pl

remove = [3852594, 4859338, 6313425, 7106245, 7174707, 7510397, 7970740, 8346342, 9157245, 9965715, \
           10124866, 10273246, 10339342, 10351085, 10514430, 10593351, 10730618, 10920273, 11234888, \
           11395018, 11772920, 12069449]

def match(KID, X):
    m, l = X.shape
    KID_success = X[0,:]
    n = len(KID)
    print n, len(KID_success)
    Xm = np.zeros((m, n))
    for i in range(n):
        l = np.where(KID_success == KID[i])[0]
        Xm[:,i] = X[:,l].reshape(Xm[:,i].shape)
    return Xm

# Target list:
data = np.genfromtxt('/Users/angusr/Desktop/astero_ages.txt').T
KIDs = data[0]
print len(KIDs)

# Figure out which KIDs are missing
#compare = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ3_output/Periods_3.txt').T
#comp = compare[0]
#for i in range(len(comp)):
#    l = np.where(KIDs == comp[i])[0]
#    if len(l) == 0:
#        print comp[i]
#raw_input('enter')

# max = 528
quarters = range(3, 17)
grid = np.ndarray((len(quarters), 4, len(KIDs)))

for quarter in quarters:
    data = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ%s_output/Periods_%s.txt' %(quarter, quarter)).T
    diff = 528 - len(data[0])
    # Remove extra targs
    for i in remove:
    	x = np.where(data[0] == i)[0]
        if len(x) > 0:
            KIDs = np.concatenate((data[0][:x], data[0][x+1:]))
            periods = np.concatenate((data[1][:x], data[1][x+1:]))
            period_errs = np.concatenate((data[2][:x], data[2][x+1:]))
            sines = np.concatenate((data[3][:x], data[3][x+1:]))
    data = np.empty((4, len(KIDs)))
    data[0,:] = KIDs
    data[1,:] = periods
    data[2,:] = period_errs
    data[3,:] = sines
    print len(data[0])
    
    # add padding to make all arrays the same length
    if diff > 0:
        padding = np.zeros(diff)
        KIDs = np.concatenate((data[0], padding))
        periods = np.concatenate((data[1], padding))
        period_errs = np.concatenate((data[2], padding))
        sines = np.concatenate((data[3], padding))
 
        data = np.empty((4, len(KIDs)))
        data[0,:] = KIDs
        data[1,:] = periods
        data[2,:] = period_errs
        data[3,:] = sines

    data_matched = match(KIDs, data)
    
    grid[(quarter-3),::] = data_matched
    #raw_input('enter')

print grid[0]
raw_input('enter')
for i in range(len(KIDs)):
    for j in range(len(quarters)):
        if grid[i] > 0. and period_errs[i] > 0.:
            pl.errorbar(quarters[3], periods[i], yerr = period_errs[i], marker = 'o', color = 'b', markersize = 5)
    else: pl.axvline(quarters[3], linewidth = 0.5, color = 'r', linestyle = '--')
    pl.ylabel('Period')
    pl.xlabel('Quarters')

    # # Harmonic lines
    # pl.axhline(acf_median, linewidth = 0.5, color = 'b')
    # pl.axhline(acf_median/2., linewidth = 0.5, color = 'b', linestyle = '--')
    # period_multiple = 0; harmonic=2
    # if acf_median != 0:
    #     while period_multiple < upper_y_lim:
    #         period_multiple = acf_median*harmonic
    #         pylab.axhline(period_multiple, linewidth = 0.5, color = 'b', linestyle = '--')
    #         harmonic += 1
    
    pl.show()
    raw_input('enter')
