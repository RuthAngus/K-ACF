import numpy as np
import matplotlib.pyplot as pl

def match(KID, X):
    n = len(KID)
    m, l = X.shape
    KID_astero = X[0,:]
    Xm = np.zeros((m, n))
    for i in range(n):
        l = np.where(KID_astero == KID[i])[0]
        Xm[:,i] = X[:,l].reshape(Xm[:,i].shape)
    return Xm

# Target list:
data = np.genfromtxt('/Users/angusr/Desktop/astero_ages.txt').T
KIDs = data[0]
print len(KIDs)

# max = 528
quarters = range(3, 17)
grid = np.ndarray((len(quarters), 4, len(KIDs)))
print np.shape(grid)

quarter = 3
for quarter in quarters:
    data = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ%s_output/Periods_%s.txt' %(quarter, quarter)).T
    diff = 528 - len(data[0])
    if diff > 0:
        padding = np.zeros(diff)
        KIDs = np.concatenate(data[0], padding)
        periods = np.concatenate(data[1], padding)
        period_errs = np.concatenate(data[2], padding)
        sines = np.concatenate(data[3], padding)

    data_matched = match(KIDs, data)
    periods = data_matched[-3,:]
    period_errs = data_matched[-2,:]
    sines = data_matched[-1,:]

    print np.shape(grid[0,::])
    print np.shape(data_matched)
    grid[0,::] = data_matched
    #grid[1,:][quarter-3] = periods
    #grid[2,:][quarter-3] = period_errs
    #grid[3,:][quarter-3] = sine
    raw_input('enter')



for i in range(len(KIDs)):
    if periods[i] > 0. and period_errs[i] > 0.:
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
