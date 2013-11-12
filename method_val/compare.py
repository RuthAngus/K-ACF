import random
import numpy as np
import spot_sim_ruth
import scipy
import glob
import ACF_star_spot
import ss_period_plots
import ss_index
import pylab as pl
import scipy.io
import pyfits
from matplotlib.pyplot import step
import scipy.optimize as optimization


plotpar = {'axes.labelsize': 16,
           'text.fontsize': 25,
           'legend.fontsize': 14,
           'xtick.labelsize': 17,
           'ytick.labelsize': 17, 
           'text.usetex': True}
pl.rcParams.update(plotpar)

def trend_removal(x, y):
    p0 = scipy.polyfit(x, (y-x), 1)
    p = scipy.poly1d(p0)
    x = true_periods
    dy = p(x)
    y = mps
    z = y - dy
    return z

def weighted_mean(periods, errors):
    # remove zeros 
    x = periods > 0
    periods = periods[x]
    errors = errors[x]
    x = errors > 0
    periods = periods[x]
    errors = errors[x]
    
    if len(errors) > 0:
        weights = 1./errors**2
        weights /= sum(weights) # Normalise
        period = np.average(periods, weights = weights) # Find weighted mean
        err = np.sqrt(sum(errors**2))/len(errors)

        return period, err
    else:
        return 0., 0.


# Load list of stars with successful period measurements
stars = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/ss_ind_quarterstest.txt')
print 'successful: ', len(stars)
    
nstars = 1000
nqs = 12

''' Create empty arrays for measured params
    mps = 2-d array containing measured periods for each quarter for each star
    mp_errs = the same for the errors
    KIDs = 2-d array indexing stars
    meanp = weighted mean period measurements for each star
    meanp_err = the weighted mean uncertainties
    medianp = median period values
    tps = true periods'''
    
# Load true periods
tps = [np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/grid/%sparams.txt' %(i+1))[7] for i in range(nstars)]

# Load measured periods
mps = [np.genfromtxt('/Users/angusr/angusr/ACF/PDCQss%s_output/results.txt' %(q+3)).T[1][1:] for q in range(nqs)]
mp_errs = [np.genfromtxt('/Users/angusr/angusr/ACF/PDCQss%s_output/results.txt' %(q+3)).T[6][1:] for q in range(nqs)]

meanp, meanp_err = zip(*[weighted_mean(mps[i], mp_errs[i]) for i in range(len(mps))])
medianp = [np.median(mps[i]) for i in range(len(mps))]
      
mps = medianp
np.savetxt('/Users/angusr/angusr/ACF/star_spot_sim/measured_vs_true2.txt', \
                  np.transpose((tps, mps, mp_errs)))

a = m_periods > 0 # Remove non-detections
# Mask removes the outliers
mask = [np.nan for i in range(len(tps[a])) if 17 < tps[a][i] < 22 and 6 < mps[a][i] < 10]
mask = np.isfinite(mask)

mask = np.zeros(len(tps[a]))
for i in range(len(tps[a])):
    if 17 < tps[a][i] < 22 and 6 < mps[a][i] < 10:
        mask[i] = np.nan
mask = np.isfinite(mask)

# Main compare plot
cols = ['#FF9933', '#339999']
pl.close(1)
pl.figure(1)
ax1 = pl.subplot2grid((4,4), (0,0), colspan=4, rowspan = 3)
x = np.arange(0,40,0.1)
pl.plot(x, 1.05*x, color = cols[0], linestyle = '--')
pl.plot(x, 0.95*x, color = cols[0], linestyle = '--')
pl.plot(x, x, color = cols[1], linestyle = '--')
pl.plot(.5*x, x, color = cols[1], linestyle = '--')
pl.plot(x, 0.4*x, color = cols[0], linestyle = '--')
pl.plot(x, 0.6*x, color = cols[0], linestyle = '--')
pl.plot(2*x, x, color = cols[1], linestyle = '--')
pl.errorbar(tps[a], mps[a], yerr = mp_errs[a], \
            fmt = 'k.', markersize = 2, capsize = 0, ecolor = '0.7' )
pl.ylabel('$\mathrm{Recovered~period~(days)}$', fontsize = 25)
pl.xlim(0, 25)
pl.ylim(0, 25)
pl.gca().set_xticklabels([])
   
# Residuals
ax2 = pl.subplot2grid((4,4), (3, 0), colspan = 4)
pl.errorbar(tps[a][mask], (-tps[a][mask] + mps[a][mask]), yerr = mp_errs[a][mask], \
            fmt = 'k.', ecolor = '0.7', markersize = 2, capsize = 0)
pl.subplots_adjust(top=0.96, bottom = 0.06)
pl.xlabel('$\mathrm{Injected~period~(days)}$', fontsize = 25)
pl.ylabel('$\mathrm{Residuals}$', fontsize = 25)
pl.xlim(0, 25)
pl.ylim(-5,2)
pl.subplots_adjust(hspace = 0., bottom = 0.1)
pl.yticks(range(-5, 2), range(-5, 2))
pl.savefig('/Users/angusr/angusr/ACF/star_spot_sim/resultstest')
pl.savefig('Compare')
pl.show()

# Removing trend - need to figure out how to find the uncertainty
z = trend_removal(tps[a][mask], mps[a][mask])

# Plot with trend removed
pl.close(2)
pl.figure(2)
x = np.arange(0,40,0.1)
pl.plot(x, 1.05*x, color = cols[0], linestyle = '--')
pl.plot(x, 0.95*x, color = cols[0], linestyle = '--')
pl.plot(x, x, color = cols[1], linestyle = '--')
pl.plot(.5*x, x, color = cols[1], linestyle = '--')
pl.plot(x, 0.6*x, color = cols[0], linestyle = '--')
pl.plot(x, 0.4*x, color = cols[0], linestyle = '--')
pl.plot(2*x, x, color = cols[1], linestyle = '--')
pl.errorbar(tps[a], (z[a]), yerr = mp_errs[a], fmt = 'k.', ecolor = '0.7', markersize = 2, capsize = 0)
b = mps == 0
pl.ylabel('Recovered Period (days)', fontsize = 25)
pl.xlabel('Injected Period (days)', fontsize = 25)
pl.xlim(0, 25)
pl.ylim(0, 25)
pl.savefig('trend_removed')
pl.show()
   
