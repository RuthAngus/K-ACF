import numpy as np
import pylab as pl
from matplotlib.pyplot import step

def weighted_mean(all_periods, all_errors):
    errors = []
    for i in range(len(all_errors)):
        errors.append(1.0/((all_errors[i])**2))
    weights = errors
    weights = weights/sum(weights) # Normalise
    period = np.average(all_periods, weights = weights) # Find weighted mean
    for i in range(len(all_errors)):
        all_errors[i] = all_errors[i]**2 # Add errors in quadrature
    err_sum = sum(all_errors)
    error = np.sqrt(err_sum)

    return period, error


# Successful stars
nstars = 1000
stars = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/ss_ind_quarterstest.txt')
print 'successful: ', len(stars)
    
which_quarter = ['ss3','ss4','ss5', 'ss6', 'ss7', 'ss8', 'ss9',\
                 'ss10', 'ss11', 'ss12', 'ss13', 'ss14']
    
all_periods = np.ndarray((nstars, 12))
all_errors = np.ndarray((nstars, 12))
all_KIDs = np.ndarray((nstars, 12))
KID = np.ndarray((nstars, 12))
periods = np.ndarray((nstars, 12))
errors = np.ndarray((nstars, 12))

temp_KID = np.array(range(1,nstars+1))
''' Reading in results (12 quarters for each star) '''
for year in range(0, len(which_quarter)):
    data = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ%s_output/results.txt' %which_quarter[year]).T
    KID.T[year] = temp_KID
    periods.T[year] = data[1][1:]
    errors.T[year] = data[6][1:]

all_periods = periods
all_errors = errors
all_KIDs = KID

''' Finding mean and median measured periods '''
m_periods = []; m_errors = []; median_periods = []
for i in range(len(all_periods)):
    if len(all_periods[i]) > 0:
        per, err = weighted_mean(all_periods[i], all_errors[i])
        m_periods.append(per); m_errors.append(err)
        median_periods.append(np.median(all_periods[i]))
    else:
        m_periods.append(0); m_errors.append(0); median_periods.append(0)


''' Finding true periods'''
true_periods = []; star_list = []; period_list = []
for i in range(nstars):
    # print 'Star = ', i+1
    # print 'Measured period = ', m_periods[i], '+/-', m_errors[i]
    # print 'Median = ', median_periods[i]
    data = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/grid/%sparams.txt' \
                                   %(i+1))
    t_period = data[7]
    # print 'True period = ', t_period
    true_periods.append(t_period); star_list.append(i)
        

m_periods = median_periods

test_periods = np.zeros(1000); n = 0
for i in range(10):
    for j in range(100):
        data = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/sim_period%s.txt' %(j+1))
        test_periods[n] = data[4]
        n += 1
        

pl.close(1)
pl.figure(1)
x = np.arange(0,40,0.1)
pl.plot(x, 1.05*x, 'r--')
pl.plot(x, 0.95*x, 'r--')
pl.plot(x, x, 'b--')
pl.plot(.5*x, x, 'b--')
pl.plot(.333333*x, x, 'b--')
pl.plot(2*x, x, 'b--')
pl.plot(3*x, x, 'b--')
pl.plot(true_periods, m_periods, 'k.')
pl.xlabel('True Period (days)')
pl.ylabel('Measured Period (days)')
pl.xlim(0, 40)
pl.ylim(0, 40)
pl.title('Measured vs true period')
# pl.savefig('/Users/angusr/angusr/ACF/star_spot_sim/resultstest')
# pl.savefig('Compare')
pl.show()
star_list = np.array(star_list)
true_periods = np.array(true_periods)
period_list = np.array(period_list)
np.savetxt('/Users/angusr/angusr/ACF/star_spot_sim/measured_vs_true.txt', \
           np.transpose((star_list, true_periods, m_periods, m_errors)))


