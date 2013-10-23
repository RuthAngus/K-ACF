import numpy as np
import pylab as pl
from matplotlib.pyplot import step

def completeness_plots(all_bins, success_bins):
    
    ''' Plots of true, measured and completeness '''    
    pl.subplot(3,1,1)
    pl.ylabel('True (number)')
    step(all_bins[1][1:], all_bins[0])
    pl.ylim(0, max(all_bins[0])+20)
    
    pl.subplot(3,1,2)
    step(success_bins[1][1:], success_bins[0])
    pl.ylabel('Success (number)')
    pl.ylim(0, max(success_bins[0])+20)
    
    pl.subplot(3,1,3)
    pl.ylabel('Completeness')
    complete = np.zeros(len(all_bins[0]))
    for i in range(len(complete)):
                   complete[i] = (float(success_bins[0][i])/float(all_bins[0][i]))*100
    print 'Completeness = ', complete
    step(all_bins[1][1:], complete)
    pl.xlabel('Log(Period)')
    pl.ylim(0, max(complete)+20)
    # pl.savefig('Completeness')
    pl.show()

    return complete


# Load measured vs true stats
data = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/measured_vs_true.txt').T
list_of_success = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/ss_ind_quarterstest.txt')
star_list = data[0]; m_periods = data[2]; true_periods = data[1]; m_err = data[3]
    
print 'Total number of simulated stars = %s' %len(star_list)
print 'Total number of measured periods = %s' %len(list_of_success)
print 'CALCULATING COMPLETENESS...'

''' Bin stars according to period '''
nbins = 20
success_period = []
m_success_period = []
m_success_period_err = []
success_tau = []
success_amp = []
taus = np.zeros(len(star_list))
Amps = np.zeros(len(star_list))

for i in range(len(star_list)):
    # Load true parameter values
    data = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/grid/%sparams.txt' %(i+1)).T
	
    kid_x = data[0]; KIC_no = str(data[1]); Amps[i] = data[2]; nspots = data[3]; ff = data[4]; \
	amp = data[5]; noise = data[6]; taus[i] = data[8] #FIXME: check which amps

# Was its period successfully measured? FIXME: check period plots!
# append periods, taus and amps to list
    x = np.where(kid_x == list_of_success)
    if len(x[0]) > 0:
	success_period.append(true_periods[x]) # FIXME: true period or the measured period here?
	m_success_period.append(m_periods[x])
	m_success_period_err.append(m_err[x])
	success_tau.append(taus[x])
	success_amp.append(Amps[x])

   
# Binning periods, taus, amps, etc
b_space = 2./float(nbins)
all_bins = np.histogram(np.log10(true_periods), bins = np.arange(0, 2., b_space))
success_bins = np.histogram(np.log10(success_period), bins = np.arange(0, 2., b_space))
all_taus = np.histogram(np.log10(taus), bins = np.arange(0, 1., 1./float(nbins)))
tau_bins = np.histogram(np.log10(success_tau), bins = np.arange(0, 1., 1./float(nbins)))
all_amps = np.histogram(np.log10(Amps), bins = np.arange(0, np.log10(15), 1./float(nbins)))
amp_bins = np.histogram(np.log10(success_amp), bins = np.arange(0, np.log10(15), 1./float(nbins)))
print 'Truth = ', all_bins[0]
print 'Successfully measured = ', success_bins[0]

pl.close(2)
pl.figure(2)
period_complete = completeness_plots(all_bins, success_bins)

raw_input('enter')

pl.close(3)
pl.figure(3)
tau_complete = completeness_plots(all_taus, tau_bins)

raw_input('enter')

pl.close(4)
pl.figure(4)
amp_complete = completeness_plots(all_amps, amp_bins)

raw_input('enter')

''' Plot of completeness for period, tau and amp '''
pl.close(5)
pl.figure(5)

pl.subplot(3,1,1)
pl.ylabel('Period')
step(all_bins[1][1:], period_complete)
pl.xlabel('Log(Period)')
pl.ylim(0, max(period_complete)+20)

pl.subplot(3,1,2)
pl.ylabel('Tau')
step(all_taus[1][1:], tau_complete)
pl.xlabel('Log(Tau)')
pl.ylim(0, max(tau_complete)+20)

pl.subplot(3,1,3)
pl.ylabel('Amplitude')
step(all_amps[1][1:], amp_complete)
pl.xlabel('Log(Amplitude)')
pl.ylim(0, max(amp_complete)+20)

#-------------------------------------------------------------------

# Reliability - figure out which stars have period measurements. Then compare their real periods
# with the measured periods. number of stars with detected periods/number of stars with
# measured periods within 20% of real.
# Stars with period measurements = list_of_success
# Their true periods = success_period
# Their measured periods = m_success_period
# Their uncertainties = m_success_period_err

''' Plots of completeness, reliability and contamination for period '''
success_period = np.array(success_period)
m_success_period = np.array(m_success_period)
m_success_period_err = np.array(m_success_period_err)

m_bins = np.histogram(np.log10(m_success_period), bins = np.arange(0, 2., b_space))

reliability = np.zeros(nbins)
contamination = np.zeros(nbins)
for i in range(nbins):
    sbin = []; mbin = []
    for j in range(len(success_period)):
        if b_space*i < np.log10(success_period[j]) < b_space*(i+1):
	    sbin.append(success_period[j])
            mbin.append(m_success_period[j])

    if len(sbin) > 0:
        reliability[i], contamination[i] = comp_rel_cont(np.array(sbin), np.array(mbin))
    else:
        reliability[i] = 0.
	contamination[i] = 0.

pl.close(6)
pl.figure(6)
pl.subplot(3,1,1)
step(np.arange(0, (2.-2./nbins), 2./nbins), period_complete)
pl.ylabel('Completeness')
pl.ylim(0, 120)
pl.subplot(3,1,2)
step(np.arange(0, 2., 2./nbins), reliability)
pl.ylabel('Reliability')
pl.ylim(0, 120)
pl.subplot(3,1,3)
step(np.arange(0, 2., 2./nbins), contamination)
pl.ylabel('Contamination')
pl.xlabel('Log(Period)')
pl.ylim(0, 120)
# pl.savefig('Completeness')


 # for j in range(len(list_of_success)):
        # if kid_x == list_of_success[j]:
		# success_period.append(true_periods[i]) # FIXME: Do I want the true period or the measured period here?
		# m_success_period.append(m_periods[i])
		# m_success_period_err.append(m_err[i])
		# success_tau.append(taus[i])
		# success_amp.append(Amps[i])
		# if np.log10(true_periods[i]) > 1.4:
			# print 'SUCCESS!', true_periods[i]


