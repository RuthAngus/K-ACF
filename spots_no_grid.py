# Wrapper for testing period measurement method
# This version does individual quarters
# This version produces star spots from the entire Kepler lc, then adds to individual quarters
# This version produces the lightcurves to be added to real lcs later
# Spots 5 does all 24 stars
# Spots 6 names stars as numbers not strings
# Spots 9 pretty much works and includes completeness
# Spots 10 doesn't use a grid - it generates random parameter values each time

#import matplotlib
#matplotlib.use('Agg')
import random
import numpy as np
import spot_sim_ruth
import scipy
import glob
import ACF_star_spot
import ss_period_plots
import ss_index
import pylab as p
import scipy.io
import pyfits
from matplotlib.pyplot import step

#======================================================
def run_master():                                      ##
    nstars = 1000
    # spot_gen(nstars)                                       ## (8 x number of periods, 3 x spot lifetimes = 24 light curves)
    # add_to_real(nstars)                           ## (inject into 10 Kepler light curves with 4 differen-t
    # ss_index.index('3')                             ## amplitudes = 960 light curves)
    # ss_index.index('4')                             ##
    # ss_index.index('5')                             ## (Index lightcurves for the ACF code)
    # ss_index.index('6')                             ##
    # ss_index.index('7')                             ##
    # ss_index.index('8')                             ##
    # ss_index.index('9')                             ##
    # ss_index.index('10')                            ##
    # ss_index.index('11')                            ##
    # ss_index.index('12')                            ##
    # ss_index.index('13')                            ##
    # ss_index.index('14')                            ##
    run_ACF(nstars)                                     ## (Calculate ACFs for all 960 lcs) 960 = 3 x 8 x 10 x 4
    recording_period_measurements(nstars)               ## (make note of the periods measured)
    period_plots(nstars)                                ## (Produce period results for each quarter)
    compare(nstars)                             
    # (Compare true vs measured periods)
    # population(stars)                                  #
    completeness()                                ##
    # random_test()
    return                                             ##
#======================================================

 

    
# Generate false light curves and saves them within star_spot_sim. The light
# curves are numbered with 2 digits. The first is corresponds to the star
# number and the second to the quarter (qs 3-6 = 1, qs 7-10 = 2, etc).
# It saves a text file containing the star's period as /ACF/star_spot_sim/sim_period*.txt

#-----------------------------------------------------------------------------------------------------------------
def spot_gen(nstars):

    period_check = np.zeros(nstars)

    print 'Generating light curves'
    ROOTDIR = '/Users/angusr/angusr/ACF/star_spot_sim'
    
    for index in range(0,nstars):
        
        # Generate random values for period that are uniform in log space.

        print 'Star %s...' %(index+1)
            
        rand_period = np.random.uniform(low = 0.0, high = 2.0)
        myperiod = 10**rand_period
        print 'Period = ', myperiod
        period_check[index] = myperiod

        # Generate random values of tau that are between 1 and 10 (uniform in log space too)

        rand_tau = np.random.uniform(low = 0.0, high = 1.0)
        mytau = 10**rand_tau
        print 'Tau = ', mytau
        
        grid = spot_sim_ruth.runSim(number_of_data_sets = 12, \
                                     myperiod = myperiod, star_name = index, spot_tau = mytau)
        grid = [grid[0], grid[1], grid[2], grid[3], myperiod, mytau]
        ''' grid contains [number of spots, filling factor, amplitude, noise, period, tau] '''

        # grid = myperiod
        # np.savetxt('/Users/angusr/angusr/ACF/star_spot_sim/tests/sim_period%s.txt' %(index+1), myperiod)
        # f = open('/Users/angusr/angusr/ACF/star_spot_sim/tests/sim_period%s.txt' %(index+1), 'w')
        # f.write('%s' %myperiod)
        # f.close()
        np.savetxt('/Users/angusr/angusr/ACF/star_spot_sim/sim_period%s.txt' %(index+1), grid)
    return 
#-------------------------------------------------------------------------------------------

def add_to_real(nstars):

    print 'Adding simulated lcs to real Kepler lcs'

    names = ['010972873', '011137075', '011502597', '011717120', '005607242', \
             '006370489', '006442183', '006531928', '006603624', '009574283']

    print 'Loading data...'
    kplr_lcs = [[],[],[],[],[],[],[],[],[],[]]
    kplr_ts = [[],[],[],[],[],[],[],[],[],[]]
    #kplr_lcs = np.ndarray([12,10])
    #kplr_ts = np.ndarray([12,10])

    for n in range(len(names)):
        
        files = glob.glob('/Users/angusr/angusr/data2/all_Qs/kplr%s-*llc.fits' % (names[n]))
        quarter_times_start = []; quarter_times_fin = []
        for quarter in range(0,12):
            time2 = []; lightcurve = []
            hdulist = pyfits.open(files[quarter])
            tbdata = hdulist[1].data #"add the first extension to tbdata"
            x = np.where(np.isfinite(tbdata['TIME'])) #Remove NANs
            time = tbdata['TIME'][x]; lc = tbdata['PDCSAP_FLUX'][x]
            x = np.where(np.isfinite(tbdata['PDCSAP_FLUX']))
            time = tbdata['TIME'][x]; lc = tbdata['PDCSAP_FLUX'][x]
            # Join quarters
            quarter_times_start.append(time[0]); quarter_times_fin.append(time[-1])
            lightcurve.extend(np.array(lc)); time2.extend(np.array(time))
            time = np.array(time2- min(time2))
            time = time + min(time2) # Add real starting value back to time array
            x = np.where(time == quarter_times_start[quarter])
            y = np.where(time == quarter_times_fin[quarter])
            start = int(x[0]); stop = int(y[0])
            time = time[start:stop] 
            lightcurve = lightcurve[start:stop] # Slice real data set into individual quarters
            lightcurve = lightcurve/np.median(lightcurve)
            lightcurve = lightcurve - np.median(lightcurve)
            kplr_lcs[n].append([]); kplr_ts[n].append([])
            kplr_lcs[n][quarter].append(lightcurve); kplr_ts[n][quarter].append(lightcurve)
            #kplr_lcs[(quarter ), n] = lightcurve
            #kplr_ts[(quarter), n] = time
    
    counter = 1

    for i in range(nstars):  

        # Choose Kepler lc at random
        sel_lc = np.random.randint(0,10)
    
        # Generate random values of tau that are between 1 and 15 (uniform in log space too)

        rand_amp = np.random.uniform(low = 0.0, high = np.log10(15))
        myamp = [10**rand_amp]
            
        print 'Star = ', i, 'Kplr lc = %s' %names[sel_lc], 'Amp = %s' %myamp, '#%s' %counter

        for quarter in range(0,12):
            print 'quarter = ', quarter+3

            '''Load simulated data'''
            mat = scipy.io.loadmat('/Users/angusr/angusr/ACF/star_spot_sim/%s/sim_%s.png.mat'\
                                   %((quarter+3), (i+1)))
            pars = mat['pars']; ts = mat['ts']
            lc2 = mat['ts']
            sim_time = lc2[0]; sim_lc =  lc2[2]

            time = kplr_ts[n][quarter][0]
            lightcurve = kplr_lcs[n][quarter][0]

            ''' Checking that the simulated and real quarters have the same number \
            of data points - the time stamps come from one particular lc, so there are mismatches!'''
            if len(time) > len(sim_time):
                for mis in range(len(time) - len(sim_time)):
                    print 'len(time) > len(sim_time)'
                    time = list(time); lightcurve = list(lightcurve)
                    time.remove(time[-1]); lightcurve.remove(lightcurve[-1])
                    time = np.array(time); lightcurve = np.array(lightcurve)
            if len(sim_time) > len(time):
                for mis in range(len(sim_time) - len(time)):
                    print 'len(sim_time) > len(time)'
                    time = list(time); lightcurve = list(lightcurve)
                    time.append(time[-1] + (time[1] - time[0])); lightcurve.append(lightcurve[-1])
                    time = np.array(time); lightcurve = np.array(lightcurve)

            '''Add simulated lc to real lc, where Amp is the scaling relation'''
            new_lc = myamp*sim_lc + lightcurve

            ts[2] = new_lc

            ''' Save back into mat file format'''
            scipy.io.savemat('/Users/angusr/angusr/ACF/star_spot_sim/%s/sim_%s.png.mat' \
                                %((quarter+3), counter), {'pars': pars, 'ts': ts})

        ''' Save parameters! '''
        grid_list = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/sim_period%s.txt' %(i+1))
        np.savetxt('/Users/angusr/angusr/ACF/star_spot_sim/grid/%sparams.txt' %counter, \
                        (np.transpose((counter, int(names[sel_lc]), myamp[0], grid_list[0], \
                                        grid_list[1], grid_list[2], grid_list[3], \
                            grid_list[4], grid_list[5]))))

        ''' Star number, KIC number, amplitude, number of spots, filling factor, amplitude, noise, period, tau '''
        counter += 1
                        

    print '%s iterations' %counter
    return


#-----------------------------------------------------------------------------------------------------------------
# Run ACF on new light curves. It takes the light curves from /ACF/star_spot_sim/3_6 etc. It saves the
# results to /ACF/PDCQss3-6 etc. It creates a file called results.txt which is then converted into a file called
# Periods_3-6.txt or similar.

#-----------------------------------------------------------------------------------------------------------------
def run_ACF(no_stars):
    test_list = ['ss3','ss4','ss5', 'ss6', 'ss7', 'ss8', 'ss9', 'ss10', 'ss11', 'ss12', 'ss13', 'ss14']
    test_list2 = ['3','4','5', '6', '7', '8', '9', '10', '11', '12', '13', '14']
    for a in range(0,len(test_list)):
        ACF_star_spot.corr_run(test_list[a], test_list2[a], number_of_stars = no_stars+1)
    return
#-----------------------------------------------------------------------------------------------------------------



# This takes the key results saved in results.txt and moves them into a file called Periods_ss3-6.txt etc

#-----------------------------------------------------------------------------------------------------------------
def recording_period_measurements(no_stars):
    print 'Recording period measurements... '
    list_of_quarters = ['ss3','ss4','ss5', 'ss6', 'ss7', 'ss8', 'ss9', 'ss10', 'ss11', 'ss12', 'ss13', 'ss14']
    for year in list_of_quarters:
        data = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ%s_output/results.txt' %year).T
        KID = np.array(range(1,no_stars+1))
        ACF_period =data[1][1:]
        period_error = data[6][1:]
        sine_period = data[2][1:]
        # KID = data[0][0]
        # #KID = range(1,no_stars+1)
        # ACF_period =data[1][0]
        # print ACF_period
        # period_error = data[6][0]
        # sine_period = data[2][0]
        print type(KID), type(ACF_period), type(period_error), type(sine_period)
        np.savetxt('/Users/angusr/angusr/ACF/PDCQ%s_output/Periods_%stest.txt' %(year,year), \
            np.transpose((KID, ACF_period, period_error, sine_period)))
        # file = open('/Users/angusr/angusr/ACF/PDCQ%s_output/Periods_%stest.txt' %(year, year), 'r')
        # for i in range(len(KID)):
            # file.write('%s   %s   %s   %s' %(str(period[x]), str(dlag_per_err[x]), str(sine_per[x])))
        # file.close()
    return
#-----------------------------------------------------------------------------------------------------------------




# This calls the period_plots_mod routine. It takes the files in PDCQss3-6 etc containing period information.
# It saves figures to /ACF/ss_year_figs and a file listing the stars with successfully measured periods
# as /ACF/star_spot_sim/ss_years.txt

#-----------------------------------------------------------------------------------------------------------------
def period_plots(no_stars):
    print 'Plotting period measurements...' 
    ss_period_plots.period_extract(ind_quarters = True, no_stars = no_stars)
    return 
#-----------------------------------------------------------------------------------------------------------------

# Need to gather period information from all over the place! This should print the measured
# period and the real period alongside each other
# Need to get it to read the stars listed in /star_spot_sim/ combine period measurements from 

#-----------------------------------------------------------------------------------------------------------------
def compare(nstars):

    # Successful stars
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
    
    m_periods = []; m_errors = []; median_periods = []
    for i in range(len(all_periods)):
        if len(all_periods[i]) > 0:
            per, err = weighted_mean(all_periods[i], all_errors[i])
            m_periods.append(per); m_errors.append(err)
            median_periods.append(np.median(all_periods[i]))
        else:
            m_periods.append(0); m_errors.append(0); median_periods.append(0)


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
        

    p.close(2)
    p.figure(2)
    x = np.arange(0,40,0.1)
    p.plot(x, x, 'b--')
    p.plot(.5*x, x, 'b--')
    p.plot(.333333*x, x, 'b--')
    p.plot(2*x, x, 'b--')
    p.plot(3*x, x, 'b--')
    p.plot(true_periods, m_periods, 'k.')
    # p.plot(true_periods, median_periods, 'r.')
    # p.errorbar(true_periods, m_periods, yerr = errors, fmt = 'k.')
    p.xlabel('True Period (days)')
    p.ylabel('Measured Period (days)')
    # p.xlim(min(periods) - 0.5, max(periods) + 0.5)
    # p.ylim(min(true_periods) - 1.5, max(true_periods) + 1.5)
    p.xlim(0, 40)
    p.ylim(0, 40)
    p.title('Measured vs true period')
    p.savefig('/Users/angusr/angusr/ACF/star_spot_sim/resultstest')
    star_list = np.array(star_list)
    true_periods = np.array(true_periods)
    period_list = np.array(period_list)
    np.savetxt('/Users/angusr/angusr/ACF/star_spot_sim/measured_vs_true.txt', \
                  np.transpose((star_list, true_periods, m_periods)))

    p.close(3)
    p.figure(3)
    p.subplot(1,2,1)
    p.plot(np.log10(true_periods), 'k.')
    p.subplot(1,2,2)
    p.plot(np.log10(test_periods), 'k.')

    return period_list
#-----------------------------------------------------------------------------------------------------------------    
# This function calculates the completeness: percentage of sample with successful period measurement,
# the reliability: fraction of objects with a given true period, detected within < 20 per cent of this value,
# the contamination: fraction of objects with measured period > 20 per cent different to the true period.
#-----------------------------------------------------------------------------------------------------------------

def completeness():

    # Load measured vs true stats
    data = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/measured_vs_true.txt').T
    list_of_success = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/ss_ind_quarterstest.txt')
    star_list = data[0]; m_periods = data[2]; true_periods = data[1]

    print 'Total number of simulated stars = %s' %len(star_list)
    print 'Total number of measured periods = %s' %len(list_of_success)
    print 'CALCULATING COMPLETENESS...'

    ''' Bin stars according to period '''
    nbins = 10
    success_period = []
    
    for i in range(len(star_list)):
        # Load true parameter values
        data = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/grid/%sparams.txt' %(i+1)).T
        kid_x = data[0]; KIC_no = str(data[1]); Amp = data[2]; nspots = data[3]; ff = data[4]; \
            amp = data[5]; noise = data[6]; true_period = data[7]; tau = data[8]

        # Was its period successfully measured?
        for j in range(len(list_of_success)):
            if kid_x == list_of_success[j]:
                 success_period.append(true_period) # FIXME: Do I want the true period or the measured period here?
                 

    all_bins = np.histogram(np.log10(true_periods), bins = nbins)
    success_bins = np.histogram(np.log10(success_period), bins = nbins)
    print 'Truth = ', all_bins[0]
    print 'Successfully measured = ', success_bins[0]

    b_space = 2./float(nbins)
    p.close(1)
    p.figure(1)
    
    p.subplot(3,1,1)
    p.ylabel('True (number)')
    step(np.arange(b_space,(2 + b_space), b_space), all_bins[0])
    #p.xlim(0,2)
    p.ylim(0, max(all_bins[0])+20)
    
    p.subplot(3,1,2)
    step(np.arange(b_space,(2 + b_space),b_space), success_bins[0])
    p.ylabel('Success (number)')
    p.ylim(0, max(success_bins[0])+20)
    
    p.subplot(3,1,3)
    p.ylabel('Completeness')
    complete = np.zeros(nbins)
    for i in range(len(complete)):
                   complete[i] = (float(success_bins[0][i])/float(all_bins[0][i]))*100
    print 'Completeness = ', complete
    step(np.arange(b_space,(2 + b_space),b_space), complete)
    p.xlabel('Log(Period)')
    p.ylim(0, max(complete)+20)
    
    

    # Reliability - figure out which stars have period measurements. Then compare their real periods
    # with the measured periods. number of stars with detected periods/number of stars with
    # measured periods within 20% of real.
    # Stars with period measurements = list_of_success
    # Their true periods = success_period
    # Their measured periods...
    measured_periods = np.zeros(len(list_of_success))
    uncert = np.zeros(len(list_of_success))
    median_periods = []

    

    all_KIDs = np.ndarray((len(list_of_success), 12))
    KID = np.ndarray((len(list_of_success), 12))
    periods = np.ndarray((len(list_of_success), 12))
    errors = np.ndarray((len(list_of_success), 12))

    ''' Reading in results (12 quarters for each star) '''
    for year in range(3,14):
        data = np.genfromtxt('/Users/angusr/angusr/ACF/PDCssQ%s_output/results.txt' %year).T
        periods.T[year] = data[1][1:]
        errors.T[year] = data[6][1:]
    
   
    for i in range(len(periods)):
        if len(periods[i]) > 0:
            per, err = weighted_mean(periods[i], errors[i])
            m_periods.append(per); m_errors.append(err)
            median_periods.append(np.median(all_periods[i]))
        else:
            m_periods.append(0); m_errors.append(0); median_periods.append(0)


    # Was its period successfully measured?
    for i in range(1, 1000):
        for j in range(len(list_of_success)):
            if i == list_of_success[j]:
                success_period.append(true_period) # FIXME: Do I want the true period or the measured period here?
                 
    
    return

#-----------------------------------------------------------------------------------------------------------------
                 # ''' Load measured period (for now just take quarter 3)'''
    #              #data = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQss3_output/results.txt').T
    #              #measured_periods.append(data[1][1:][kid_x-1])
    #              #errors.append(data[6][1:][kid_x-1])
    #              for k in range(3, 14):
    #                  data = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQss%s_output/results.txt' %k).T
    #                  measured_periods.append(data[1][1:][kid_x-1])
    #                  errors.append(data[6][1:][kid_x-1])
    #              m_periods, m_errors = weighted_mean(measured_periods, errors)
    # print 'm_periods', m_periods





# Plots histograms of true period and measured period (shows interesting offset?!?)
#-----------------------------------------------------------------------------------------------------------------

def population(nstars):

    # Load measured vs true stats
    data = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/measured_vs_true.txt')
    x = np.isfinite(data[2])
    period_list = data[2][x]

    print 'Total number of simulated stars = %s' %nstars
    print 'Total number of measured periods = %s' %len(period_list)

    print 'CALCULATING COMPLETENESS...'

    ''' Bin stars according to period '''
    nbins = 10
    #pbins = np.zeros((nbins,1))
    #bins = np.zeros(nbins)
    #true_bins = np.zeros(nbins)
    bins = np.ones(nbins)
    true_bins = np.ones(nbins)
    b_space = 2./float(nbins)
    bdries = []
    for k in range(1, nbins):
            bdries.append(k**(b_space*k))
    bdries.append(100.)

    measured_periods = []
    true_periods = []

    for i in range(1, nstars):
        # Load simulation data
        data = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/grid/%sparams.txt' %(i+1)).T
        kid_x = data[0]; KIC_no = str(data[1]); Amp = data[2]; nspots = data[3]; ff = data[4]; \
            amp = data[5]; noise = data[6]; true_period = data[7]; tau = data[8]

        #print 'Star=', kid_x, 'KIC=', KIC_no, 'Amp=', Amp, 'Nspots=', nspots, \
        #'ff=', ff, 'amp2=', amp, 'noise=', noise, 'tau=', tau
        

        ''' Load measured period '''
        m_period = []; m_err = []
        for j in range(3, 14):
            data = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQss%s_output/results.txt' %j).T
            periods = data[1][1:][i]
            errors = data[6][1:][i]
            
            m_period.append(periods)
            m_err.append(errors)
        m_period = np.array(m_period); m_err = np.array(m_err)
        x = np.where(m_period > 0.)
        m_period = m_period[x]
        m_err = m_err[x]
        #print m_period

        if len(m_period) > 1:
            period, error = weighted_mean(m_period, m_err)

            # For now just take the quarter 3 value!
            period = m_period[0]; error = m_err[0]
            measured_periods.append(period)
            true_periods.append(true_period)
            #print 'Measured period = ', period, '+/-',  error
        
            for m in range(len(bdries)-1):
                if bdries[m] < period < bdries[m+1]:
                    bins[m] += 1
                if bdries[m] < true_period < bdries[m+1]:
                    true_bins[m] += 1

            print 'Star = ', kid_x, 'TRUE PERIOD = ', true_period, 'Measured period = ', period

    print bins
    print true_bins
    #print (bins/true_bins)*100

    ''' Population study'''
    p.close(1)
    p.figure(1)
    p.subplot(2,1,1)
    #p.hist(true_periods)
    p.ylabel('True (number)')
    step(np.log10(bdries), true_bins)
    p.ylim(0, max(true_bins)+1)
    p.subplot(2,1,2)
    #p.hist(measured_periods)#, np.log10(bdries))
    step(np.log10(bdries), bins)
    p.ylabel('Measured (number)')
    p.ylim(0, max(bins)+1)
    #p.subplot(3,1,3)
    #step(np.log10(bdries), (bins/true_bins)*100)
    #p.ylabel('Completeness (%)')
    p.xlabel('Log(Period)')
    #p.ylim(0, max((bins/true_bins)*100 + 10))

    

    return

#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
def weighted_mean(all_periods, all_errors):
    #for i in range(len(all_errors)):
    #    if all_errors[i] == 0.0:
    #        all_errors[i] = 0.01
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

#-----------------------------------------------------------------------------------------------------------------


def random_test():
    rand_periods = np.zeros(1000)
    periods = np.zeros(1000)
    for i in range(1000):
        rand_periods[i] = np.random.uniform(low = 0.0, high = 2.0)
        periods[i] = 10**rand_periods[i]

    true_periods = np.zeros(1000)
    for i in range(1000):
        data = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/tests/sim_period%s.txt' %(i+1))
        true_periods[i] = data


    p.close(4)
    p.figure(4)
    p.subplot(3,1,1)
    p.plot(rand_periods, 'k.')
    p.subplot(3,1,2)
    p.plot(periods, 'k.')
    p.subplot(3,1,3)
    p.plot(np.log10(true_periods) ,'k.')
    # p.plot(range(100,200), np.log10(true_periods), 'k.')
    # p.plot(range(200,300), np.log10(true_periods), 'k.')
    # p.plot(range(300,400), np.log10(true_periods), 'k.')
    # p.plot(range(400,500), np.log10(true_periods), 'k.')

    ''' Plotting as close to original periods as I can'''

    p.close(10)
    p.figure(10)
    p.subplot(1,2,1)
    orig_periods = np.zeros(100)
    for i in range(100):
        data = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/sim_period%s.txt' %(i+1)).T
        p.axhline(np.log10(data[4]), color = 'k')
    p.subplot(1,2,2)
    for i in range(100):
        data = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/grid/%sparams.txt' %(i+1)).T
        p.axhline(np.log10(data[7]), color = 'k')
    
   
    return
    
