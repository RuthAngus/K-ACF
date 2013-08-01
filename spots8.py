# Wrapper for testing period measurement method
# This version does individual quarters
# This version produces star spots from the entire Kepler lc, then adds to individual quarters
# This version produces the lightcurves to be added to real lcs later
# Spots 5 does all 24 stars
# Spots 6 names stars as numbers not strings

import matplotlib
matplotlib.use('Agg')
import random
import numpy as np
import spot_sim_ruth5
import scipy
import glob
import ACF_star_spot8
import ss_period_plots3
import ss_index2
import pylab as p
import scipy.io
import pyfits
from matplotlib.pyplot import step

# This is the master function that calls everything below.

#======================================================
def run_master():                                    ##
    # spot_gen(8, [5,2,1])                             ## (8 x number of periods, 3 x spot lifetimes = 24 light curves)
    # add_to_real(10, 4, 24)                           ## (inject into 10 Kepler light curves with 4 differen-t
    # ss_index2.index('3')                             ## amplitudes = 960 light curves)
    # ss_index2.index('4')                             ##
    # ss_index2.index('5')                             ## (Index lightcurves for the ACF code)
    # ss_index2.index('6')                             ##
    # ss_index2.index('7')                             ##
    # ss_index2.index('8')                             ##
    # ss_index2.index('9')                             ##
    # ss_index2.index('10')                            ##
    # ss_index2.index('11')                            ##
    # ss_index2.index('12')                            ##
    # ss_index2.index('13')                            ##
    # ss_index2.index('14')                            ##
    # ##run_ACF(10)                                     ## (Calculate ACF for 1 lc) 1 = 1 x 1 x 1 x 1
    # run_ACF(960)                                     ## (Calculate ACFs for all 960 lcs) 960 = 3 x 8 x 10 x 4
    # ##recording_period_measurements(10)               ##   
    # recording_period_measurements(960)               ## (make note of the periods measured)
    # ##period_plots(10)                                ##
    # period_plots(960)                                ## (Produce period results for each quarter)
    # compare()                                        ## (Compare true vs measured periods)
    #population(960)                              ##
    completeness(960)
    return                                           ##
#======================================================

 

    
# Generate false light curves and saves them within star_spot_sim. The light
# curves are numbered with 2 digits. The first is corresponds to the star
# number and the second to the quarter (qs 3-6 = 1, qs 7-10 = 2, etc).
# It saves a text file containing the star's period as /ACF/star_spot_sim/sim_period*.txt

#-----------------------------------------------------------------------------------------------------------------
def spot_gen(no_periods, taus):

    grid_list = []
    print 'Generating light curves'
    ROOTDIR = '/Users/angusr/angusr/ACF/star_spot_sim'

    
    n = 0
    for i in range(0, len(taus)):
        
        # Generate light curves with 3 values of tau
        tau = taus[i]
        for index in range(0,no_periods):
            # Create lightcurves for 24 stars: 3 spot lifetimes and 8 periods.
            # Generate random values for period that are uniform in log space.

            #print 'Star %s...' %(index+1)
            print 'Star %s...' %(n+1)
            
            ''' Try generating numbers from a distribution that is uniform in log space '''
            #myperiod = [(np.sqrt(random.gauss(3., .5))**2)]
            random_gen = np.random.uniform(low = 0.0, high = 2.0)
            myperiod = [10**random_gen]
            save_period_in_txt_file = float(myperiod[0])


            # These are lengths of time that correspond to quarters 3-6, 7-10
            # and 11-14 respectively in days.
            durations = [363.68, 347.27, 367.30]

            # Create 3 lightcurves, one for each year of data,
            # for a star with a given period.
            #spot_sim_ruth_mod.runSim(number_of_data_sets = 3, myperiod = myperiod, star_name = index)
            grid = spot_sim_ruth5.runSim(number_of_data_sets = 12, \
                                         myperiod = myperiod, star_name = n, spot_tau = tau)
            grid = [grid[0], grid[1], grid[2], grid[3], save_period_in_txt_file, tau]
            ''' grid contains [number of spots, filling factor, amplitude, noise, period, tau] '''
            np.savetxt('/Users/angusr/angusr/ACF/star_spot_sim/sim_period%stest.txt' %(n+1), grid)
            n += 1
            grid_list.append(grid)
    return grid_list
#-------------------------------------------------------------------------------------------

def add_to_real(number_of_kplr_lcs, number_of_amps, nstars):

    print 'Adding simulated lcs to real Kepler lcs'
    
    counter = 1
    
    ''' Setting number of Kepler lcs '''
    if number_of_kplr_lcs == 10:
        names = ['010972873', '011137075', '011502597', '011717120', '005607242', \
             '006370489', '006442183', '006531928', '006603624', '009574283']
    elif number_of_kplr_lcs == 1:
        names = ['010972873']
        
    ''' Setting number of amplitudes '''
    if number_of_amps == 4:
        Amps = [1.0, 5.0, 10.0, 15.0]
    elif number_of_amps == 1:
        Amps = [1.0]

    for n in range(0, len(names)):  # 10
        print 'Kepler lc = ', names[n]
        files = glob.glob('/Users/angusr/angusr/data2/all_Qs/kplr%s-*llc.fits' % (names[n]))
        quarter_times_start = []; quarter_times_fin = []
        print len(files)

        for a in range(0, len(Amps)): # 4 (40)
            print 'Amplitude = ', Amps[a]
            
            for star in range(1, nstars+1): # 24 (960)
            #for star in range(1, 4):
                print 'Star = ', star, 'Kplr lc = %s' %names[n], 'Amp = %s' %Amps[a], '#%s' %counter

                for quarter in range(0,12):
                    print 'quarter = ', quarter+3

                    '''Load real data'''
                    time2 = []; lightcurve = []
                    hdulist = pyfits.open(files[quarter])
                    tbdata = hdulist[1].data #"add the first extension to tbdata"
                    x = np.where(np.isfinite(tbdata['TIME'])) #Remove NANs
                    time = tbdata['TIME'][x]; lc = tbdata['PDCSAP_FLUX'][x]
                    x = np.where(np.isfinite(tbdata['PDCSAP_FLUX']))
                    time = tbdata['TIME'][x]; lc = tbdata['PDCSAP_FLUX'][x]
                    quarter_times_start.append(time[0]); quarter_times_fin.append(time[-1])
                    lightcurve = list(lightcurve)
                    lightcurve.extend(np.array(lc)); time2.extend(np.array(time))
                    time = np.array(time2- min(time2))
                    addit = min(time2)
                    time = time + addit # Add real starting value back to time array
                    #print 'start = ', quarter_times_start[quarter]
                    #print 'stop = ', quarter_times_fin[quarter]
                    x = np.where(time == quarter_times_start[quarter])
                    y = np.where(time == quarter_times_fin[quarter])
                    start = int(x[0]); stop = int(y[0])
                    time = time[start:stop] 
                    lightcurve = lightcurve[start:stop] # Slice real data set into individual quarters
                    lightcurve = lightcurve/np.median(lightcurve)
                    lightcurve = lightcurve - np.median(lightcurve)

                
                    '''Load simulated data'''
                    mat = scipy.io.loadmat('/Users/angusr/angusr/ACF/star_spot_sim/%s/sim_%s.png.mat'\
                                           %((quarter+3), star))
                    pars = mat['pars']; ts = mat['ts']
                    lc2 = mat['ts']
                    sim_time = lc2[0]; sim_lc =  lc2[2]

                    '''Plotting'''
                    # p.close(1)
                    # p.figure(1)
                    # p.subplot(3,1,1)
                    # p.plot(sim_time, sim_lc, 'r.')
                    # p.subplot(3,1,2)
                    # p.plot(time, lightcurve, 'b.')
                    # p.subplot(3,1,3)

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
                    new_lc = Amps[a]*sim_lc + lightcurve

                    ts[2] = new_lc

                    ''' Save back into mat file format'''
                    scipy.io.savemat('/Users/angusr/angusr/ACF/star_spot_sim/%s/sim_%s.png.mat' \
                                        %((quarter+3), counter), {'pars': pars, 'ts': ts})

                ''' Save parameters! '''
                grid_list = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/sim_period%stest.txt' %star)
                print 'GRID_LIST', grid_list
                np.savetxt('/Users/angusr/angusr/ACF/star_spot_sim/grid/%sparams.txt' %counter, \
                              (np.transpose((counter, int(names[n]), Amps[a], grid_list[0], \
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
    # Running ACF code on false light curves
    test_list = ['ss3','ss4','ss5', 'ss6', 'ss7', 'ss8', 'ss9', 'ss10', 'ss11', 'ss12', 'ss13', 'ss14']
    test_list2 = ['3','4','5', '6', '7', '8', '9', '10', '11', '12', '13', '14']
    for a in range(0,len(test_list)):
        ACF_star_spot8.corr_run(test_list[a], test_list2[a], number_of_stars = no_stars+1)
    return
#-----------------------------------------------------------------------------------------------------------------



# This takes the key results saved in results.txt and moves them into a file called Periods_ss3-6.txt etc

#-----------------------------------------------------------------------------------------------------------------
def recording_period_measurements(no_stars):
    print 'Recording period measurements... '
    list_of_quarters = ['ss3','ss4','ss5', 'ss6', 'ss7', 'ss8', 'ss9', 'ss10', 'ss11', 'ss12', 'ss13', 'ss14']
    for year in list_of_quarters:
        data = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ%s_output/results.txt' %year).T
        KID = range(1,no_stars+1)
        ACF_period =data[1][1:]
        period_error = data[6][1:]
        sine_period = data[2][1:]
        # KID = data[0][0]
        # #KID = range(1,no_stars+1)
        # ACF_period =data[1][0]
        # print ACF_period
        # period_error = data[6][0]
        # sine_period = data[2][0]
        np.savetxt('/Users/angusr/angusr/ACF/PDCQ%s_output/Periods_%stest.txt' %(year,year), \
            np.transpose((KID, ACF_period, period_error, sine_period))) 
    return
#-----------------------------------------------------------------------------------------------------------------




# This calls the period_plots_mod routine. It takes the files in PDCQss3-6 etc containing period information.
# It saves figures to /ACF/ss_year_figs and a file listing the stars with successfully measured periods
# as /ACF/star_spot_sim/ss_years.txt

#-----------------------------------------------------------------------------------------------------------------
def period_plots(no_stars):
    print 'Plotting period measurements...' 
    ss_period_plots3.period_extract(ind_quarters = True, no_stars = no_stars)
    return 
#-----------------------------------------------------------------------------------------------------------------




# Need to gather period information from all over the place! This should print the measured
# period and the real period alongside each other
# Need to get it to read the stars listed in /star_spot_sim/ combine period measurements from 

#-----------------------------------------------------------------------------------------------------------------
def compare():

    ''' Reading in real periods '''
    stars = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/ss_ind_quarterstest.txt')
    
    which_quarter = ['ss3','ss4','ss5', 'ss6', 'ss7', 'ss8', 'ss9',\
                     'ss10', 'ss11', 'ss12', 'ss13', 'ss14']
    
    all_periods = []
    all_errors = []
    for star in range(0,len(stars)):
        all_periods.append([])
        all_errors.append([])

    ''' Reading in results '''
    for year in range(0, len(which_quarter)):
        # data = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ%s_output/Periods_%s.txt' \
                                # %(which_quarter[year], which_quarter[year])).T
        data = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ%s_output/Periods_%stest.txt'\
                                %(which_quarter[year], which_quarter[year])).T
        KID = data[0]
        periods = data[1]
        errors = data[2]
        for j in range(0,len(KID)):
            for star in range(0, len(stars)):
                if KID[j] == stars[star]:
                    all_periods[star].append(periods[j])
                    all_errors[star].append(errors[j])
   

    '''Find the weighted mean. Weights are 1/errors**2'''

    ''' Find weights '''
    errors = all_errors
    for i in range(0, len(all_errors)):
        for j in range(0, len(all_errors[i])):
            all_errors[i][j] = 1.0/((all_errors[i][j])**2)
    weights = all_errors

    ''' Normalise weights '''
    for i in range(0, len(weights)):
        sum_per_star = sum(weights[i])
        weights[i] = weights[i]/sum_per_star
        
    periods = np.zeros(len(stars))

    ''' find the weighted mean and add errors in quadrature '''
    for i in range(0, len(all_periods)):
        periods[i] = np.average(all_periods[i], weights = weights[i])
        errors[i] = (errors[i])**2
        err_sum = sum(errors[i])
        errors[i] = np.sqrt(err_sum)

    ''' Print and plot results '''
    p.close(2)
    p.figure(2)

    true_periods = []; star_list = []; period_list = []
    for i in range(0, len(periods)):
        print 'Star = ', stars[i]
        
        print 'Measured period = ', periods[i], '+/-', errors[i]
        trueperiod = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/sim_period%stest.txt' \
                                  %(int(stars[i])))
        period = trueperiod[4]
        #p.plot(i+1, period, 'b.')
        #p.legend()
        print 'True period = ', period
        true_periods.append(period)
        star_list.append(stars[i])
        period_list.append(periods[i])

    x = np.arange(0,40,0.1)
    p.plot(x, x, 'b--')
    p.plot(.5*x, x, 'b--')
    p.plot(.333333*x, x, 'b--')
    p.plot(2*x, x, 'b--')
    p.plot(3*x, x, 'b--')
    p.errorbar(periods, true_periods, yerr = errors, fmt = 'r.')
    p.xlabel('True Period (days)')
    p.ylabel('Measured Period (days)')
    #p.xlim(min(periods) - 0.5, max(periods) + 0.5)
    #p.ylim(min(true_periods) - 1.5, max(true_periods) + 1.5)
    p.xlim(0, 40)
    p.ylim(0, 40)
    p.title('Measured vs true period')
    p.savefig('/Users/angusr/angusr/ACF/star_spot_sim/resultstest')
    np.savetxt('/Users/angusr/angusr/ACF/star_spot_sim/measured_vs_true.txt', \
                  (star_list, true_periods, period_list))

    return period_list
#-----------------------------------------------------------------------------------------------------------------    
# This function calculates the completeness: percentage of sample with successful period measurement,
# the reliability: fraction of objects with a given true period, detected within < 20 per cent of this value,
# the contamination: fraction of objects with measured period > 20 per cent different to the true period.
#-----------------------------------------------------------------------------------------------------------------

def completeness(nstars):

    # Load measured vs true stats
    data = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/measured_vs_true.txt')
    x = np.isfinite(data[2])
    period_list = data[2][x]

    print 'Total number of simulated stars = %s' %nstars
    print 'Total number of measured periods = %s' %len(period_list)
    print 'CALCULATING COMPLETENESS...'

    ''' Bin stars according to period '''
    nbins = 10
    true_bins = np.ones(nbins)
    success_bins = numpy.ones(nbins)
    b_space = 2./float(nbins)
    bdries = []
    for k in range(1, nbins):
            bdries.append(k**(b_space*k))
    bdries.append(100.)

    true_periods = []
    success_period = []
    measured_periods = []
    errors = []

    nstars = 24
    for i in range(1, nstars):
        # Load data
        data = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/grid/%sparams.txt' %(i+1)).T
        kid_x = data[0]; KIC_no = str(data[1]); Amp = data[2]; nspots = data[3]; ff = data[4]; \
            amp = data[5]; noise = data[6]; true_period = data[7]; tau = data[8]

        # bin
        for m in range(len(bdries)-1):
            if bdries[m] < true_period < bdries[m+1]:
                true_bins[m] += 1

        # Was its period successfully measured?
        list_of_success = np.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/ss_ind_quarterstest.txt')
        for j in range(len(list_of_success)):
            if kid_x == list_of_success[j]:
                 success_period.append(true_period) # Do I want the true period or the measured period here?

                 
                 ''' Load measured period (for now just take quarter 3)'''
                 #data = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQss3_output/results.txt').T
                 #measured_periods.append(data[1][1:][kid_x-1])
                 #errors.append(data[6][1:][kid_x-1])
                 for k in range(3, 14):
                     data = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQss%s_output/results.txt' %k).T
                     measured_periods.append(data[1][1:][kid_x-1])
                     errors.append(data[6][1:][kid_x-1])
                 m_periods, m_errors = weighted_mean(measured_periods, errors)
    print 'm_periods', m_periods

    # re-bin
    for m in range(len(bdries)-1):
        for n in range(0, len(success_period)):
            if bdries[m] < success_period[n] < bdries[m+1]:
                success_bins[m] += 1    

    print 'Truth = ', true_bins
    print 'Successfully measured = ', success_bins

    p.close(1)
    p.figure(1)
    p.ylabel('True (number)')
    step(np.log10(bdries), (success_bins/true_bins)*100)
    p.ylim(0,110)

    print 'CALCULATING RELIABILITY...'
    ''' Load measured period '''
    m_period = []; m_err = []
    for j in range(3, 14):
        data = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQss%s_output/results.txt' %j).T
        periods = data[1][1:][i]
        errors = data[6][1:][i]
    

    return

#-----------------------------------------------------------------------------------------------------------------
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

    print all_periods, all_errors
    
    for i in range(len(all_errors)):
        if all_errors[i] == 0.0:
            all_errors[i] = 0.01
    errors = all_errors
    for i in range(len(all_errors)):
        all_errors[i] = 1.0/((all_errors[i])**2)
    weights = all_errors
    weights = weights/sum(weights) # Normalise
    period = np.average(all_periods, weights = weights) # Find weighted mean
    for i in range(len(errors)):
        errors[i] = errors[i]**2 # Add errors in quadrature
    err_sum = sum(errors)
    error = np.sqrt(err_sum)

    print period, error
    return period, error

#-----------------------------------------------------------------------------------------------------------------
