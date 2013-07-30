# Wrapper for testing period measurement method
# This version does individual quarters
# This version produces star spots from the entire Kepler lc, then adds to individual quarters
# This version produces the lightcurves to be added to real lcs later
# Spots 5 does all 24 stars
# Spots 6 names stars as numbers not strings

#import matplotlib
#matplotlib.use('Agg')
import random
import numpy
import spot_sim_ruth5
import scipy
import glob
import ACF_star_spot7
import ss_period_plots3
import ss_index2
import pylab
import scipy.io
import pyfits

# This is the master function that calls everything below.

#======================================================
def run_master():
    spot_gen(10,[5])
    # # ####spot_gen(8, [5,2,1])           ## (8 x number of periods, 3 x spot lifetimes = 24 light curves) 
    add_to_real(1, 1)
    # ####add_to_real(10, 4)                       ## (inject into 10 Kepler light curves with 4 differen-t
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
    # run_ACF(10)                                     ## (Calculate ACF for 1 lc) 1 = 1 x 1 x 1 x 1
    ####run_ACF(960)                              ## (Calculate ACFs for all 960 lcs) 960 = 3 x 8 x 10 x 4
    recording_period_measurements(10)
    ####recording_period_measurements(960)               ## (make note of the periods measured)
    period_plots(10)
    ####period_plots(960)                                ## (Produce period results for each
    # quarter)
    compare()                                        ## (Compare true vs measured periods)
    return                                             ##
#======================================================

 

    
# Generate false light curves and saves them within star_spot_sim. The light
# curves are numbered with 2 digits. The first is corresponds to the star
# number and the second to the quarter (qs 3-6 = 1, qs 7-10 = 2, etc).
# It saves a text file containing the star's period as /ACF/star_spot_sim/sim_period*.txt

#-----------------------------------------------------------------------------------------------------------------
def spot_gen(no_periods, taus):
    
    print 'Generating light curves'
    ROOTDIR = '/Users/angusr/angusr/ACF/star_spot_sim'

    
    n = 0
    for i in range(0, len(taus)):
        # Generate light curves with 3 values of tau
        tau = taus[i]
        for index in range(0,no_periods):
            # Create lightcurves for 24 stars: 3 spot lifetimes and 8 periods.
            # Generate random values for period that are uniform in log space.
        
            print 'Star %s...' %(index+1)
            
            ''' Try generating numbers from a distribution that is uniform in log space '''
            #myperiod = [(numpy.sqrt(random.gauss(3., .5))**2)]
            random_gen = numpy.random.uniform(low = 0.0, high = 2.0)
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
            grid = [grid[0], grid[1], grid[2], grid[3], save_period_in_txt_file]
            ''' grid contains [number of spots, filling factor, amplitude, noise, period] '''
            numpy.savetxt('/Users/angusr/angusr/ACF/star_spot_sim/sim_period%stest.txt' %n, grid)
            n += 1
    return
#-------------------------------------------------------------------------------------------

def add_to_real(number_of_kplr_lcs, number_of_amps):

    counter = 1

    ''' Setting number of Kepler lcs '''
    if number_of_kplr_lcs == 10:
        names = ['010972873', '011137075', '011502597', '011717120', '005607242', \
             '006370489', '006442183', '006531928', '006603624', '009574283']
    elif number_of_kplr_lcs == 1:
        names = ['010972873']

    for n in range(0, len(names)):
        print 'Kepler lc = ', names[n]
        files = glob.glob('/Users/angusr/angusr/data2/all_Qs/kplr%s-*llc.fits' % (names[n]))
        quarter_times_start = []; quarter_times_fin = []
        print len(files)


        ''' Setting number of amplitudes '''
        if number_of_amps == 4:
            Amps = [1.0, 5.0, 10.0, 15.0]
        elif number_of_amps == 1:
            Amps = [1.0]
            
        for a in range(0, len(Amps)):
            print 'Amplitude = ', Amps[a]
            
            for star in range(1, number_of_kplr_lcs+1):
            #for star in range(1, 4):
                print 'Star = ', star

                for quarter in range(0,12):
                    print 'quarter = ', quarter

                    '''Load real data'''
                    time2 = []; lightcurve = []
                    hdulist = pyfits.open(files[quarter])
                    tbdata = hdulist[1].data #"add the first extension to tbdata"
                    x = numpy.where(numpy.isfinite(tbdata['TIME'])) #Remove NANs
                    time = tbdata['TIME'][x]; lc = tbdata['PDCSAP_FLUX'][x]
                    x = numpy.where(numpy.isfinite(tbdata['PDCSAP_FLUX']))
                    time = tbdata['TIME'][x]; lc = tbdata['PDCSAP_FLUX'][x]
                    quarter_times_start.append(time[0]); quarter_times_fin.append(time[-1])
                    lightcurve = list(lightcurve)
                    lightcurve.extend(numpy.array(lc)); time2.extend(numpy.array(time))
                    time = numpy.array(time2- min(time2))
                    addit = min(time2)
                    time = time + addit # Add real starting value back to time array
                    #print 'start = ', quarter_times_start[quarter]
                    #print 'stop = ', quarter_times_fin[quarter]
                    x = numpy.where(time == quarter_times_start[quarter])
                    y = numpy.where(time == quarter_times_fin[quarter])
                    start = int(x[0]); stop = int(y[0])
                    time = time[start:stop] 
                    lightcurve = lightcurve[start:stop] # Slice real data set into individual quarters
                    lightcurve = lightcurve/numpy.median(lightcurve)
                    lightcurve = lightcurve - numpy.median(lightcurve)

                

                    '''Load simulated data'''
                    print quarter+3, star
                    mat = scipy.io.loadmat('/Users/angusr/angusr/ACF/star_spot_sim/%s/sim_%s.png.mat'\
                    #mat = scipy.io.loadmat('/Users/angusr/angusr/ACF/star_spot_sim/%s/orig_sim_%s.png.mat'\
                    #mat = scipy.io.loadmat('/Users/angusr/angusr/ACF/star_spot_sim/%s/testtesttest%s.png.mat'\
                                           %((quarter+3), star))
                    pars = mat['pars']; ts = mat['ts']
                    lc2 = mat['ts']
                    sim_time = lc2[0]; sim_lc =  lc2[2]

                    '''Plotting'''
                    # pylab.close(1)
                    # pylab.figure(1)
                    # pylab.subplot(3,1,1)
                    # pylab.plot(sim_time, sim_lc, 'r.')
                    # pylab.subplot(3,1,2)
                    # pylab.plot(time, lightcurve, 'b.')
                    # pylab.subplot(3,1,3)

                    ''' Checking that the simulated and real quarters have the same number \
                    of data points - the time stamps come from one particular lc, so there are mismatches!'''
                    if len(time) > len(sim_time):
                        for mis in range(len(time) - len(sim_time)):
                            print 'len(time) > len(sim_time)'
                            time = list(time); lightcurve = list(lightcurve)
                            time.remove(time[-1]); lightcurve.remove(lightcurve[-1])
                            time = numpy.array(time); lightcurve = numpy.array(lightcurve)
                    if len(sim_time) > len(time):
                        for mis in range(len(sim_time) - len(time)):
                            print 'len(sim_time) > len(time)'
                            time = list(time); lightcurve = list(lightcurve)
                            time.append(time[-1] + (time[1] - time[0])); lightcurve.append(lightcurve[-1])
                            time = numpy.array(time); lightcurve = numpy.array(lightcurve)

                    '''Add simulated lc to real lc, where Amp is the scaling relation'''
                    new_lc = Amps[a]*sim_lc + lightcurve

                    # pylab.plot(time, new_lc, 'g.')

                    ts[2] = new_lc

                    ''' Save back into mat file format'''
                    scipy.io.savemat('/Users/angusr/angusr/ACF/star_spot_sim/%s/sim_%s.png.mat' \
                                        %((quarter+3), counter), {'pars': pars, 'ts': ts})
                    #raw_input('enter')
                numpy.savetxt('/Users/angusr/angusr/ACF/star_spot_sim/grid/KIC_amp_%s.txt' %counter, \
                                      (numpy.transpose((counter, int(names[n]), Amps[a]))))
                counter += 1
                        

    print '%s iterations' %counter
    return



# Run ACF on new light curves. It takes the light curves from /ACF/star_spot_sim/3_6 etc. It saves the
# results to /ACF/PDCQss3-6 etc. It creates a file called results.txt which is then converted into a file called
# Periods_3-6.txt or similar.

#-----------------------------------------------------------------------------------------------------------------
def run_ACF(no_stars):
    # Running ACF code on false light curves
    test_list = ['ss3','ss4','ss5', 'ss6', 'ss7', 'ss8', 'ss9', 'ss10', 'ss11', 'ss12', 'ss13', 'ss14']
    test_list2 = ['3','4','5', '6', '7', '8', '9', '10', '11', '12', '13', '14']
    for a in range(0,len(test_list)):
        ACF_star_spot7.corr_run(test_list[a], test_list2[a], number_of_stars = no_stars+1)
    return
#-----------------------------------------------------------------------------------------------------------------



# This just takes the key results saved in results.txt and moves them into a file called Periods_ss3-6.txt etc

#-----------------------------------------------------------------------------------------------------------------
def recording_period_measurements(no_stars):
    print 'Recording period measurements... '
    list_of_quarters = ['ss3','ss4','ss5', 'ss6', 'ss7', 'ss8', 'ss9', 'ss10', 'ss11', 'ss12', 'ss13', 'ss14']
    for year in list_of_quarters:
        data = numpy.genfromtxt('/Users/angusr/angusr/ACF/PDCQ%s_output/results.txt' %year).T
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
        numpy.savetxt('/Users/angusr/angusr/ACF/PDCQ%s_output/Periods_%stest.txt' %(year,year), \
            numpy.transpose((KID, ACF_period, period_error, sine_period))) 
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
    stars = numpy.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/ss_ind_quarterstest.txt')
    
    which_quarter = ['ss3','ss4','ss5', 'ss6', 'ss7', 'ss8', 'ss9',\
                     'ss10', 'ss11', 'ss12', 'ss13', 'ss14']
    
    all_periods = []
    all_errors = []
    for star in range(0,len(stars)):
        all_periods.append([])
        all_errors.append([])

    ''' Reading in results '''
    for year in range(0, len(which_quarter)):
        # data = numpy.genfromtxt('/Users/angusr/angusr/ACF/PDCQ%s_output/Periods_%s.txt' \
                                # %(which_quarter[year], which_quarter[year])).T
        data = numpy.genfromtxt('/Users/angusr/angusr/ACF/PDCQ%s_output/Periods_%stest.txt'\
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
        
    periods = numpy.zeros(len(stars))

    ''' find the weighted mean and add errors in quadrature '''
    for i in range(0, len(all_periods)):
        periods[i] = numpy.average(all_periods[i], weights = weights[i])
        errors[i] = (errors[i])**2
        err_sum = sum(errors[i])
        errors[i] = numpy.sqrt(err_sum)

    ''' Print and plot results '''
    pylab.close(2)
    pylab.figure(2)

    true_periods = []; star_list = []; period_list = []
    for i in range(0, len(periods)):
        print 'Star = ', stars[i]
        
        print 'Measured period = ', periods[i], '+/-', errors[i]
        trueperiod = numpy.genfromtxt('/Users/angusr/angusr/ACF/star_spot_sim/sim_period%stest.txt' \
                                  %(int(stars[i]-1)))
        period = trueperiod[4]
        #pylab.plot(i+1, period, 'b.')
        #pylab.legend()
        print 'True period = ', period
        true_periods.append(period)
        star_list.append(stars[i])
        period_list.append(periods[i])

    x = numpy.arange(0,10,0.1)
    pylab.plot(x,x,'b--')
    pylab.errorbar(periods, true_periods, yerr = errors, fmt = 'r.')
    pylab.xlabel('True Period')
    pylab.ylabel('Measured Period')
    pylab.xlim(min(periods) - 0.5, max(periods) + 0.5)
    pylab.ylim(min(true_periods) - 1.5, max(true_periods) + 1.5)
    pylab.title('Measured vs true period')
    pylab.savefig('/Users/angusr/angusr/ACF/star_spot_sim/resultstest')
    print star_list
    print true_periods
    print period_list
    numpy.savetxt('/Users/angusr/angusr/ACF/star_spot_sim/measured_vs_true.txt', \
                  (star_list, true_periods, period_list))
   
    return
#-----------------------------------------------------------------------------------------------------------------    













