# This is The version of the code that runs over all the subsets (or just those defined in read_in.py)
import scipy
from numpy.random import normal
import matplotlib.image as mpimg
import random
import numpy
import atpy
import pylab
import copy
import glob
from sets import Set
import collections


from scipy import signal
import KOI_tools_b12 as kt
import filter
import gls
import mpfit

#if rpy not available, comment out this..
# from rpy2 import robjects as ro
# r = ro.r
# import rpy2.robjects.numpy2ri
# rpy2.robjects.numpy2ri.activate()
# and set no_rpy to True...
no_rpy = True


dir = '/Users/angusr/angusr/ACF'
gap_days = 0.02043365  # assume for long cadence
jump_arr = scipy.array([131.51139, 169.51883, 169.75000, 182.00000, 200.31000, 231.00000, 246.19000, 256.00000, 260.22354, 281.00000, 291.00000, 322.00000, 352.37648, 373.23000, 384.00000, 398.00000, 443.48992, 475.50000, 504.00000, 539.44868, 567.00000, 599.00000, 630.17387, 661.00000, 691.00000, 711.20000, 735.36319, 762.00000, 808.51558, 845.00000, 874.50000, 906.84469, 937.00000, 970.00000, 1001.20718, 1032.50000, 1063.50000 ,1071.00000, 1093.60000])

qmin = 1
qmax = 16


###############################################################################################################

def corr_run(quarter_sel, tr_out = False, tr_type = None):
    ''' Reads index and light curves and runs the ACF calculation modules
        tr_out = True if transits are to be removed
        tr_type = EB if EB not planet'''
    #read in index and create list of unique KIDs
    index = atpy.Table('%s/index_Q%s.ipac' %(dir,quarter_sel), type = 'ipac')
    index = index.where(index.long_cadence == 1)
    
    #index = index.rows(scipy.r_[0:10:1])
    print 'Index includes %s files...'% len(index)
    id_list = scipy.array(list(Set(index.keplerid)))
    print '...and %d individual stars' % len(id_list)
    id_list.sort()
    print quarter_sel

    # Sorting out quarters for the table
    if 2 < quarter_sel < 17:
        t = index.where(index.quarter == int(quarter_sel))
    elif quarter_sel == '0-4':
        t = index.where(index.quarter == 1)
    elif quarter_sel == '0-9':
        t = index.where(index.quarter == 1)
    elif quarter_sel == '5-9':
        t = index.where(index.quarter == 5)
    elif quarter_sel == '2-5':
        t = index.where(index.quarter == 2)
    elif quarter_sel == '3-6':
        t = index.where(index.quarter == 3)
    elif quarter_sel == '4-9':
        t = index.where(index.quarter == 4)
    elif quarter_sel == '7-10':
        t = index.where(index.quarter == 7)
    elif quarter_sel == '11-14':
        t = index.where(index.quarter == 11)
    t.sort('keplerid')
    print t
    if len(t) != len(id_list):
        print len(t), len(id_list)
        print 'Error: table length != KID list length'
        return
    #raw_input('enter')


    # Create empty arrays
    acf_peak_per = scipy.ones(len(id_list)) * -9999.0
    sine_per = scipy.ones(len(id_list)) * -9999.0
    sine_height = scipy.ones(len(id_list)) * -9999.0
    med_dlag_per = scipy.ones(len(id_list)) * -9999.0
    dlag_per_err = scipy.ones(len(id_list)) * -9999.0
    h1 = scipy.ones(len(id_list)) * -9999.0
    hlocgrad = scipy.ones(len(id_list)) * -9999.0
    hloc_grad_scatter = scipy.ones(len(id_list)) * -9999.0
    width_grad = scipy.ones(len(id_list)) * -9999.0
    width_grad_scatter = scipy.ones(len(id_list)) * -9999.0
    w1 = scipy.ones(len(id_list)) * -9999.0
    lh1 = scipy.ones(len(id_list)) * -9999.0
    num_of_peaks = scipy.ones(len(id_list)) * -9999.0
    harmonic_det = scipy.ones(len(id_list)) * -9999.0
    amp_all = scipy.ones(len(id_list)) * -9999.0
    amp_per = scipy.ones(len(id_list)) * -9999.0
    period = scipy.ones(len(id_list)) * -9999.0

    # loop over KIDs
    for x in scipy.arange(len(id_list)):
        kid_x = int(id_list[x])
        print 'x = %d/%d, kid = %d' %(x, len(id_list), kid_x)

        # Load light curve (load_lc returns normalised, quarter-joined LC)
        print 'Loading LC...'
        condition, lc_tab, qt_max, tablen = load_lc(quarter_sel, kid_x, index, tr_out, tr_type)
        if condition == False: continue
        print 'qt_max', qt_max

        print lc_tab, len(lc_tab)
        #raw_input('enter')

        print 'Time interval =', lc_tab.time[0]-lc_tab.time[1]
        print min(lc_tab.time)

        print 'Plotting raw flux'
        pylab.clf()
        pylab.plot(lc_tab.time,lc_tab.flux)
        pylab.title('Raw flux')
        #raw_input('enter')

        # main ACF figure (raw, corr, amp, acf) 
        pylab.figure(1,(12, 9))
        pylab.clf()
        pylab.subplot(4,1,1)
        pylab.title('ID: %s' %(kid_x), fontsize = 16)
        for j in scipy.arange(tablen):
            if j % 2 == 0:
                pylab.axvspan(qt_max[j], qt_max[j+1], facecolor = 'k', alpha=0.1)
        pylab.plot(lc_tab.time, lc_tab.flux_pdc, 'k-')
        for k in scipy.arange(len(jump_arr)):
            pylab.axvline(jump_arr[k], ls = '--', c = 'b')
        pylab.xlim(lc_tab.time.min(), lc_tab.time.max())
        pylab.ylim(min(lc_tab.flux_pdc[scipy.isfinite(lc_tab.flux_pdc) == True]), \
                   max(lc_tab.flux_pdc[scipy.isfinite(lc_tab.flux_pdc) == True]))
        pylab.ylabel('Raw Flux')
        pylab.subplot(4,1,2)
        pylab.plot(lc_tab.time, lc_tab.flux, 'k-')
        pylab.xlim(lc_tab.time.min(), lc_tab.time.max())
        pylab.ylim(lc_tab.flux.min(), lc_tab.flux.max())
        pylab.ylabel('Norm Flux')

        pylab.figure(2,(12, 9))
        pylab.clf()
        pylab.subplot(3,1,1)
        pylab.title('ID: %s' %(kid_x), fontsize = 16)
        for j in scipy.arange(tablen):
            if j % 2 == 0:
                pylab.axvspan(qt_max[j], qt_max[j+1], facecolor = 'k', alpha=0.1)
        pylab.plot(lc_tab.time, lc_tab.flux_pdc, 'k-')
        for k in scipy.arange(len(jump_arr)):
            pylab.axvline(jump_arr[k], ls = '--', c = 'b')
        pylab.xlim(lc_tab.time.min(), lc_tab.time.max())
        pylab.ylim(min(lc_tab.flux_pdc[scipy.isfinite(lc_tab.flux_pdc) == True]), \
                   max(lc_tab.flux_pdc[scipy.isfinite(lc_tab.flux_pdc) == True]))
        pylab.ylabel('Raw Flux')
        pylab.subplot(3,1,2)
        pylab.plot(lc_tab.time, lc_tab.flux, 'k-')
        pylab.xlim(lc_tab.time.min(), lc_tab.time.max())
        pylab.ylim(lc_tab.flux.min(), lc_tab.flux.max())
        pylab.ylabel('Norm Flux')
        ax = pylab.gca()
        pylab.text(0.415, -0.15, 'Time (days)', transform = ax.transAxes)
        pylab.text(0.415, -1.4, 'Period (days)', transform = ax.transAxes)

        # max period searched for is len(flux) / 2
        max_psearch_len = len(lc_tab.flux) / 2.0
        
        # Calculate ACF
        print 'Calculating ACF...'

        acf_tab, acf_per_pos, acf_per_height, acf_per_err, locheight, asym,  = \
            acf_calc(quarter_sel, time = lc_tab.time, flux = lc_tab.flux, interval = gap_days, kid = kid_x, max_psearch_len = max_psearch_len)

        pgram_tab, sine_per[x], sine_height[x] = \
            pgram_calc(quarter_sel, time = lc_tab.time, flux = lc_tab.flux, interval = gap_days, kid = kid_x, max_psearch_len = max_psearch_len)


        pylab.figure(1)
        pylab.subplot(4,1,4)
        pylab.plot(acf_tab.lags_days, acf_tab.acf_smooth, 'k-')
        pylab.axhline(0, ls = '--', c = 'k')
        for i in scipy.arange(len(acf_per_pos)):
            pylab.axvline(acf_per_pos[i], ls = '--', c = 'b')
        pylab.ylabel('ACF')
        pylab.xlim(0, acf_tab.lags_days.max())

        pylab.figure(2)
        pylab.subplot(3,1,3)
        pylab.plot(pgram_tab.period, pgram_tab.pgram, 'k-')
        pylab.axvline(sine_per[x], c = 'r', ls = '--')
        pylab.ylim(0, 1.1*max(pgram_tab.pgram))
        pylab.axhline(0.3, c = 'g', ls = '--')
        pylab.xlim(0, pgram_tab.period.max())
        pylab.ylabel('Periodogram')
        pylab.savefig('%s/PDCQ%s_output/plots_pgram/%s_pgram.png' % (dir, quarter_sel, kid_x))

        pylab.figure(12)
        pylab.clf()
        pylab.plot(acf_tab.lags_days, acf_tab.acf_smooth, 'k-')
        for m in scipy.arange(len(acf_per_pos)):
            pylab.axvline(acf_per_pos[m], ls = '--', c = 'r')
        

        # plot and calculate acf peak statistics
        if acf_per_pos[0] != -9999:
            med_dlag_per[x], dlag_per_err[x], acf_peak_per[x], h1[x], w1[x], lh1[x], \
                hlocgrad[x], hloc_grad_scatter[x], width_grad[x], width_grad_scatter[x], \
                num_of_peaks[x], harmonic_det[x], sel_peaks, one_peak_only, peak_ratio =\
                plot_stats(quarter_sel, lc_tab.time, lc_tab.flux, kid_x, acf_per_pos, acf_per_height, acf_per_err, locheight, asym)

 # # plot and calculate acf peak statistics
 #        if acf_per_pos[0] != -9999:
 #            med_dlag_per[x], dlag_per_err[x], acf_peak_per[x], h1[x], w1[x], lh1[x], \
 #                hlocgrad[x], hloc_grad_scatter[x], width_grad[x], width_grad_scatter[x], \
 #                num_of_peaks[x], harmonic_det[x], sel_peaks, one_peak_only =\
 #                plot_stats(quarter_sel, lc_tab.time, lc_tab.flux, kid_x, acf_per_pos, acf_per_height, acf_per_err, locheight, asym)


            #if med_dlag_per[x] >0: period[x] = med_dlag_per[x]
            #else: period[x] = acf_peak_per[x]

            if max(acf_per_height[:2]) > 0:
                if med_dlag_per[x] >0: period[x] = med_dlag_per[x] #################################
                else: period[x] = acf_peak_per[x]
            else:
                print '!!!!!! PEAK HEIGHT < 0 !!!!!!'
                #period[x] = -9999.0
                period[x] = 0.0
            if one_peak_only == 1:
                print '!!!!!! ONLY ONE PEAK FOUND < 0.1 !!!!!!'
                #period[x] = -9999.0
                period[x] = 0.0
            if lh1[x] < 0.1:
                print '!!!!!! LOCAL PEAK HEIGHT < 0.1 !!!!!!'
                #period[x] = -9999.0
                period[x] = 0.0

            print 'PEAK RATIO = ', peak_ratio
            

            # plot period lines on full plot
            pylab.figure(1)
            pylab.subplot(4,1,4)
            if med_dlag_per[x] >0:
                for n in scipy.arange(len(sel_peaks)):
                    pylab.axvline(sel_peaks[n], ls = '--', c = 'r')
            pylab.axvline(period[x], ls = '--', c = 'k')
            if med_dlag_per[x] > 0:
                pylab.axvspan(med_dlag_per[x]-dlag_per_err[x], med_dlag_per[x]+dlag_per_err[x],\
                              facecolor = 'k', alpha=0.2)

            # variability stats
            print 'calculating var for P_med...'
            amp_all[x], amp_per[x], per_cent, var_arr_real = \
                calc_var(quarter_sel, kid = kid_x, time_in = lc_tab.time, flux = lc_tab.flux, period = acf_peak_per[x])

            pylab.figure(1)
            pylab.subplot(4,1,3)
            pylab.plot(per_cent, var_arr_real, 'k.')
            pylab.axhline(amp_per[x], ls = '--', c = 'b')
            pylab.xlim(lc_tab.time.min(),lc_tab.time.max())
            pylab.ylim(var_arr_real.min(), var_arr_real.max())
            ax = pylab.gca()
            pylab.text(0.415, -0.15, 'Time (days)', transform = ax.transAxes)
            pylab.text(0.415, -1.4, 'Period (days)', transform = ax.transAxes)
            pylab.ylabel('Amplitudes')
            pylab.savefig('%s/PDCQ%s_output/plots_acf/%s_full.png' % (dir, quarter_sel, kid_x))

            # make zoomed in plot to show periodicity
            pylab.figure(5,(12, 9))
            pylab.clf()
            pylab.subplot(3,1,1)
            pylab.title('Zoom of ID: %s, showing peak period and av lag period' %kid_x)
            maxpts = 40.0
            if scipy.floor(lc_tab.time.max() / acf_peak_per[x]) < maxpts:
                maxpts = float(scipy.floor(lc_tab.time.max() / acf_peak_per[x]))
            inc = lc_tab.time - lc_tab.time.min() <= (maxpts*acf_peak_per[x])
            pylab.plot(lc_tab.time[inc], lc_tab.flux[inc], 'k-')
            for i in scipy.arange(int(maxpts)):
                pylab.axvline(lc_tab.time.min() + (i)*(acf_peak_per[x]), ls = '--', c = 'b')        
            pylab.xlim(lc_tab.time[inc].min(), lc_tab.time[inc].max())
            pylab.ylabel('Flux (with Peak Pos)')
            pylab.subplot(3,1,2)
            pylab.plot(lc_tab.time[inc], lc_tab.flux[inc], 'k-')
            for i in scipy.arange(int(maxpts)):
                pylab.axvline(lc_tab.time.min() + (i)*(med_dlag_per[x]), ls = '--', c = 'r')        
            pylab.xlim(lc_tab.time[inc].min(), lc_tab.time[inc].max())
            pylab.ylabel('Flux (with Med Pos)')
            pylab.xlabel('Time (days)')
            pylab.subplot(3,1,3)
            pylab.subplots_adjust(top=0.96, bottom = 0.06)
            print acf_tab.lags_days[acf_tab.lags_days <= (maxpts*acf_peak_per[x])]
	    print len(acf_tab.lags_days[acf_tab.lags_days <= (maxpts*acf_peak_per[x])])
	    print acf_tab.acf_smooth[acf_tab.lags_days <= (maxpts*acf_peak_per[x])]
	    print len(acf_tab.acf_smooth[acf_tab.lags_days <= (maxpts*acf_peak_per[x])])
            pylab.plot(acf_tab.lags_days[acf_tab.lags_days <= (maxpts*acf_peak_per[x])], \
                       acf_tab.acf_smooth[acf_tab.lags_days <= (maxpts*acf_peak_per[x])])
            pylab.ylabel('ACF')
            pylab.xlabel('Lag (days)')
            pylab.axvline(acf_peak_per[x], ls = '--', c = 'b')
            pylab.axvline(med_dlag_per[x], ls = '--', c = 'r')
            if acf_tab.lags_days.max() > 10 * acf_peak_per[x]:
                pylab.xlim(0,10 * acf_peak_per[x])
            else: pylab.xlim(0,max(acf_tab.lags_days))
            pylab.savefig('%s/PDCQ%s_output/plots_z/%s_zoom.png' % (dir, quarter_sel, kid_x))

            print '**************************', 'KID = ', kid_x, 'PEAK HEIGHT = ', max(acf_per_height[:2]), 'LOCAL PEAK HEIGHT = ', lh1[x]

        # save file of stats for individual star
        t_stats = atpy.Table()
        t_stats.add_column('acf_per_pos', acf_per_pos)
        t_stats.add_column('acf_per_height', acf_per_height)
        t_stats.add_column('acf_per_err', acf_per_err)
        t_stats.add_column('asym', asym)
        t_stats.add_column('locheight', locheight)
        savefilen = '%s/PDCQ%s_output/dat_files/%s_stats.txt' % (dir, quarter_sel, kid_x)
        atpy.asciitables.write_ascii(t_stats, savefilen, delimiter = ' ', overwrite = True)

        print 'PERIOD = ', period[x]
        
        #raw_input('enter')

    # save stats that are 1 per kid
    t.remove_columns(['filename'])
    t.add_column('period', period) #period
    t.add_column('sine_per', sine_per) #sine period
    t.add_column('sine_height', sine_height)
    t.add_column('acf_peak_per', acf_peak_per)
    t.add_column('med_dlag_per', med_dlag_per)
    t.add_column('dlag_per_err', dlag_per_err) #error
    t.add_column('h1', h1)
    t.add_column('w1', w1)
    t.add_column('lh1', lh1)
    t.add_column('hlocgrad', hlocgrad)
    t.add_column('hloc_grad_scatter', hloc_grad_scatter)
    t.add_column('width_grad', width_grad)
    t.add_column('width_grad_scatter', width_grad_scatter)
    t.add_column('num_of_peaks', num_of_peaks)
    t.add_column('harmonic_det', harmonic_det)
    t.add_column('amp_all', amp_all)
    t.add_column('amp_per', amp_per)
    savefilen = '%s/PDCQ%s_output/results.txt' % (dir,quarter_sel)
    atpy.asciitables.write_ascii(t, savefilen, delimiter = ' ', overwrite = True)
    
    return

###############################################################################################################

def load_lc(quarter_load, kid_x = None, index = None, tr_out = False, tr_type = False):
    # read in tableset of LCs for all Qs, normalise and quarter join
    dir_lc = '/Users/angusr/angusr/data2/Q%s_public' %quarter_load + '/Light_Curves'
    print 'Reading in LCs...'
    tset, found = kt.GetLc(id = kid_x, dir = dir_lc, tr_out = tr_out, \
                           filename = index.filename[(index.keplerid == kid_x) * (index.quarter <= qmax) * (index.quarter >= qmin)])
    if found == 0:
        print '****not found****'
        return False, atpy.Table(), 0, 0
    
    # find number of Quarters and max(time) per quarter
    if tr_out == False: lc = tset[0]
    else: lc = tset[1]
    quarters = scipy.array(list(set(lc.Q)))
    tablen = len(quarters)
    qt_max = scipy.zeros(tablen)
    for z in scipy.arange(len(quarters)):
        qt_max[z] = max(lc.TIME[lc.Q == quarters[z]])
    qt_max = scipy.append(0,qt_max)
    finite = scipy.where((scipy.isfinite(lc.TIME) == True) * (scipy.isfinite(lc.SAP_FLUX) == True) * \
                         (scipy.isfinite(lc.PDCSAP_FLUX) == True))

    # keep unchanged set for table
    flux_base = lc.SAP_FLUX[finite]
    q_arr_base = lc.Q[finite]
    for qt in quarters:
        flux_base[q_arr_base == qt] = (flux_base[q_arr_base == qt] / scipy.median(flux_base[q_arr_base == qt])) - 1.0

    time = lc.TIME
    flux_raw = lc.SAP_FLUX
    flux = lc.PDCSAP_FLUX
    
    if tr_out == True:
        print 'Removing eclipses...'    
        # remove in-eclipse regions
        phase, inTr = kt.TransitPhase(tset)
        npl, nobs = inTr.shape
        for ipl in scipy.arange(npl):
            trlist = inTr[ipl,:].astype(bool)
            time[trlist] = scipy.nan
            flux_raw[trlist] = scipy.nan
            flux[trlist] = scipy.nan

    if tr_type == 'EB':
        print 'Removing EB transits...'
        phase, inTr = EBTransitPhase(tset, kid_x)
        npl, nobs = inTr.shape
        for ipl in scipy.arange(npl):
            trlist = inTr[ipl,:].astype(bool)
            time[trlist] = scipy.nan
            flux_raw[trlist] = scipy.nan
            flux[trlist] = scipy.nan
        pylab.plot(time, flux, 'r.')

    finite_tr = scipy.where((scipy.isfinite(time) == True) * \
                            (scipy.isfinite(flux) == True) * (scipy.isfinite(flux_raw) == True) * (scipy.isfinite(lc.Q) == True))
    time = time[finite_tr]
    flux = flux[finite_tr]
    flux_raw = flux_raw[finite_tr]
    quarter_arr = lc.Q[finite_tr]

    # median normalise each quarter
    for qt in quarters:
        flux[quarter_arr == qt] = (flux[quarter_arr == qt] / scipy.median(flux[quarter_arr == qt])) - 1.0
        flux_raw[quarter_arr == qt] = (flux_raw[quarter_arr == qt] / scipy.median(flux_raw[quarter_arr == qt])) - 1.0
    
    # fit 1d polynomial and remove >5sig outliers 
    coeff = scipy.polyfit(time,flux, 1) 
    p = scipy.poly1d(coeff) 
    med, sig = filter.medsig(flux - p(time))
    sigclip = abs(flux - p(time)) <= 5*sig
    time = time[sigclip]
    flux = flux[sigclip]
    flux_base = flux_base[sigclip]

    # to get noise properties
    flux_filt = filter.filt1d(flux, 10, 5)
    med_noise, sig_noise = filter.medsig(flux-flux_filt)

    # find gaps greater than 1 (1.1) data points and linear interp with noise
    diff1 = time[1:] - time[:-1]
    diff1 = scipy.append(diff1, gap_days)
    gap_find = diff1 > 1.1*gap_days
    gap_times = time[gap_find]
    time_index = scipy.r_[0:len(time):1]

    fill_arr_t = scipy.zeros(1)
    fill_arr_f = scipy.zeros(1)
    fill_arr_nan = scipy.zeros(1)

    print 'Filling gaps...'
    for m in scipy.arange(len(gap_times)):
        time_start = time[time_index[gap_find][m]]
        flux_start = flux[time_index[gap_find][m]]
        time_end = time[time_index[gap_find][m] + 1]
        flux_end = flux[time_index[gap_find][m] + 1]
        span =  time_end - time_start
        if span < 2.0*gap_days:
            fill = scipy.array([time_start, time_start+gap_days, time_end])
        else: fill = scipy.r_[time_start: time_end: gap_days]
        fill = fill[1:-1]
        if time_end - fill.max() > 1.1*gap_days: fill = scipy.append(fill, fill.max()+gap_days)

        f = scipy.interpolate.interp1d([time_start, time_end], [flux_start, flux_end])
        gap_new = f(fill)
        if sig_noise > 0: gap_new += normal(0,sig_noise,len(fill))

        fill_arr_t = scipy.append(fill_arr_t, fill)
        fill_arr_f = scipy.append(fill_arr_f, gap_new)
        fill_arr_nan = scipy.append(fill_arr_nan, scipy.ones(len(fill))*scipy.nan)

    # combine time and flux arrays with their filled sections
    fill_arr_t = fill_arr_t[1:]
    fill_arr_f = fill_arr_f[1:]
    fill_arr_nan = fill_arr_nan[1:]
    fill_arr_t = scipy.append(fill_arr_t, time)
    fill_arr_f = scipy.append(fill_arr_f, flux)
    fill_arr_nan = scipy.append(fill_arr_nan, flux_base)

    if len(fill_arr_t) == 0:
        print '*************** empty time array ***************'
        return False, atpy.Table(), 0, 0
    
    # put in table and sort
    tf = atpy.Table()
    tf.add_column('time', fill_arr_t)
    tf.add_column('flux', fill_arr_f)
    tf.add_column('flux_pdc', fill_arr_nan)
    tf.sort('time')

    return True, tf, qt_max, tablen

###############################################################################################################

def acf_calc(quarter_acf, time, flux, interval, kid, max_psearch_len):
    ''' Calculate ACF, calls error calc function'''
    
    # ACF calculation in pylab, close fig when finished
    pylab.figure(50)
    pylab.clf()
    lags, acf, lines, axis = pylab.acorr(flux, maxlags = max_psearch_len)
    pylab.close(50)

    #convolve smoothing window with Gaussian kernel
    gauss_func = lambda x,sig: 1./numpy.sqrt(2*numpy.pi*sig**2) * numpy.exp(-0.5*(x**2)/(sig**2)) #define a Gaussian
    conv_func = gauss_func(numpy.arange(-28,28,1.),9.) #create the smoothing kernel
    acf_smooth = numpy.convolve(acf,conv_func,mode='same') #and convolve
    lenlag = len(lags)
    lags = lags[int(lenlag/2.0):lenlag][:-1] * interval
    acf = acf[int(lenlag/2.0): lenlag][0:-1]
    acf_smooth = acf_smooth[int(lenlag/2.0): lenlag][1:]


    print 'plotting raw acf'
    pylab.clf()
    pylab.plot(acf)
    #pylab.xlim(numpy.median(acf), 1000)
    #pylab.xlim(8900,9200)
    #pylab.ylim(0.95, 1.0)
    pylab.title('Raw acf')
    #raw_input('enter')

    # find max using usmoothed acf (for plot only)
    max_ind_us, max_val_us = extrema(acf, max = True, min = False)

    # find max/min using smoothed acf
    max_ind_s, max_val_s = extrema(acf_smooth, max = True, min = False)
    min_ind_s, min_val_s = extrema(acf_smooth, max = False, min = True)
    maxmin_ind_s, maxmin_val_s = extrema(acf_smooth, max = True, min = True)

    if len(max_ind_s) > 0 and len(min_ind_s) > 0:
        # ensure no duplicate peaks are detected
        t_max_s = atpy.Table()
        t_max_s.add_column('ind', max_ind_s)
        t_max_s.add_column('val', max_val_s)
        t_min_s = atpy.Table()
        t_min_s.add_column('ind', min_ind_s)
        t_min_s.add_column('val', min_val_s)
        t_maxmin_s = atpy.Table()
        t_maxmin_s.add_column('ind', maxmin_ind_s)
        t_maxmin_s.add_column('val', maxmin_val_s)

        ma_i = collections.Counter(t_max_s.ind)
        dup_arr = [i for i in ma_i if ma_i[i]>1]
        if len(dup_arr) > 0:
            for j in scipy.arange(len(dup_arr)):
                tin = t_max_s.where(t_max_s.ind != dup_arr[j])
                tout = t_max_s.where(t_max_s.ind == dup_arr[j])
                tout = tout.rows([0])
                tin.append(tout)
            t_max_s = copy.deepcopy(tin)

        ma_i = collections.Counter(t_min_s.ind)
        dup_arr = [i for i in ma_i if ma_i[i]>1]
        if len(dup_arr) > 0:
            for j in scipy.arange(len(dup_arr)):
                tin = t_min_s.where(t_min_s.ind != dup_arr[j])
                tout = t_min_s.where(t_min_s.ind == dup_arr[j])
                tout = tout.rows([0])
                tin.append(tout)
            t_min_s = copy.deepcopy(tin)

        ma_i = collections.Counter(t_maxmin_s.ind)
        dup_arr = [i for i in ma_i if ma_i[i]>1]
        if len(dup_arr) > 0:
            for j in scipy.arange(len(dup_arr)):
                tin = t_maxmin_s.where(t_maxmin_s.ind != dup_arr[j])
                tout = t_maxmin_s.where(t_maxmin_s.ind == dup_arr[j])
                tout = tout.rows([0])
                tin.append(tout)
            t_maxmin_s = copy.deepcopy(tin)

        t_max_s.sort('ind')
        t_min_s.sort('ind')
        t_maxmin_s.sort('ind')
    
        # relate max inds to lags
        maxnum = len(t_max_s.ind)
        acf_per_pos = lags[t_max_s.ind]
        acf_per_height = acf[t_max_s.ind]
    
        print 'Calculating errors and asymmetries...'
        # Calculate peak widths, asymmetries etc
        acf_per_err, locheight, asym= \
            calc_err(quarter_acf, kid = kid, lags = lags, acf = acf, inds = t_maxmin_s.ind, vals = t_maxmin_s.val, maxnum = maxnum)

    else:
        acf_per_pos = scipy.array([-9999])
        acf_per_height = scipy.array([-9999])
        acf_per_err = scipy.array([-9999])
        locheight = scipy.array([-9999])
        asym = scipy.array([-9999])
        
    # save corrected LC and ACF 
    t_lc = atpy.Table()
    t_lc.add_column('time', time)
    t_lc.add_column('flux', flux)
    t_lc.write('%s/PDCQ%s_output/dat_files/%s_LC.ipac' % (dir, quarter_acf, kid), overwrite = True)

    t_acf = atpy.Table()
    t_acf.add_column('lags_days', lags)
    t_acf.add_column('acf', acf)
    t_acf.add_column('acf_smooth', acf_smooth)
    t_acf.write('%s/PDCQ%s_output/dat_files/%s_ACF.ipac' % (dir, quarter_acf, kid), overwrite = True)

    pylab.figure(6,(10, 5.5))
    pylab.clf()
    pylab.plot(t_acf.lags_days, t_acf.acf, 'k-')
    pylab.plot(t_acf.lags_days, t_acf.acf_smooth, 'r-')
    for j in scipy.arange(len(max_ind_us)):
        pylab.axvline(t_acf.lags_days[max_ind_us[j]], ls = '--', c = 'k', lw=1)
    for i in scipy.arange(len(acf_per_pos)):
        pylab.axvline(acf_per_pos[i], ls = '--', c = 'r', lw = 2)
    if t_acf.lags_days.max() > 10 * acf_per_pos[0]:
        pylab.xlim(0,10 * acf_per_pos[0])
    else: pylab.xlim(0,max(t_acf.lags_days))
    pylab.xlabel('Period (days)')
    pylab.ylabel('ACF')
    pylab.title('ID: %s, Smoothed \& Unsmoothed ACF' %(kid))
    pylab.savefig('%s/PDCQ%s_output/plots_sm/%s_sm.png' % (dir, quarter_acf, kid))

    return  t_acf, acf_per_pos, acf_per_height, acf_per_err, locheight, asym

###############################################################################################################

def pgram_calc(quarter_pgram, time, flux, interval, kid, max_psearch_len):
    ''' Calculate Sine Fitting Periodogram'''
    # Calculate sine lsq periodogram        
    pmin = 0.1
    pmax = interval * max_psearch_len
    nf = 1000
    sinout = gls.sinefit(time, flux, err = None, pmin = pmin, pmax = pmax, nper = nf, \
                             doplot = False, return_periodogram = True)

    sine_per = sinout[0]
    sine_height = max(sinout[5])

    # save periodogram
    t_pg = atpy.Table()
    t_pg.add_column('period', sinout[4])
    t_pg.add_column('pgram', sinout[5])
    t_pg.write('%s/PDCQ%s_output/dat_files/%s_PGRAM.ipac' % (dir, quarter_pgram, kid), overwrite = True)
    
    return t_pg, sine_per, sine_height

###############################################################################################################

def calc_err(quarter_calc, kid, lags, acf, inds, vals, maxnum):
    ''' Calculate peak widths, heights and asymmetries '''
    if len(inds) == 0: return -9999, -9999, -9999
    acf_per_err = scipy.ones(maxnum) * -9999
    asym = scipy.ones(maxnum) * -9999
    mean_height = scipy.ones(maxnum) * -9999

    # Asymmetry plot
    pylab.figure(4,(10, 5.5))
    pylab.clf()
    pylab.title('ID: %d, Asymmetries of 1st 10 peaks' %kid)
    pylab.plot(lags, acf, 'b-')

    acf_ind = scipy.r_[0:len(acf):1]
    num = len(vals)
    if maxnum*2 > num: maxnum -= 1    
    # loop through maxima, assuming 1st index is for minima
    for i in scipy.arange(maxnum):
        # find values and indices of centre left and right
        centre_v = vals[2*i+1]
        centre_i = inds[2*i+1]
        pylab.axvline(lags[centre_i], ls = '--', c = 'r')

        # select value half way between max and min to calc width and asymmetry
        if 2*i + 2 >= num:
            # if it goes beyond limit of lags
            acf_per_err[i] = -9999
            mean_height[i] = -9999
            asym[i] = -9999
        else:
            left_v = vals[2*i]
            left_i = inds[2*i]
            right_v = vals[2*i+2]
            right_i = inds[2*i+2]
            
            sect_left_acf = acf[left_i:centre_i+1]
            sect_left_ind = acf_ind[left_i:centre_i+1]
            sect_right_acf = acf[centre_i:right_i+1]
            sect_right_ind = acf_ind[centre_i:right_i+1]

            height_r = centre_v - right_v
            height_l = centre_v - left_v

            # calc height from peak down 0.5 * mean of side heights
            mean_height[i] = 0.5*(height_l + height_r)
            mid_height = centre_v - 0.5*mean_height[i]
            if mid_height <= min(sect_left_acf): mid_height = min(sect_left_acf)
            if mid_height <= min(sect_right_acf): mid_height = min(sect_right_acf)

            sect_left_acf_r = sect_left_acf[::-1]
            sect_left_ind_r = sect_left_ind[::-1]
            for j in scipy.arange(len(sect_left_acf)):
                if sect_left_acf_r[j] <= mid_height:
                    if j == 0:
                        lag_mid_left = lags[sect_left_ind_r[j]]
                        break
                    else:    
                        pt1 = sect_left_acf_r[j-1]
                        pt2 = sect_left_acf_r[j]
                        lag1 = lags[sect_left_ind_r[j-1]]
                        lag2 = lags[sect_left_ind_r[j]]
                        if pt1 < pt2:
                            ptarr = scipy.array([pt1, pt2])
                            lagarr = scipy.array([lag1, lag2])
                        else:
                            ptarr = scipy.array([pt2, pt1])
                            lagarr = scipy.array([lag2, lag1])   
                        f_left = scipy.interpolate.interp1d(ptarr, lagarr)
                        lag_mid_left = f_left(mid_height)
                        break

            for j in scipy.arange(len(sect_right_acf)):
                if sect_right_acf[j] <= mid_height:
                    if j == 0:
                        lag_mid_right = lags[sect_right_ind[j]]
                        break
                    else:   
                        pt1 = sect_right_acf[j-1]
                        pt2 = sect_right_acf[j]
                        lag1 = lags[sect_right_ind[j-1]]
                        lag2 = lags[sect_right_ind[j]]
                        if pt1 < pt2:
                            ptarr = scipy.array([pt1, pt2])
                            lagarr = scipy.array([lag1, lag2])
                        else:
                            ptarr = scipy.array([pt2, pt1])
                            lagarr = scipy.array([lag2, lag1])
                        f_right = scipy.interpolate.interp1d(ptarr, lagarr)
                        lag_mid_right = f_right(mid_height)
                        break
        
            pos_l = lag_mid_right - lags[centre_i] 
            pos_r = lags[centre_i] - lag_mid_left
            asym[i] = pos_r / pos_l
            acf_per_err[i] = pos_l + pos_r
            if asym[i] <= 0:
                acf_per_err[i] = -9999
                mean_height[i] = -9999
                asym[i] = -9999

        if asym[i] != -9999:    
            pylab.plot([lag_mid_right, lag_mid_left], [mid_height, mid_height], 'm-')
            pylab.axvline(lag_mid_left, ls = '--', c = 'g')
            pylab.axvline(lag_mid_right, ls = '--', c = 'g')

    if 10 * lags[inds[0]] < max(lags):
        pylab.xlim(0, 10 * lags[inds[1]])
    pylab.xlabel('Lags (days)')
    pylab.ylabel('ACF')
    pylab.savefig('%s/PDCQ%s_output/plots_asym/%s_asym.png' % (dir, quarter_calc, kid))
    
    return acf_per_err, mean_height, asym

###############################################################################################################

def plot_stats(quarter_stats, time, flux, kid_x, acf_per_pos_in, acf_per_height_in, acf_per_err_in, locheight_in, asym_in):
    ''' Plot and calculate statistics of peaks '''
    acf_per_pos_in = acf_per_pos_in[asym_in != -9999]
    acf_per_height_in = acf_per_height_in[asym_in != -9999]
    acf_per_err_in = acf_per_err_in[asym_in != -9999]
    locheight_in = locheight_in[asym_in != -9999]
    asym_in = asym_in[asym_in != -9999]

    if len(acf_per_pos_in) == 0: return -9999, -9999, -9999, -9999, -9999, -9999,\
        -9999, -9999, -9999, -9999, 0, -9999, -9999, 1, 0.0

    x = 10  #number of periods used for calc
    hdet = 0  #start with 0 harmonic, set to 1 if 1/2P is 1st peak  
    # deal with cases where 1/2 P is 1st peak
    if len(acf_per_pos_in) >= 2:
        one_peak_only = 0
        ind = scipy.r_[1:len(acf_per_pos_in)+1:1]
        if locheight_in[1] > locheight_in[0]:
            print '1/2 P found (1st)'
            hdet = 1  # mark harmonic found
            pk1 = acf_per_pos_in[1]
            acf_per_pos_in = acf_per_pos_in[1:]
            acf_per_height_in = acf_per_height_in[1:]
            acf_per_err_in = acf_per_err_in[1:]
            locheight_in = locheight_in[1:]
            asym_in = asym_in[1:]
            '''if 1 == 1:
            pknumin = int(raw_input('pk in: '))
            if pknumin > 0: hdet = 1  # mark harmonic found
            pk1 = acf_per_pos_in[pknumin]
            acf_per_pos_in = acf_per_pos_in[pknumin:]
            acf_per_height_in = acf_per_height_in[pknumin:]
            acf_per_err_in = acf_per_err_in[pknumin:]
            locheight_in = locheight_in[pknumin:]
            asym_in = asym_in[pknumin:]'''
        else:
            print 'no harmonic found or harmonic not 1st peak'
            pk1 = acf_per_pos_in[0]
    else:
        print 'Only 1 peak'
        one_peak_only = 1
        pk1 = acf_per_pos_in[0] 


    # select only peaks which are ~multiples of 1st peak (within phase 0.2)
    acf_per_pos_0_test = scipy.append(0,acf_per_pos_in)
    ind_keep = scipy.r_[0:len(acf_per_pos_0_test):1]
    fin = False
    while fin == False:
        delta_lag_test = acf_per_pos_0_test[1:] - acf_per_pos_0_test[:-1]
        delta_lag_test = scipy.append(0, delta_lag_test)
        phase = ((delta_lag_test % pk1) / pk1)
        phase[phase > 0.5] -= 1.0
        excl = abs(phase) > 0.2
        ind_temp = scipy.r_[0:len(delta_lag_test):1]
        if len(phase[excl]) == 0: break
        else:
            ind_rem = ind_temp[excl][0]
            rem = acf_per_pos_0_test[ind_rem]
            ind_keep = ind_keep[acf_per_pos_0_test != rem]
            acf_per_pos_0_test = acf_per_pos_0_test[acf_per_pos_0_test != rem]
            
    ind_keep = ind_keep[1:] - 1
    keep_pos = acf_per_pos_in[ind_keep]
    keep_pos = scipy.append(0, keep_pos)
    delta_keep_pos = keep_pos[1:] - keep_pos[:-1]
    # remove very small delta lags points (de-noise peak detections)
    if len(ind_keep) > 1:
        ind_keep = ind_keep[delta_keep_pos > 0.3*delta_keep_pos[0]]
        delta_keep_pos = delta_keep_pos[delta_keep_pos > 0.3*delta_keep_pos[0]]

    if len(ind_keep) > 1:
        ind_gap = ind_keep[delta_keep_pos > 2.2*delta_keep_pos[0]]
        if len(ind_gap) != 0:
            ind_keep = ind_keep[ind_keep < ind_gap[0]]
            delta_keep_pos = delta_keep_pos[ind_keep < ind_gap[0]]

    # limit to x lags for plot and calc
    if len(acf_per_pos_in[ind_keep]) < x: x = len(acf_per_pos_in[ind_keep])
        
    acf_per_pos = acf_per_pos_in[ind_keep][:x]
    acf_per_height = acf_per_height_in[ind_keep][:x]
    acf_per_err = acf_per_err_in[ind_keep][:x]
    asym = asym_in[ind_keep][:x]
    locheight = locheight_in[ind_keep][:x]
    print '%d Peaks kept' %len(acf_per_pos)


    if len(acf_per_pos) == 1:
        return -9999, 0.0, pk1, acf_per_height[0], acf_per_err[0], locheight[0], -9999, -9999, -9999, -9999, 1, hdet, -9999, 1, 0.0

    ''' Delta Lag '''
    acf_per_pos_0 = scipy.append(0,acf_per_pos)
    delta_lag = acf_per_pos_0[1:] - acf_per_pos_0[:-1]
    av_delt = scipy.median(delta_lag)
    delt_mad = 1.483*scipy.median(abs(delta_lag - av_delt)) # calc MAD
    delt_mad = delt_mad / scipy.sqrt(float(len(delta_lag)-1.0)) # err = MAD / sqrt(n-1)
    med_per = av_delt
    mad_per_err = delt_mad
    
    pylab.figure(3,(12, 9))
    pylab.clf()
    pylab.subplot(3,2,1)
    pylab.plot(acf_per_pos, delta_lag, 'bo')
    pylab.axhline(av_delt, ls = '--', c = 'k')
    pylab.axhspan(av_delt-delt_mad, av_delt+delt_mad, facecolor = 'k', alpha=0.2)
    pylab.xlim(0, 1.1*acf_per_pos.max())
    pylab.ylim(0.99*delta_lag.min(), 1.01*delta_lag.max())
    pylab.ylabel('Delta Lag (days)')
    # use 1st selected peak
    pylab.plot(pk1, delta_lag[acf_per_pos == pk1][0], 'ro')

    ''' Abs Height '''
    pylab.subplot(3,2,2)
    pylab.plot(acf_per_pos, acf_per_height, 'bo')
    pylab.xlim(0, 1.1*acf_per_pos.max())
    if acf_per_height.min() >= 0 and acf_per_height.max() >= 0: pylab.ylim(0.9*acf_per_height.min(), 1.1*(acf_per_height.max()))
    elif acf_per_height.min() < 0 and acf_per_height.max() >= 0: pylab.ylim(1.1*acf_per_height.min(), 1.1*(acf_per_height.max()))
    elif acf_per_height.max() < 0: pylab.ylim(1.1*acf_per_height.min(), 0.9*(acf_per_height.max()))
    pylab.plot(pk1, acf_per_height[acf_per_pos == pk1][0], 'ro')
    pylab.axhline(acf_per_height[acf_per_pos == pk1][0], ls = '--', c = 'k')
    ax = pylab.gca()
    pylab.text(0.42, 0.9, 'H1=%.3f' %(acf_per_height[acf_per_pos == pk1][0]), transform = ax.transAxes)
    pylab.ylabel('Abs Height')

    ''' Local Height '''
    pylab.figure(3)
    pylab.subplot(3,2,3)    
    pylab.plot(acf_per_pos, locheight, 'bo')
    pylab.xlim(0, 1.1*acf_per_pos.max())
    if min(locheight) >= 0 and max(locheight) >= 0:
        pylab.ylim(0.99*min(locheight), 1.01*max(locheight))
        ymax = 1.01*max(locheight)
    elif min(locheight) < 0 and max(locheight) >= 0:
        pylab.ylim(1.01*min(locheight), 1.01*max(locheight))
        ymax = 1.01*max(locheight)
    elif max(locheight) < 0:
        pylab.ylim(1.01*min(locheight), 0.99*max(locheight))
        ymax = 0.99*max(locheight)
    pylab.plot(pk1, locheight[acf_per_pos == pk1][0], 'ro')
    pylab.ylabel('Local Height')

    if len(acf_per_pos) > 2 and no_rpy == False:
        print 'Fitting line to Local Height...'
        st_line = lambda p, x: p[0] + p[1] * x
        '''st_line_err = lambda p, x, y, fjac: [0, (y - st_line(p, x)), None]
        fa = {'x': acf_per_pos, 'y': locheight}
        p = [scipy.median(locheight), 0.0]
        m = mpfit.mpfit(st_line_err, p, functkw = fa, quiet = True)
        p = m.params
        errs = m.perror
        h_grad = p[1]
        h_timescale = 1.0 / h_grad
        h_grad_err = errs[1]'''
        r.assign('xdf1', acf_per_pos)
        r.assign('ydf1', locheight)
        p1 = r('''
        xdf <- c(xdf1)
        ydf <- c(ydf1)
        library(quantreg)
        rqmodel1 <- rq(ydf~xdf)
        plot(xdf,ydf)
        abline(rqmodel1,col=3)
        sumr = summary(rqmodel1, se = 'boot')
        grad <- rqmodel1[['coefficients']][['xdf']]
        interc <- rqmodel1[['coefficients']][['(Intercept)']]
        err_grad <- sumr[[3]][[3]]
        err_interc <- sumr[[3]][[4]]
        resids = resid(rqmodel1)
        output <- c(resids, grad, interc, err_grad, err_interc)
        ''')
        res =scipy.array(p1[0:len(acf_per_pos)])
        pqr = scipy.array(p1[len(acf_per_pos):])
        h_grad = pqr[0]
        h_timescale = 1.0 / h_grad
        h_grad_err = pqr[3]
        h_grad_scatter = sum((st_line([pqr[1],pqr[0]], acf_per_pos) - locheight) ** 2) / scipy.sqrt(float(len(acf_per_pos)))
        pylab.plot(acf_per_pos, st_line([pqr[1],pqr[0]], acf_per_pos), 'r-')
        ax = pylab.gca()
        pylab.text(0.3, 0.9,'m=%.5f, TS=%.2f' %(h_grad, abs(h_timescale)), transform = ax.transAxes)
    else:
        h_grad = -9999
        h_timescale = -9999
        h_grad_err = -9999
        h_grad_scatter = -9999
 
    ''' Width '''
    pylab.subplot(3,2,4)
    pylab.plot(acf_per_pos, acf_per_err, 'bo')
    pylab.xlim(0, 1.1*acf_per_pos.max())
    pylab.ylim(0.99*acf_per_err.min(), 1.01*acf_per_err.max())
    ymax = 1.01*acf_per_err.max()
    pylab.plot(pk1, acf_per_err[acf_per_pos == pk1][0], 'ro')
    pylab.ylabel('Width (days)')
    pylab.xlabel('Lag (days)')

    if len(acf_per_pos) > 2 and no_rpy == False:
        print 'Fitting line to Width...'
        r.assign('xdf1', acf_per_pos)
        r.assign('ydf1', acf_per_err)
        p1 = r('''
        xdf <- c(xdf1)
        ydf <- c(ydf1)
        library(quantreg)
        rqmodel1 <- rq(ydf~xdf)
        plot(xdf,ydf)
        abline(rqmodel1,col=3)
        sumr = summary(rqmodel1, se = 'boot')
        grad <- rqmodel1[['coefficients']][['xdf']]
        interc <- rqmodel1[['coefficients']][['(Intercept)']]
        err_grad <- sumr[[3]][[3]]
        err_interc <- sumr[[3]][[4]]
        resids = resid(rqmodel1)
        output <- c(resids, grad, interc, err_grad, err_interc)
        ''')
        res =scipy.array(p1[0:len(acf_per_pos)])
        pqr = scipy.array(p1[len(acf_per_pos):])
        w_grad = pqr[0]
        w_timescale = 1.0 / w_grad
        w_grad_err = pqr[3]
        w_grad_scatter = sum((st_line([pqr[1],pqr[0]], acf_per_pos) - acf_per_err) ** 2) / scipy.sqrt(float(len(acf_per_pos)))
        pylab.plot(acf_per_pos, st_line([pqr[1],pqr[0]], acf_per_pos), 'r-')
        ax = pylab.gca()
        pylab.text(0.3, 0.9, 'm=%.5f, TS=%.2f' %(w_grad, abs(w_timescale)), transform = ax.transAxes)
    else:
        w_grad = -9999
        w_timescale = -9999
        w_grad_err = -9999
        w_grad_scatter = -9999
   
    pylab.subplot(3,2,5)
    pylab.plot(acf_per_pos, asym , 'bo')
    pylab.xlim(0, 1.1*acf_per_pos.max())
    pylab.plot(pk1, asym[acf_per_pos == pk1][0], 'ro')
    pylab.ylim(0.99*min(asym), 1.01*max(asym))
    pylab.ylabel('Asymmetry')
    pylab.xlabel('Lag (days)')
    pylab.suptitle('ID: %s, P\_dlag = %.3fd +/-%.3fd, P\_pk = %.3fd' \
                   %(kid_x, med_per, mad_per_err, pk1), fontsize = 16)
    pylab.savefig('%s/PDCQ%s_output/plots_stats/%s_stats.png' % (dir, quarter_stats, kid_x))

    peak_ratio = len(acf_per_pos)/len(acf_per_pos_in)

    return med_per, mad_per_err, pk1, acf_per_height[acf_per_pos == pk1][0], acf_per_err[acf_per_pos == pk1][0],\
        locheight[acf_per_pos == pk1][0], h_grad, h_grad_scatter, w_grad, w_grad_scatter, len(acf_per_pos), hdet, acf_per_pos, one_peak_only, peak_ratio

###############################################################################################################

def calc_var(quarter_var, kid = None, time_in = None, flux = None, period = None):
    ''' calculate 5-95th percentile of flux (i.e. amplitude) whole LC and in each period block '''

    # Of whole LC...
    sort_flc= sorted(flux)
    fi_ind_flc = int(len(sort_flc) * 0.05)
    nifi_ind_flc = int(len(sort_flc) * 0.95)
    amp_flc = sort_flc[nifi_ind_flc] - sort_flc[fi_ind_flc]

    # Of each block...
    if period > 0:
        num = int(scipy.floor(len(time_in) / period))
        var_arr = scipy.zeros(num) + scipy.nan
        per_cent = scipy.zeros(num) + scipy.nan
        for i in scipy.arange(num-1):
            t_block = time_in[(time_in-time_in.min() >= i*period) * (time_in-time_in.min() < (i+1)*period)]
            f_block = flux[(time_in-time_in.min() >= i*period) * (time_in-time_in.min() < (i+1)*period)]
            if len(t_block) < 2: continue
            sort_f_block = sorted(f_block)
            fi_ind = int(len(sort_f_block) * 0.05)
            nifi_ind = int(len(sort_f_block) * 0.95)
            range_f = sort_f_block[nifi_ind] - sort_f_block[fi_ind]
            var_arr[i] = range_f
            per_cent[i] = scipy.mean(t_block)

        var_arr_real = var_arr[scipy.isfinite(var_arr) == True]
        per_cent = per_cent[scipy.isfinite(var_arr) == True]
        var_med = scipy.median(var_arr_real)

        
        savefilen = '%s/PDCQ%s_output/dat_files/%s_amps.txt' % (dir, quarter_var, kid)
        savefile = open(savefilen, 'w')
        savefile.write('per_centre   amp\n')
        for z in scipy.arange(len(per_cent)):    
            savefile.write('%f   %f\n' %(per_cent[z], var_arr_real[z]))
        savefile.close()
    else:
        var_med = -9999
        per_cent = scipy.array([-9999])
        var_arr_real = scipy.array([-9999])
    
    return amp_flc, var_med, per_cent, var_arr_real

###############################################################################################################

def EBTransitPhase(tset, kid_x):
    ebp = atpy.Table('%s/eb_pars.txt' %dir, type = 'ascii')
    
    lc = tset.tables[1]   
    time = lc.TIME
    flux = lc.PDCSAP_FLUX
    nobs = len(time)
    lg = scipy.isfinite(time)
    pylab.figure(52)
    pylab.clf()
    pylab.plot(time[lg], flux[lg])
    npl = 2
    phase = scipy.zeros((npl, nobs))
    inTr = scipy.zeros((npl, nobs), 'int')
    period = ebp.P[ebp.KID == kid_x]
    for ipl in scipy.arange(npl):
        if ipl == 0: t0 = ebp.Ep1[ebp.KID == kid_x]
        if ipl == 1: t0 = ebp.Ep2[ebp.KID == kid_x]
        if ipl == 0: dur = ebp.Dur1[ebp.KID == kid_x]
        if ipl == 1: dur = ebp.Dur2[ebp.KID == kid_x]
        dur /= period
        counter = 0
        while (time[lg] - t0).min() < 0:
            t0 -= period
            counter += 1
            if counter > 1000: break
        ph = ((time - t0) % period) / period
        ph[ph < -0.5] += 1
        ph[ph > 0.5] -= 1
        phase[ipl,:] = ph
        inTr[ipl,:] = (abs(ph) <= dur/1.5)
    return phase, inTr

###############################################################################################################

def extrema(x, max = True, min = True, strict = False, withend = False):
	"""
	This function will index the extrema of a given array x.
	Options:
		max		If true, will index maxima
		min		If true, will index minima
		strict		If true, will not index changes to zero gradient
		withend	If true, always include x[0] and x[-1]
	This function will return a tuple of extrema indexies and values
	"""
	# This is the gradient
	from numpy import zeros
	dx = zeros(len(x))
	from numpy import diff
	dx[1:] = diff(x)
	dx[0] = dx[1]
	
	# Clean up the gradient in order to pick out any change of sign
	from numpy import sign
	dx = sign(dx)
	
	# define the threshold for whether to pick out changes to zero gradient
	threshold = 0
	if strict:
		threshold = 1
		
	# Second order diff to pick out the spikes
	d2x = diff(dx)
	
	if max and min:
		d2x = abs(d2x)
	elif max:
		d2x = -d2x
	
	# Take care of the two ends
	if withend:
		d2x[0] = 2
		d2x[-1] = 2
	
	# Sift out the list of extremas
	from numpy import nonzero
	ind = nonzero(d2x > threshold)[0]
	return ind, x[ind]



#for q in range(3, 17):
#    corr_run(q)
