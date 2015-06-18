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
from scipy import optimize
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import signal

#if rpy not available, comment out this..
#from rpy2 import robjects as ro
#r = ro.r
#import rpy2.robjects.numpy2ri
#rpy2.robjects.numpy2ri.activate()
# and set no_rpy to True...
no_rpy = True

from misc import KOI_tools_b12 as kt
from misc import filter
from misc import gls
from misc import mpfit
from misc import peakdetect as pd

# Set up constants and directories
# max and min quarters to include
qmin = 3
qmax = 14

dir = '/home/gil-pc/amy/tmp/set_exmis%d/Q%s_%s' %(qmin, qmax)
gap_days = 0.02043365  # assume for long cadence
# known jumps and discontinuities 
jump_arr = scipy.array([131.51139, 169.51883, 169.75000, 182.00000, 200.31000, 231.00000, 246.19000, 256.00000, 260.22354, 281.00000, 291.00000, 322.00000, 352.37648, 373.23000, 384.00000, 398.00000, 443.48992, 475.50000, 504.00000, 539.44868, 567.00000, 599.00000, 630.17387, 661.00000, 691.00000, 711.20000, 735.36319, 762.00000, 808.51558, 845.00000, 874.50000, 906.84469, 937.00000, 970.00000, 1001.20718, 1032.50000, 1063.50000 ,1071.00000, 1093.60000, 1098.00000, 1117.00000, 1121.40000, 1126.40000, 1154.40000, 1161.00000, 1182.00000, 1215.40000, 1245.30000, 1269.00000, 1273.20000, 1290.00000, 1306.00000, 1337.00000])


###############################################################################################################

def corr_run(tr_out = False, tr_type = None, doplot = False):
    ''' Reads index and light curves and runs the ACF calculation modules
        tr_out = True if transits are to be removed
        tr_type = EB if EB not planet'''
    #read in index and create list of unique KIDs
    index = atpy.Table('/home/orbit/amy/Documents/tmp/FGK/indices/missing_8q.ipac', type = 'ipac')  #index with KIDs and filenames
    index.sort('keplerid')
    
    print 'Index includes %s files...'% len(index)
    id_list = scipy.array(list(Set(index.keplerid)))
    print '...and %d individual stars' % len(id_list)
    id_list.sort()

    # get one row of main table for each KID (currently index has 1 row per quarter)
    for z in scipy.arange(len(id_list)):
        if z == 0:
            t = index.where(index.keplerid == id_list[z]).rows([0])
        else: t.append(index.where(index.keplerid == id_list[z]).rows([0]))

    t.sort('keplerid')
    if len(t) != len(id_list):
        print 'Error: table length != KID list length'
        return

    print 'len(id_list)', len(id_list)
    raw_input('enter')

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
    period_err_arr = scipy.ones(len(id_list)) * -9999.0

    # loop over KIDs
    for x in scipy.arange(len(id_list)):
        kid_x = int(id_list[x])
        print 'x = %d/%d, kid = %d' %(x, len(id_list), kid_x)
        print 'set%d, qmin%d, qmax%d' %(fgkset, qmin, qmax)

        # Load light curve (load_lc returns normalised, quarter-joined LC)
        print 'Loading LC...'
        condition, lc_tab, qt_max, tablen, trange = load_lc(kid_x, index, tr_out, tr_type)
        if condition == False: continue

        if doplot == True:
            # main ACF figure (raw, corr, amp, acf)
            pylab.figure(1,(12, 9))
            pylab.clf()
            pylab.subplots_adjust(bottom=0.06, top =0.97, left = 0.09, right=0.98, hspace = 0.3)
            pylab.subplot(4,1,1)
            pylab.title('ID: %s' %(kid_x), fontsize = 16)
            for j in scipy.arange(tablen):
                if j % 2 == 0:
                    pylab.axvspan(qt_max[j], qt_max[j+1], facecolor = 'k', alpha=0.1)
            pylab.plot(lc_tab.time, lc_tab.flux_raw, 'k-')
            for k in scipy.arange(len(jump_arr)):
                pylab.axvline(jump_arr[k], ls = '--', c = 'b')
            pylab.xlim(lc_tab.time.min(), lc_tab.time.max())
            pylab.ylim(min(lc_tab.flux_raw[scipy.isfinite(lc_tab.flux_raw) == True]), \
                       max(lc_tab.flux_raw[scipy.isfinite(lc_tab.flux_raw) == True]))
            pylab.ylabel('Raw Flux')
            pylab.subplot(4,1,2)
            toplot_flux = copy.copy(lc_tab.flux)
            for n in scipy.arange(len(toplot_flux)):
                if toplot_flux[n] == 0.0: toplot_flux[n] = scipy.nan
            pylab.plot(lc_tab.time, toplot_flux, 'k-')
            pylab.xlim(lc_tab.time.min(), lc_tab.time.max())
            pylab.ylim(lc_tab.flux.min(), lc_tab.flux.max())
            pylab.ylabel('Norm Flux')

            pylab.figure(2,(12, 9))
            pylab.clf()
            pylab.subplot(3,1,1)
            pylab.subplots_adjust(bottom=0.06, top =0.96, left = 0.09, right=0.98, hspace = 0.25)
            pylab.title('ID: %s' %(kid_x), fontsize = 16)
            for j in scipy.arange(tablen):
                if j % 2 == 0:
                    pylab.axvspan(qt_max[j], qt_max[j+1], facecolor = 'k', alpha=0.1)
            pylab.plot(lc_tab.time, lc_tab.flux_raw, 'k-')
            for k in scipy.arange(len(jump_arr)):
                pylab.axvline(jump_arr[k], ls = '--', c = 'b')
            pylab.xlim(lc_tab.time.min(), lc_tab.time.max())
            pylab.ylim(min(lc_tab.flux_raw[scipy.isfinite(lc_tab.flux_raw) == True]), \
                       max(lc_tab.flux_raw[scipy.isfinite(lc_tab.flux_raw) == True]))
            pylab.ylabel('Raw Flux')
            pylab.subplot(3,1,2)

            toplot_flux = copy.copy(lc_tab.flux)
            for n in scipy.arange(len(toplot_flux)):
                if toplot_flux[n] == 0.0: toplot_flux[n] = scipy.nan

            pylab.plot(lc_tab.time, toplot_flux, 'k-')
            pylab.xlim(lc_tab.time.min(), lc_tab.time.max())
            pylab.ylim(lc_tab.flux.min(), lc_tab.flux.max())
            pylab.ylabel('Norm Flux')
            pylab.xlabel('Time (days)')

        # max period searched for is len(flux) / 2
        max_psearch_len = len(lc_tab.flux) / 2.0
        
        # Calculate ACF
        print 'Calculating ACF...'
        acf_tab, acf_per_pos, acf_per_height, acf_per_err, locheight, asym,  = \
            acf_calc(time = lc_tab.time, flux = lc_tab.flux, interval = gap_days, kid = kid_x, max_psearch_len = max_psearch_len)

        # Calculate periodogram
        print 'Calculating periodogram...'
        pgram_tab, sine_per[x], sine_height[x] = \
            pgram_calc(time = lc_tab.time, flux = lc_tab.flux, interval = gap_days, kid = kid_x, max_psearch_len = max_psearch_len)

        if doplot == True:
            pylab.figure(1)
            pylab.subplot(4,1,4)
            pylab.plot(acf_tab.lags_days, acf_tab.acf, 'k-')
            pylab.axhline(0, ls = '--', c = 'k')
            for i in scipy.arange(len(acf_per_pos)):
                pylab.axvline(acf_per_pos[i], ls = '--', c = 'b')
            pylab.ylabel('ACF')
            pylab.xlabel('Period (days)')
            pylab.xlim(0, acf_tab.lags_days.max())

            pylab.figure(2)
            pylab.subplot(3,1,3)
            pylab.plot(pgram_tab.period, pgram_tab.pgram, 'k-')
            pylab.axvline(sine_per[x], c = 'r', ls = '--')
            pylab.ylim(0, 1.1*max(pgram_tab.pgram))
            pylab.axhline(0.3, c = 'g', ls = '--')
            pylab.xlim(0, pgram_tab.period.max())
            pylab.ylabel('Periodogram')
            pylab.xlabel('Period (days)')
            pylab.savefig('%s/ACF_output/plots_pgram/%s_pgram.png' % (dir, kid_x))
    
        # plot and calculate acf peak statistics
        if acf_per_pos[0] != -9999:
            med_dlag_per[x], dlag_per_err[x], acf_peak_per[x], h1[x], w1[x], lh1[x], \
                hlocgrad[x], hloc_grad_scatter[x], width_grad[x], width_grad_scatter[x], \
                num_of_peaks[x], harmonic_det[x], sel_peaks =\
                plot_stats(lc_tab.time, lc_tab.flux, kid_x, acf_per_pos, acf_per_height, acf_per_err, locheight, asym, doplot)

            if med_dlag_per[x] >0:
                period[x] = med_dlag_per[x]
                period_err_arr[x] = dlag_per_err[x]
            else:
                period[x] = acf_peak_per[x]
                period_err_arr[x] = w1[x] / 2.3548   # sigma = fwhm / 2.3548

            if doplot == True:
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
            print 'Calculating variability stats...'
            amp_all[x], amp_per[x], per_cent, var_arr_real = \
                calc_var(kid = kid_x, time_in = trange.time, flux = trange.flux, period = acf_peak_per[x])

            if doplot == True:
                pylab.figure(1)
                pylab.subplot(4,1,3)
                pylab.plot(per_cent, var_arr_real, 'k.')
                pylab.axhline(amp_per[x], ls = '--', c = 'b')
                pylab.xlim(lc_tab.time.min(),lc_tab.time.max())
                pylab.ylim(var_arr_real.min(), var_arr_real.max())
                pylab.xlabel('Time (days)')
                pylab.ylabel('Amplitudes')
                pylab.savefig('%s/ACF_output/plots_acf/%s_full.png' % (dir, kid_x))

                # make zoomed in plot to show periodicity
                pylab.figure(5,(12, 9))
                pylab.clf()
                pylab.subplot(3,1,1)
                pylab.title('Zoom of ID: %s, showing peak period and av lag period' %kid_x)
                maxpts = 40.0
                if scipy.floor(lc_tab.time.max() / acf_peak_per[x]) < maxpts:
                    maxpts = float(scipy.floor(lc_tab.time.max() / acf_peak_per[x]))
                inc = lc_tab.time - lc_tab.time.min() <= (maxpts*acf_peak_per[x])
                toplot_flux = copy.copy(lc_tab.flux[inc])
                for n in scipy.arange(len(toplot_flux)):
                    if toplot_flux[n] == 0.0: toplot_flux[n] = scipy.nan
                pylab.plot(lc_tab.time[inc], toplot_flux, 'k-')
                for i in scipy.arange(int(maxpts)):
                    pylab.axvline(lc_tab.time.min() + (i)*(acf_peak_per[x]), ls = '--', c = 'b')        
                pylab.xlim(lc_tab.time[inc].min(), lc_tab.time[inc].max())
                pylab.ylabel('Flux (with Peak Pos)')
                pylab.subplot(3,1,2)

                toplot_flux = copy.copy(lc_tab.flux[inc])
                for n in scipy.arange(len(toplot_flux)):
                    if toplot_flux[n] == 0.0: toplot_flux[n] = scipy.nan

                pylab.plot(lc_tab.time[inc], toplot_flux, 'k-')
                for i in scipy.arange(int(maxpts)):
                    pylab.axvline(lc_tab.time.min() + (i)*(med_dlag_per[x]), ls = '--', c = 'r')        
                pylab.xlim(lc_tab.time[inc].min(), lc_tab.time[inc].max())
                pylab.ylabel('Flux (with Med Pos)')
                pylab.xlabel('Time (days)')
                pylab.subplot(3,1,3)
                pylab.subplots_adjust(top=0.96, bottom = 0.06)
                pylab.plot(acf_tab.lags_days[acf_tab.lags_days <= (maxpts*acf_peak_per[x])], \
                           acf_tab.acf[acf_tab.lags_days <= (maxpts*acf_peak_per[x])])
                pylab.ylabel('ACF')
                pylab.xlabel('Lag (days)')
                pylab.axvline(acf_peak_per[x], ls = '--', c = 'b')
                pylab.axvline(med_dlag_per[x], ls = '--', c = 'r')
                if acf_tab.lags_days.max() > 10 * acf_peak_per[x]:
                    pylab.xlim(0,10 * acf_peak_per[x])
                else: pylab.xlim(0,max(acf_tab.lags_days))
                pylab.savefig('%s/ACF_output/plots_z/%s_zoom.png' % (dir, kid_x))

        # save file of stats for individual star
        t_stats = atpy.Table()
        t_stats.add_column('acf_per_pos', acf_per_pos)
        t_stats.add_column('acf_per_height', acf_per_height)
        t_stats.add_column('acf_per_err', acf_per_err)
        t_stats.add_column('asym', asym)
        t_stats.add_column('locheight', locheight)
        savefilen = '%s/ACF_output/dat_files/%s_stats.txt' % (dir, kid_x)
        atpy.asciitables.write_ascii(t_stats, savefilen, delimiter = ' ', overwrite = True)
        
        #raw_input('enter')

    # save stats that are 1 per kid
    t.remove_columns(['filename'])
    t.add_column('period', period)
    t.add_column('period_err', period_err_arr)
    t.add_column('sine_per', sine_per)
    t.add_column('sine_height', sine_height)
    t.add_column('acf_peak_per', acf_peak_per)
    t.add_column('med_dlag_per', med_dlag_per)
    t.add_column('dlag_per_err', dlag_per_err)
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
    savefilen = '%s/ACF_output/corr_output_pt3.txt' % (dir)
    atpy.asciitables.write_ascii(t, savefilen, delimiter = ' ', overwrite = True)
    
    return

###############################################################################################################

def load_lc(kid_x = None, index = None, tr_out = False, tr_type = False):
    # read in tableset of LCs for all Qs, normalise and quarter join
    dir_lc = dir + '/Light_Curves'
    print 'Reading in LCs...'

    filename = index.filename[(index.keplerid == kid_x) * (index.quarter <= qmax) * (index.quarter >= qmin)]
    tset, found = kt.GetLc(id = kid_x, dir = dir_lc, tr_out = tr_out, \
                           filename = filename)
    if found == 0:
        print '****not found****'
        return False, atpy.Table(), 0, 0, atpy.Table()
    
    # find number of Quarters and max(time) per quarter
    if tr_out == False: lc = tset[0]
    else: lc = tset[1]

    quarters = scipy.array(list(set(lc.Q)))
    quarters = scipy.sort(quarters)
    print 'quarters: ', quarters
    tablen = len(quarters)
    qt_max = scipy.zeros(tablen)
    qt_min = scipy.zeros(tablen)
    for z in scipy.arange(len(quarters)):
        lc_time_new = lc.TIME[scipy.isfinite(lc.TIME) == True]
        qt_max[z] = max(lc.TIME[lc.Q == quarters[z]])
        qt_min[z] = min(lc.TIME[lc.Q == quarters[z]])
    qt_max = scipy.append(0,qt_max)
    qt_min = scipy.append(0,qt_min)
    qt_max = scipy.append(qt_max,5000)
    qt_min = scipy.append(qt_min,5000)
    
    finite_base = scipy.where((scipy.isfinite(lc.TIME) == True) * (scipy.isfinite(lc.SAP_FLUX) == True))
    
    # keep unchanged set for table
    time = lc.TIME
    time_base = copy.copy(time[finite_base])
    flux_base = lc.SAP_FLUX[finite_base]
    q_arr_base = lc.Q[finite_base]
    flux = lc.PDCSAP_FLUX

    if tr_out == True:
        print 'Removing eclipses...'    
        # remove in-eclipse regions
        phase, inTr = kt.TransitPhase(tset)
        npl, nobs = inTr.shape
        for ipl in scipy.arange(npl):
            trlist = inTr[ipl,:].astype(bool)
            flux[trlist] = scipy.nan
            
    if tr_type == 'EB':
        print 'Removing EB transits...'
        phase, inTr = EBTransitPhase(tset, kid_x)
        npl, nobs = inTr.shape
        for ipl in scipy.arange(npl):
            trlist = inTr[ipl,:].astype(bool)
            flux[trlist] = scipy.nan
        pylab.plot(time, flux, 'r.')

    finite = scipy.where((scipy.isfinite(time) == True) * (scipy.isfinite(flux) == True) * (scipy.isfinite(lc.Q) == True))

    time = time[finite]
    flux = lc.PDCSAP_FLUX[finite]
    quarter_arr = lc.Q[finite]

    # median normalise each quarter
    for qt in quarters:
        flux[quarter_arr == qt] = (flux[quarter_arr == qt] / scipy.median(flux[quarter_arr == qt])) - 1.0
        flux_base[q_arr_base == qt] = (flux_base[q_arr_base == qt] / scipy.median(flux_base[q_arr_base == qt])) - 1.0
    
    # fit 1d polynomial and remove >5sig outliers 
    coeff = scipy.polyfit(time,flux, 1) 
    p = scipy.poly1d(coeff) 
    med, sig = filter.medsig(flux - p(time))
    sigclip = abs(flux - p(time)) <= 5*sig
    time = time[sigclip]
    flux = flux[sigclip]

    trange = atpy.Table()
    trange.add_column('time', time)
    trange.add_column('flux', flux)
    trange.sort('time')

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
        # fill temporarily with NAN
        gap_new = fill*scipy.nan
        fill_arr_t = scipy.append(fill_arr_t, fill)
        fill_arr_f = scipy.append(fill_arr_f, gap_new)

    # combine time and flux arrays with their filled sections
    fill_arr_t = fill_arr_t[1:]
    fill_arr_f = fill_arr_f[1:]
    fill_arr_t = scipy.append(fill_arr_t, time)
    fill_arr_f = scipy.append(fill_arr_f, flux)

    if len(fill_arr_t) == 0:
        print '*************** empty time array ***************'
        return False, atpy.Table(), 0, 0

    finite = scipy.where((scipy.isfinite(fill_arr_t) == True) * (scipy.isfinite(fill_arr_f) == True))
    fill_arr_t = fill_arr_t[finite]
    fill_arr_f = fill_arr_f[finite]
    
    # put in table and sort
    tf1 = atpy.Table()
    tf1.add_column('time', fill_arr_t)
    tf1.add_column('flux', fill_arr_f)
    tf1.sort('time')

    tf2 = atpy.Table()
    tf2.add_column('time_base', time_base)
    tf2.add_column('flux_base', flux_base)
    tf2.sort('time_base')

    # map to exactly even time sampling  
    tnew = scipy.r_[tf1.time.min():tf1.time.max():gap_days]
    f_interp = scipy.interpolate.interp1d(tf1.time, tf1.flux)
    f_interp_raw = scipy.interpolate.interp1d(tf2.time_base, tf2.flux_base)
    fnew = f_interp(tnew)
    fnew_raw = f_interp_raw(tnew)
    
    # refil with zeros
    for x in scipy.arange(len(quarters)+1):
        sel = (tnew > qt_max[x]) * (tnew < qt_min[x+1])
        fnew[sel] = scipy.zeros(len(fnew[sel]))

    tf = atpy.Table()
    tf.add_column('time', tnew)
    tf.add_column('flux', fnew)
    tf.add_column('flux_raw', fnew_raw)

    return True, tf, qt_max, tablen, trange

###############################################################################################################

def acf_calc(time, flux, interval, kid, max_psearch_len):
    ''' Calculate ACF, calls error calc function'''
    
    # ACF calculation in pylab, close fig when finished
    pylab.figure(50)
    pylab.clf()
    lags, acf, lines, axis = pylab.acorr(flux, maxlags = max_psearch_len)
    pylab.close(50)

    lenlag = len(lags)
    lags = lags[int(lenlag/2.0):lenlag][:-1] * interval
    acf = acf[int(lenlag/2.0): lenlag][0:-1]
    ptout = pd.peakdetect(acf, lookahead = 100*interval, delta = 0.02)

    if len(ptout[0]) == 0: 
        t_acf = atpy.Table()
        t_acf.add_column('lags_days', lags)
        t_acf.add_column('acf', acf)
        t_acf.write('%s/ACF_output/dat_files/%s_ACF.ipac' % (dir, kid), overwrite = True)
        acf_per_pos = scipy.array([-9999])
        acf_per_height = scipy.array([-9999])
        acf_per_err = scipy.array([-9999])
        locheight = scipy.array([-9999])
        asym = scipy.array([-9999])
        return t_acf, acf_per_pos, acf_per_height, acf_per_err, locheight, asym
    
    max_inds = scipy.array(scipy.array(ptout[0])[:,0], dtype = 'int')
    max_height = scipy.array(ptout[0])[:,1]
    min_inds = scipy.array(ptout[1])[:,0]
    min_height = scipy.array(ptout[1])[:,1]

    maxnum = len(max_inds)
    if maxnum > 0:
        t = atpy.Table()
        t.add_column('ind', scipy.append(max_inds, min_inds))
        t.add_column('val', scipy.append(max_height, min_height))
        t.sort('ind')

        acf_per_pos = lags[max_inds]
        acf_per_height = max_height
        print 'Calculating errors and asymmetries...'
        # Calculate peak widths, asymmetries etc
        acf_per_err, locheight, asym= \
            calc_err(kid = kid, lags = lags, acf = acf, inds = t.ind, vals = t.val, maxnum = maxnum)

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
    t_lc.write('%s/ACF_output/dat_files/%s_LC.ipac' % (dir, kid), overwrite = True)

    t_acf = atpy.Table()
    t_acf.add_column('lags_days', lags)
    t_acf.add_column('acf', acf)
    t_acf.write('%s/ACF_output/dat_files/%s_ACF.ipac' % (dir, kid), overwrite = True)

    return  t_acf, acf_per_pos, acf_per_height, acf_per_err, locheight, asym

###############################################################################################################

def pgram_calc(time, flux, interval, kid, max_psearch_len):
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
    t_pg.write('%s/ACF_output/dat_files/%s_PGRAM.ipac' % (dir, kid), overwrite = True)
    
    return t_pg, sine_per, sine_height

###############################################################################################################

def calc_err(kid, lags, acf, inds, vals, maxnum):
    ''' Calculate peak widths, heights and asymmetries '''
    if len(inds) == 0: return -9999, -9999, -9999
    acf_per_err = scipy.ones(maxnum) * -9999
    asym = scipy.ones(maxnum) * -9999
    mean_height = scipy.ones(maxnum) * -9999

    acf_ind = scipy.r_[0:len(acf):1]
    num = len(vals)
    if maxnum*2 > num: maxnum -= 1    
    # loop through maxima, assuming 1st index is for minima
    for i in scipy.arange(maxnum):
        # find values and indices of centre left and right
        centre_v = vals[2*i+1]
        centre_i = inds[2*i+1]

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
    
    return acf_per_err, mean_height, asym

###############################################################################################################

def plot_stats(time, flux, kid_x, acf_per_pos_in, acf_per_height_in, acf_per_err_in, locheight_in, asym_in, doplot):
    ''' Plot and calculate statistics of peaks '''
    acf_per_pos_in = acf_per_pos_in[asym_in != -9999]
    acf_per_height_in = acf_per_height_in[asym_in != -9999]
    acf_per_err_in = acf_per_err_in[asym_in != -9999]
    locheight_in = locheight_in[asym_in != -9999]
    asym_in = asym_in[asym_in != -9999]

    if len(acf_per_pos_in) == 0: return -9999, -9999, -9999, -9999, -9999, -9999,\
        -9999, -9999, -9999, -9999, 0, -9999, -9999

    x = 4  #number of peaks used for calculating P
    hdet = 0  #start with 0 harmonic, set to 1 if 1/2P is 1st peak  
    # deal with cases where 1st peak is 1/2 real period
    if len(acf_per_pos_in) >= 3:
        ind = scipy.r_[1:len(acf_per_pos_in)+1:1]
        if (locheight_in[1] > locheight_in[0]) and (locheight_in[1] > locheight_in[2]):
            print '1/2 P found (1st)'
            hdet = 1  # mark harmonic found
            pk1 = acf_per_pos_in[1]
            acf_per_pos_in = acf_per_pos_in[1:]
            acf_per_height_in = acf_per_height_in[1:]
            acf_per_err_in = acf_per_err_in[1:]
            locheight_in = locheight_in[1:]
            asym_in = asym_in[1:]
        else:
            print 'no harmonic found or harmonic not 1st peak'
            pk1 = acf_per_pos_in[0]

    elif len(acf_per_pos_in) == 2:
        ind = scipy.r_[1:len(acf_per_pos_in)+1:1]
        if (locheight_in[1] > locheight_in[0]):
            print '1/2 P found (1st)'
            hdet = 1  # mark harmonic found
            pk1 = acf_per_pos_in[1]
            acf_per_pos_in = acf_per_pos_in[1:]
            acf_per_height_in = acf_per_height_in[1:]
            acf_per_err_in = acf_per_err_in[1:]
            locheight_in = locheight_in[1:]
            asym_in = asym_in[1:]
        else:
            print 'no harmonic found or harmonic not 1st peak'
            pk1 = acf_per_pos_in[0]

    else:
        print 'Only 1 peak'
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
        ind_gap = ind_keep[delta_keep_pos > 1.1*delta_keep_pos[0]]   
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

    if len(acf_per_pos) == 1 or len(acf_per_pos) == 2:
        return -9999, 0.0, pk1, acf_per_height[0], acf_per_err[0], locheight[0], -9999, -9999, -9999, -9999, 1, hdet, -9999

    # Fit line to peak positions
    
    fitfunc = lambda p, x: (p[0]*x) + 0.0
    errfunc = lambda p, x, y, err, fjac = None: [0, (fitfunc(p, x) - y) / err, None]
    # add 0.1-1 error scaling for height and same for width
    err_h = 1.1 - (locheight**2 / locheight.max()**2)
    fa = {'x': scipy.r_[0:len(acf_per_pos):1]+1, 'y': acf_per_pos, 'err': err_h}
    limited = scipy.array([0])
    limits = scipy.array([0.0])
    p0 = scipy.array([acf_per_pos[0]])

    m = mpfit.mpfit(errfunc, p0, functkw = fa, quiet = True)
    p1 = m.params
    perror = m.perror
    std_err  = scipy.sqrt( (scipy.sum((acf_per_pos - fitfunc(p1, scipy.r_[0:len(acf_per_pos):1]+1))**2) / (len(acf_per_pos)-2)) \
                           / scipy.sum( ((scipy.r_[0:len(acf_per_pos):1]+1) - scipy.mean(scipy.r_[0:len(acf_per_pos):1]+1))**2))
    med_per = p1
    mad_per_err = std_err

    if doplot == True:
        pylab.figure(51)
        pylab.clf()
        pylab.errorbar(scipy.r_[0:len(acf_per_pos):1]+1, acf_per_pos, yerr = acf_per_pos[0]*err_h, fmt = '.')
        pylab.plot(scipy.r_[0:len(acf_per_pos):1]+1, acf_per_pos, 'bo')
        pylab.xlim(0, 1.1*(len(acf_per_pos)))
        pylab.ylim(0, 1.03*(acf_per_pos.max()+acf_per_pos[0]*err_h.max()))
        pylab.plot(scipy.r_[0:len(acf_per_pos)+1:1], fitfunc(p1, scipy.r_[0:len(acf_per_pos)+1:1]), 'r--')
        pylab.xlabel('Peak')
        pylab.ylabel('Period (days)')
        pylab.title('ID: %s, P_grad = %.3f +/-%.3f days' %(kid_x, med_per, mad_per_err), fontsize = 16)
        pylab.savefig('%s/ACF_output/plots_stats/%s_perstats.png' % (dir, kid_x))
   
        pylab.figure(3,(12, 9))
        pylab.clf()
        pylab.subplot(3,2,1)
        pylab.errorbar(scipy.r_[0:len(acf_per_pos):1]+1, acf_per_pos, yerr = acf_per_pos[0]*err_h, fmt = '.')
        pylab.plot(scipy.r_[0:len(acf_per_pos):1]+1, acf_per_pos, 'bo')
        pylab.xlim(0, 1.1*(1+len(acf_per_pos)))
        pylab.ylim(0, 1.1*(acf_per_pos.max()+acf_per_pos[0]*err_h.max()))
        pylab.plot(scipy.r_[0:len(acf_per_pos)+1:1], fitfunc(p1, scipy.r_[0:len(acf_per_pos)+1:1]), 'r--')
   
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
        if doplot == True:
            pylab.plot(acf_per_pos, st_line([pqr[1],pqr[0]], acf_per_pos), 'r-')
            ax = pylab.gca()
            pylab.text(0.3, 0.9,'m=%.5f, TS=%.2f' %(h_grad, abs(h_timescale)), transform = ax.transAxes)
    else:
        h_grad = -9999
        h_timescale = -9999
        h_grad_err = -9999
        h_grad_scatter = -9999

    if doplot == True:
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
        if doplot == True:
            pylab.plot(acf_per_pos, st_line([pqr[1],pqr[0]], acf_per_pos), 'r-')
            ax = pylab.gca()
            pylab.text(0.3, 0.9, 'm=%.5f, TS=%.2f' %(w_grad, abs(w_timescale)), transform = ax.transAxes)
    else:
        w_grad = -9999
        w_timescale = -9999
        w_grad_err = -9999
        w_grad_scatter = -9999

    if doplot == True:
        pylab.subplot(3,2,5)
        pylab.plot(acf_per_pos, asym , 'bo')
        pylab.xlim(0, 1.1*acf_per_pos.max())
        pylab.plot(pk1, asym[acf_per_pos == pk1][0], 'ro')
        pylab.ylim(0.99*min(asym), 1.01*max(asym))
        pylab.ylabel('Asymmetry')
        pylab.xlabel('Lag (days)')
        pylab.suptitle('ID: %s, P_dlag = %.3f +/-%.3f d, P_pk = %.3fd' \
                       %(kid_x, med_per, mad_per_err, pk1), fontsize = 16)
        pylab.savefig('%s/ACF_output/plots_stats/%s_stats.png' % (dir, kid_x))

    return med_per, mad_per_err, pk1, acf_per_height[acf_per_pos == pk1][0], acf_per_err[acf_per_pos == pk1][0],\
        locheight[acf_per_pos == pk1][0], h_grad, h_grad_scatter, w_grad, w_grad_scatter, len(acf_per_pos), hdet, acf_per_pos

###############################################################################################################

def calc_var(kid = None, time_in = None, flux = None, period = None):
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

        var_arr_real = var_arr[(scipy.isfinite(var_arr) == True) * (var_arr != 0.0)]
        per_cent = per_cent[(scipy.isfinite(var_arr) == True) * (var_arr != 0.0)]
        var_med = scipy.median(var_arr_real)

        savefilen = '%s/ACF_output/dat_files/%s_amps.txt' % (dir, kid)
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
