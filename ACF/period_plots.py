# This script takes the Periods.txt files from various directories and compiles lists of targets and various measurements of period. It then makes a plot for each star with all the measurements of period. It also tells you whether you have made a reliable period measurement for each star. Additional - It will plot either all_data, individual quarters or yearly data. 

import numpy as np
import scipy
import pylab
import atpy
mps = []
ids = []
nerrs = []

def period_extract(all_data = False, ind_quarters = False, years = False):

    final_KID_list = []
    # Load data
    data = np.genfromtxt('/Users/angusr/Desktop/astero_ages.txt').T
    all_targs = data[0]
    
    # Load results
    data3 = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ3_output/Periods_3.txt').T
    data4 = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ4_output/Periods_4.txt').T  
    data5 = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ5_output/Periods_5.txt').T
    data6 = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ6_output/Periods_6.txt').T
    data7 = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ7_output/Periods_7.txt').T
    data8 = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ8_output/Periods_8.txt').T
    data9 = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ9_output/Periods_9.txt').T
    data10 = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ10_output/Periods_10.txt').T
    data11 = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ11_output/Periods_11.txt').T
    data12 = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ12_output/Periods_12.txt').T
    data13 = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ13_output/Periods_13.txt').T
    data14 = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ14_output/Periods_14.txt').T
    data15 = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ15_output/Periods_15.txt').T
    data16 = np.genfromtxt('/Users/angusr/angusr/ACF/PDCQ16_output/Periods_16.txt').T
    KID3 = data3[0]
    periods3 = data3[1]
    errors3 = data3[2]
    sine3 = data3[3]
    KID4 = data4[0]
    periods4 = data4[1]
    errors4 = data4[2]
    sine4 = data4[3]
    KID5 = data5[0]
    periods5 = data5[1]
    errors5 = data5[2]
    sine5 = data5[3]
    KID6 = data6[0]
    periods6 = data6[1]
    errors6 = data6[2]
    sine6 = data6[3]
    KID7 = data7[0]
    periods7 = data7[1]
    errors7 = data7[2]
    sine7 = data7[3]
    KID8 = data8[0]
    periods8 = data8[1]
    errors8 = data8[2]
    sine8 = data8[3]
    KID9 = data9[0]
    periods9 = data9[1]
    errors9 = data9[2]
    sine9 = data9[3]
    KID10 = data10[0]
    periods10 = data10[1]
    errors10 = data10[2]
    sine10 = data10[3]
    KID11 = data11[0]
    periods11 = data11[1]
    errors11 = data11[2]
    sine11 = data11[3]
    KID12 = data12[0]
    periods12 = data12[1]
    errors12 = data12[2]
    sine12 = data12[3]
    KID13 = data13[0]
    periods13 = data13[1]
    errors13 = data13[2]
    sine13 = data13[3]
    KID14 = data14[0]
    periods14 = data14[1]
    errors14 = data14[2]
    sine14 = data14[3]
    KID15 = data15[0]
    periods15 = data15[1]
    errors15 = data15[2]
    sine15 = data15[3]
    KID16 = data16[0]
    periods16 = data16[1]
    errors16 = data16[2]
    sine16 = data16[3]                                       
    print 'Found %d Targets in quarter 3' %len(KID3)
    print 'Found %d Targets in quarter 4' %len(KID4)
    print 'Found %d Targets in quarter 5' %len(KID5)
    print 'Found %d Targets in quarter 6' %len(KID6)
    print 'Found %d Targets in quarter 7' %len(KID7)
    print 'Found %d Targets in quarter 8' %len(KID8)
    print 'Found %d Targets in quarter 9' %len(KID9)
    print 'Found %d Targets in quarter 10' %len(KID10)
    print 'Found %d Targets in quarter 11' %len(KID11)
    print 'Found %d Targets in quarter 12' %len(KID12)
    print 'Found %d Targets in quarter 13' %len(KID13)
    print 'Found %d Targets in quarter 14' %len(KID14)
    print 'Found %d Targets in quarter 15' %len(KID15)
    print 'Found %d Targets in quarter 16' %len(KID16)

        
    KID_method = [KID3, KID4, KID5, KID6, KID7, KID8, KID9, KID10, KID11, KID12, \
                  KID13, KID14, KID15, KID16]
    period_method = [periods3, periods4, periods5, periods6, periods7, periods8, \
                     periods9, periods10, periods11, periods12, periods13, periods14, periods15, periods16]
    error_method = [errors3, errors4, errors5, errors6, errors7, errors8, \
                    errors9 , errors10, errors11, errors12, errors13, errors14, errors15, errors16]
    sine_method = [sine3, sine4, sine5, sine6, sine7, sine8, sine9, sine10, \
                   sine11, sine12, sine13, sine14, sine15, sine16]  
    KID_meth_string = ['Q3', 'Q4', 'Q5', 'Q6', 'Q7', 'Q8', 'Q9', 'Q10', 'Q11', 'Q12', 'Q13', 'Q14', 'Q15', 'Q16']
    xtick_values = [1,2,3,4,5,6,7,8,9,10,11,12, 13, 14]
    xtick_labels = ['3','4','5','6','7','8','9','10','11','12','13','14', '15', '16']
    fig_dir = 'ind_qs_figs'

    # Loop over all targets
    all_KIDs = []
    all_periods = []
    all_errors = []
    all_sines = []
    for i in all_targs:
        #Assemble master list of all targets with periods, errors and sines. 
        KIDs, periods, errors, sines = find_periods(i, KID_method, period_method, error_method, sine_method)
        all_KIDs.append(KIDs)
        all_periods.append(periods)
        all_errors.append(errors)
        all_sines.append(sines)

        # Find medians
        acf_median = np.median(periods)
        acf_errs = np.median(errors)
        sine_median = np.median(sines)

        # Find rms about median
        acf_rms = np.sqrt((sum((periods - acf_median)**2))/np.float(len(periods)))
        sine_rms = np.sqrt((sum((sines - sine_median)**2))/np.float(len(sines)))
        diff = sine_median - acf_median

        # PLOTTING
        KID1s = range(1,len(KID_method)+1)
        fake_errors = np.zeros(len(KID_method))
        pylab.figure(1)
        pylab.clf()

        # Only plot if a period measurement was made. 
        for j in range(0,len(KID_method)):
            if periods[j] > 0. and errors[j] > 0.:
                pylab.errorbar(KID1s[j], periods[j], yerr = errors[j], marker = 'o', color = 'b', markersize = 5)
            else: pylab.axvline(j+1, linewidth = 0.5, color = 'r', linestyle = '--')
            if sines[j] > 0.:
                red_herring = 0
                #pylab.errorbar(KID1s[j], sines[j], yerr = fake_errors[j], marker = 'o', color = 'r', markersize = 5)
        pylab.ylabel('Period')
        pylab.xlabel('Quarters')
        #if max(sines)>max(periods):
        #    upper_y_lim = max(sines) + 5
        #else:
        upper_y_lim = max(periods) + 5
        pylab.xlim(0, len(KID_meth_string)+1)
        pylab.ylim(0,upper_y_lim)
        #pylab.axhline(sine_median, linewidth = 0.5, color = 'r')

        '''Adding harmonic lines'''
        pylab.axhline(acf_median, linewidth = 0.5, color = 'b')
        pylab.axhline(acf_median/2., linewidth = 0.5, color = 'b', linestyle = '--')
        period_multiple = 0; harmonic=2
        if acf_median != 0:
            while period_multiple < upper_y_lim:
                period_multiple = acf_median*harmonic
                pylab.axhline(period_multiple, linewidth = 0.5, color = 'b', linestyle = '--')
                harmonic += 1
        
        #if 0<sine_median<1000 and sine_median != acf_median:
            #red_herring = 0
            #pylab.text(10.1, sine_median, 'Sine', color = 'r')
        #if 0<acf_median<1000:
            #pylab.text(10.1, acf_median, 'ACF', color = 'b')
            
        KID_string = np.int(i)
        KID_string = np.str(KID_string)
        pylab.title('%s' %KID_string)
        
        '''Making sure that all measurements lie within 15% of the harmonic lines'''
        number_of_measurements = 0
        number_of_good_measurements = 0
        margin = acf_median/7.5
        for p in range(0,len(periods)):
            if periods[p] > 0.0 and errors[p] > 0:
                number_of_measurements += 1
                #if acf_median - acf_median/margin < periods[p] < acf_median + acf_median/margin:
                if acf_median - margin < periods[p] < acf_median + margin:
                    number_of_good_measurements += 1
                #elif acf_median*2.0 - acf_median*2.0/margin < periods[p] < acf_median*2.0 + acf_median*2.0/margin:
                elif acf_median*2.0 - margin < periods[p] < acf_median*2.0 + margin:
                    number_of_good_measurements += 1
                elif acf_median*3.0 - margin < periods[p] < acf_median*3.0 + margin:
                    number_of_good_measurements += 1
                #elif acf_median*0.5 - acf_median*0.5/margin < periods[p] < acf_median*0.5 + acf_median*0.5/margin:
                elif acf_median*0.5 - margin < periods[p] < acf_median*0.5 + margin:
                    number_of_good_measurements += 1
                ##elif  sine_median - sine_median/(margin*2) < periods[p] < sine_median + sine_median/(margin*2):
                ##number_of_good_measurements += 1
        #pylab.axhspan(acf_median - acf_median/margin, acf_median + acf_median/margin, facecolor='b', alpha=0.1)
        #pylab.axhspan(acf_median*2.0 - acf_median*2.0/margin, acf_median*2.0 + acf_median*2.0/margin, facecolor='b', alpha=0.1)
        #pylab.axhspan(acf_median*0.5 - acf_median*0.5/margin, acf_median*0.5 + acf_median*0.5/margin, facecolor='b', alpha=0.1)
        ##pylab.axhspan(sine_median - sine_median/(margin*2), sine_median + sine_median/(margin*2), facecolor='r', alpha=0.1)
        pylab.axhspan(acf_median - margin, acf_median + margin, facecolor='b', alpha=0.1)
        pylab.axhspan(acf_median*2.0 - margin, acf_median*2.0 + margin, facecolor='b', alpha=0.1)
        pylab.axhspan(acf_median*3.0 - margin, acf_median*3.0 + margin, facecolor='b', alpha=0.1)
        pylab.axhspan(acf_median*0.5 - margin, acf_median*0.5 + margin, facecolor='b', alpha=0.1)

        '''Requires that periods are measured in at least 2/3'''
        if years == True:
            relax = 1.0
        else:
            relax = 2./3.
        if number_of_measurements < int((len(KID_meth_string))*relax): 
            #pylab.text(0,0, 'MISSING MEASUREMENTS, require %s/%s' %( int((len(KID_meth_string))*2./3.), len(KID_meth_string)))
            pylab.axvspan(0, 120, facecolor='r', alpha=0.07)
            selected = False
            #require that 2/3 measurements are good!
        else:
            bonus = number_of_measurements - int((len(KID_meth_string))*relax) + 1
            outliers = number_of_measurements - number_of_good_measurements
            if outliers > bonus:
                pylab.axvspan(0, 20, facecolor='r', alpha=0.07)
                #pylab.text(0,0, 'INCONSISTENT MEASUREMENTS, require < %s outliers' %(outliers - bonus + 3))
                selected = False
            else:
                final_KID_list.append(KID_string)
                print acf_median, i
                mps.append(acf_median)
                ids.append(int(i))
                nerrs.append(acf_errs)
                selected = True
        
          
        pylab.xticks(xtick_values, xtick_labels)
        print KID_string
        #if years == True and selected == True:
            #print max(periods) - min(periods)
        pylab.savefig('/Users/angusr/angusr/ACF/%s/%s.pdf' %(fig_dir, KID_string))
        #raw_input('enter')

    print final_KID_list
    print '%s targets selected' %len(final_KID_list)
    if ind_quarters == True:
        txt_tit = 'ind_quarters'
    elif years == True:
        txt_tit = 'years'
    elif all_data == True:
        txt_tit = 'all_data'
    # np.savetxt('/Users/angusr/Desktop/results.txt', final_KID_list)
    np.savetxt('/Users/angusr/Desktop/results2.txt', (ids, mps, nerrs))
    return


def find_periods(KID, list_of_KID_methods, list_of_period_methods, list_of_error_methods, list_of_sine_methods):
    KID1 = KID
    KID_method = list_of_KID_methods
    period_method = list_of_period_methods
    error_method = list_of_error_methods
    sine_method = list_of_sine_methods

    KIDs = np.zeros(len(KID_method))
    periods = np.zeros(len(KID_method))
    errors = np.zeros(len(KID_method))
    sines = np.zeros(len(KID_method))
    n = 0; m = 0
    for KID_list in KID_method:
        for star in range(0,len(KID_list)):
            if KID_list[star] == KID1:
                KIDs[n] = KID_method[m][star]
                periods[n] = period_method[m][star]
                errors[n] = error_method[m][star]
                sines[n] = sine_method[m][star]
                n += 1
        m += 1

    return KIDs, periods, errors, sines

period_extract(ind_quarters = True)
