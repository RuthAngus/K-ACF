# This script takes the Periods.txt files from various directories and compiles lists of targets and various measurements of period. It then makes a plot for each star with all the measurements of period. It also tells you whether you have made a reliable period measurement for each star. Additional - It will plot either all_data, individual quarters or yearly data. This is the false lc version. You should need the period finder because there are no missing stars

#This version does the above for the fake star spot data.

# This version is for individual quarters


import numpy
import scipy
import pylab
import atpy

def period_extract(all_data = False, ind_quarters = False, years = False, no_stars = 24):

    final_KID_list = []
    ''' importing file of all KIDs'''
    all_targs = range(0,no_stars)
    #all_targs = ['01','02','03','04','05','06','07','08','09','10', ]
    #all_targs = numpy.genfromtxt('/Users/angusr/Documents/rotation/all_targets.txt').T

    ''' Extracting data from text files.'''
    data3 = numpy.genfromtxt('/Users/angusr/angusr/ACF/PDCQss3_output/Periods_ss3test.txt').T
    data4 = numpy.genfromtxt('/Users/angusr/angusr/ACF/PDCQss4_output/Periods_ss4test.txt').T    
    data5 = numpy.genfromtxt('/Users/angusr/angusr/ACF/PDCQss5_output/Periods_ss5test.txt').T
    data6 = numpy.genfromtxt('/Users/angusr/angusr/ACF/PDCQss6_output/Periods_ss6test.txt').T
    data7 = numpy.genfromtxt('/Users/angusr/angusr/ACF/PDCQss7_output/Periods_ss7test.txt').T
    data8 = numpy.genfromtxt('/Users/angusr/angusr/ACF/PDCQss8_output/Periods_ss8test.txt').T
    data9 = numpy.genfromtxt('/Users/angusr/angusr/ACF/PDCQss9_output/Periods_ss9test.txt').T
    data10 = numpy.genfromtxt('/Users/angusr/angusr/ACF/PDCQss10_output/Periods_ss10test.txt').T
    data11 = numpy.genfromtxt('/Users/angusr/angusr/ACF/PDCQss11_output/Periods_ss11test.txt').T
    data12 = numpy.genfromtxt('/Users/angusr/angusr/ACF/PDCQss12_output/Periods_ss12test.txt').T
    data13 = numpy.genfromtxt('/Users/angusr/angusr/ACF/PDCQss13_output/Periods_ss13test.txt').T
    data14 = numpy.genfromtxt('/Users/angusr/angusr/ACF/PDCQss14_output/Periods_ss14test.txt').T
    # data_q36 = numpy.genfromtxt('/Users/angusr/angusr/ACF/PDCQss3-6_output/Periods_ss3-6test.txt').T
    # data_q710 = numpy.genfromtxt('/Users/angusr/angusr/ACF/PDCQss7-10_output/Periods_ss7-10test.txt').T
    # data_q1114 = numpy.genfromtxt('/Users/angusr/angusr/ACF/PDCQss11-14_output/Periods_ss11-14test.txt').T

    KID3 = data3[0]; periods3 = data3[1]; errors3 = data3[2]; sine3 = data3[3]
    KID4 = data4[0]; periods4 = data4[1]; errors4 = data4[2]; sine4 = data4[3]
    KID5 = data5[0]; periods5 = data5[1]; errors5 = data5[2]; sine5 = data5[3]
    KID6 = data6[0]; periods6 = data6[1]; errors6 = data6[2]; sine6 = data6[3]
    KID7 = data7[0]; periods7 = data7[1]; errors7 = data7[2]; sine7 = data7[3]
    KID8 = data8[0]; periods8 = data8[1]; errors8 = data8[2]; sine8 = data8[3]
    KID9 = data9[0]; periods9 = data9[1]; errors9 = data9[2]; sine9 = data9[3]
    KID10 = data10[0]; periods10 = data10[1]; errors10 = data10[2]; sine10 = data10[3]
    KID11 = data11[0]; periods11 = data11[1]; errors11 = data11[2]; sine11 = data11[3]
    KID12 = data12[0]; periods12 = data12[1]; errors12 = data12[2]; sine12 = data12[3]
    KID13 = data13[0]; periods13 = data13[1]; errors13 = data13[2]; sine13 = data13[3]
    KID14 = data14[0]; periods14 = data14[1]; errors14 = data14[2]; sine14 = data14[3]
    # KIDq710 = data_q710[0]; periodsq710 = data_q710[1]; errorsq710 = data_q710[2]; sineq710 = data_q710[3]
    # KIDq1114 = data_q1114[0]; periodsq1114 = data_q1114[1]; errorsq1114 = data_q1114[2]; sineq1114 = data_q1114[3]
    # KIDq36 = data_q36[0]; periodsq36 = data_q36[1]; errorsq36 = data_q36[2]; sineq36 = data_q36[3]

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
    # print 'Found %d Targets in quarters 3-6' %len(KIDq36)
    # print 'Found %d Targets in quarters 7-10' %len(KIDq710)
    # print 'Found %d Targets in quarters 11-14' %len(KIDq1114)

    if all_data == True:
        KID_method = [KID3, KID4, KID5, KID6, KID7, KID8, KID9, KID10, KID11, KID12, KID13, KID14, KIDq04, KIDq25, KIDq59, KIDq36, KIDq710, KIDq1114, KIDq09]
        period_method = [periods3, periods4, periods5, periods6, periods7, periods8, periods9, periods10, periods11, periods12, periods13, periods14, periodsq04, periodsq25, periodsq59, periodsq36, periodsq710, periodsq1114, periodsq09]
        error_method = [errors3, errors4, errors5, errors6, errors7, errors8, errors9 , errors10, errors11, errors12, errors13, errors14, errorsq04, errorsq25, errorsq59, errorsq36, errorsq710, errorsq1114, errorsq09]
        sine_method = [sine3, sine4, sine5, sine6, sine7, sine8, sine9, sine10, sine11, sine12, sine13, sine14, sineq04, sineq25, sineq59, sineq36, sineq710, sineq1114, sineq09]  
        KID_meth_string = ['Q3', 'Q4', 'Q5', 'Q6', 'Q7', 'Q8', 'Q9', 'Q10', 'Q11', 'Q12', 'Q13', 'Q14', 'Qs 1-4', 'Qs 2-5', 'Qs 5-9', 'Qs 3-6', 'Qs 7-10', 'Qs 11-14', 'Qs 1-8']
        xtick_values = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
        xtick_labels = ['3','4','5','6','7','8','9','10','11','12','13','14','1-4','2-5', '3-6', '5-9','7-10','11-14','1-8']
        fig_dir = 'ss_all_data_figs'
        
    elif ind_quarters == True:
        KID_method = [KID3, KID4, KID5, KID6, KID7, KID8, KID9, KID10, KID11, KID12, KID13, KID14]
        period_method = [periods3, periods4, periods5, periods6, periods7, periods8, periods9, periods10, periods11, periods12, periods13, periods14]
        error_method = [errors3, errors4, errors5, errors6, errors7, errors8, errors9 , errors10, errors11, errors12, errors13, errors14]
        sine_method = [sine3, sine4, sine5, sine6, sine7, sine8, sine9, sine10, sine11, sine12, sine13, sine14]  
        KID_meth_string = ['Q3', 'Q4', 'Q5', 'Q6', 'Q7', 'Q8', 'Q9', 'Q10', 'Q11', 'Q12', 'Q13', 'Q14']
        xtick_values = [1,2,3,4,5,6,7,8,9,10,11,12]
        xtick_labels = ['3','4','5','6','7','8','9','10','11','12','13','14']
        fig_dir = 'ss_ind_qs_figs'

    elif years == True:
        KID_method = [KIDq36, KIDq710, KIDq1114]
        period_method = [periodsq36, periodsq710, periodsq1114]
        error_method = [errorsq36, errorsq710, errorsq1114]
        sine_method = [sineq36, sineq710, sineq1114]  
        KID_meth_string = ['Qs 3-6', 'Qs 7-10', 'Qs 11-14']
        xtick_values = [1,2,3]
        xtick_labels = ['3-6','7-10','11-14']
        fig_dir = 'ss_year_figs'


    '''Loop over all targets'''
    for i in range(0,len(all_targs)):
        
        if years == True:
            periods = [periodsq36[i], periodsq710[i], periodsq1114[i]]
            errors = [errorsq36[i], errorsq710[i], errorsq1114[i]]

        elif ind_quarters == True:
            periods = [periods3[i], periods4[i], periods5[i], periods6[i], periods7[i], periods8[i], periods9[i], periods10[i], periods11[i], \
                       periods12[i], periods13[i], periods14[i]]
            errors = [errors3[i], errors4[i], errors5[i], errors6[i], errors7[i], errors8[i], errors9[i], errors10[i], errors11[i], \
                      errors12[i], errors13[i], errors14[i]]
        
        ''' Find medians '''
        acf_median = numpy.median(periods)
        
        '''PLOTTING'''
        KID1s = range(1,len(KID_method)+1)
        pylab.figure(1)
        pylab.clf()

        '''only plot if a period measurement was made. '''
        for j in range(0,len(KID_method)):
            if periods[j] > 0. and errors[j] > 0.:
                pylab.errorbar(KID1s[j], periods[j], yerr = errors[j], marker = 'o', color = 'b', markersize = 5)
            else: pylab.axvline(j+1, linewidth = 0.5, color = 'r', linestyle = '--')
        pylab.ylabel('Period')
        upper_y_lim = max(periods) + 5
        pylab.xlim(0, len(KID_meth_string)+1)
        pylab.ylim(0,upper_y_lim)
        
        '''Adding harmonic lines'''
        pylab.axhline(acf_median, linewidth = 0.5, color = 'b')
        pylab.axhline(acf_median/2., linewidth = 0.5, color = 'b', linestyle = '--')
        period_multiple = 0; harmonic=2
        
        if acf_median > 0:
            while period_multiple < upper_y_lim:
                period_multiple = acf_median*harmonic
                pylab.axhline(period_multiple, linewidth = 0.5, color = 'b', linestyle = '--')
                harmonic += 1
            
        KID_string = numpy.int(i)
        KID_string = numpy.str(KID_string)
        pylab.title('%s' %KID_string)

        
        '''Making sure that all measurements lie within 15% of the harmonic lines'''
        number_of_measurements = 0
        number_of_good_measurements = 0
        margin = acf_median/7.5
        for p in range(0,len(periods)):
            if periods[p] > 0.0 and errors[p] > 0:
                number_of_measurements += 1
                if acf_median - margin < periods[p] < acf_median + margin:
                    number_of_good_measurements += 1
                elif acf_median*2.0 - margin < periods[p] < acf_median*2.0 + margin:
                    number_of_good_measurements += 1
                elif acf_median*3.0 - margin < periods[p] < acf_median*3.0 + margin:
                    number_of_good_measurements += 1
                elif acf_median*0.5 - margin < periods[p] < acf_median*0.5 + margin:
                    number_of_good_measurements += 1
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
            pylab.text(0,0, 'MISSING MEASUREMENTS, require %s/%s' %( int((len(KID_meth_string))*2./3.), len(KID_meth_string)))
            pylab.axvspan(0, 120, facecolor='r', alpha=0.07)
            selected = False
            #require that 2/3 measurements are good!
        else:
            bonus = number_of_measurements - int((len(KID_meth_string))*relax) + 1
            outliers = number_of_measurements - number_of_good_measurements
            if outliers > bonus:
                pylab.axvspan(0, 20, facecolor='r', alpha=0.07)
                pylab.text(0,0, 'INCONSISTENT MEASUREMENTS, require < %s outliers' %(outliers - bonus + 3))
                selected = False
            else:
                final_KID_list.append((numpy.float(KID_string))+1)
                selected = True
        
          
        pylab.xticks(xtick_values, xtick_labels)
        #print KID_string
        #if years == True and selected == True:
            #print max(periods) - min(periods)
        pylab.savefig('/Users/angusr/angusr/ACF/%s/%s.pdf' %(fig_dir, KID_string))
        
        '''The difference between newfigs2 and newfigs3 is that newfigs3 has data that has had the ACF code run on it more recently newfigs5 has all data, all_data_figs has all data, ind_qs_figs will have individual quarters only and year_figs will have just long term data'''

    
    print final_KID_list
    print '%s targets selected' %len(final_KID_list)
    if ind_quarters == True:
        txt_tit = 'ss_ind_quarterstest'
    elif years == True:
        txt_tit = 'ss_yearstest'
    elif all_data == True:
        txt_tit = 'ss_all_datatest'
    numpy.savetxt('/Users/angusr/angusr/ACF/star_spot_sim/%s.txt' %txt_tit, final_KID_list)
    return
