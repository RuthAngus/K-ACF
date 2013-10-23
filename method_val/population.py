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
