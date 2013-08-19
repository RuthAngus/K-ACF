# This is the version that reads in Kepler data, uses its time sampling and adds star spots to it.
# It does individual quarters
# This version produces star spots from the entire Kepler lc, then adds to individual quarters


import scipy
import scipy.io
import pylab
import numpy
import glob
import pyfits

ROOTDIR = '/Users/angusr/angusr/ACF/star_spot_sim'

#tau = 1.5
def mklc(nspot = 200, incl = (scipy.pi)*5./12., amp = 0.01, \
         tau = 30.5, diffrot = 0.0, \
         dur = 20.0, samp = 0.01, noise = 0.001, doplot = False, myperiod = [10.0], Amplitude = 1.0, quarter = 3):

         # spots.py calls runSim which calls mklc
    
    ''' This is a simplified version of the class-based routines in
    spot_model.py. It generates a light curves for dark, point like
    spots with no limb-darkening. 

    Parameters:
    nspot = desired number of spots present on star at any
            one time
    amp = desired light curve amplitude
    tau = characteristic spot life-time
    dur = light curve duration
    samp = time-sampling
    diffrot = fractional difference between equatorial and polar
              rotation period
    doplot: set to True to produce plots
    (unit of time is equatorial rotation period)'''

    # IMPORT KEPLER LIGHTCURVE

    name = '006370489'
    time2 = []; lightcurve = []
    quarter_times2 = []
    quarter_times_start = []
    quarter_times_fin = []
    files = glob.glob('/Users/angusr/angusr/data2/all_Qs/kplr%s-*llc.fits' % (name))
    #print len(files)
    for f in range(0,12):
        hdulist = pyfits.open(files[f])
        tbdata = hdulist[1].data #"add the first extension to tbdata"
        x = numpy.where(numpy.isfinite(tbdata['TIME']))
        time = tbdata['TIME'][x]
        lc = tbdata['PDCSAP_FLUX'][x]
        x = numpy.where(numpy.isfinite(tbdata['PDCSAP_FLUX']))
        time = tbdata['TIME'][x]
        lc = tbdata['PDCSAP_FLUX'][x]
        quarter_times2.append(tbdata['TIME'])
        quarter_times_start.append(time[0])
        quarter_times_fin.append(time[-1])
        lightcurve.extend(numpy.array(lc))
        time2.extend(numpy.array(time))

    
    
    # myperiod = float(myperiod[0])
    print 'Period = ', myperiod
    dur = (max(time2) - min(time2))

    # SET UP THE SPOTS
    print 'Setting up the spots...'

    # (crude estimate of) total number of spots needed during entire
    # time-series
    nspot_tot = int(nspot * dur / 2 / tau)
    # uniform distribution of spot longitudes
    lon = scipy.rand(nspot_tot) * 2 * scipy.pi 
   # distribution of spot latitudes uniform in sin(latitude)
    lat = scipy.arcsin(scipy.rand(nspot_tot)) 
    # spot rotation rate optionally depends on latitude



    
    period = ((scipy.sin(lat) - 0.5) * diffrot + 1.0 )*myperiod

    
    period0 = scipy.ones(nspot_tot)*myperiod
    # all spots have the same maximum area
    # (crude estimate of) filling factor needed per spot
    ff = amp / scipy.sqrt(nspot)
    scale_fac = 1
    amax = scipy.ones(nspot_tot) * ff * scale_fac
    # all spots have the evolution timescale
    decay = scipy.ones(nspot_tot) * tau
    # uniform distribution of spot peak times 
    # start well before and end well after time-series limits (to
    # avoid edge effects)
    extra = 3 * decay.max()
    pk = scipy.rand(nspot_tot) * (dur + 2 * extra) - extra

    # COMPUTE THE LIGHT CURVE
    print 'Computing the light curve...'
    

    time = numpy.array(time2- min(time2))
    addit = min(time2)

    
    npt = len(time)
    area_tot = scipy.zeros(npt)
    dF_tot = scipy.zeros(npt)
    dF_tot0 = scipy.zeros(npt)

   
    
    
    
    
    # add up the contributions of individual spots
    for i in range(nspot_tot):
        # Spot area
        if (pk[i] == 0) + (decay[i] == 0):
            area = scipy.ones(npt) * amax[i]
        else:
            area = amax[i] * \
                scipy.exp(-(time - pk[i])**2 / 2. / decay[i]**2)
        area_tot += area
        # Fore-shortening 
        phase = 2 * scipy.pi * time / period[i] + lon[i]
        phase0 = 2 * scipy.pi * time / period0[i] + lon[i]
        mu = scipy.cos(incl) * scipy.sin(lat[i]) + \
            scipy.sin(incl) * scipy.cos(lat[i]) * scipy.cos(phase)
        mu0 = scipy.cos(incl) * scipy.sin(lat[i]) + \
            scipy.sin(incl) * scipy.cos(lat[i]) * scipy.cos(phase0)
        mu[mu < 0] = 0.0
        mu0[mu0 < 0] = 0.0
        # Flux
        dF_tot -= area * mu
        dF_tot0 -= area * mu0

   
    
    print 'Adding noise...'
    # ADD NOISE
    noi = pylab.normal(0, 1, npt) * noise
    dF_tot += noi
    dF_tot0 += noi
        
    amp_eff = dF_tot.max()-dF_tot.min()
    nspot_eff = area_tot / scale_fac / ff

    if doplot == True:
        print 'Used %d spots in total over %d rotation periods.' \
            % (nspot_tot, dur)
        print 'Mean filling factor of individual spots was %.4f.' \
            % ff
        print 'Desired amplitude was %.4f, actual amplitude was %.4f.' \
            % (amp, amp_eff)
        print 'Desired number of spots at any one time was %d.' % nspot
        print 'Actual number of spots was %d (min), %d (average), %d (max).' \

        pylab.close(1)
        pylab.figure(1, (10, 4))
        xoff = 0.1
        xoffr = 0.05
        yoff = 0.13
        yoffr = 0.07
        xwi = (1.0 - xoff - xoffr)
        ywi = (1.0 - yoff - yoffr) / 2
        xll = xoff 
        yll = 1.0 - yoffr - ywi
        ax1 = pylab.axes([xll, yll, xwi, ywi])
        pylab.plot(time, area_tot * 100, 'k-')
        #pylab.title('i=%.1f <Ns>=%d A=%.4f Fs=%.4f tau=%.2fP sig=%.4f diff=%.1f' \
        #            % (incl * 180 / scipy.pi, nspot_eff.mean(), amp_eff, ff, tau, noise, diffrot))
        #pylab.title('i=%.1f <Ns>=%d tau=%.2fP Period = %d days' \
        #            % (incl * 180 / scipy.pi, nspot_eff.mean(), tau, myperiod)) 
        pylab.ylabel('fill. fac. (%)')
        yll = 1.0 - yoffr - 2 * ywi
        axc = pylab.axes([xll, yll, xwi, ywi], sharex = ax1)
        pylab.plot(time, 1 + dF_tot, 'k-', label = 'diff rot')
        pylab.plot(time, 1 + dF_tot0 - amp_eff, 'r-', label = 'no diff rot')
        pylab.ylabel('rel. flux')
        pylab.xlabel('time (days)')
        pylab.xlim(time.min(), time.max())
        pylab.legend()
        #pylab.savefig('/Users/angusr/angusr/ACF/star_spot_sim', )

    
    

    

    # Split into quarters
        
    print 'quarter = ', quarter
    time = time + addit
    #print 'start = ', quarter_times_start[quarter]
    #print 'stop = ', quarter_times_fin[quarter]
    print quarter_times_start
    print quarter_times_fin
    x = numpy.where(time == quarter_times_start[quarter])
    y = numpy.where(time == quarter_times_fin[quarter])

    start = int(x[0])
    stop = int(y[0])
    
    time = time[start:stop]
    area_tot = area_tot[start:stop]
    dF_tot = dF_tot[start:stop]
    dF_tot0 = dF_tot0[start:stop]
    lightcurve = lightcurve[start:stop]

    # Add normalised Kepler lc to simulated lc
    lightcurve = lightcurve/numpy.median(lightcurve)
    lightcurve = lightcurve - numpy.median(lightcurve)
    dF_tot = dF_tot - numpy.median(dF_tot)
    dF_tot0 = dF_tot0 - numpy.median(dF_tot0)
    #dF_tot = (dF_tot)*numpy.median(lightcurve) + lightcurve
    #dF_tot0 = (dF_tot0)*numpy.median(lightcurve) + lightcurve
    # Normalise
    #dF_tot = dF_tot/numpy.median(dF_tot)
    #dF_tot0 = dF_tot0/numpy.median(dF_tot0)

    npt = len(time)

    res0 = scipy.array([nspot_eff.mean(), ff, amp_eff, noise])
    res1 = scipy.zeros((4, npt))

    pylab.close(2)
    pylab.figure(2)
    pylab.subplot(3,1,1)
    pylab.plot(time, dF_tot, 'r.')
    pylab.ylim(min(dF_tot), max(dF_tot))
    pylab.subplot(3,1,2)
    pylab.ylim(min(lightcurve),max(lightcurve))
    pylab.plot(time, lightcurve, 'b.')
    pylab.subplot(3,1,3)
    pylab.ylim(min(lightcurve),max(lightcurve))
    dF_tot2 = dF_tot + lightcurve*Amplitude
    pylab.plot(time, dF_tot, 'g.')
    pylab.ylim(min(dF_tot), max(dF_tot))
    
    res1[0,:] = time
    res1[1,:] = area_tot
    res1[2,:] = dF_tot
    res1[3,:] = dF_tot0

    
    print 'Done'
    return res0, res1

def runSim(amp = 0.01, N = 1, doplot = True, dur = 50, number_of_data_sets = 3, myperiod = [10.0], star_name = 1, spot_tau = 1, Amplitude = 1.0):
    tau = spot_tau
    diffrot = scipy.zeros(N) + 0.5
    years = ['3-6', '7-10', '11-14'] 
    for year in years:
        for i in range(N): # One lightcurve per set of parameters
            pars, ts = mklc(tau, diffrot[i], doplot = doplot, \
                            dur = 50, myperiod = myperiod, quarter = quarter)
            if doplot == True:
                pylab.savefig('%s/%s/sim_%04dtest.png' % (ROOTDIR, (year), (star_name+1)))
                #scipy.io.savemat('%s/%s/orig_sim_%s.png' % (ROOTDIR, (q+3), (star_name+1)), \
                #scipy.io.savemat('%s/%s/testtesttest%s.png' % (ROOTDIR, (q+3), (star_name+1)), \
                #                 {'pars': pars, 'ts': ts})
                scipy.io.savemat('%s/%s/sim_%s.png' % (ROOTDIR, (year), (star_name+1)), \
                                 {'pars': pars, 'ts': ts})
            
    return pars 

# pars = [number of spots, filling factor, amplitude, noise], ts = [time, area_tot, dF_tot, dF_tot0] 

