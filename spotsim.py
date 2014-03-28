import scipy
import scipy.io
import pylab
import numpy as np
import glob
import pyfits
import kepler_quarters
from load_data import load

def mklc(x, nspot, incl, amp, tau, diffrot, dur, samp, noise, myperiod, \
        Amplitude, doplot = False):

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

    print 'Period = ', myperiod
    dur = (max(x) - min(x))

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

    time = np.array(x - min(x))
    addit = min(x)

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

    # Normalise
    dF_tot = dF_tot/np.median(dF_tot) - 1.
    dF_tot0 = dF_tot0/np.median(dF_tot0) -1.

    res0 = scipy.array([nspot_eff.mean(), ff, amp_eff, noise])
    res1 = scipy.zeros((4, npt))

    res1[0,:] = time
    res1[1,:] = area_tot
    res1[2,:] = dF_tot
    res1[3,:] = dF_tot0

    print 'Done'
    return res0, res1

if __name__ == "__main__":

    datadir = "/Users/angusr/angusr/data2/all_Qs"
    savedir = "/Users/angusr/angusr/injections"

    nspot = 200
    incl = np.pi*5./12.
    amp = 0.01
    diffrot = 0.5
    dur = 50.
    samp = 0.01
    noise = 0.001
    Amplitude = 1.

    n = 2
    ps = 10**(np.random.rand(n)*2.) # FIXME: yes Dan, these should be np.exp
    taus = 10**np.random.rand(n)
    myamp = 10**(np.random.rand(n)*1.18)
    names = ['010972873', '011137075', '011502597', '011717120', '005607242', \
             '006370489', '006442183', '006531928', '006603624', '009574283']

    for i in range(n):
        for KID in names:
            lc_files = np.array(glob.glob("%s/kplr%s*"%(datadir, KID)))
            for q, lc_file in enumerate(lc_files):

                # load data
                x, y, yerr = load(lc_file)

                # generate simulated lcs
                doplot = True
                pars, ts = mklc(x, nspot, incl, amp, taus[i], diffrot, dur, samp, noise, \
                        ps[i], Amplitude, doplot)

                if doplot == True:
                    pylab.savefig('%s/%s_%s'%(savedir, (i+1), q))

                # save raw simulations
                scipy.io.savemat('%s/sim_%s_%s_raw'%(savedir, (i+1), q), \
                        {'pars': pars, 'ts': ts})

                # add to real Kepler lcs and save
                nspots, ff, amp, noise = pars
                # ts = [time, area covered in spots, flux, flux w diff rot]
                time, area, flux, diff_flux = ts
                ts[2] = y + ts[2]*myamp[i]
                lc = np.ndarray((len(ts)+1, len(ts[0])))
                lc[1:,:] = ts
                lc[-1:,:] = yerr
                scipy.io.savemat('%s/sim_%s_%s'%(savedir, (i+1), q), \
                        {'pars': pars, 'ts': lc})

                # save file of parameters
                np.savetxt('%s/pars_%s_%s.txt'%(savedir, (i+1), q), \
                        np.transpose((float(KID), ps[i], taus[i], myamp[i], nspots, \
                        ff, amp, noise)))

    # all quarters
    for i in range(n):
        for KID in names:

            # load data
            lc_files = np.array(glob.glob("%s/kplr%s*"%(datadir, KID)))
            x, y, yerr = load(lc_files[0])
            for q in range(1, len(lc_files)):
                x = np.concatenate((x, (load(lc_files[q])[0])))
                y = np.concatenate((y, (load(lc_files[q])[1])))
                yerr = np.concatenate((yerr, (load(lc_files[q])[2])))

            # generate simulated lcs
            doplot = True
            pars, ts = mklc(x, nspot, incl, amp, taus[i], diffrot, dur, samp, noise, \
                    ps[i], Amplitude, doplot)

            if doplot == True:
                pylab.savefig('%s/%s_%s'%(savedir, (i+1), 'all'))

            # save raw simulations
            scipy.io.savemat('%s/sim_%s_%s_raw'%(savedir, (i+1), 'all'), \
                    {'pars': pars, 'ts': ts})

            # add to real Kepler lcs and save
            nspots, ff, amp, noise = pars
            # ts = [time, area covered in spots, flux, flux w diff rot]
            time, area, flux, diff_flux = ts
            ts[2] = y + ts[2]*myamp[i]
            lc = np.ndarray((len(ts)+1, len(ts[0])))
            lc[1:,:] = ts
            lc[-1:,:] = yerr
            scipy.io.savemat('%s/sim_%s_%s'%(savedir, (i+1), 'all'), \
                    {'pars': pars, 'ts': lc})

            # save file of parameters
            np.savetxt('%s/pars_%s_%s.txt'%(savedir, (i+1), 'all'), \
                    np.transpose((float(KID), ps[i], taus[i], myamp[i], nspots, \
                    ff, amp, noise)))

    # years
    nyrs = nqs = 4
    for i in range(n):
        for KID in names:
            for yr in range(nyrs):
                # load data
                lc_files = np.array(glob.glob("%s/kplr%s*"%(datadir, KID)))
                x, y, yerr = load(lc_files[yr*4])
                for q in range(yr*nrs, (yr*nyrs)+(nqs-1)):
                    x = np.concatenate((x, (load(lc_files[q+1])[0])))
                    y = np.concatenate((y, (load(lc_files[q+1])[1])))
                    yerr = np.concatenate((yerr, (load(lc_files[q+1])[2])))

                # generate simulated lcs
                doplot = True
                pars, ts = mklc(x, nspot, incl, amp, taus[i], diffrot, dur, samp, noise, \
                        ps[i], Amplitude, doplot)

                if doplot == True:
                    pylab.savefig('%s/%s_yr%s'%(savedir, (i+1), (yr+1)))

                # save raw simulations
                scipy.io.savemat('%s/sim_%s_yr%s_raw'%(savedir, (i+1), (yr+1)), \
                        {'pars': pars, 'ts': ts})

                # add to real Kepler lcs and save
                nspots, ff, amp, noise = pars
                # ts = [time, area covered in spots, flux, flux w diff rot]
                time, area, flux, diff_flux = ts
                ts[2] = y + ts[2]*myamp[i]
                lc = np.ndarray((len(ts)+1, len(ts[0])))
                lc[1:,:] = ts
                lc[-1:,:] = yerr
                scipy.io.savemat('%s/sim_%s_yr%s'%(savedir, (i+1), (yr+1)), \
                        {'pars': pars, 'ts': lc})

                # save file of parameters
                np.savetxt('%s/pars_%s_yr%s.txt'%(savedir, (i+1), (yr+1)), \
                        np.transpose((float(KID), ps[i], taus[i], myamp[i], nspots, \
                        ff, amp, noise)))
