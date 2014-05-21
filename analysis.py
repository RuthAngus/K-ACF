# Read in and assess all period measurement data

import numpy as np
import matplotlib.pyplot as pl
import glob
from match import assemble
from mpl_toolkits.mplot3d import Axes3D

# # take a list of measurements (quarterly or yearly)
# # and decide whether they're allowed
# def allowance(pars, k, ps, perr, all_KIDs, all_ps, all_perrs, rmsyn):
#
#     m, f, l, d = pars
#
#     # find median and rms of quarter by quarter data
#     med = np.median(ps)
#     rms = np.sqrt(np.sum((ps-med)**2))
#
#     # only periods lying within fraction of median for each star
#     # or 1/2* and 2* the median
#     mr = m*med
#     harmonics = [.5, 1, 2.]
#     allowed, allowed_err = [], []
#     rms = np.zeros(len(harmonics))
#
#     # find the set of 'allowed' period measurements
#     for j, h in enumerate(harmonics):
#         a = (ps < h*(med+mr)) * (ps > h*(med-mr))
#         allowed.append(ps[a])
#         allowed_err.append(perr[a])
# #             rms[j] = np.sqrt(np.sum(ps[a] - med)**2))
#
#     allowed = [a for ll in allowed for a in ll]
#
#     # calculate the rms in each harmonic region
#     n = len(harmonics)
#     rms = np.zeros(n)
#     b = (med-(.75*med) > ps)
#     if len(ps[b]) > 0:
#         rms[0] = (np.sqrt((1./n)*np.sum((ps[b]-(med*harmonics[0]))**2)))
#     b = (ps < med+(1.5*med)) * (med-(.75*med) < ps)
#     if len(ps[b]) > 0:
#         rms[1] = (np.sqrt((1./n)*np.sum((ps[b]-(med*harmonics[1]))**2)))
#     b = (med+(1.5*med) < ps)
#     if len(ps[b]) > 0:
#         rms[2] = (np.sqrt((1./n)*np.sum((ps[b]-(med*harmonics[2]))**2)))
#     RMS = np.sqrt(np.sum(rms**2))
#
#     # find the fraction of accepted periods for each star
#     fraction = 0
#     if len(allowed) > 0:
#         fraction = float(len(allowed))/float(len(ps))
#
#     # only accept star if fraction is greater than f and
#     # there are more than l total period measurements
#     # then save the all-quarter period
#     if fraction > f and len(ps) > l: # and np.sqrt((med-all_ps[i])**2)<d*med:
#         accept = k
#         accept_p = all_ps
#         if RMS < all_perrs:
#             accept_perr = all_perrs
#         else:
#             if rmsyn == True: accept_perr = RMS
#             else: accept_perr = all_perrs
#
#     return accept, accept_p, accept_perr

def analyse(pars, DIR, loadfile, globfile, targlist, rmsyn):

    m, f, l, d = pars
    # load 'all' data
    data = np.genfromtxt(loadfile).T
    all_KIDs = data[0]
    all_ps = data[1]
    all_perrs = data[2]

    # load target list
    KIDs = targlist

    accept = []; accept_p = []; accept_perr = []
    # select results files
    lcfiles = glob.glob(globfile)
    for i, lc in enumerate(lcfiles):
        qdata = np.genfromtxt(lc).T
        ind_KIDs = qdata[0]
        k = ind_KIDs[0]
        ps = qdata[1]
        perr = qdata[2]

        # get rid of periods that are zero
        a = ps > 0
        ps = ps[a]
        perr = perr[a]

#         accept, accept_p, accept_perr = allowance(pars, k, ps, perr, all_KIDs[i], \
#                 all_ps[i], all_perrs[i], rmsyn=rmsyn)

        # find median and rms of quarter by quarter data
        med = np.median(ps)
#         rms = np.sqrt(np.mean(ps**2))
        rms = np.sqrt(np.sum((ps-med)**2))
#
        # only periods lying within fraction of median for each star
        # or 1/2* and 2* the median
        mr = m*med
        harmonics = [.5, 1, 2.]
        allowed, allowed_err = [], []
        rms = np.zeros(len(harmonics))
#
        # find the set of 'allowed' period measurements
        for j, h in enumerate(harmonics):
            a = (ps < h*(med+mr)) * (ps > h*(med-mr))
            allowed.append(ps[a])
            allowed_err.append(perr[a])
#             rms[j] = np.sqrt(np.sum(ps[a] - med)**2))
#
        allowed = [a for ll in allowed for a in ll]
#
        # calculate the rms in each harmonic region
        n = len(harmonics)
        rms = np.zeros(n)
#
        b = (med-(.75*med) > ps)
        if len(ps[b]) > 0:
            rms[0] = (np.sqrt((1./n)*np.sum((ps[b]-(med*harmonics[0]))**2)))
        b = (ps < med+(1.5*med)) * (med-(.75*med) < ps)
        if len(ps[b]) > 0:
            rms[1] = (np.sqrt((1./n)*np.sum((ps[b]-(med*harmonics[1]))**2)))
        b = (med+(1.5*med) < ps)
        if len(ps[b]) > 0:
            rms[2] = (np.sqrt((1./n)*np.sum((ps[b]-(med*harmonics[2]))**2)))
#
        RMS = np.sqrt(np.sum(rms**2))
#
#         # find the fraction of accepted periods for each star
#         if len(ps) > 0:
#             fraction = float(len(ps[a]))/float(len(ps))
#         else: fraction = 0
#
        # find the fraction of accepted periods for each star
        fraction = 0
        if len(allowed) > 0:
            fraction = float(len(allowed))/float(len(ps))
#
        # only accept star if fraction is greater than f and
        # there are more than l total period measurements
        # then save the all-quarter period
        if fraction > f and len(ps) > l: # and np.sqrt((med-all_ps[i])**2)<d*med:
            accept.append(k)
            accept_p.append(all_ps[i])
            if RMS < all_perrs[i]:
                accept_perr.append(all_perrs[i])
            else:
                if rmsyn == True: accept_perr.append(RMS)
                else: accept_perr.append(all_perrs[i])

    accept = np.array(accept)
    accept_p = np.array(accept_p)
    accept_perr = np.array(accept_perr)

#     a = accept_perr > 10.
#     for i in range(len(accept_p[a])):
#         print accept_p[a][i], accept_perr[a][i]

    # save results
#     np.savetxt("/Users/angusr/angusr/%s/periods.txt" %DIR, \
#             np.transpose((accept, accept_p, accept_perr)))

    return accept, accept_p, accept_perr

if __name__ == "__main__":

    # tuning params
    m = .2 # periods must lie within what fraction of the median?
    f = .7 # what fraction of total measurements required to be consistent?
    l = 0. # what number of quarters must be available?
    d = 100. # within what fraction of the median must the residual p measurement lie?

#     DIR = "injections"
#     DIR = "Suz_simulations"
    DIR = "ACF2"

    pars = [m, f, l, d]
    loadfile = "/Users/angusr/angusr/%s/results.txt" %DIR
    globfile = "/Users/angusr/angusr/%s/*_results.txt" %DIR

#     targlist = np.array(range(1, 1001))
#     targlist = np.array(range(800))
    targlist = np.genfromtxt("/Users/angusr/Python/Gyro/data/astero_targets.txt").T
    print len(targlist), 'targets in total'

    # run analysis on quarter by quarter
    KID, period, p_err = analyse(pars, DIR, loadfile, globfile, targlist, rmsyn=False)
    print len(KID), 'periods measured by us'

    # run analysis on year by year
    globfile = "/Users/angusr/angusr/%s/*yr_result.txt" %DIR
    pars = [.05, .9, 3., 100.]
    KIDyr, periodyr, p_erryr = \
            analyse(pars, DIR, loadfile, globfile, targlist, rmsyn=False)
    print len(KIDyr), 'yr periods measured by us'

    new = []; newp = []; newperr = []
    for i, kid in enumerate(KIDyr):
        if len(KID[KID==kid]) == 0:
            new.append(KIDyr[i])
            newp.append(periodyr[i])
            newperr.append(p_erryr[i])

    KID = np.concatenate((KID, np.array(new)))
    period = np.concatenate((period, np.array(newp)))
    p_err = np.concatenate((p_err, np.array(newperr)))

    data = assemble(KID, period, p_err).T
    KID = data[0]
    period = data[1]
    temp = data[3]
    age = data[13]
    print len(period), 'stars with periods'

#     # 3d plot
#     a = temp > 0
#     print len(temp[a]), 'with non-zero temps'
#     pl.clf()
#     fig = pl.figure()
#     ax = fig.gca(projection='3d')
#     ax.scatter(period[a], temp[a], age[a], c = 'b', marker = 'o')
#     ax.set_xlabel('Rotational period (days)')
#     ax.set_ylabel('bv')
#     ax.set_zlabel('Age (Gyr)')
#     pl.show()
