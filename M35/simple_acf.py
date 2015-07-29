import numpy as np
import matplotlib.pyplot as plt
import emcee
import glob

# dan's acf function
def dan_acf(x, axis=0, fast=False):
    """
    Estimate the autocorrelation function of a time series using the FFT.
    :param x:
        The time series. If multidimensional, set the time axis using the
        ``axis`` keyword argument and the function will be computed for every
        other axis.
    :param axis: (optional)
        The time axis of ``x``. Assumed to be the first axis if not specified.
    :param fast: (optional)
        If ``True``, only use the largest ``2^n`` entries for efficiency.
        (default: False)
    """
    x = np.atleast_1d(x)
    m = [slice(None), ] * len(x.shape)

    # For computational efficiency, crop the chain to the largest power of
    # two if requested.
    if fast:
        n = int(2**np.floor(np.log2(x.shape[axis])))
        m[axis] = slice(0, n)
        x = x
    else:
        n = x.shape[axis]

    # Compute the FFT and then (from that) the auto-correlation function.
    f = np.fft.fft(x-np.mean(x, axis=axis), n=2*n, axis=axis)
    m[axis] = slice(0, n)
    acf = np.fft.ifft(f * np.conjugate(f), axis=axis)[m].real
    m[axis] = 0
    return acf / acf[m]

def process(x, y):

    # fit straight line
    gap_days = 0.02043365
    AT = np.vstack((x, np.ones_like(x)))
    ATA = np.dot(AT, AT.T)
    m, b = np.linalg.solve(ATA, np.dot(AT, y))
    y = y - (m*x + b)

    # perform acf
    acf = dan_acf(y)

    # smooth with Gaussian kernel convolution
    Gaussian = lambda x, sig: 1./(2*np.pi*sig**.5) * \
                 np.exp(-0.5*(x**2)/(sig**2))
    #create the smoothing kernel
    conv_func = Gaussian(np.arange(-28,28,1.), 9.)
    acf_smooth = np.convolve(acf, conv_func, mode='same')

    # create 'lags' array
    lags = np.arange(len(acf))*gap_days

#     # trim acf
#     trim = 20
#     acf_smooth, lags = acf_smooth[trim:], lags[trim:]

    # find all the peaks
    peaks = np.array([i for i in range(1, len(lags)-1)
                     if acf_smooth[i-1] < acf_smooth[i] and
                     acf_smooth[i+1] < acf_smooth[i]])

    # find the lag of highest correlation
    m = acf_smooth ==  max(acf_smooth[peaks])
    highest = lags[m]

    # find the first and second peaks
    if acf_smooth[peaks[0]] > acf_smooth[peaks[1]]:
        period = lags[peaks[0]]
    else: period = lags[peaks[1]]
    print period

    plt.clf()
    plt.subplot(2, 1, 1)
    plt.plot(x, y, "k.")
    plt.subplot(2, 1, 2)
    for i in peaks:
        plt.axvline(lags[i], color="r")
    plt.axvline(highest, color="g")
    plt.axvline(period, color="k")
    plt.plot(lags, acf_smooth)
    plt.show()

if __name__ == "__main__":

    DIR = "."  # edit me!
    fnames = glob.glob("%s/*.dat" % DIR)

    for i, fname in enumerate(fnames[1:]):
        id = fname.split("/")[-1].split("_")[0]  # edit me!
        x, y, _, _ = np.genfromtxt(fname, skip_header=1).T
        yerr = np.ones_like(y) * 1e-5  # FIXME

        process(x, y)
