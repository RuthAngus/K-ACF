import numpy as np
from run_ACF import year_acf, all_acf
import glob
import scipy
import Kepler_ACF

def load_data(lc_file):
    data = np.genfromtxt(lc_file).T
    x = data[0]
    y = data[1]
    yerr = y*1e-5
    return x, y, yerr

if __name__ == "__main__":

        savedir = "/Users/angusr/angusr/Suz_simulations"
        loaddir = "/Users/angusr/angusr/Suz_simulations/final"

        lc_files = np.array(glob.glob("%s/lightcurve_*.txt" %loaddir))

        for kid in range(len(lc_files)):
            kid = str(kid)
            l = len(kid)
            if l == 1:
                lc_file = "%s/lightcurve_000%s.txt" %(loaddir, kid)
            elif l == 2:
                lc_file = "%s/lightcurve_00%s.txt" %(loaddir, kid)
            elif l == 3:
                lc_file = "%s/lightcurve_0%s.txt" %(loaddir, kid)
            elif l == 4:
                lc_file = "%s/lightcurve_%s.txt" %(loaddir, kid)

            time, flux, flux_err = load_data(lc_file)
            print lc_file
            if kid == '9':
                raw_input('enter')
            Kepler_ACF.corr_run(time, flux, flux_err, kid, 'sim', savedir)
