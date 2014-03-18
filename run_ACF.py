# This is the wrapper that calls KeplerACF.py

import numpy as np
import glob
from Kepler_ACF import corr_run

id_list = ['3544595','10124866']
lc_file = np.array(glob.glob('/Users/angusr/angusr/Kepler/%s/kplr*_llc.fits' %id_list[0]))

corr_run(id_list[0], lc_file[0], 0)

# for i, file in enumerate(lc_file):
#     corr_run(id_list[0], lc_file, quarter)
