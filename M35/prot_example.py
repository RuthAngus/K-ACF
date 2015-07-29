import numpy as np
import matplotlib.pyplot as plt
from Kepler_ACF import corr_run
import glob

DIR = "."  # edit me!
fnames = glob.glob("%s/*.dat" % DIR)

for i, fname in enumerate(fnames):
    id = fname.split("/")[-1].split("_")[0]  # edit me!
    x, y, _, _ = np.genfromtxt(fname, skip_header=1).T
    yerr = np.ones_like(y) * 1e-5  # FIXME
    corr_run(x, y, yerr, id, "acf")
