import numpy as np
import matplotlib.pyplot as plt
from Kepler_ACF import corr_run

x, y, _, _ = np.genfromtxt("J06071008+2404374_xy_v2_ap2.0.dat",
                           skip_header=1).T

plt.clf()
plt.plot(x, y, "k.")
plt.savefig("lc")
