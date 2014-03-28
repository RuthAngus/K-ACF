Kepler_ACF.py takes time, flux and flux_err and runs the acf.

run_ACF.py is the wrapper that runs Kepler_ACF.

savedata.py saves the results in a more useful format.

analysis.py analyses the rotation period results.

spotsim.py does the star spot simulations.

corot_ACF - This will return an ACF for a CoRoT lc (or any ascii file!) At the moment it follows the - 'if the ACF doesn't have a significant initial peak, or if there are no repeating peaks', return period = 0.0 - principle. To switch this off comment out lines 177, 182, 187
