spots8 - wrapper.  

spot_sim_ruth5 - generates simulated light curves using spot model with a range of spot lifetimes, number of spots, etc.

ss_index2 - indexes simulated light curves for the ACF code.

ACF_star_spot8 - calculates ACFs for the simulated lightcurves. 

ss_period_plots3 - candidate selection process. A star will only be selected if a reliable and consisted period is measured in two thirds or more of the Kepler quarters. 

corot_ACF - This will return an ACF for a CoRoT lc (or any ascii file!) At the moment it follows the - 'if the ACF doesn't have a significant initial peak, or if there are no repeating peaks', return period = 0.0 - principle. To switch this off comment out lines 177, 182, 187

spots_no_grid - calculates completeness, reliability and contamination for 1000 simulated light curves with different, randomly selected rotation periods, spot lifetimes and amplitudes. 
