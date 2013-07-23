# This code must be run before ACF.py or similar. it has been used to index fits files from# multiple quarters# ss_index2 indexes individual quartersimport globimport numpy as npimport osimport os.pathimport pyfitsimport stringimport atpyimport timedef index(quarters):    """Produce ASCII index of all Kepler light curve files"""    # Locate all N2 light curve files (but not ones with REBIN suffix)    lcfile = np.array(glob.glob('/Users/angusr/angusr/ACF/star_spot_sim/%s/sim_*.png.mat' %quarters))     #lcfile = np.array(glob.glob('/Users/angusr/angusr/ACF/star_spot_sim/%s/testtesttest*.png.mat' %quarters))    Nfiles = len(lcfile)    print 'Found %d LC files' % Nfiles    kid = np.zeros(Nfiles, 'int32')    tmid = np.zeros(Nfiles, 'int32')    pq = np.zeros(Nfiles, 'int32')    aq = np.zeros(Nfiles, 'int32')    galaxy = np.zeros(Nfiles, 'int32')    variable = np.zeros(Nfiles, 'int32')    ra = np.zeros(Nfiles, 'float32') + np.nan    dec = np.zeros(Nfiles, 'float32') + np.nan    pmra = np.zeros(Nfiles, 'float32') + np.nan    pmdec = np.zeros(Nfiles, 'float32') + np.nan    gmag = np.zeros(Nfiles, 'float32') + np.nan    rmag = np.zeros(Nfiles, 'float32') + np.nan    imag = np.zeros(Nfiles, 'float32') + np.nan    zmag = np.zeros(Nfiles, 'float32') + np.nan    kepmag = np.zeros(Nfiles, 'float32') + np.nan    jmag = np.zeros(Nfiles, 'float32') + np.nan    hmag = np.zeros(Nfiles, 'float32') + np.nan    kmag = np.zeros(Nfiles, 'float32') + np.nan    gk = np.zeros(Nfiles, 'float32') + np.nan    teff = np.zeros(Nfiles, 'float32') + np.nan    logg = np.zeros(Nfiles, 'float32') + np.nan    feh = np.zeros(Nfiles, 'float32') + np.nan    ebv = np.zeros(Nfiles, 'float32') + np.nan    av = np.zeros(Nfiles, 'float32') + np.nan    rad = np.zeros(Nfiles, 'float32') + np.nan    glat = np.zeros(Nfiles, 'float32') + np.nan    glon = np.zeros(Nfiles, 'float32') + np.nan    pmtotal = np.zeros(Nfiles, 'float32') + np.nan    mtime = np.zeros(Nfiles, '15S')    module = np.zeros(Nfiles, 'int32')    output = np.zeros(Nfiles, 'int32')    quarter = np.zeros(Nfiles, 'int32')    long_cadence =np.zeros(Nfiles, 'int32')       ''' This is where you tell the table what the stars are called. A 1, 2, 3 or 4 digit \    name should be valid '''    for ss in range(0,len(lcfile)):        kid[ss] = str(ss + 1)        # if quarters < 10:        #     if lcfile[ss][47:48] == 'a':        #         kid[ss] = lcfile[ss][46:47]        #     elif lcfile[ss][48:49] == 'a':        #         kid[ss] = lcfile[ss][46:48]        #         #----------------------------        #     if lcfile[ss][47:48] == '.':        #         kid[ss] = lcfile[ss][46:47]        #     elif lcfile[ss][48:49] == '.':        #         kid[ss] = lcfile[ss][46:48]        # elif quarters > 9:        #     if lcfile[ss][48:49] == 'a':        #         kid[ss] = lcfile[ss][47:48]        #     elif lcfile[ss][49:50] == 'a':        #         kid[ss] = lcfile[ss][47:49]        #         #----------------------------        #     if lcfile[ss][48:49] == '.':        #         kid[ss] = lcfile[ss][47:48]        #     elif lcfile[ss][49:50] == '.':        #         kid[ss] = lcfile[ss][47:49]        # elif quarters == '3-6':        #     if lcfile[ss][48:49] == '.':        #         kid[ss] = lcfile[ss][47:48]        #     elif lcfile[ss][49:50] == '.':        #         kid[ss] = lcfile[ss][47:49]        #     elif lcfile[ss][50:51] == '.':        #         kid[ss] = lcfile[ss][47:50]        #     else:        #         kid[ss] = lcfile[ss][47:51]        # elif quarters == '7-10':        #     if lcfile[ss][49:50] == '.':        #         kid[ss] = lcfile[ss][48:49]        #     elif lcfile[ss][50:51] == '.':        #         kid[ss] = lcfile[ss][48:50]        #     elif lcfile[ss][51:52] == '.':        #         kid[ss] = lcfile[ss][48:51]        #     else:        #         kid[ss] = lcfile[ss][48:52]        # elif quarters == '11-14':        #     if lcfile[ss][50:51] == '.':        #         kid[ss] = lcfile[ss][49:50]        #     elif lcfile[ss][51:52] == '.':        #         kid[ss] = lcfile[ss][49:51]        #     elif lcfile[ss][52:53] == '.':        #         kid[ss] = lcfile[ss][49:52]        #     else:        #         kid[ss] = lcfile[ss][49:53]    print 'kid = ', kid    for i in np.arange(Nfiles):        print i, '/', Nfiles# Extract file name        print lcfile[i]    #Nobj = Nfiles    Nobj = np.size(np.unique(kid))    print 'Nobj = ',Nobj    print 'Found %d unique objects' % Nobj    # Store in ATpy table    t = atpy.Table()    t.add_column('keplerid', kid, description = 'Unique Kepler ID')    t.add_column('tmid', tmid, description = 'Unique 2MASS catalog ID')    t.add_column('ra', ra, description = 'RA (degrees, J2000)')    t.add_column('dec', dec, description = 'Dec (degrees, J2000)')    t.add_column('pmra', pmra, description = 'RA proper motion (as/yr)')    t.add_column('pmdec', pmdec, description = 'Dec proper motion (as/yr)')    t.add_column('gmag', gmag, description = 'g-band magnitude')    t.add_column('rmag', rmag, description = 'r-band magnitude')    t.add_column('imag', imag, description = 'i-band magnitude')    t.add_column('zmag', zmag, description = 'z-band magnitude')    t.add_column('jmag', jmag, description = '2MASS j-band magnitude')    t.add_column('hmag', hmag, description = '2MASS h-band magnitude')    t.add_column('kmag', kmag, description = '2MASS k-band magnitude')    t.add_column('kepmag', kepmag, description = 'Kepler-band magnitude')    t.add_column('gkcol', gk, description = 'g-k colour (mag)')    t.add_column('teff', teff, description = 'Effective temperature (K)')    t.add_column('logg', logg, description = 'Log10 surface gravity (cm/sec2)')    t.add_column('feh', feh, description = 'Log10 Fe/H metallicity (log(Fe/H))')    t.add_column('ebv', ebv, description = 'E(B-V) reddening (mag)')    t.add_column('av', av, description = 'A_V extinction (mag)')    t.add_column('rad', rad, description = 'Radius (solar radii)')    #t.add_column('pq', pq, description = 'Photometric quality indicator')    t.add_column('aq', aq, description = 'Astrophysical quality indicator')    t.add_column('galaxy', galaxy, description = 'Star/Galaxy indicator')    t.add_column('variable', variable, description = 'Constant/Variable indicator')    t.add_column('glat', glat, description = 'Galactic latitude')    t.add_column('glon', glon, description = 'Galactic longitude')    t.add_column('pmtotal', pmtotal, description = 'Total proper motion')    t.add_column('mtime', mtime, description = 'last modification time')    t.add_column('filename', lcfile, description = 'light curve filename')    t.add_column('module', module, description = 'CCD Module (1-25)')    t.add_column('output', output, description = 'CCD Output (1-4)')    t.add_column('quarter', quarter, description = 'Quarter number')    t.add_column('long_cadence', long_cadence, description = 'Long cadence flag (0-1)')    #t = t.where(t.long_cadence == 1)    print 'len', len(t)    print type(kid)        savefile = '/Users/angusr/angusr/ACF/index_Qss%stest.ipac' %quarters    t.describe()    print 'Saving index in %s' % savefile    t.write(savefile, overwrite = True)    print len(lcfile)    return