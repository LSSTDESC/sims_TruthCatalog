import os
import numpy as np
import h5py
import sqlite3
import sys
from datetime import datetime as dt

from multiprocessing import Process, Lock, BoundedSemaphore

from GCR import GCRQuery
import GCRCatalogs

import lsst.sims.photUtils as sims_photUtils
from lsst.sims.catUtils.dust import EBVbase

#import desc.sims_truthcatalog.conversion_utils as conversion_utils

MAX_PARALLEL= 10    # to be tuned..  May depend on chunk size
KNL_FACTOR = 8      # KNL is about 8 times slower than Cori. Adjust 

_global_cosmoDC2_data = {}

def _logit(log_lock, log_path, to_write):
    if log_lock is not None:
        log_lock.acquire()            # timeout?
    with open(log_path, 'a') as out_file:
        out_file.write(to_write)
    if log_lock is not None:
        log_lock.release()
def _good_indices(galaxies, bad_gs):
    '''
    Return a list of indices for the good galaxies
    Parameters:

    galaxies:     list of galaxy ids, monotone increasing
    bad_gs:       list of galaxy ids of galaxies deemed bad. Monotone increasing

    Return:       list of indices referring to those elements in galaxies
                  which are not in bad_gs
    '''
    bad_ix = 0
    good_ixes = []
    
    for gal_ix in range(0, len(galaxies) ):
        while (bad_gs[bad_ix] < galaxies[gal_ix]):
            bad_ix += 1
            if bad_ix == len(bad_gs) :    # no more bad galaxies
                good_ixes += [g for g in range(gal_ix, len(galaxies))]
                return good_ixes
        if galaxies[gal_ix] < bad_gs[bad_ix] :  # a good one
            good_ixes.append(gal_ix)
    return good_ixes

def _write_sqlite(dbfile, galaxy_ids, ra, dec, redshift, flux_by_band_MW,
                  flux_by_band_noMW, good_ixes):
    from desc.sims_truthcatalog.conversion_utils import write_column_descriptions
    with sqlite3.connect(dbfile) as conn:
        #conversion_utils.write_column_descriptions(conn)
        write_column_descriptions(conn)
        cursor = conn.cursor()

        cmd = '''CREATE TABLE IF NOT EXISTS truth_summary
        (id BIGINT, host_galaxy BIGINT, ra DOUBLE, dec DOUBLE,
        redshift FLOAT, is_variable INT, is_pointsource INT,
        flux_u FLOAT, flux_g FLOAT, flux_r FLOAT, 
        flux_i FLOAT, flux_z FLOAT, flux_y FLOAT,
        flux_u_noMW FLOAT, flux_g_noMW FLOAT, flux_r_noMW FLOAT, 
        flux_i_noMW FLOAT, flux_z_noMW FLOAT, flux_y_noMW FLOAT)'''
        cursor.execute(cmd)
        conn.commit()
        #print("Created table if not exists truth_summary")
        values = ((int(galaxy_ids[i_obj]),int(-1),
                   ra[i_obj],dec[i_obj],
                   redshift[i_obj], 0, 0,
                   flux_by_band_MW['u'][i_obj], flux_by_band_MW['g'][i_obj],
                   flux_by_band_MW['r'][i_obj], flux_by_band_MW['i'][i_obj],
                   flux_by_band_MW['z'][i_obj], flux_by_band_MW['y'][i_obj],
                   flux_by_band_noMW['u'][i_obj], flux_by_band_noMW['g'][i_obj],
                   flux_by_band_noMW['r'][i_obj], flux_by_band_noMW['i'][i_obj],
                   flux_by_band_noMW['z'][i_obj], flux_by_band_noMW['y'][i_obj])
                  for i_obj in good_ixes)

        cursor.executemany('''INSERT INTO truth_summary
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''',
                           values)
        conn.commit()

    #sys.stdout.flush()

def _process_chunk(db_lock, log_lock, sema, sed_fit_name, cosmoDC2_data,
                   first_gal, self_dict, bad_gals):
    """
    Do all chunk-specific work:  compute table contents for a
    collection of galaxies and write to db

    Parameters
    ----------
    db_lock          Used to avoid conflicts writing to sqlite output
    log_lock         Used to avoid conflicts writing to per-healpixel log
    sema             A semaphore. Release when done
    sed_fit_name     File where sed fits for this healpixel are
    cosmoDC2_data    Values from cosmoDC2 for this healpixel, keyed by
                     column name
    first_gal        index of first galaxy in our chunk (in sed fit list)
    self_dict        Random useful values stored in GalaxyTruthWriter
    bad_gals         List of galaxy ids, monotone increasing, to be
                     skipped
    """

    dry = self_dict['dry']
    chunk_size = self_dict['chunk_size']
    dbfile = self_dict['dbfile']
    logfile = self_dict['logfile']

    if dry:
        _logit(log_lock, logfile,
               '_process_chunk invoke for first_gal {}, chunk size {}'.format(first_gal, chunk_size))
        if sema is None:
            return
        sema.release()

        #exit(0)
        return
    
    lsst_bp_dict = self_dict['lsst_bp_dict']    
    galaxy_ids = []
    ra = []
    dec = []
    redshift = []
    ebv_vals = None
    ebv_vals_init = False    # does this belong somewhere else?
    ccm_w = None
    total_gals = self_dict['total_gals']


    chunk_start = first_gal
    chunk_end = min(first_gal + chunk_size, total_gals)
    with h5py.File(sed_fit_name, 'r') as sed_fit_file:

        sed_names = sed_fit_file['sed_names'][()]
        sed_names = [s.decode() for s in sed_names]  # becse stored as bytes

        gals_this_chunk = chunk_end - chunk_start
        subset = slice(chunk_start, chunk_end)
        galaxy_ids = sed_fit_file['galaxy_id'][()][subset]
        to_log = 'Start with galaxy #{}, id={}\n# galaxies for _process_chunk: {}\n'.format(first_gal, galaxy_ids[0], len(galaxy_ids))
        _logit(log_lock, logfile, to_log)

        # get the cross-match between the sed fit and cosmoDC2
        cosmo_len = len(cosmoDC2_data['galaxy_id'])

        crossmatch_dex = np.searchsorted(cosmoDC2_data['galaxy_id'],
                                         galaxy_ids)
        np.testing.assert_array_equal(
            galaxy_ids, cosmoDC2_data['galaxy_id'][crossmatch_dex])

        ra = sed_fit_file['ra'][()][subset]
        dec = sed_fit_file['dec'][()][subset]
        np.testing.assert_array_equal(ra,
                                      cosmoDC2_data['ra'][crossmatch_dex])
        np.testing.assert_array_equal(dec,
                                      cosmoDC2_data['dec'][crossmatch_dex])

        good_ixes = _good_indices(galaxy_ids.tolist(), bad_gals[0])
        if (len(good_ixes) == 0):
            if sema is not None:
                sema.release()
            return
        else:
            _logit(log_lock, logfile,
                   'Found {} good indices for chunk starting with {}\n'.format(len(good_ixes), chunk_start))
        flux_by_band_MW = {}
        flux_by_band_noMW = {}

        # Calculate E(B-V) for dust extinction in Milky Way along relevant
        # lines of sight
        band_print="Processing band {}, first gal {}, time {}\n"
        if not ebv_vals_init:
            equatorial_coords = np.array([np.radians(ra), np.radians(dec)])
            ebv_model = EBVbase()
            ebv_vals = ebv_model.calculateEbv(equatorialCoordinates=equatorial_coords,
                                              interp=True)
            ebv_vals_init = True

        for i_bp, bp in enumerate('ugrizy'):
            if (i_bp == 0 or i_bp == 5):
                _logit(log_lock, logfile,
                       band_print.format(bp, first_gal, dt.now()))
            fluxes_noMW = {}
            fluxes = {}
            for component in ['disk', 'bulge']:
                fluxes_noMW[component] = np.zeros(gals_this_chunk, dtype=float)
                fluxes[component] = np.zeros(gals_this_chunk, dtype=float)

            for component in ['disk', 'bulge']:
                #print("   Processing component ", component)
                sed_arr = sed_fit_file['%s_sed' % component][()][subset]
                av_arr = sed_fit_file['%s_av' % component][()][subset]
                rv_arr = sed_fit_file['%s_rv' % component][()][subset]
                mn_arr = sed_fit_file['%s_magnorm' % component][()][i_bp,:][subset]
                z_arr = cosmoDC2_data['redshift'][crossmatch_dex]
                gii = 0
                done = False
                for i_gal, (s_dex, mn, av,
                            rv, zz, ebv) in enumerate(zip(sed_arr, mn_arr,
                                                          av_arr, rv_arr,
                                                          z_arr, ebv_vals)):
                    if done: break
                    while good_ixes[gii] < i_gal :
                        gii += 1
                        if gii == len(good_ixes):   # ran out of good ones
                            done = True
                            break
                    if done: break
                    if good_ixes[gii] > i_gal :  # skipped over it; it's bad
                        continue
                    # Leave space for it in the arrays, but values
                    # for all the fluxes will be left at 0

                    # read in the SED file from the library
                    sed_file_name = os.path.join(self_dict['sed_lib_dir'],
                                                 sed_names[s_dex])
                    sed = sims_photUtils.Sed()
                    sed.readSED_flambda(sed_file_name)

                    # find and apply normalizing flux
                    fnorm = sims_photUtils.getImsimFluxNorm(sed, mn)
                    sed.multiplyFluxNorm(fnorm)

                    # add internal dust
                    if ccm_w is None or not np.array_equal(sed.wavelen, ccm_w):
                        ccm_w = np.copy(sed.wavelen)
                        a_x, b_x = sed.setupCCM_ab()
                    sed.addDust(a_x, b_x, A_v=av, R_v=rv)

                    # apply redshift
                    sed.redshiftSED(zz, dimming=True)

                    # flux, in Janskys, without Milky Way dust extinction
                    f_noMW = sed.calcFlux(lsst_bp_dict[bp])

                    # apply Milky Way dust
                    # (cannot reuse a_x, b_x because wavelength grid changed
                    # when we called redshiftSED)
                    a_x_mw, b_x_mw = sed.setupCCM_ab()
                    sed.addDust(a_x_mw, b_x_mw, R_v=3.1, ebv=ebv)

                    f_MW = sed.calcFlux(lsst_bp_dict[bp])

                    fluxes_noMW[component][i_gal] = f_noMW
                    fluxes[component][i_gal] = f_MW
                if (component == 'disk') and (bp == 'r'):
                    redshift = z_arr

            total_fluxes = fluxes_noMW['disk'] + fluxes_noMW['bulge']
            total_fluxes_MW = fluxes['disk'] + fluxes['bulge']

            dummy_sed = sims_photUtils.Sed()

            # add magnification due to weak lensing
            kappa = cosmoDC2_data['convergence'][crossmatch_dex]
            gamma_sq = (cosmoDC2_data['shear_1'][crossmatch_dex]**2
                            + cosmoDC2_data['shear_2'][crossmatch_dex]**2)
            magnification = 1.0/((1.0-kappa)**2-gamma_sq)
            magnified_fluxes = magnification*total_fluxes
            magnified_fluxes_MW = magnification*total_fluxes_MW
            flux_by_band_noMW[bp] = magnified_fluxes
            flux_by_band_MW[bp] = magnified_fluxes_MW

    #  Open connection to sqlite db and write
    #print('Time before db write is {}, first gal={}'.format(dt.now(), first_gal))
    #sys.stdout.flush()
    if not db_lock.acquire(timeout=120.0):
        _logit(log_lock, logfile,
               "Failed to acquire db lock, first gal=", first_gal)
        if sema is None:
            return
        sema.release()
        exit(1)

    try:
        _write_sqlite(dbfile, galaxy_ids, ra, dec, redshift,
                      flux_by_band_MW, flux_by_band_noMW, good_ixes)
        db_lock.release()
        if sema is not None:
            sema.release()

        _logit(log_lock, logfile,
               'Time after db write: {}, first_gal={}\n'.format(dt.now(),first_gal))
        exit(0)
    except Exception as ex:
        db_lock.release()
        if sema is not None:
            sema.release()
        raise(ex)
        
class GalaxyTruthWriter(object):
    '''
    Writes truth catalog for static galaxies to sqlite3 file for a 
    single healpixel
    '''
    def __init__(self, output_dir, hpid, sed_fit_dir, mag_cut, chunk_size,
                 start=0, nchunk=None, parallel=10, dry=False, call=False,
                 knl=False):
        #self.sed_fit_dir = os.path.join(sed_fit_dir,
        #                                'DC2/cosmoDC2_v1.1.4/sedLookup')
        self.sed_fit_dir = sed_fit_dir
        assert os.path.isdir(self.sed_fit_dir)
        self.sed_fit_name = os.path.join(self.sed_fit_dir,
                                         'sed_fit_%d.h5' % hpid)
        assert os.path.isfile(self.sed_fit_name)
        assert os.path.isdir(output_dir)
        self.dbfile = os.path.join(output_dir,
                                   'truth_summary_hp{}.sqlite3'.format(hpid))
        self.logfile = os.path.join(output_dir,
                                    'truth_summary_log_hp{}.txt'.format(hpid))
        logfile = self.logfile
        self.hpid = hpid
        self.mag_cut = mag_cut
        self.chunk_size=chunk_size
        self.start = start
        self.nchunk = nchunk
        self.dry = dry
        self.parallel = parallel
        self.call = call
        self.knl = knl

        logstring='GalaxyTruthWriter invoked with arguments\n  output-dir={}\n'
        logstring += 'healpixel={}\nsed_fit_dir={}\nmag_cut={}\nchunk_size={}\nchunk={}\nparallel={}\ndry\n'
        _logit(None, logfile, logstring.format(output_dir, hpid, sed_fit_dir,
                                               mag_cut,chunk_size, nchunk,
                                               parallel,dry))

        
    def do_indices(dbfile):
        with sqlite3.connect(dbfile) as conn:
            cursor = conn.cursor()    
            cmd = '''CREATE INDEX IF NOT EXISTS  ra_dec ON truth_summary (ra,dec)'''
            cursor.execute(cmd)
            cmd = '''CREATE UNIQUE INDEX IF NOT EXISTS gal_id ON truth_summary (id)'''
            cursor.execute(cmd)
            conn.commit
            #print("created indexes")

    def write(self):
        """
        Do all the 'shared' work (get cosmoDC2 information, figure
        out how to chunk,  etc.)
        """
        _logit(None, self.logfile,
               'Enter GalaxyTruthWriter.write, time {}\n'.format(dt.now()))
        db_lock = Lock()
        log_lock = Lock()
        self.log_lock = log_lock
        sema = BoundedSemaphore(self.parallel)
        self_dict = {}
        self_dict['dry'] = self.dry
        self_dict['chunk_size'] = self.chunk_size
        self_dict['dbfile'] = self.dbfile
        self_dict['logfile'] = self.logfile

        bad_gals = []
        cosmoDC2_data = {}

        # read in LSST bandpasses
        lsst_bp_dict=sims_photUtils.BandpassDict.loadTotalBandpassesFromFiles()
        
        if self.dry:
            self.total_gals = 17e6
            bad_gals = list()
            self.cosmoDC2_data = {}
        else:
            cat = GCRCatalogs.load_catalog('cosmoDC2_v1.1.4_image')

            _logit(log_lock, self.logfile,"cosmoDC2 catalog loaded")
            # get galaxy_id and redshift for crossmatching with SED fit files;
            # we will also get the magnitudes that should be reproduced
            # by our synthetic photometry (hp_query makes sure we only load
            # the healpixel we are interested in)
            hp_query = GCRQuery('healpix_pixel==%d' % self.hpid)
            cosmoDC2_data = cat.get_quantities(
                ['galaxy_id', 'redshift', 'ra', 'dec',
                 'mag_true_u_lsst', 'mag_true_g_lsst', 'mag_true_r_lsst',
                 'mag_true_i_lsst', 'mag_true_z_lsst', 'mag_true_y_lsst',
                 'mag_u_lsst', 'mag_g_lsst', 'mag_r_lsst',
                 'mag_i_lsst', 'mag_z_lsst', 'mag_y_lsst',
                 'shear_1', 'shear_2', 'convergence'],
                native_filters=[hp_query])


            # make sure cosmoDC2_data is sorted by galaxy_id
            sorted_dex = np.argsort(cosmoDC2_data['galaxy_id'])
            for colname in cosmoDC2_data.keys():
                cosmoDC2_data[colname] = cosmoDC2_data[colname][sorted_dex]
            #self.cosmoDC2_data = cosmoDC2_data
            _global_cosmoDC2_data = cosmoDC2_data
            
            bad_ixes = np.where(cosmoDC2_data['mag_r_lsst'][()]> self.mag_cut)
            bad_gals = list([cosmoDC2_data['galaxy_id'][()][i] for i in bad_ixes])
          
            self_dict['lsst_bp_dict'] = lsst_bp_dict

            # the parent directory for the SED library
            self_dict['sed_lib_dir'] = os.environ['SIMS_SED_LIBRARY_DIR']
        
            # Find total number of galaxies in the healpixel
            with h5py.File(self.sed_fit_name, 'r') as sed_fit_file:
                self.total_gals = len(sed_fit_file['galaxy_id'][()])

        max_chunk = int(np.ceil(self.total_gals/float(self.chunk_size)))

        if self.nchunk != None:
            max_chunk = min(max_chunk, self.nchunk)

        _logit(log_lock, self.logfile, 'max_chunk is {}\n'.format(max_chunk))

        max_parallel = min(self.parallel, max_chunk)

        self_dict['total_gals'] = self.total_gals

        first = self.start
        starts = tuple(first + i*self.chunk_size for i in range(max_chunk))
        remaining = max_chunk

        # For Haswell allow 6000 seconds for 50k chunk size, which is
        # conservative.   Multiply by 8 or so for knl
        arch_scale = 0.12
        if self.knl: arch_scale = KNL_FACTOR * arch_scale
        to = np.ceil(self.chunk_size * arch_scale) + 20
        _logit(log_lock, self.logfile,
               "Start forking at {}. Semaphore timeout value is: {}\n".format(dt.now(), to))

        # May eliminate this, allong with --call option
        if self.call:
            for i in range(max_chunk):
                _process_chunk(db_lock, Log_lock, None, self.sed_fit_name,
                               _global_cosmoDC2_data, starts[i], self_dict,
                               bad_gals)
            return

        # Use semaphore to keep from starting too many processes
        p_list = []
        for i in range(max_chunk):
            if (not sema.acquire(timeout=to)):
                _logit(log_lock, self.logfile,
                       "Unable to obtain process slot before timeout")
                for pp in p_list:
                    if pp.is_alive():
                        pp.terminate()
                    pp.join()
                exit(1)
            p = Process(target=_process_chunk,name='chunk_{}'.format(i),
                        args=(db_lock, log_lock, sema, self.sed_fit_name,
                              cosmoDC2_data,
                              starts[i],self_dict, bad_gals))
            p.start()
            p_list.append(p)

        # Now wait for the last ones to complete.   Should be able to
        # acquire the original semaphore count
        for i in range(self.parallel):
            if (not sema.acquire(timeout=to)):
                _logit(log_lock, self.logfile,
                       "Final batch of processes did not complete in time")
                for pp in p_list:
                    if pp.is_alive():
                        pp.terminate()
                    pp.join()
                exit(1)

        if not self.dry:
            # Create indexes
            #  Maybe first check that everything got processed?
            # Or even do this with a separate script or at least class
            # method?
            pass
            #do_indices(self.dbfile)
 
        _logit(log_lock, self.logfile, 'done at time {}\n'.format(dt.now()))

