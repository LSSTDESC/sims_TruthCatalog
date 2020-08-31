"""
Module to write truth tables for lensed AGNs in DC2 Run3.0i.
"""
import sys
from collections import namedtuple
import logging
import sqlite3
import numpy as np
import pandas as pd
from lsst.sims.photUtils import PhotometricParameters
from lsst.sims.utils import angularSeparation
from .synthetic_photometry import find_sed_file, SyntheticPhotometry
from .agn_truth import agn_mag_norms


__all__ = ['write_lensed_agn_truth_summary',
           'write_lensed_agn_variability_truth']


logging.basicConfig(format="%(asctime)s %(name)s: %(message)s",
                    stream=sys.stdout)


def write_lensed_agn_truth_summary(lensed_agn_truth_cat, outfile,
                                   bands='ugrizy', verbose=False):
    """
    Write the truth_summary table for the lensed AGNs.

    Parameters
    ----------
    lensed_agn_truth_cat: str
        The sqlite3 file containing the model parameters for the lensed AGNs.
    outfile: str
        Filename of the output sqlite3 file.
    bands: list-like ['ugrizy']
        LSST bands.
    verbose: bool [False]
        Verbosity flag.  This flag has no effect in this function.
    """
    logger = logging.getLogger('write_lensed_agn_truth_summary')
    if verbose:
        logger.setLevel(logging.INFO)

    logger.info('processing %s', lensed_agn_truth_cat)
    table_name = 'truth_summary'
    create_table_sql = f'''CREATE TABLE IF NOT EXISTS {table_name}
        (id TEXT, host_galaxy BIGINT, ra DOUBLE, dec DOUBLE,
        redshift FLOAT, is_variable INT, is_pointsource INT,
        flux_u FLOAT, flux_g FLOAT, flux_r FLOAT,
        flux_i FLOAT, flux_z FLOAT, flux_y FLOAT,
        flux_u_noMW FLOAT, flux_g_noMW FLOAT, flux_r_noMW FLOAT,
        flux_i_noMW FLOAT, flux_z_noMW FLOAT, flux_y_noMW FLOAT)'''
    sed_file = find_sed_file('agnSED/agn.spec.gz')
    with sqlite3.connect(lensed_agn_truth_cat) as conn, \
         sqlite3.connect(outfile) as output:
        # Create the output table if it does not exist.
        output.cursor().execute(create_table_sql)
        output.commit()
        # Query for the columns containing the model info for each SN
        # from the lensed_sne table.
        query = '''select unique_id, ra, dec, redshift, magnorm,
                   magnification, av_mw, rv_mw from lensed_agn'''
        cursor = conn.execute(query)
        # Loop over each object and write its fluxes for all relevant
        # visits to the output table.
        values = []
        for unique_id, ra, dec, z, magnorm, magnification, gAv, gRv in cursor:
            synth_phot = SyntheticPhotometry(sed_file, magnorm, z)
            row = [unique_id, -1, ra, dec, z, 1, 1]
            fluxes_noMW = {_: synth_phot.calcFlux(_)*magnification
                           for _ in bands}
            synth_phot.add_dust(gAv, gRv, 'Galactic')
            for band in bands:
                row.append(synth_phot.calcFlux(band)*magnification)
            for band in bands:
                row.append(fluxes_noMW[band])
        output.cursor().executemany(f'''insert into {table_name} values
                                    (?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                                     ?, ?, ?, ?, ?, ?, ?, ?, ?)''', values)


def write_lensed_agn_variability_truth(opsim_db_file, lensed_agn_truth_cat,
                                       outfile, fp_radius=2.05, bands='ugrizy',
                                       start_mjd=59580.,
                                       end_mjd=61395.,
                                       agn_walk_start_date=58350.,
                                       verbose=False,
                                       num_objects=None):
    """
    Write the lensed AGN fluxes to the lensed_agn_variabilty_truth table.

    Parameters
    ----------
    opsim_db_file: str
        OpSim db file.  This will be the minion 1016 db file that was
        modified by DESC for DC2.
    lensed_agn_truth_cat: str
        The sqlite3 file containing the model parameters for the lensed AGNs.
    outfile: str
        Filename of the output sqlite3 file.
    fp_radius: float [2.05]
        Radius in degrees of the smallest acceptance cone containing the
        LSST focalplane projected onto the sky.
    bands: list-like or string ['ugrizy']
        The LSST bands.
    start_mjd: float [59580.]
        Starting MJD for opsim db query.  The default is the nominal starting
        date for DC2.
    end_mjd: float [61395.]
        Ending MJD for opsim db query. The default is the end of year 5.
        Starting date for visits to process.
    agn_walk_start_date: float [58350.]
        The starting MJD for AGN light curves.  This starts 240 days
        earlier than the nominal minion 1016 start date to accommodate
        AGN time delays.  It must match the value set at https://github.com/LSSTDESC/SLSprinkler/blob/v1.0.0/scripts/dc2/dc2_utils/variability.py#L22
    verbose: bool [False]
        Verbosity flag.
    num_objects: int [None]
        Maximum number of objects from the input catalog to process.
        If None, then process all objects.
    mjd_max: float [61395]
    """
    logger = logging.getLogger('write_lensed_agn_variability_truth')
    if verbose:
        logger.setLevel(logging.INFO)
    logger.info('processing %s', lensed_agn_truth_cat)
    table_name = 'lensed_agn_variability_truth'
    create_table_sql = f'''create table if not exists {table_name}
                           (id TEXT, obsHistID INT, MJD FLOAT, bandpass TEXT,
                            delta_flux FLOAT, num_photons FLOAT)'''
    sed_file = find_sed_file('agnSED/agn.spec.gz')
    # Read the opsim db data into a dataframe.
    with sqlite3.connect(opsim_db_file) as conn:
        opsim_df = pd.read_sql(
            f'''select obsHistID, descDitheredRA, descDitheredDec,
                expMJD, filter from Summary where expMJD >= {start_mjd}
                and expMJD < {end_mjd}''', conn)
    opsim_df['ra'] = np.degrees(opsim_df['descDitheredRA'])
    opsim_df['dec'] = np.degrees(opsim_df['descDitheredDec'])

    # Loop over objects in the truth_cat containing the model
    # parameters and write the fluxes to the output table for the
    # relevant visits.
    par_table_name = 'lensed_agn'
    colnames = (['unique_id', 'ra', 'dec', 'redshift', 't_delay', 'magnorm',
                 'magnification', 'seed']
                + [f'agn_tau_{_}' for _ in bands]
                + [f'agn_sf_{_}' for _ in bands]
                + ['av_mw', 'rv_mw'])
    query = f'select {",".join(colnames)} from {par_table_name}'
    AGNParams = namedtuple('AGNParams', colnames)
    with sqlite3.connect(lensed_agn_truth_cat) as conn, \
         sqlite3.connect(outfile) as output:
        # Create the output table if it does not exist.
        output.cursor().execute(create_table_sql)
        output.commit()
        # Get the number of AGN entries in the input db file.
        cursor = conn.execute(f'select count(*) from {par_table_name}')
        nobjs = list(cursor)[0][0]
        if num_objects is None:
            num_objects = nobjs
        # Query for the columns containing the model info for each AGN
        # from the lensed_agn table.
        cursor = conn.execute(query)
        # Loop over each object and write its fluxes for all relevant
        # visits to the output table.
        for i, row in zip(range(num_objects), cursor):
            pars = AGNParams(*row)._asdict()
            logger.info('%s  %d  %d', pars['unique_id'], i, num_objects)
            decmin, decmax = pars['dec'] - fp_radius, pars['dec'] + fp_radius
            df = pd.DataFrame(opsim_df.query(f'{decmin} <= dec <= {decmax}'))
            df['ang_sep'] = angularSeparation(df['ra'].to_numpy(),
                                              df['dec'].to_numpy(),
                                              pars['ra'], pars['dec'])
            df = df.query(f'ang_sep <= {fp_radius}')
            if len(df) == 0:
                continue
            # Compute baseline fluxes in each band.
            synth_phot0 = SyntheticPhotometry(sed_file, pars['magnorm'],
                                              redshift=pars['redshift'],
                                              iAv=0, gAv=pars['av_mw'],
                                              gRv=pars['rv_mw'])
            flux0 = {_: synth_phot0.calcFlux(_) for _ in bands}
            values = []
            for band in bands:
                phot_params = PhotometricParameters(nexp=1, exptime=30, gain=1,
                                                    bandpass=band)
                bp = synth_phot0.bp_dict[band]
                my_df = df.query(f'filter == "{band}"')
                mjds = my_df['expMJD'].to_numpy() - pars['t_delay']
                mag_norms = (agn_mag_norms(mjds,
                                           pars['redshift'],
                                           pars[f'agn_tau_{band}'],
                                           pars[f'agn_sf_{band}'],
                                           pars['seed'],
                                           start_date=agn_walk_start_date)
                             + pars['magnorm'])
                for visit, mjd, mag_norm in zip(my_df['obsHistID'],
                                                my_df['expMJD'], mag_norms):
                    synth_phot = SyntheticPhotometry(sed_file, mag_norm,
                                                     redshift=pars['redshift'],
                                                     iAv=0, gAv=pars['av_mw'],
                                                     gRv=pars['rv_mw'])
                    delta_flux = (pars['magnification']*
                                  (synth_phot.calcFlux(band) - flux0[band]))
                    num_photons = (pars['magnification']
                                   *synth_phot.sed.calcADU(bp, phot_params))
                    values.append((pars['unique_id'], visit, mjd, band,
                                   delta_flux, num_photons))
            output.cursor().executemany(f'''insert into {table_name} values
                                            (?, ?, ?, ?, ?, ?)''', values)
