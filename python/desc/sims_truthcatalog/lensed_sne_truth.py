"""
Module to write truth tables for lensed SNe in DC2 Run3.0i.
"""
import sys
import logging
import sqlite3
import numpy as np
import pandas as pd
from lsst.sims.photUtils import PhotometricParameters
from lsst.sims.utils import angularSeparation
from .sne_truth import SNSynthPhotFactory


__all__ = ['write_lensed_sn_truth_summary', 'write_lensed_sn_variability_truth',
           'write_lensed_sn_light_curves']


logging.basicConfig(format="%(asctime)s %(name)s: %(message)s",
                    stream=sys.stdout)


def write_lensed_sn_truth_summary(lensed_sne_truth_cat, outfile, verbose=False):
    """
    Write the truth_summary table for the lensed SNe.

    Parameters
    ----------
    lensed_sne_truth_cat: str
        The sqlite3 file containing the model parameters for the lensed SNe.
    outfile: str
        Filename of the output sqlite3 file.
    verbose: bool [False]
        Verbosity flag.  No effect in this function.
    """
    logger = logging.getLogger('write_lensed_sn_truth_summary')
    if verbose:
        logger.setLevel(logging.INFO)

    logger.info('processing %s', lensed_sne_truth_cat)
    table_name = 'truth_summary'
    create_table_sql = f'''CREATE TABLE IF NOT EXISTS {table_name}
        (id TEXT, host_galaxy BIGINT, ra DOUBLE, dec DOUBLE,
        redshift FLOAT, is_variable INT, is_pointsource INT,
        flux_u FLOAT, flux_g FLOAT, flux_r FLOAT,
        flux_i FLOAT, flux_z FLOAT, flux_y FLOAT,
        flux_u_noMW FLOAT, flux_g_noMW FLOAT, flux_r_noMW FLOAT,
        flux_i_noMW FLOAT, flux_z_noMW FLOAT, flux_y_noMW FLOAT)'''
    with sqlite3.connect(lensed_sne_truth_cat) as conn, \
         sqlite3.connect(outfile) as output:
        # Create the output table if it does not exist.
        output.cursor().execute(create_table_sql)
        output.commit()
        # Query for the columns containing the model info for each SN
        # from the lensed_sne table.
        query = '''select unique_id, ra, dec, redshift from lensed_sne'''
        cursor = conn.execute(query)
        # Loop over each object and write its fluxes for all relevant
        # visits to the output table.
        values = []
        for unique_id, ra, dec, z in cursor:
            values.append((unique_id, -1, ra, dec, z, 1, 1,
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,))
        output.cursor().executemany(f'''insert into {table_name} values
                                    (?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                                     ?, ?, ?, ?, ?, ?, ?, ?, ?)''', values)
        output.commit()


def write_lensed_sn_light_curves(row, output, opsim_df, fp_radius, table_name,
                                 start_mjd=59580., end_mjd=61405):
    """
    Write the light curve values to the variability table.
    """
    unique_id, ra, dec, z, t_delay, magnification, t0, x0, x1, c = row
    sp_factory = SNSynthPhotFactory(z=z, t0=t0, x0=x0, x1=x1, c=c,
                                    snra=ra, sndec=dec)
    # Find the opsim db entries corresponding to the time span
    # when the SN is active and which are within fp_radius
    # degrees of the SN position.
    tmin, tmax = sp_factory.mintime() + t_delay, sp_factory.maxtime() + t_delay
    tmin = max(start_mjd, tmin)
    tmax = min(end_mjd, tmax)
    dmin, dmax = dec - fp_radius, dec + fp_radius
    df = pd.DataFrame(opsim_df.query(f'{tmin} <= expMJD <= {tmax} and '
                                     f'{dmin} <= dec <= {dmax}'))
    if len(df) == 0:
        return
    df['ang_sep'] = angularSeparation(df['ra'].to_numpy(), df['dec'].to_numpy(),
                                      ra, dec)
    df = df.query(f'ang_sep <= {fp_radius}')

    # Insert the rows into the variability truth table.
    values = []
    for visit, band, mjd in zip(df['obsHistID'], df['filter'], df['expMJD']):
        phot_params = PhotometricParameters(nexp=1, exptime=30, gain=1,
                                            bandpass=band)
        synth_phot = sp_factory.create(mjd - t_delay)
        flux = magnification*synth_phot.calcFlux(band)
        bp = synth_phot.bp_dict[band]
        num_photons = magnification*synth_phot.sed.calcADU(bp, phot_params)
        values.append((unique_id, visit, mjd, band, flux, num_photons))
    output.cursor().executemany(f'''insert into {table_name} values
                                    (?,?,?,?,?,?)''', values)
    output.commit()


def write_lensed_sn_variability_truth(opsim_db_file, lensed_sne_truth_cat,
                                      outfile, fp_radius=2.05, verbose=False,
                                      num_objects=None, start_mjd=59580.,
                                      end_mjd=61405.):
    """
    Write the lensed SNe fluxes to the lensed_sn_variabilty_truth table.

    Parameters
    ----------
    opsim_db_file: str
        OpSim db file.  This will be the minion 1016 db file that was
        modified by DESC for DC2.
    lensed_sne_truth_cat: str
        The sqlite3 file containing the model parameters for the lensed SNe.
    outfile: str
        Filename of the output sqlite3 file.
    fp_radius: float [2.05]
        Radius in degrees of the smallest acceptance cone containing the
        LSST focalplane projected onto the sky.
    verbose: bool [False]
        Verbosity flag.
    num_objects: int [None]
        Maximum number of objects from the input catalog to process.
        If None, then process all objects.
    """
    logger = logging.getLogger('write_lensed_sn_variability_truth')
    if verbose:
        logger.setLevel(logging.INFO)

    logger.info('processing %s', lensed_sne_truth_cat)
    table_name = 'lensed_sn_variability_truth'
    create_table_sql = f'''create table if not exists {table_name}
                           (id TEXT, obsHistID INT, MJD FLOAT, bandpass TEXT,
                            delta_flux FLOAT, num_photons FLOAT)'''

    # Read the opsim db data into a dataframe.
    with sqlite3.connect(opsim_db_file) as conn:
        opsim_df = pd.read_sql(
            '''select obsHistID, descDitheredRA, descDitheredDec,
               expMJD, filter from Summary''', conn)
    opsim_df['ra'] = np.degrees(opsim_df['descDitheredRA'])
    opsim_df['dec'] = np.degrees(opsim_df['descDitheredDec'])

    par_table_name = 'lensed_sne'
    # Loop over objects in the truth_cat containing the model
    # parameters and write the fluxes to the output table for the
    # relevant visits.
    with sqlite3.connect(lensed_sne_truth_cat) as conn, \
         sqlite3.connect(outfile) as output:
        # Create the output table if it does not exist.
        output.cursor().execute(create_table_sql)
        output.commit()
        # Get the number of SN entries in the input db file.
        cursor = conn.execute(f'select count(*) from {par_table_name}')
        nobjs = list(cursor)[0][0]
        if num_objects is None:
            num_objects = nobjs
        # Query for the columns containing the model parameters for
        # each SN from the lensed_sne table.
        query = f'''select unique_id, ra, dec, redshift, t_delay,
                    magnification, t0, x0, x1, c from {par_table_name}'''
        cursor = conn.execute(query)
        # Loop over each object and write its fluxes for all relevant
        # visits to the output table.
        for i, row in zip(range(num_objects), cursor):
            logger.info('%d  %d', i, num_objects)
            write_lensed_sn_light_curves(row, output, opsim_df, fp_radius,
                                         table_name, start_mjd=start_mjd,
                                         end_mjd=end_mjd)
