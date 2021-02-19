"""
Module to compute truth catalog values for lensed hosts in DC2 Run3.0i.
"""
import os
import sys
import glob
from collections import defaultdict
import logging
import sqlite3
from astropy.io import fits
import pandas as pd
from lsst.sims.photUtils import PhotometricParameters
from .synthetic_photometry import find_sed_file, SyntheticPhotometry
from .sqlite_utils import write_column_descriptions


__all__ = ['write_lensed_host_truth']


logging.basicConfig(format="%(asctime)s %(name)s: %(message)s",
                    stream=sys.stdout)


def get_mag_norms(host_type, component, image_dir, bands='ugrizy'):
    """
    Get magnorm values from FITS stamp file headers.

    Parameters
    ----------
    host_type: str
        'agn' or 'sne'
    component: str
        'bulge' or 'disk'
    image_dir: str
        Directory containing the FITS stamps.
    bands: str or list-like ['ugrizy']
        Bands for which to return magnorms.

    Returns
    -------
    dict, keyed by 'LENS_ID', of dict of magnorm values keyed by band.
    """
    pattern = os.path.join(image_dir, f'{host_type}_lensed_{component}s',
                           '*.fits')
    files = glob.glob(pattern)
    mag_norms = dict()
    for item in files:
        with fits.open(item) as hdus:
            header = hdus[0].header
            mag_norms[header['LENS_ID']] \
                = {band: header[f'MAGNORM{band.upper()}'] for band in bands}
    return mag_norms


def get_lensed_host_fluxes(host_truth_db_file, image_dir, bands='ugrizy',
                           components=('bulge', 'disk'),
                           host_types=('agn', 'sne'), verbose=False):
    """
    Loop over entries in `agn_hosts` and `sne_hosts` tables in
    the host_truth_db_file and compute fluxes (with and without MW
    extinction) in each band.  Return dicts of fluxes and coordinates
    keyed by object ids.

    Parameters
    ----------
    host_truth_db_file: str
        File containing model parameters for lensed host of AGNs and SNe.
    image_dir: str
        Directory containing the FITS stamps.
    bands: str or list-like ['ugrizy']
        Bands for which to return magnorms.
    components: list-like [('bulge', 'disk')]
        Galaxy components of lensed hosts.
    host_types: list-like [('agn', 'sne')]
        Types of hosted objects.
    verbose: bool [False]
        Verbose flag.

    Returns
    -------
    (dict, dict, dict): dicts of fluxes with MW extinction, w/out MW
        extinction, and a dict of (ra, dec, redshift) tuples, all keyed
        by object id.
    """
    logger = logging.getLogger('get_lensed_host_fluxes')
    if verbose:
        logger.setLevel(logging.INFO)

    logger.info('processing %s', host_truth_db_file)
    band_fluxes = lambda: {band:0 for band in bands}
    fluxes = defaultdict(band_fluxes)
    fluxes_noMW = defaultdict(band_fluxes)
    num_photons = defaultdict(band_fluxes)
    coords = dict()
    mag_norms = dict()
    with sqlite3.connect(host_truth_db_file) as conn:
        for host_type in host_types:
            df = pd.read_sql(f'select * from {host_type}_hosts', conn)
            for component in components:
                mag_norms[component] = get_mag_norms(host_type, component,
                                                     image_dir)
            for iloc in range(len(df)):
                logger.info('%s  %d  %d', host_type, iloc, len(df))
                row = df.iloc[iloc]
                ra = row['ra_lens']
                dec = row['dec_lens']
                redshift = row['redshift']
                unique_id = str(row['unique_id'])
                coords[unique_id] = [ra, dec, redshift]
                gAv = row['av_mw']
                gRv = row['rv_mw']
                for component in components:
                    if unique_id not in mag_norms[component]:
                        continue
                    sed_file = find_sed_file(
                        row[f'sed_{component}_host'].lstrip('b').strip("'"))
                    iAv = row[f'av_internal_{component}']
                    iRv = row[f'rv_internal_{component}']
                    for band in bands:
                        mag_norm = mag_norms[component][unique_id][band]
                        synth_phot = SyntheticPhotometry(sed_file, mag_norm,
                                                         redshift=redshift,
                                                         iAv=iAv, iRv=iRv)
                        fluxes_noMW[unique_id][band] \
                            += synth_phot.calcFlux(band)
                        synth_phot.add_dust(gAv, gRv, 'Galactic')
                        fluxes[unique_id][band] += synth_phot.calcFlux(band)
                        bp = synth_phot.bp_dict[band]
                        photpars = PhotometricParameters(nexp=1, exptime=30,
                                                         gain=1, bandpass=band)
                        num_photons[unique_id][band] \
                            += synth_phot.sed.calcADU(bp, photpars)
    return dict(fluxes), dict(fluxes_noMW), dict(num_photons), coords


def write_lensed_host_truth(host_truth_db_file, image_dir, outfile,
                            bands='ugrizy', verbose=False):
    """
    Write the truth_summary fluxes for the lensed hosts.

    Parameters
    ----------
    host_truth_db_file: str
        File containing model parameters for lensed host of AGNs and SNe.
    image_dir: str
        Directory containing the FITS stamps.
    outfile: str
        Filename of output sqlite3 file.
    bands: str or list-like ['ugrizy']
        Bands for which to return magnorms.
    verbose: bool [False]
        Verbose flag.
    """
    cmd = '''CREATE TABLE IF NOT EXISTS truth_summary
        (id TEXT, host_galaxy BIGINT, ra DOUBLE, dec DOUBLE,
        redshift FLOAT, is_variable INT, is_pointsource INT,
        flux_u FLOAT, flux_g FLOAT, flux_r FLOAT,
        flux_i FLOAT, flux_z FLOAT, flux_y FLOAT,
        flux_u_noMW FLOAT, flux_g_noMW FLOAT, flux_r_noMW FLOAT,
        flux_i_noMW FLOAT, flux_z_noMW FLOAT, flux_y_noMW FLOAT,
        num_photons_u FLOAT, num_photons_g FLOAT, num_photons_r FLOAT,
        num_photons_i FLOAT, num_photons_z FLOAT, num_photons_y FLOAT)'''
    fluxes, fluxes_noMW, num_photons, coords \
        = get_lensed_host_fluxes(host_truth_db_file, image_dir, verbose=verbose)
    host_galaxy = -1
    is_variable = 0
    is_pointsource = 0
    ids = sorted(list(set(fluxes.keys()).intersection(fluxes_noMW.keys())))
    values = []
    for unique_id in ids:
        values.append([unique_id, host_galaxy] + coords[unique_id]
                      + [is_variable, is_pointsource]
                      + [fluxes[unique_id][_] for _ in bands]
                      + [fluxes_noMW[unique_id][_] for _ in bands]
                      + [num_photons[unique_id][_] for _ in bands])
    if os.path.isfile(outfile):
        raise OSError(f'{outfile} already exists.')
    with sqlite3.connect(outfile) as conn:
        write_column_descriptions(conn)
        cursor = conn.cursor()
        cursor.execute(cmd)
        conn.commit()
        cursor.executemany('INSERT INTO truth_summary '
                           'VALUES (?,?,?,?,?,?,?,'
                                   '?,?,?,?,?,?,'
                                   '?,?,?,?,?,?,'
                                   '?,?,?,?,?,?)', values)
        conn.commit()
        # index to speed up location searches
        cursor.execute('create index radec_ix on truth_summary(ra, dec)')
        conn.commit()
