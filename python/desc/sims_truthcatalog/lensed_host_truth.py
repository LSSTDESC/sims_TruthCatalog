"""
Module to compute truth catalog values for lensed hosts in DC2 Run3.0i.
"""
import os
import glob
from collections import defaultdict
import sqlite3
from astropy.io import fits
import numpy as np
import pandas as pd
from .synthetic_photometry import find_sed_file, SyntheticPhotometry
from .write_sqlite import write_sqlite


__all__ = ['write_lensed_host_truth']


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
                           host_types=('agn', 'sne')):
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

    Returns
    -------
    (dict, dict, dict): dicts of fluxes with MW extinction, w/out MW
        extinction, and a dict of (ra, dec, redshift) tuples, all keyed
        by object id.
    """
    band_fluxes = lambda: {band:0 for band in bands}
    fluxes = defaultdict(band_fluxes)
    fluxes_noMW = defaultdict(band_fluxes)
    coords = dict()
    mag_norms = dict()
    with sqlite3.connect(host_truth_db_file) as conn:
        for host_type in host_types:
            df = pd.read_sql(f'select * from {host_type}_hosts', conn)
            for component in components:
                mag_norms[component] = get_mag_norms(host_type, component,
                                                     image_dir)
            for iloc in range(len(df)):
                row = df.iloc[iloc]
                ra = row['ra_lens']
                dec = row['dec_lens']
                redshift = row['redshift']
                unique_id = str(row['unique_id'])
                coords[unique_id] = (ra, dec, redshift)
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
    return dict(fluxes), dict(fluxes_noMW), coords


def write_lensed_host_truth(host_truth_db_file, image_dir, outfile,
                            bands='ugrizy'):
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
    """
    fluxes, fluxes_noMW, coords = get_lensed_host_fluxes(host_truth_db_file,
                                                         image_dir)
    if os.path.isfile(outfile):
        raise OSError(f'{outfile} already exists.')
    ids = sorted(list(set(fluxes.keys()).intersection(fluxes_noMW.keys())))
    ras, decs, redshifts = [], [], []
    flux_by_band_MW = {_: [] for _ in bands}
    flux_by_band_noMW = {_: [] for _ in bands}
    for unique_id in ids:
        ra, dec, redshift = coords[unique_id]
        ra.append(ra)
        dec.append(dec)
        redshifts.append(redshift)
        for band in bands:
            flux_by_band_MW[band].append(fluxes[unique_id][band])
            flux_by_band_noMW[band].append(fluxes_noMW[unique_id][band])
    galaxy_ids = -1*np.ones(len(ids))
    is_variable = np.zeros(len(ids))
    is_pointsource = np.zeros(len(ids))
    good_ixes = range(len(ids))
    write_sqlite(outfile, ids, galaxy_ids, ras, decs, redshifts,
                 is_variable, is_pointsource, flux_by_band_MW,
                 flux_by_band_noMW, good_ixes)
