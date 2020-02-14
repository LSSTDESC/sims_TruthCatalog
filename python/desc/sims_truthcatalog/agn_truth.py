"""
Module to write truth catalogs for AGNs using the AGNs parameters db as input.
"""
import os
import math
from collections import defaultdict
import sqlite3
import numpy as np
import pandas as pd
from lsst.sims.utils import angularSeparation
from .synthetic_photometry import SyntheticPhotometry
from .write_sqlite import write_sqlite


__all__ = ['AGNTruthWriter', 'agn_light_curve']


def agn_light_curve(mjds, redshift, tau, sf, seed, start_date=58580.,
                    dt_frac=1e-2):
    """
    Return the delta fluxes wrt the infinite-time average flux for the
    provided AGN light curve parameters.

    Parameters
    ----------
    mjds: list-like
        Times at which to evaluate the light curve delta flux values in MJD.
    redshift: float
        Redshift of the AGN, used to account for time-dilation between
        rest frame and observer frame.
    tau: float
        Variability time scale in days.
    sf: float
        Structure function parameter, i.e., asymptotic rms variability on
        long time scales.
    seed: int
        Random number seed.
    start_date: float [58580.]
        Start date for light curve generation in MJD.
    dt_frac: float [1e-2]
        Fraction of tau to use as the time step for generating the
        light curve.

    Returns
    -------
    np.array of delta magnitude values.

    Notes
    -----
    This implements the model described in Macleod et al. 2010, ApJ, 721, 1014.
    https://arxiv.org/abs/1004.0276
    """
    duration_rest_frame = (max(mjds) - start_date)/(1. + redshift)
    dt = tau*dt_frac
    nbins = int(math.ceil(duration_rest_frame/dt)) + 1
    time_indexes = np.round((mjds - start_date)/(1. + redshift)/dt).astype(int)
    time_index_map = defaultdict(list)
    ct_index = 0
    for i, t_index in enumerate(time_indexes):
        time_index_map[t_index].append(i)

    rng = np.random.RandomState(seed)
    es = rng.normal(0, 1., nbins)*np.sqrt(dt_frac)
    dx2 = 0
    x1 = 0
    x2 = 0
    d_m_out = np.zeros(len(mjds))
    for it in range(nbins):
        dx1 = dx2
        dx2 = -dx1*dt_frac + sf*es[it] + dx1
        x1 = x2
        x2 += dt
        if it in time_index_map:
            for it_out in time_index_map[it]:
                local_end = (mjds[it_out] - start_date)/(1. + redshift)
                dm_val = (local_end*(dx1 - dx2) + dx2*x1 - dx1*x2)/(x1 - x2)
                d_m_out[it_out] = dm_val
    return d_m_out

class AGNTruthWriter:
    '''
    Write Summary and Variable truth tables for unlensed AGNs.
    '''
    def __init__(self, outfile, agn_db_file):
        '''
        Parameters
        ----------
        outfile: str
            Name of the sqlite3 file to contain the truth tables.
        agn_db_file: str
            The sqlite3 file containing the AGN model parameters.
        '''
        self.outfile = outfile
        if os.path.isfile(outfile):
            raise OSError(f'{outfile} already exists.')
        if not os.path.isfile(agn_db_file):
            raise FileNotFoundError(f'{agn_db_file} not found.')
        with sqlite3.connect(agn_db_file) as conn:
            self.agn_db = pd.read_sql('select * from agn_params', conn)

    def write(self):
        pass
