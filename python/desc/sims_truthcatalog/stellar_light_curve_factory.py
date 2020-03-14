import os
import sys
import sqlite3
import numpy as np
import pandas as pd
from lsst.sims.utils import defaultSpecMap
import lsst.sims.photUtils as sims_photUtils
import lsst.sims.catUtils.mixins.VariabilityMixin as variability


__all__ = ['StellarLightCurveFactory']


class VariabilityGenerator(variability.StellarVariabilityModels,
                           variability.MLTflaringMixin,
                           variability.ParametrizedLightCurveMixin):
    """
    This is a class that mimics the API of an InstanceCatalog. This makes it
    easier for us to generate light curves from sources with differing
    variability models.
    """

    def __init__(self, chunk):
        """
        chunk = (simobjid, hpid, sedFilename, magNorm, ebv, varParamStr,
                 parallax, ra, dec)
        """

        # initialize some member variables expected by the
        # stellar variability models
        self.photParams = sims_photUtils.PhotometricParameters(nexp=1,
                                                               exptime=30)
        self.lsstBandpassDict \
            = sims_photUtils.BandpassDict.loadTotalBandpassesFromFiles()
        self._actually_calculated_columns = []
        for bp in 'ugrizy':
            self._actually_calculated_columns.append('lsst_%s' % bp)

        sed_dir = os.environ['SIMS_SED_LIBRARY_DIR']
        self.quiescent_mags = {}
        for bp in 'ugrizy':
            self.quiescent_mags[bp] = np.zeros(len(chunk), dtype=float)

        self.parallax = np.zeros(len(chunk), dtype=float)
        self.ebv = np.zeros(len(chunk), dtype=float)
        self.varParamStr = np.empty(len(chunk), dtype=(str,300))
        self.simobjid = np.empty(len(chunk), dtype=int)
        ccm_wav = None

        # find the quiescent magnitudes of the stars
        for i_star, star in enumerate(chunk):
            self.simobjid[i_star] = int(star[0])
            full_file_name = os.path.join(sed_dir,
                                          defaultSpecMap[star[2]])
            spec = sims_photUtils.Sed()
            spec.readSED_flambda(full_file_name)
            fnorm = sims_photUtils.getImsimFluxNorm(spec, float(star[3]))
            spec.multiplyFluxNorm(fnorm)
            if ccm_wav is None or not np.array_equal(spec.wavelen, ccm_wav):
                ccm_wav = np.copy(spec.wavelen)
                a_x, b_x = spec.setupCCM_ab()
            ebv_val = float(star[4])
            spec.addDust(a_x, b_x, ebv=ebv_val, R_v=3.1)
            mag_list = self.lsstBandpassDict.magListForSed(spec)
            for i_bp, bp in enumerate('ugrizy'):
                self.quiescent_mags[bp][i_star] = mag_list[i_bp]
            self.ebv[i_star] = ebv_val
            self.varParamStr[i_star] = star[5]
            self.parallax[i_star] = float(star[6])

        self.parallax *= (np.pi/648000000.0) # convert from mas to radians

    def column_by_name(self, name):
        if name.startswith('quiescent'):
            bp = name[-1]
            return np.copy(self.quiescent_mags[bp])
        elif name == 'ebv':
            return np.copy(self.ebv)
        elif name=='parallax':
            return np.copy(self.parallax)
        elif name=='simobjid':
            return np.copy(self.simobjid)
        raise RuntimeError("\n\nCannot get column %s\n\n" % name)


class StellarLightCurveFactory:
    def __init__(self, stars_db_file=None, mjds=None):
        if stars_db_file is None:
            stars_db_file = ('/global/projecta/projectdirs/lsst/'
                             'groups/SSim/DC2/dc2_stellar_healpixel.db')
        self.conn = sqlite3.connect(stars_db_file)

        # Make sure that the Kepler light curves are loaded into memory
        variability.ParametrizedLightCurveMixin()\
                   .load_parametrized_light_curves()
        self.t0 = 59580.
        self.mjds = np.arange(0, 365*5, 5) + self.t0 if mjds is None else mjds

    def __del__(self):
        self.conn.close()

    def create(self, simobjid):
        bands = 'ugrizy'
        query = f"""SELECT simobjid, -1, sedFilename, magNorm, ebv,
                    varParamStr, parallax, ra, decl FROM stars where
                    simobjid={simobjid} limit 1"""
        chunk = list(self.conn.execute(query))
        foo = VariabilityGenerator(chunk)
        dmags = foo.applyVariability(foo.varParamStr, expmjd=self.mjds)
        # Transpose the columns so that the ordering is star, filter, mjd
        dmags = dmags.transpose([1, 0, 2])
        return {band: _ for band, _ in zip(bands, dmags[0])}
