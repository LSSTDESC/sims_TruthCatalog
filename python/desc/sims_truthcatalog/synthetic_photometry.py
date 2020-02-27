"""
Interface to lsst.sims.photUtils code to perform synthetic photometry,
i.e., calculate object fluxes in LSST standard bandpasses.
"""
import os
import copy
import numpy as np
from lsst.sims.catUtils.dust import EBVbase
import lsst.sims.photUtils as sims_photUtils
from lsst.sims.utils import defaultSpecMap


__all__ = ['SyntheticPhotometry', 'find_sed_file']


def find_sed_file(sed_file):
    """
    Return the full path to the SED file assuming it is in the
    lsst_sims SED library.
    """
    try:
        spec_file = defaultSpecMap[sed_file]
    except KeyError:
        spec_file = sed_file
    full_path = os.path.join(os.environ['SIMS_SED_LIBRARY_DIR'], spec_file)
    if not os.path.isfile(full_path):
        raise FileNotFoundError(full_path)
    return full_path


class CCMmodel:
    """
    Helper class to cache a(x) and b(x) arrays evaluated on wavelength
    grids for intrinsic and Galactic extinction calculations.
    """
    def __init__(self):
        self.wavelen = dict()
        self.a = dict()
        self.b = dict()
        for ext_type in ('intrinsic', 'Galactic'):
            self.wavelen[ext_type] = None

    def add_dust(self, sed_obj, Av, Rv, ext_type):
        """
        Add dust reddening to the SED object.

        Parameters
        ----------
        sed_obj: lsst.sims.photUtils.Sed
            SED object to which to add the reddening.
        Av: float
            Extinction coefficient.
        Rv: float
            Extinction coefficient.
        ext_type: str
            Extinction type: 'intrinsic' or 'Galactic'
        """
        if (self.wavelen[ext_type] is None or
            not np.array_equal(sed_obj.wavelen, self.wavelen[ext_type])):
            self.a[ext_type], self.b[ext_type] = sed_obj.setupCCM_ab()
            self.wavelen[ext_type] = copy.deepcopy(sed_obj.wavelen)
        sed_obj.addDust(self.a[ext_type], self.b[ext_type],
                        A_v=Av, R_v=Rv)


class SyntheticPhotometry:
    """
    Class to provide interface to lsst.sims.photUtils code.
    """
    # The dust_models dict is a shared class-level resource among
    # SyntheticPhotometry instances to take advantage of the caching of
    # the a(x) and b(x) arrays used by the dust models.
    dust_models = dict(ccm=CCMmodel())
    lsst_bp_dict = sims_photUtils.BandpassDict.loadTotalBandpassesFromFiles()
    ebv_model = EBVbase()
    def __init__(self, sed_file, mag_norm, redshift=0, iAv=0, iRv=3.1,
                 gAv=0, gRv=3.1, dust_model_name='ccm', bp_dict=None):
        """
        Parameters
        ----------
        sed_file: str
            File containing the unnormalized SED as columns of wavelength
            in nm and flux-density as flambda.  If None, then don't create
            an Sed object.
        mag_norm: float
            Monochromatic magnitude of the object at 500nm.  This provides
            the normalization of the SED.
        redshift: float [0]
            The redshift of the object.
        iAv: float [0]
            Reference extinction parameter for internal reddening.
            iAv = 0 corresponds to no extinction.
        iRv: float [3.1]
            Extinction ratio.  iRv = 3.1 corresponds to a nominal Milky Way
            extinction law assuming the CCM model.
        gAv: float [0]
            Reference extinction parameter for Milky Way reddening.
        gRv: float [3.1]
            Galactic extinction ratio.  gRv = 3.1 corresponds to a nominal
            Milky Way extinction law assuming the CCM model.
        dust_model_name: str ['ccm']
            Name of the dust model to use in the shared dust_model dict.
            'ccm' corresponds to the model from Cardelli, Clayton, & Mathis
            1989 ApJ, 345, 245C.
        bp_dict: dict [None]
            Dictionary of bandpass objects.  If None, then the LSST total
            throughputs will be used.
        """
        self.sed_file = sed_file
        self.mag_norm = mag_norm
        self.redshift = redshift
        self.iAv = iAv
        self.iRv = iRv
        self.gAv = gAv
        self.gRv = gRv
        self.dust_model_name = dust_model_name
        if bp_dict is None:
            self.bp_dict = self.lsst_bp_dict
        if sed_file is not None:
            self._create_sed()

    @staticmethod
    def create_from_sed(sed, redshift):
        synth_phot = SyntheticPhotometry(None, None)
        synth_phot.sed = sed
        if redshift != 0:
            synth_phot.sed.redshiftSED(redshift, dimming=True)
            synth_phot.sed.resampleSED(wavelen_match
                                       =synth_phot.bp_dict.wavelenMatch)
        return synth_phot

    @staticmethod
    def get_MW_AvRv(ra, dec, Rv=3.1):
        eq_coord = np.array([[np.radians(ra)], [np.radians(dec)]])
        ebv = SyntheticPhotometry.ebv_model\
                                 .calculateEbv(equatorialCoordinates=eq_coord,
                                               interp=True)
        Av = Rv*ebv
        return Av, Rv

    def add_MW_dust(self, ra, dec, Rv=3.1):
        Av, _ = self.get_MW_AvRv(ra, dec, Rv=Rv)
        self.add_dust(Av, Rv, 'Galactic')
        return Av, Rv

    def add_dust(self, Av, Rv, component):
        """
        Apply dust to the SED for the specified component.
        """
        if component not in ('intrinsic', 'Galactic'):
            raise RuntimeError(f'Unrecognized dust component: {component}')
        self.dust_models[self.dust_model_name].add_dust(self.sed, Av, Rv,
                                                        component)

    def _create_sed(self):
        """
        Function to create the lsst.sims.photUtils Sed object.
        """
        self.sed = sims_photUtils.Sed()
        self.sed.readSED_flambda(self.sed_file)
        fnorm = sims_photUtils.getImsimFluxNorm(self.sed, self.mag_norm)
        self.sed.multiplyFluxNorm(fnorm)

        if self.iAv != 0:
            self.add_dust(self.sed, self.iAv, self.iRv, 'intrinsic')
        if self.redshift != 0:
            self.sed.redshiftSED(self.redshift, dimming=True)
        self.sed.resampleSED(wavelen_match=self.bp_dict.wavelenMatch)
        if self.gAv != 0:
            self.add_dust(self.gAv, self.gRv, 'Galactic')

    def calcFlux(self, band):
        """
        Calculate the flux in the desired band.

        Parameters
        ----------
        band: str
            `ugrizy` band.

        Returns
        -------
        Flux in the band in nanojanskys.
        """
        return self.sed.calcFlux(self.bp_dict[band])*1e9
