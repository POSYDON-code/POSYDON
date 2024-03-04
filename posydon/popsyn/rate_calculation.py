
__author__ = ['Simone Bavera <Simone.Bavera@unige.ch>',
              'Max Briel <max.briel@unige.ch>',]

from posydon.utils.constants import Zsun

DEFAULT_MODEL = {
    'delta_t' : 100, # Myr
    'SFR' : 'IllustrisTNG',
    'sigma_SFR' : None,
    'Z_max' : 5.,
    'select_one_met' : False,
    'dlogZ' : None, # e.g, [np.log10(0.0142/2),np.log10(0.0142*2)]
    'Zsun' : Zsun,
    'compute_GRB_properties' : False,
    'GRB_beaming' : 1., # e.g., 0.5, 'Goldstein+15'
    'GRB_efficiency' : 0., # e.g., 0.01
    'E_GRB_iso_min' : 0., # e.g., 1e51 erg 
}