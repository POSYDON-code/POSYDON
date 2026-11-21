"""Explosion criteria of the core-collapse mechanisms.

For a star classified as a CCSN, or a PPISN or PISN candidate, the criterion
of the mechanism decides whether the explosion succeeds ("NS") or fails
("BH", direct collapse), or leaves the decision to the remnant mass (None).
"""

from abc import ABC, abstractmethod

from sklearn.neighbors import KNeighborsClassifier

from posydon.binary_evol.SN.engine_data import (
    NS_STATE,
    CO_core_at_He_depletion,
    Couch20_table,
    Patton20_CO_core_models,
    Sukhbold16_table,
)
from posydon.utils.posydonwarning import Pwarn

SN_MECHANISMS = [
    "Fryer+12-rapid",
    "Fryer+12-delayed",
    "direct",
    "Sukhbold+16-engine",
    "Patton&Sukhbold20-engine",
    "Couch+20-engine",
    "Maltsev+25-engine",
]

N_NEIGHBORS = 1

# Deprecated mechanisms, with the settings that reproduce them. They still
# work, but raise a DeprecationWarning and might be removed in the future.
DEPRECATED_MECHANISMS = {
    "direct_he_core": {"mechanism": "direct",
                       "conserve_hydrogen_envelope": False},
}


def replace_deprecated_mechanism(config):
    """Return config with a deprecated mechanism replaced by its equivalent.

    Configurations with any other mechanism are returned unchanged.
    """
    option = config['mechanism']
    if option not in DEPRECATED_MECHANISMS:
        return config
    replacement = DEPRECATED_MECHANISMS[option]
    settings = ", ".join(f"{key}={value!r}"
                         for key, value in replacement.items())
    conflicts = {key: config[key] for key, value in replacement.items()
                 if key != 'mechanism' and config[key] != value}
    if conflicts:
        raise ValueError(
            f"mechanism={option!r} is equivalent to {settings}, so it "
            f"cannot be combined with {conflicts}."
        )
    Pwarn(f"mechanism={option!r} is deprecated and will be removed in a "
          f"future release. Use {settings} instead.", "DeprecationWarning")
    new_config = config.copy()
    new_config.update(replacement)
    return new_config


# Mechanisms that are not available, with the settings used instead.
UNAVAILABLE_MECHANISMS = {
    "Maltsev+25-engine": {"mechanism": "Patton&Sukhbold20-engine",
                          "engine": "N20"},
}


def replace_unavailable_mechanism(config):
    """Return config with an unavailable mechanism replaced by its fallback.

    Configurations with any other mechanism are returned unchanged.
    """
    option = config['mechanism']
    if option not in UNAVAILABLE_MECHANISMS:
        return config
    replacement = UNAVAILABLE_MECHANISMS[option]
    Pwarn(f"mechanism={option!r} (engine={config['engine']!r}) is not "
          f"available in this version of POSYDON. Using "
          f"mechanism={replacement['mechanism']!r} with "
          f"engine={replacement['engine']!r} instead.",
          "UnsupportedModelWarning")
    new_config = config.copy()
    new_config.update(replacement)
    return new_config


# Possible outcomes of the explosion criteria:
#   "NS" : successful explosion, creates a neutron star
#          (Note that this does not mean NSs don't experience fallback)
#   "BH" : failed explosion, direct collapse to a black hole
#   None : no criterion; the remnant mass decides NS vs BH
class CCSN_check_base(ABC):
    """Base class of the explosion criteria.

    Calling it returns "NS", "BH" or None (see the outcomes above);
    subclasses implement _outcome.
    """
    def __init__(self,**kwargs):
        self.mechanism = kwargs['mechanism']
        self.verbose = kwargs['verbose']

    def __call__(self,star):
        """Return "NS", "BH" or None for a star classified as CCSN."""
        outcome = self._outcome(star)
        if self.verbose:
            print(f"{type(self).__name__}: explosion outcome "
                  f"{outcome or 'None (the remnant mass decides)'}")
        return outcome

    @abstractmethod
    def _outcome(self,star):
        """Return the explosion outcome of the star."""
        raise NotImplementedError

    def __repr__(self):
        return f"{type(self).__name__}(mechanism={self.mechanism!r})"


class CCSN_check_mass_based(CCSN_check_base):
    """
    NS vs BH is decided directly by the remnant mass prescription.

    Used by mechanisms where the remnant mass prescription sets the outcome on
    its own (Fryer+12, direct collapse). The remnant is a BH if its
    gravitational mass is at or above max_NS_mass, and a NS otherwise.
    """
    def _outcome(self,star):
        return None


class CCSN_check_from_tables(CCSN_check_base):
    """
    A general class for CCSN explosion outcomes that are taken from data
    tables with pre-computed core-collapse models. Currently applied to the
    mechanisms of Sukhbold et al. 2016 and Couch et al. 2020.

    Includes a kNN classifier (where n_neighbors = 1 is the default value to
    better capture the intrinsic nonlinearity of the problem), and predicts
    the outcome as "NS" for successful SNe, or "BH" for failed ones (direct
    collapse).
    """
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.engine = kwargs['engine']
        self.n_neighbors = N_NEIGHBORS
        table = self._table()
        self._classifier = KNeighborsClassifier(self.n_neighbors,
                                                weights="distance")
        self._classifier.fit(table[["He_c_mass"]].to_numpy(),
                             table["stellar_state"].to_numpy())
        if self.verbose:
            print(f"{type(self).__name__}: kNN classifier fitted on "
                  f"{len(table)} models (engine {self.engine!r}, "
                  f"n_neighbors={self.n_neighbors})")

    @abstractmethod
    def _table(self):
        """Return the models, with columns 'He_c_mass' and 'stellar_state'."""
        raise NotImplementedError

    def _outcome(self,star):
        state = self._classifier.predict([[star.he_core_mass]])[0]
        if self.verbose:
            print(f"{type(self).__name__}: He core "
                  f"{star.he_core_mass:.3f} Msun")
        return "NS" if state == NS_STATE else "BH"

    def __repr__(self):
        return (f"{type(self).__name__}(mechanism={self.mechanism!r}, "
                f"engine={self.engine!r})")


class CCSN_check_Sukhbold16(CCSN_check_from_tables):
    """
    Explosion outcome from the tables of Sukhbold et al. (2016), ApJ, 821, 38.

    The most massive model (He core 14.61 Msun) forms a BH for every engine, so
    every star with He core mass above that limit does too.
    """
    def _table(self):
        return Sukhbold16_table(self.engine)


class CCSN_check_Couch20(CCSN_check_from_tables):
    """
    Explosion outcome from the models of Couch, Warren & O'Connor (2020),
    ApJ, 890, 127, for the turbulence strength given by the option engine.

    For turbulence strengths of 1.2 and above the heaviest model (He core
    14.68 Msun) explodes, so every star with He core mass above that limit
    does too. For 1.3 and 1.4, all models form NSs.
    """
    def _table(self):
        return Couch20_table(self.engine)


# Ertl et al. (2016): (k1, k2) of each engine calibration
ERTL16_PARAMETERS = {
    "N20": (0.194, 0.0580),
    "S19.8": (0.274, 0.0470),
    "W15": (0.225, 0.0495),
    "W18": (0.283, 0.0430),
    "W20": (0.284, 0.0393),
    "Ertl2020": (0.182, 0.0608),
}


def explodes_Ertl16(M4, mu4, k1, k2):
    """Ertl et al. (2016): the explosion fails if k1·(M4·mu4) + k2 < mu4."""
    return not (k1 * (mu4 * M4) + k2 < mu4)


class CCSN_check_Patton20(CCSN_check_base):
    """
    Explosion outcome from the CO-core models of Patton & Sukhbold (2020).

    The star's CO core at He depletion (carbon mass fraction and CO core
    mass) is mapped onto Patton & Sukhbold's models, and an explodability
    criterion (given by the subclass) judges their pre-supernova structure.
    CO cores below the models (<= 2.5 Msun) explode; those above them
    (>= 10 Msun) collapse directly into a BH. M4 and mu4 are recorded on
    the star.
    """
    M_CO_min = 2.5
    M_CO_max = 10.0
    quantities = ("M4", "mu4")

    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.approx_at_he_depletion = kwargs['approx_at_he_depletion']
        self._Patton20 = Patton20_CO_core_models(self.quantities)
        if self.verbose:
            print(f"{type(self).__name__}: {self._Patton20!r} fitted on the "
                  "Patton & Sukhbold (2020) models")

    def _outcome(self,star):
        X_C, M_CO = CO_core_at_He_depletion(star, self.approx_at_he_depletion)
        values = self._Patton20(X_C, M_CO)
        star.M4, star.mu4 = values["M4"], values["mu4"]
        if self.verbose:
            print(f"{type(self).__name__}: CO core at He depletion "
                  f"X_C={X_C:.3f}, M_CO={M_CO:.3f} Msun: "
                  + ", ".join(f"{q}={v:.4f}" for q, v in values.items()))
        if M_CO <= self.M_CO_min:
            return "NS"
        if M_CO >= self.M_CO_max:
            return "BH"
        return "NS" if self._explodes(values) else "BH"

    @abstractmethod
    def _explodes(self,values):
        """True if the explosion of a model with these quantities succeeds."""
        raise NotImplementedError


class CCSN_check_Patton20_Ertl16(CCSN_check_Patton20):
    """
    Patton & Sukhbold (2020) CO-core models with the criterion of Ertl et al.
    (2016), using k1 and k2 of the engine calibration given by the option
    engine.
    """
    def __init__(self,**kwargs):
        self.engine = kwargs['engine']
        if self.engine not in ERTL16_PARAMETERS:
            raise ValueError(
                f"Invalid engine {self.engine!r} for the mechanism "
                f"'Patton&Sukhbold20-engine'. Valid engines are: "
                f"{', '.join(ERTL16_PARAMETERS)}."
            )
        self.k1, self.k2 = ERTL16_PARAMETERS[self.engine]
        super().__init__(**kwargs)

    def _explodes(self,values):
        return explodes_Ertl16(values["M4"], values["mu4"], self.k1, self.k2)

    def __repr__(self):
        return (f"{type(self).__name__}(mechanism={self.mechanism!r}, "
                f"engine={self.engine!r})")


CCSN_EXPLOSION_CRITERIA = {
    "Fryer+12-rapid": CCSN_check_mass_based,
    "Fryer+12-delayed": CCSN_check_mass_based,
    "direct": CCSN_check_mass_based,
    "Sukhbold+16-engine": CCSN_check_Sukhbold16,
    "Couch+20-engine": CCSN_check_Couch20,
    "Patton&Sukhbold20-engine": CCSN_check_Patton20_Ertl16,
}


def get_CCSN_prescription(**kwargs):
    """Return the explosion criterion of the mechanism kwargs['mechanism']."""
    option = kwargs['mechanism']
    if option in CCSN_EXPLOSION_CRITERIA:
        return CCSN_EXPLOSION_CRITERIA[option](**kwargs)
    if option in SN_MECHANISMS:
        raise NotImplementedError(
            f"Mechanism {option!r} is not available in this version of "
            "POSYDON."
        )
    raise ValueError(
        f"Invalid option {option!r} given to mechanism (i.e. the CCSN "
        f"prescription). Valid options are: {', '.join(SN_MECHANISMS)}."
    )
