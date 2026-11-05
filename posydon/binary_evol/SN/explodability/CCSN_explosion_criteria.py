from abc import ABC, abstractmethod

from sklearn.neighbors import KNeighborsClassifier

from posydon.binary_evol.SN.engine_data import SUKHBOLD16_NS, Sukhbold16_table
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

# Possible outcomes of the explosion criteria:
#   "NS"          : successful explosion, creates a neutron star
#                   (Note that this does not mean NSs don't experience fallback)
#   "BH_fallback" : successful explosion, black hole formed by fallback
#   "BH_direct"   : failed explosion, direct collapse to a black hole
#   None          : no criterion; the remnant mass decides NS vs BH
EXPLOSION_OUTCOMES = ("NS", "BH_fallback", "BH_direct", None)


class CCSN_check_base(ABC):
    def __init__(self,**kwargs):
        self.mechanism = kwargs['mechanism']
        self.verbose = kwargs['verbose']

    @abstractmethod
    def __call__(self,star):
        """Return one of EXPLOSION_OUTCOMES for a star classified as CCSN."""
        raise NotImplementedError

    def __repr__(self):
        return f"{type(self).__name__}(mechanism={self.mechanism!r})"


class CCSN_check_mass_based(CCSN_check_base):
    """
    NS vs BH is decided directly by the remnant mass prescription.

    Used by mechanisms where the remnant mass prescription sets the outcome on
    its own (Fryer+12, direct collapse). The remnant is a BH if its mass
    exceeds max_NS_mass, and a NS otherwise.
    """
    def __call__(self,star):
        return None


class CCSN_check_Sukhbold16(CCSN_check_base):
    """
    Explosion outcome from the tables of Sukhbold et al. (2016), ApJ, 821, 38.

    A k-nearest-neighbours classifier, trained on the He core masses of the
    models of the chosen engine, predicts "NS" for a successful explosion or
    "BH_direct" for a failed one. With n_neighbors = 1 the star gets the
    outcome of the model with the closest He core mass. The most massive
    model (He core 14.61 Msun) forms a BH for every engine, so every star
    above the tables does too.
    """
    n_neighbors = 1

    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.engine = kwargs['engine']
        table = Sukhbold16_table(self.engine)
        self._classifier = KNeighborsClassifier(self.n_neighbors,
                                                weights="distance")
        self._classifier.fit(table[["He_c_mass"]].to_numpy(),
                             table["stellar_state"].to_numpy())

    def __call__(self,star):
        state = self._classifier.predict([[star.he_core_mass]])[0]
        return "NS" if state == SUKHBOLD16_NS else "BH_direct"

    def __repr__(self):
        return (f"{type(self).__name__}(mechanism={self.mechanism!r}, "
                f"engine={self.engine!r})")


CCSN_EXPLOSION_CRITERIA = {
    "Fryer+12-rapid": CCSN_check_mass_based,
    "Fryer+12-delayed": CCSN_check_mass_based,
    "direct": CCSN_check_mass_based,
    "Sukhbold+16-engine": CCSN_check_Sukhbold16,
}


def get_CCSN_prescription(**kwargs):
    option = kwargs['mechanism']
    if option in CCSN_EXPLOSION_CRITERIA:
        return CCSN_EXPLOSION_CRITERIA[option](**kwargs)
    if option in SN_MECHANISMS:
        raise NotImplementedError(
            f"The explosion criterion for mechanism {option!r} has not been "
            "ported to the new step_SN structure yet."
        )
    raise ValueError(
        f"Invalid option {option!r} given to mechanism (i.e. the CCSN "
        f"prescription). Valid options are: {', '.join(SN_MECHANISMS)}."
    )
