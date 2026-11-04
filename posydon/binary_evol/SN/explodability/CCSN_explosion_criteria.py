from abc import ABC, abstractmethod

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
# work, but raise a DeprecationWarning and will be removed in the future.
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

# Possible outcomes of an explosion criterion:
#   "NS"          : successful explosion leaving a neutron star
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
    """No explosion criterion: NS vs BH is decided by the remnant mass.

    Used by mechanisms whose remnant mass prescription sets the outcome on
    its own (Fryer+12, direct collapse): the remnant is a BH if its mass
    exceeds max_NS_mass, and a NS otherwise.
    """
    def __call__(self,star):
        return None


CCSN_EXPLOSION_CRITERIA = {
    "Fryer+12-rapid": CCSN_check_mass_based,
    "Fryer+12-delayed": CCSN_check_mass_based,
    "direct": CCSN_check_mass_based,
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
