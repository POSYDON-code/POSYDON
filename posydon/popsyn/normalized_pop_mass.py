"""Compute the underlying stellar population mass for a given simulation."""


__authors__ = [
    "Devina Misra <devina.misra@unige.ch>",
    "Simone Bavera <Simone.Bavera@unige.ch>",
    "Konstantinos Kovlakas <Konstantinos.Kovlakas@unige.ch>",
]


from posydon.popsyn import independent_sample
from scipy.integrate import quad


def initial_total_underlying_mass(df=None, **kwargs):
    """Compute the initial total mass of the population.

    Parameters
    ----------
    df : DataFrame
        Data frame of a population class. If nothing is provided the code will
        sample initial conditions.
    primary_mass_min : type
        Description of parameter `primary_mass_min`.
    primary_mass_max : type
        Description of parameter `primary_mass_max`.
    primary_mass_scheme : type
        Description of parameter `primary_mass_scheme`.
    secondary_mass_scheme : type
        Description of parameter `secondary_mass_scheme`.
    **kwargs : type
        Description of parameter `**kwargs`.

    Returns
    -------
    type
        Description of returned object.

    """
    # ENHANCEMENT: the code should assume default values of the POSYDON sampler
    # if not provided in kwargs.

    if df is None:
        initial_ZAMS_mass = independent_sample.generate_independent_samples(
            **kwargs)
        initial_ZAMS_mass_1 = initial_ZAMS_mass[2]
        initial_ZAMS_mass_2 = initial_ZAMS_mass[3]
        initial_ZAMS_TOTAL_mass = (sum(initial_ZAMS_mass_1)
                                   + sum(initial_ZAMS_mass_2))
    elif isinstance(df, float):
        initial_ZAMS_TOTAL_mass = df
    else:
        sel = df['event'] == 'ZAMS'
        initial_ZAMS_TOTAL_mass = sum(df['S1_mass'][sel]+df['S2_mass'][sel])
    #Eirini's changes to the IMF to include the 
    def imf_part_1(m, m_min, alpha1):
        return (m/m_min)**-alpha1

    def imf_part_2(m, m_1, m_min, alpha1, alpha2):
        return ((m_1/m_min)**-alpha1) * ((m/m_1)**-alpha2)

    def imf_part_3(m, m_1, m_2, m_min, alpha1, alpha2, alpha3):
        return ((m_1/m_min)**-alpha1)*((m_2/m_1)**-alpha2)*((m/m_2)**-alpha3)

    def mean_mass_of_binary(f0, f_bin, m_1, m_2, m_min, m_max,
                            alpha1, alpha2, alpha3):
        a1, a2, a3 = alpha1, alpha2, alpha3
        mean_mass = f0 * (1 + (f_bin/2)) * (
            (1/(2-a1)) * ((m_1**(2-a1) - m_min**(2-a1))/(m_min**-a1))
            + ((1/(2-a2))
               * (m_1/m_min)**-a1 * ((m_2**(2-a2) - m_1**(2-a2))/(m_1**-a2)))
            + ((1/(2-a3)) * (m_1/m_min)**-a1 * (m_2/m_1)**-a2
               * ((m_max**(2-a3) - m_2**(2-a3))/(m_2**-a3))))
        return mean_mass
    #Eirini's suggestion is to built a universal mean_mass_of_binary that can be used both for model mass and the pop mass.
    #Also if posydon goes to smaller masses this way will be flexible to that too.
    def mean_mass_binary(f0, f_bin, m_1, m_2, m_min, m_max,
                            alpha1, alpha2, alpha3):
        if m_min <= 0.08:
            #Include the all the imf parts
            a1, a2, a3 = 1, 1, 1
        elif 0.08 < m_min <= 0.5:
            m1 = m_min
            a1,a2,a3 = 0, 1, 1
            alpha1 = 0
        else: 
            a1,a2,a3 = 0, 0, 1
            alpha1,alpha2 = 0,0 
            m_2 = m_min 
        mean_mass = f0 * (1 + (f_bin/2))* (1 + (f_bin/2))*(quad(lambda x: a1*imf_part_1(x, m_min, a1)*x , m_min, m_1)[0]
        + quad(lambda x: a2*imf_part_2(x,m_1,m_min,alpha1,alpha2)*x, m_1, m_2)[0]
        + quad(lambda x: a3*imf_part_3(x, m_1, m_2, m_min, alpha1, alpha2, alpha3)*x, m_2, m_max)[0])
        
        return mean_mass   
    #Normalization parameter.
    def f_0(a1,a2,a3,alpha1,alpha2,alpha3):
        f0 = 1/(quad(a1*imf_part_1, m_min, m_1, args=(m_min, alpha1))[0]
            + quad(a2*imf_part_2, m_1, m_2, args=(m_1, m_min, alpha1, alpha2))[0]
            + quad(a3*imf_part_3, m_2, m_max, args=(m_1, m_2, m_min, alpha1,alpha2, alpha3))[0])
        return f0
    
    def mean_mass_simulated(alpha3, m_a, m_b):
        a = alpha3
        mean_mass_sim = (3/2) * ((1-a) / (2-a)) * (
            (m_b**(2-a)-m_a**(2-a))/(m_b**(1-a)-m_a**(1-a)))
        return mean_mass_sim

    def fraction_simulated(f_bin, m_1, m_2, m_min, m_max,
                           alpha1, alpha2, alpha3):
        a1, a2, a3 = alpha1, alpha2, alpha3
        f_model = f_bin * f0 * (1/(1-a3)) * (m_1/m_min)**(-a1) * (m_2/m_1)**(
            -a2) * ((m_b**(1-a3)-m_a**(1-a3)) / (m_2**-a3))
        return f_model

    # Kroupa P., 2001, MNRAS, 322, 231
    if (kwargs['primary_mass_scheme'] == 'Kroupa2001'
            and kwargs['secondary_mass_scheme'] == 'flat_mass_ratio'):
        m_1 = 0.08
        m_2 = 0.5
        alpha1 = 0.3
        alpha2 = 1.3
        alpha3 = 2.3

    # Salpeter E. E., 1955, ApJ, 121, 161
    elif (kwargs['primary_mass_scheme'] == 'Salpeter'
            and kwargs['secondary_mass_scheme'] == 'flat_mass_ratio'):
        m_1 = 0.08
        m_2 = 0.5
        alpha1 = 2.35
        alpha2 = 2.35
        alpha3 = 2.35
    else:
        raise ValueError("Scheme not included yet")

    f_bin_pop = kwargs['binary_fraction_const']
    f_bin = 0.7
    m_min = 0.01
    m_max = 200.0
    m_a = kwargs['primary_mass_min']
    m_b = kwargs['primary_mass_max']

    f0 = 1/(quad(imf_part_1, m_min, m_1, args=(m_min, alpha1))[0]
            + quad(imf_part_2, m_1, m_2, args=(m_1, m_min, alpha1, alpha2))[0]
            + quad(imf_part_3, m_2, m_max, args=(m_1, m_2, m_min, alpha1,
                                                 alpha2, alpha3))[0])

    f_corr = (fraction_simulated(f_bin, m_1, m_2, m_min,
                                 m_max, alpha1, alpha2, alpha3)
              * mean_mass_simulated(alpha3, m_a, m_b)
              / mean_mass_of_binary(f0, f_bin, m_1, m_2, m_min, m_max,
                                    alpha1, alpha2, alpha3))

    f_model = fraction_simulated(f_bin, m_1, m_2, m_min, m_max,
                                 alpha1, alpha2, alpha3)

    return initial_ZAMS_TOTAL_mass / f_corr, f_corr, f_model
