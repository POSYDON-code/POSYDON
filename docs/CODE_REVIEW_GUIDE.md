# POSYDON scientific code review guide

## 1. Purpose, scope, and reference snapshot

Use this guide to review changes that could affect the correctness,
reproducibility, or interpretation of POSYDON results. It is a set of review
instructions, not a claim that the current implementation has the example
defects described below and not a certification of scientific correctness.

The guide was prepared on **20 November 2026**, using representative source,
configuration, test, and workflow files from
[`POSYDON-code/POSYDON`, commit `11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1`](https://github.com/POSYDON-code/POSYDON/tree/11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1)
(the `main` head dated 19 November 2026), together with the Fragos et al.
(2023) and Andrews et al. (2025) instrument papers. Selected primary
references are listed in Section 14. This was a source/documentation inspection;
the preparation of this guide did not run POSYDON, retrain interpolators, or
reproduce the papers' numerical results.

For each future PR, inspect its actual base/head revisions, configuration,
data/model versions, and more-specific instructions. Paths below are a map of
the inspected snapshot; follow moved implementations rather than assuming old
paths or defaults are permanent. The instrument papers establish published
methodology and assumptions; later deliberate scientific improvements are
allowed and need their own rationale and validation.

### Review priorities

Prioritize incorrect scientific outcomes, biased population estimates,
unintended model changes, numerical failures, invalid state transitions,
corrupted data, and compatibility failures. Include ordinary software defects
when they affect these outcomes. Avoid style commentary already handled by
formatters and CI, speculative physical objections without a relevant regime,
and requests for extensive experiments unrelated to the diff.

## 2. Review procedure

1. **Establish the change.** Read the diff and stated objective. Identify the
   physical or statistical result that could change, the affected configurations,
   and whether the PR intends a scientific change or an implementation-only
   change. Inspect existing relevant tests and CI results.
2. **Trace the calculation.** Follow input construction, units, transforms,
   branches, callers, event handling, saved fields, and downstream consumers.
   Check the affected component and population channel, not only a local line.
3. **Write down the scientific contract.** Identify the prescription, variable
   definitions, input domain, assumptions, valid missing states, and intended
   output. Distinguish model fidelity from numerical accuracy and physical
   uncertainty. Use Section 4 for evidence-dependent literature checks.
4. **Try to falsify a specific concern.** Check dimensions, an appropriate
   limit, a boundary case, or a minimal counterexample. Compare with the base
   revision when attributing a regression. Use independent reasoning rather
   than duplicating the same expression in a test.
5. **Assess downstream impact.** Determine whether the behavior affects an
   accessible configuration, which systems are affected, and how their fates,
   properties, weights, or rates change. Trace existing validation or sanitizing
   code that might already handle the case.
6. **Report demonstrated defects and evidence gaps separately.** Give a concise
   actionable finding only when the defect and consequence are sufficiently
   supported. Identify targeted additional validation for unresolved material
   questions. If nothing is established, do not manufacture a finding.

A missing test is normally a validation gap, not proof of an incorrect result.
A changed scientific baseline may be intentional; a matching baseline may
preserve an old error. Evaluate both against the intended science.

## 3. Module map and cross-cutting conventions

### Where to look

| Changed area | Main paths in the inspected snapshot | Relevant sections |
| --- | --- | --- |
| State flow and object histories | `posydon/binary_evol/flow_chart.py`, `binarystar.py`, `singlestar.py`, `simulationproperties.py` | 5, 11 |
| Detailed-grid evolution and track matching | `posydon/binary_evol/MESA/step_mesa.py`, `posydon/binary_evol/track_match.py` | 5, 8, 9 |
| Detached evolution, tides, winds, braking, radiation | `posydon/binary_evol/DT/step_detached.py`; `DT/tides/`, `DT/winds/`, `DT/magnetic_braking/`, `DT/gravitational_radiation/`; `DT/double_CO.py` | 5, 6, 11 |
| Common envelope and merger products | `posydon/binary_evol/CE/step_CEE.py`, `posydon/binary_evol/DT/step_merged.py` | 5, 7 |
| Supernovae, compact remnants and spins | `posydon/binary_evol/SN/step_SN.py`, `SN/profile_collapse.py`, `posydon/grids/SN_MODELS.py` | 5, 7, 9 |
| Shared physical calculations and thresholds | `posydon/utils/common_functions.py`, `constants.py`, `limits_thresholds.py` | 5–7, 11 |
| Interpolation and classification | `posydon/interpolation/IF_interpolation.py`, `interpolation.py`, `data_scaling.py`, `constraints.py`, `eep.py`, `profile_interpolation.py` | 8, 9, 11 |
| Grid assembly and post-processing | `posydon/grids/psygrid.py`, `io.py`, `lazy_hdf.py`, `scrubbing.py`, `downsampling.py`, `post_processing.py`, `termination_flags.py` | 8, 9 |
| Population generation and sampling | `posydon/popsyn/binarypopulation.py`, `independent_sample.py`, `sample_from_file.py`, `distributions.py`, `IMFs.py`, `Moes_distributions.py`, `norm_pop.py` | 10, 11 |
| Cosmic populations, rates and selections | `posydon/popsyn/synthetic_population.py`, `rate_calculation.py`, `star_formation_history.py`, `selection_effects.py`, `transient_select_funcs.py`, `GRB.py` | 10 |
| MESA inputs, defaults and pipelines | `grid_params/`, including the `POSYDON-MESA-INLISTS` submodule; associated grid/pipeline code | 5, 9, 12 |
| Tests and scientific comparisons | `posydon/unit_tests/`, `.github/workflows/`, `dev-tools/script_data/src/` | 12 |

### Conventions to establish before judging a formula

| Quantity | Required check |
| --- | --- |
| Mass and radius | Distinguish solar units from cgs, total from core/envelope mass, and baryonic from gravitational remnant mass. |
| Time | Distinguish absolute age, elapsed step time, delay time and cosmic time; verify years, Myr, Gyr and seconds at each interface. |
| Orbit | Determine whether separation denotes semimajor axis, instantaneous distance, or periastron; identify which mass ratio is used. |
| Spin | Distinguish angular velocity, angular momentum, specific angular momentum, dimensionless spin magnitude and signed projection. Check radians/second versus radians/year. |
| Metallicity | Distinguish absolute mass fraction `Z` from `Z/Zsun`, and the solar convention of the data or prescription. |
| Logarithms | Identify base 10 versus natural logarithm, normalization before taking the logarithm, inverse transforms, and zero/negative-value conventions. |
| Rates and fractions | Distinguish signed stellar mass loss, positive transfer magnitude, accreted mass and escaping mass; determine which quantity a fraction divides. |
| Probability | State whether a PDF is per linear variable, natural log, or log10 variable, and whether it is conditional or truncated. |
| Missing data | Establish the meaning of `None`, `NaN`, infinities and sentinels for the specific field and state. Undefined/not-applicable values are not physical zero. |

At the reference snapshot, `constants.py` supplies cgs constants and solar
conversions. `inspiral_timescale_from_separation` returns **Myr**, whereas binary
evolution time is generally expressed in **years**. Detached integration uses
angular velocities in radians/year internally and converts stored stellar
quantities at interfaces. These are concrete interfaces to inspect, not
permission to infer a unit solely from a variable name or stale comment.

## 4. Scientific and literature verification

### When a lookup is warranted

Consult primary sources when the PR introduces or changes a physical formula,
fitting coefficient, threshold, calibration, approximation or validity domain;
claims to fix an existing scientific error; or exposes uncertainty about the
meaning of a prescription. Also check sources when papers, documentation,
comments and executable code disagree in a way that matters to the result.

A mechanical rename or equivalent refactor does not require an unrelated
literature survey. Read the few sources needed for the changed scientific
claim. The objective is to verify that claim, not to recommend every newer
model in the field.

### Procedure for a defensible comparison

1. Find the citation already attached to the affected implementation and
   establish the selected runtime prescription. Consult the relevant POSYDON
   paper section from Section 14 for context.
2. Open the original prescription paper, its applicable version and any
   relevant correction when the available tools permit. Search by author,
   title, DOI or arXiv identifier. Use abstracts and search snippets to locate
   sources, not to verify an equation they do not show.
3. Record the exact quantity and assumptions in the source: definitions of
   mass, metallicity, age, angular variables, normalization, coordinate system,
   fitted range, evolutionary epoch and approximations. Map code variables to
   those quantities before comparing coefficients or factors.
4. Check the algebra and numerical realization independently. For piecewise
   fits, inspect every affected branch and boundary. If HTML/OCR loses a sign,
   subscript or exponent, consult the PDF or source rather than guessing.
   Resolve discrepancies between the original prescription and a later
   adaptation explicitly.
5. Cite the source actually read, with a verified equation, section or page
   when possible. Explain the mismatch and its consequence in the applicable
   regime. Do not attach an equation number from memory or cite an abstract as
   evidence of a coefficient.

Use these evidence distinctions in the review record:

- **Verified:** the relevant text/equation was inspected and compared with the
  implementation under matching assumptions.
- **Partially verified:** the source identity or general prescription was
  established, but the decisive equation, data table or version was unavailable.
- **Unverified:** an external comparison remains necessary; identify exactly
  what is missing and what evidence would resolve it.

An instructions file does not grant internet access, ADS access or journal
subscriptions. If external access is unavailable, use repository documentation,
available reference material and independent analytic checks, then state the
remaining limitation where the output interface permits. An inaccessible
source alone is not a defect, and an unresolved suspicion is not a confirmed
scientific error. Avoid uploading code or unpublished data to external services
to obtain a reference; literature queries normally need only public bibliographic
information and the physical topic.

### Model choice versus implementation error

Keep four cases distinct: a formula implemented incorrectly; a previously
selected model changed unintentionally; a justified change to the model; and an
acknowledged approximation or scientific uncertainty. CE efficiencies, natal
kick prescriptions, explosion engines, stellar winds and mixing assumptions
are often choices with substantial uncertainty. Require a clear rationale and
appropriate validation for changes, without declaring a choice incorrect merely
because another publication adopts a different one. Conversely, published text
does not excuse a dimensional inconsistency or a demonstrable transcription
error.

For new prescriptions, prefer durable documentation containing the source,
variable mapping, units, domain, assumptions and one independent validation
case. This preserves useful evidence for future reviews even without browsing.

## 5. Stellar and binary evolution contracts

### Component identity and state transitions

Trace initial component identity separately from the current primary,
secondary, donor and accretor. Changes involving `flip_stars`, reverse mass
transfer or symmetric grid use must update the appropriate component, histories,
event labels, spin/kick quantities and component-specific binary fields. A
double swap should restore the applicable state; evaluate that property only
for fields included in the operation's declared contract.

Check that changes to a stellar-state threshold or termination flag select the
intended next step. Verify both component orderings and states at boundaries:
initial RLOF, contact, detachment, collapse, CE ejection, merger, disruption,
massless remnant and maximum evolution time. Distinguish physical termination
from a numerical failure. A failed track must not silently become a valid
merger or zero-mass survivor.

Histories must remain aligned across the binary and both stars. Allow repeated
physical timestamps for instantaneous events while preserving event order and
avoiding repeated application of the same event. A track-relative clock must
not reset the binary's accumulated age. Check near-simultaneous events on both
stars and the maximum-time boundary.

### Physical consistency and conservation

Use mass, angular-momentum and energy accounting appropriate to the model.
Include systemic winds, ejecta, neutrinos, radiation, accretion losses and
spin/orbit exchange where present. Do not demand conservation of the bound
binary's mass or orbital angular momentum in an explicitly open system.

For consistently nested core definitions, check
`0 <= M_CO_core <= M_He_core <= M_star`, allowing physical equality and justified numerical
tolerance. Do not impose this ordering on unrelated core definitions. Check
core radii and convective-region geometry against their definitions. Individual
composition fractions lie in `[0, 1]`; normalization to one applies to a complete
composition, not an arbitrary subset of stored species.

Check radiative closure using emergent luminosity and photospheric quantities,
`L = 4*pi*R^2*sigma_SB*T_eff^4`, with logs and solar normalization decoded.
Nuclear, neutrino and surface luminosities need not be identical during
evolution. Interpolation constraints should preserve the appropriate relation
without silently imposing a physically inapplicable one.

For mass transfer, establish the signs and meanings of each rate and efficiency.
Trace retained mass, expelled mass, and the specific angular momentum carried
away. Check zero-transfer and conservative limits where applicable, donor
switches, Eddington-limited accretion, and updates to stellar composition and
rotation. Do not assume a universal bound on an effective CE efficiency or
monotonic loss of mass by the initially more massive star.

### MESA-dependent physics

A change in winds, overshooting, mixing, thermohaline treatment, rotation,
accretion energy, equation of state or stopping conditions may change the
underlying model grid. Establish the MESA revision **and custom patches**,
inlist layering, submodule revision, rerun options, grid/data release and
post-processing configuration. A revision label alone does not identify all
physics. Distinguish an intended alternative grid from accidentally mixed
assumptions. Do not interpret a population-level parameter change as modifying
physics already fixed in stored MESA tracks.

## 6. Orbital dynamics, detached evolution and analytic checks

Use the following checks only within their stated assumptions. In the equations,
`m1` and `m2` are component masses, `M = m1 + m2`, `a` is semimajor axis, `e` is
eccentricity, `G` is the gravitational constant, and `c` is the speed of light.
All quantities in a calculation must use a consistent unit system.

### Kepler and Roche geometry

For a Newtonian two-body orbit,

$$
P^2 = \frac{4\pi^2 a^3}{G(m_1+m_2)}.
$$

Check period/separation round trips and mass-exchange symmetry. Choose a
tolerance consistent with the repository's constants and conversion precision;
do not demand machine precision when independent conversion constants differ.

For the Eggleton approximation [R3], define `q = M_donor/M_accretor` for the
donor's lobe:

$$
\frac{R_{\rm L}}{a} =
\frac{0.49q^{2/3}}{0.6q^{2/3}+\ln(1+q^{1/3})}.
$$

Check the correct component, logarithm base and input separation convention.
The standard circular, synchronous interpretation must not silently become a
general eccentric/asynchronous prescription. If the code uses a periastron
approximation, evaluate that declared approximation explicitly. Swapping the
donor changes `q` to its reciprocal; it does not leave that donor's lobe unchanged.

### Gravitational radiation

For isolated, constant-mass, initially circular point masses evolving under
leading-order quadrupole radiation [R4], define

$$
\beta = \frac{64}{5}\frac{G^3m_1m_2(m_1+m_2)}{c^5},\qquad
t_{\rm merge} = \frac{a_0^4}{4\beta}
=\frac{5c^5a_0^4}{256G^3m_1m_2(m_1+m_2)}.
$$

Check mass symmetry, `a0^4` scaling and the circular limit of eccentric
calculations. If stopping at a finite separation `a_f`, the circular elapsed
time is `(a0^4 - a_f^4)/(4*beta)`, not the formal time to zero separation.
Compare the ODE solver and standalone delay-time helper after converting their
time units and reconciling stop conditions.

In `DT/double_CO.py` at the reference snapshot, the Peters equations are
integrated using `s = -ln(a/a0)`, `l = ln(e)` and dimensionless time. Its scale
`a0^4/beta` is **four times** the circular time to zero separation. A variable
called a characteristic timescale need not equal the merger time. For changes
to these transformations, derive the chain rule, integration direction, initial
conditions and conversion back to physical time; check dense output and event
status remapping as well as the RHS.

For valid bound weak-field orbits with GW radiation alone, leading-order `a`
decreases, `de/dt` is zero at `e=0` and negative for `0<e<1`. Check small-e
stability, high-e behavior, overflow/underflow, and the distinction between
coupled eccentric evolution and an approximation at fixed eccentricity.

The separate detached radiation prescription in
`DT/gravitational_radiation/default_gravrad.py` cites Junker & Schäfer [R5] and
includes higher-order terms. Compare it with Peters in the appropriate
weak-field limit, not by requiring exact equality at every finite separation.
Match the source's orbital-element convention and approximation order.

### Tides, winds and mass transfer

For a circular orbit with aligned, synchronized spins, fixed stellar structure
and other torques switched off, tidal derivatives should vanish under the corresponding
equilibrium-tide prescription [R6]. Outside that regime, tides can transfer
spin angular momentum into the orbit; universal shrinkage or eccentricity
damping is not an appropriate test. Inspect the chosen convective/radiative
regime, moment of inertia, spin units, and stellar evolution terms.

For slow isotropic wind loss with no accretion in the relevant adiabatic
two-body approximation, `a*M` remains constant and eccentricity is unchanged.
Do not apply that result to impulsive mass loss or a different angular-momentum
loss prescription.

For conservative transfer between circular point masses, with fixed total mass
and orbital angular momentum and no spin exchange,

$$
\frac{a_f}{a_i} =
\left(\frac{m_{1,i}m_{2,i}}{m_{1,f}m_{2,f}}\right)^2.
$$

Use this as a conditional limit of an appropriate transfer calculation, not
as a constraint on general nonconservative grid evolution.

## 7. Common envelope, supernovae and compact remnants

### Common envelope

Establish the selected CE prescription, donor identity, core definition,
envelope mass, lambda calculation, internal/recombination energy terms and
post-ejection treatment. For a single envelope and an unchanged companion,
one conventional positive-binding-energy form [R7] is

$$
E_{\rm bind} = \alpha_{\rm CE}
\left(\frac{G M_{\rm core}M_{\rm comp}}{2a_f}
-\frac{G M_{\rm donor}M_{\rm comp}}{2a_i}\right).
$$

Here `M_comp` is the unchanged companion mass, `E_bind` is the positive energy
required for ejection, and `alpha_CE` is the selected effective efficiency.
POSYDON's implementation may store binding
energy with the opposite sign; translate the convention before judging the
algebra. A double CE requires the appropriate envelope budgets and final
component masses. Keep a core definition and its associated lambda/mass/radius
consistent throughout the calculation.

Check both post-CE Roche lobes, the ejection/merger decision, and any subsequent
stable transfer or wind-loss phase. Do not equate the immediate energy-formalism
separation with the final separation after additional phases. Check behavior
for negligible envelopes, unavailable profiles, unbound envelopes and selected
fallbacks. An effective `alpha_CE > 1` or a different core boundary is not alone
a programming error. Missing lambda must follow a documented path rather than
accidentally erasing the envelope budget.

### Collapse and pair instability

Trace the selected engine through configuration, preprocessing, serialized
grid quantities and runtime collapse. Do not mix a He-core-based threshold with
a CO-core-based fit, or a core property at He depletion with one at collapse.
Check piecewise boundaries just below, at, and above each affected threshold;
intentional discontinuities and nonmonotonic explodability are possible.

Distinguish pre-SN stellar mass, baryonic remnant mass, gravitational mass,
fallback, ejecta, neutrino losses, and disk-accreted/radiated mass. Verify that
losses are counted once and the physical budget closes under the chosen
prescription. Check WD, ECSN, iron-core collapse, PPISN, PISN and direct-collapse
paths that the change can reach. Complete disruption must propagate consistently
to a massless-remnant state and binary evolution.

Compare interpolated collapse outputs, profile-based collapse and core-mass
fallback only when they use compatible assumptions; differences among these
paths are not automatically defects. Changing an engine or its parameters
requires matching stored post-processing products or a correctly selected
runtime calculation. Check SN model names, parameter metadata and absent-model
handling.

For a BH, `chi = c*J/(G*M^2)` is dimensionless, with physical magnitude bounded
by one in the applicable Kerr description. The same proxy for an extended
progenitor is not subject to that bound. Check enclosed mass, profile cell
orientation and widths, angular-momentum integration, disk losses, and the
difference between spin magnitude, tilt and signed projected spin.

### Natal kicks and orbital response

Verify the exploding component, pre-/post-event masses, orbital phase,
right-handed frame, vector rotations and velocity units. Compare energy and
angular momentum with the resulting orbit. Different conventions can produce
equivalent formulas; map them explicitly before reporting a sign error.

Useful conditional tests [R8, R9] include:

- **Zero kick and zero instantaneous mass loss:** orbital elements remain
  unchanged. Zero kick alone is insufficient when mass is lost.
- **Initially circular Newtonian orbit, instantaneous isotropic loss, no kick
  or accretion:** let `f = DeltaM/M_tot,initial`, where the denominator is the
  total pre-SN binary mass. For `0 <= f < 1/2`,
  `a_f/a_i = (1-f)/(1-2*f)` and `e_f = f/(1-f)`; at `f=1/2` the relative orbit
  is parabolic, and larger loss unbinds it. This is not a general eccentric-orbit
  formula.
- **Isotropic kick directions:** `cos(theta)` is uniform on `[-1,1]` and
  azimuth is uniform on `[0,2*pi)` in the relevant axis convention. A uniform
  polar angle is not isotropic.
- **Specified kick distribution:** distinguish a Maxwellian's one-dimensional
  dispersion from its speed moments, and check fallback/mass normalization
  exactly once. Test deterministic supplied kicks separately from sampling.

Treat unbound solutions as disruptions rather than clipping eccentricity into
a bound interval. Check systemic velocity, spin-orbit tilt and post-SN
histories, including the second collapse. Statistical sampling tests require
appropriate uncertainty; a small random sample need not match an analytic
distribution exactly.

### Accretion luminosities and observables

When changing X-ray, GRB or related selections, identify the accretion-rate
definition, radiative efficiency, compact-object type, Eddington normalization,
bandpass, beaming convention and luminosity units. Trace physical luminosity
versus isotropic-equivalent or selected observable luminosity. Check piecewise
sub-/super-Eddington branches and zero-accretion behavior without imposing a
hard luminosity cap on a prescription that explicitly permits a different
observable behavior. Consult the specific source used by that prescription.

## 8. Interpolation, classification and validation data

### Coordinates, classes and validity

Trace feature order, `M2` versus `q`, mass-ratio definition, period/separation,
metallicity, log transforms, scales and inverse transforms from training to
inference. Check class-specific preprocessing and the chosen metallicity model.
The published v2 scheme uses separate interpolators/classifiers by metallicity;
a deliberate new scheme requires its own validation rather than automatic
rejection [P2, Section III].

Preserve the distinction between valid numeric tracks, classifier-only
outcomes, physical terminal outcomes, failed tracks, missing histories and
unavailable fields. At the reference snapshot, `initial_MT` can be classified
without providing a valid numerical interpolation, and there are explicit NaN
and sentinel paths. Do not convert these into plausible physical values.

Check exact-grid points, hull boundaries, holes, rare classes, equal-mass cases
and out-of-domain behavior. `LinInterpolator` currently has a nearest-neighbor
fallback outside the interpolation hull: evaluate whether the documented
fallback and warning remain correct rather than assuming every extrapolation
attempt is a defect. Do not average through distinct outcomes or require
smoothness across physical bifurcations. Changing labels or a validity mask can
alter an entire population channel.

### Physical constraints and profiles

Inspect both raw predictions and their sanitization. Correct constraint order,
field selection, component identity and units matter for Kepler closure,
Stefan–Boltzmann closure, nuclear-luminosity sums, composition, mass/radius
relations and convective-region geometry. Luminosities represented as logarithms
must be summed in linear space. Clipping or renormalization must not conceal a
systematic error or silently bias rare outcomes.

For profile interpolation and EEP operations, check phase alignment,
mass/radius coordinates, orientation, endpoint/core boundaries, missing profiles
and integral quantities. A smooth interpolated curve is insufficient evidence
that its mass, composition and angular-momentum integrals are consistent.

### Training, tuning and independent evaluation

Determine what every grid actually does; its name alone does not establish
whether it is a training, validation or test sample. Track reuse through
scaling, imputation, feature selection, hyperparameter tuning, classifier
selection and interpolation fitting. Correlated samples from the same physical
track must not create a misleading independence claim.

Refitting on training plus validation data can be appropriate after model
selection if evaluation remains independent. A production refit including a
former test grid is also possible, but the earlier test score is not independent
validation of that refitted model. Increasing training data does not guarantee
monotonic improvement of every estimator or physical channel.

For affected algorithms, examine class-dependent errors, confusion matrices,
rare-class performance, metallicity dependence and errors near boundaries.
Check the denominator of relative errors near zero. Compare physical outputs
and downstream fates as well as aggregate statistical metrics. Require only
the scale of evaluation justified by the PR.

Check model serialization, feature metadata, training provenance, dependency
compatibility and regeneration requirements when the training or inference
contract changes. Existing pickles do not automatically adopt new preprocessing.

## 9. Grid processing, data provenance and storage

Preserve run identity through scrubbing, downsampling, joining, reruns,
post-processing and model training. Matching row counts is insufficient: verify
MESA-directory/run identifiers and ordering. Resolve duplicates and replacement
priority explicitly. A rerun under changed mixing or energy treatment is a
potential model change, even when intended to improve convergence.

For histories, distinguish MESA backsteps/retries from physical evolution.
Scrubbing must keep the accepted sequence and align the binary and both
stellar histories by their appropriate model identifiers. Sorting a single
array by age can misalign all other fields. If age is shifted at RLO onset,
apply the same convention to all affected histories and preserve total age
where needed. For profiles, use their spatial/mass coordinate rather than
assuming the independent coordinate is time.

Downsampling must preserve endpoints and the physical episodes relevant to
consumers: transfer onset/termination, short-lived transitions, depletion,
collapse, extrema and threshold crossings. Check the meaning of relative
error near zero, normalized coordinates, and whether independent downsampling
can break synchronization. A lower point count is not sufficient validation.

Check history/profile/final-value consistency, single-star versus binary flags,
and compact-object handling in post-processing. State thresholds must agree
with the intended upstream MESA definitions or deliberately document a changed
classification. Missing profiles and failed runs require explicit treatment.

At the reference snapshot, SN outputs can be stored separately under
`/grid/SN_MODELS/{SN_MODEL_NAME}` with parameter metadata; consumers include the
interpolation layer. Schema changes must account for both these datasets and
the associated readers, writers and model products. Check field types and
string widths, old-file handling, lazy-file lifetime, and HDF row/ID alignment.

Record code revision, MESA revision/patches and inlist gitlink, grid/data release,
metallicity convention, post-processing settings, trained model identity and
relevant numerical dependencies. Check that a schema or defaults change does
not silently mix old and new products. Do not regenerate scientific reference
data merely to make a failing comparison disappear.

`grid_params/POSYDON-MESA-INLISTS` is a separate repository represented by a
submodule. A parent-repository review can inspect a changed gitlink, but this
guide in POSYDON does not automatically review independent PRs in the inlists
repository. Those reviews need their own repository configuration and applicable
scientific instructions.

## 10. Sampling, population normalization and cosmological rates

### Probability measures and represented mass

Establish the sampling/proposal distribution, target distribution, supported
domain and probability measure for every variable. Inspect IMF normalization,
conditional mass-ratio bounds, primary/secondary mass limits, periods or
separations, eccentricities, single-star/binary fractions and companion masses.
Check truncated distributions in the intended measure.

A change of variable uses the corresponding Jacobian. For fixed total mass,
Kepler's law gives `dlog(a)/dlog(P) = 2/3`. In the reference implementation,
several `.pdf()` methods accept a linear argument but return density per
`log10` variable, and the period-conversion code uses this logarithmic factor.
Do not replace it with a linear-variable Jacobian without changing the complete
contract and its consumers.

Trace model weights, simulated mass and mean-mass corrections together.
Distinguish a system count, event count, stellar mass and weighted contribution.
Check normalization when filtering failures, splitting/merging batches or
changing binary fractions. Reweighting cannot reconstruct target support that
the simulation never sampled; replacing infinite weights with zero can conceal
that problem.

Verify that transient rows map to the correct initial systems and that one
system can contribute multiple legitimate events. Preserve IDs rather than
relying on incidental DataFrame ordering. Doubling a sample and its represented
mass should leave the appropriately normalized expectation unchanged; this
test does not apply to every raw count or unnormalized sum.

### Cosmic time, metallicity and observed rates

Track birth time, delay time and event time, with the same cosmology and units.
Use the correct sign when converting lookback time to cosmic age. Exclude or
handle future events consistently. Check interpolation/extrapolation in
redshift, time and metallicity, bin edges, and finite bin widths.

Distinguish `Z` from `Z/Zsun`, linear- from log-metallicity distributions, and
lognormal parameters from arithmetic means. Inspect the configured truncation
and renormalization policy; weights covering a restricted metallicity interval
do not always sum to one. Do not silently change the population being represented.

Distinguish source-frame event-rate density from observer-frame event counts.
Trace star-formation history, mass normalization, metallicity weights, time
bins, comoving volume, redshift time dilation and detection probability through
the **whole expression**. Apply conversions between Mpc/Gpc and yr/Myr/Gyr
exactly once.

In the inspected code, `get_shell_comoving_volume` includes a `1/(1+z)` factor,
despite its abbreviated name. The rate calculation's conventions must be
examined before alleging a missing or duplicate time-dilation factor. Do not
infer physical dimensions solely from a helper's name.

For selection functions, check probability bounds where the output really is
a probability, source versus redshifted masses, luminosity distance versus
comoving distance, orientations, detector assumptions, calibrated domain and
intrinsic versus observed distributions. A weight or rate need not lie in
`[0,1]` just because a detection probability does.

Useful checks include zero rate for zero event/SFH weight, correct limiting
selection probability, preservation of total weights across compatible batch
operations, and consistent normalized results when the same sample is
replicated with its normalization. Use a small analytically controlled example
before assuming a large population comparison is necessary.

## 11. Numerical reliability and reproducibility

For numerical changes, inspect cancellation, overflow/underflow, near-zero
denominators, nonfinite values, log-domain limits, dtype/broadcasting behavior
and loss of precision at long absolute times. Tolerances need units and scales;
an arbitrary global epsilon can change a physical classification. Test both
sides of reachable thresholds and explain any intentionally discontinuous
behavior.

Check integrator success, terminal events, event direction and ordering, dense
output, interpolation support, maximum time, and transformations between
independent variables. A solver returning an array is not evidence it reached
the intended event. Silent clipping, `nan_to_num`, ignored warnings or broad
exception handling can turn a failed physical calculation into an apparently
valid population member; distinguish documented safeguards from suppression.

Use an independent reference, dimensional derivation or convergence check when
the changed numerical method could bias a physical result. A regression test
that simply reproduces the changed formula does not validate that formula.
Do not demand exact floating-point equality for trajectories or stochastic
results whose contract allows numerical/statistical tolerance.

Check random-number generator ownership, seed propagation and independent
streams across processes, jobs and metallicities. Require repeatability for
the documented execution configuration. Do not assume bitwise identity across
changed worker counts unless the implementation explicitly guarantees it.
For base/head comparisons, use comparable seeds and explain when altered draw
order prevents a simple paired comparison.

Inspect shared mutable defaults, cache keys, model/data identifiers, lazy
resources and output ordering where they can couple independent calculations.
Numerical dependency changes can affect persisted interpolation models and
reproducibility; evaluate compatibility rather than casually requesting upgrades.

## 12. Verification with the actual POSYDON infrastructure

### Unit and focused checks

Use the current `pyproject.toml`, dependency constraints and
`.github/workflows/continuous_integration.yml`. At the reference snapshot,
CI uses Python 3.11, installs the test extra, and runs `pytest` with:

```bash
export PATH_TO_POSYDON=./
export PATH_TO_POSYDON_DATA=./posydon/unit_tests/_data/
export MESA_DIR=./
pytest
```

These commands are for a prepared checkout with the required fixtures and
dependencies. `config.py` appends `POSYDON_data` to the supplied data root;
do not append that directory twice. The unit-CI `MESA_DIR` placeholder does
not constitute an installed MESA environment.

For a focused diagnostic, one example is:

```bash
pytest -o addopts='' posydon/unit_tests/utils/test_common_functions.py
```

This intentionally omits repository-wide default pytest options for that
diagnostic; report it as a targeted test, not as a passing full CI gate. Use
the actual affected test path. Default coverage settings cover selected
modules, so their percentage is not proof of complete scientific coverage.
There is no dedicated `posydon/unit_tests/interpolation/` directory in the
reference tree; the utility-interpolator tests are a different target.

Do not mark a calculation verified when its tests were skipped because grids,
profiles, serialized models or compiled dependencies were unavailable. State
the exact unverified behavior and the smallest practical next check.

### Scientific baseline comparison

POSYDON already contains
`.github/workflows/baseline-validation-on-comment.yml`. At the reference
snapshot, an authorized user's `/run-baseline-check` comment can trigger it
for a same-repository PR. It uses a self-hosted runner, available grids, pinned
base/head revisions and defined binary cases. The guide is not authorization
for a reviewing agent to trigger that privileged workflow itself.

Relevant scripts are:

- `dev-tools/generate_baseline.sh`
- `dev-tools/validate_binaries.sh`
- `dev-tools/run_test_suite.sh`
- `dev-tools/script_data/src/binary_test_cases.py`
- `dev-tools/script_data/src/binaries_suite.py`
- `dev-tools/script_data/src/compare_runs.py`

Read the actual scripts and current workflow for invocation details; some
paths in the dev-tools README were stale at inspection. Inspect structural,
categorical, numerical, and warning/error differences. **Scientific differences
can be reported while the workflow succeeds.** A green status therefore does
not establish unchanged scientific behavior.

Compare equivalent configurations, data releases and stochastic settings, and
identify intended changes explicitly. Do not update accepted baselines as a
substitute for explaining differences. Formula changes benefit from an
independent analytic/reference check in addition to these regressions.

### Appropriate effort

Use the smallest set of checks that addresses the concrete risk. Routine review
does not require a full MESA grid rebuild, all-paper reproduction, large data
downloads, population retraining or production-scale runs. If such work is
needed to validate a material change, describe the experiment, expected
discriminating result, and required inputs for the maintainers. Keep available
CI and scientific regression evidence separate from inference.

## 13. Findings, severity and output

### Actionable findings

An effective finding contains the changed location; conditions under which the
problem occurs; what the implementation does; the intended behavior and
supporting evidence; and the resulting physical, statistical or software
consequence. Cite a primary source for source-dependent scientific claims,
and give a minimal example or analytic limit when practical. Place the comment
on the smallest useful changed range and avoid duplicate downstream symptoms
of the same cause.

A useful prose pattern is:

> Under [specific configuration/domain], [changed operation] produces [incorrect
> result] because [reason]. The relevant [verified relation/reference or test]
> requires [expected behavior], so this changes [scientific consequence].
> Reproduce with [minimal inputs/check]; [bounded correction, if established].

Do not fill placeholders with guessed evidence. Treat suspected defects that
depend on an unknown assumption as questions or evidence gaps, not definitive
findings. An existing pattern outside the diff needs a demonstrated connection
to the PR before being presented as its regression.

### Severity

| Priority | Interpretation for this project |
| --- | --- |
| P0 | An exceptional, broadly catastrophic defect requiring immediate attention, supported by clear evidence; do not use this merely because the code is scientific. |
| P1 | A demonstrated serious error in an accessible, important path, such as systematic corruption of evolutionary outcomes, materially biased scientific rates, or major data loss. State affected conditions and reach. |
| P2 | A real but more limited correctness issue or moderate regression. Do not upgrade it solely to obtain a comment from a restricted interface. |
| P3 / suggestion | Minor improvements or optional clarification; normally omit from automatic review unless specifically requested. |

As documented when this guide was prepared, the built-in GitHub integration
posts P0/P1 findings. These instructions do not promise delivery of P2 findings
or override that filter. A separate, broader review can apply the same guide
when maintainers need moderate issues or an explicit scientific review report.

Where the review surface supports a summary, state the affected scientific
scope, tests actually run, references actually consulted, and material checks
that remain unresolved. Follow the surface's required output format; do not
invent an inline bug merely to communicate an infrastructure limitation.
Silence or an empty findings list is not a statement that the entire codebase
or every relevant paper has been validated.

These instructions guide analysis. They do not establish an enforceable
scientific approval gate or replace maintainer judgment, tests, or required
reviews.

## 14. Sources and reference map

### Instrument papers

- **[P1] Fragos et al. (2023).** *POSYDON: A General-purpose Population
  Synthesis Code with Detailed Binary-evolution Simulations*, ApJS 264, 45.
  [DOI](https://doi.org/10.3847/1538-4365/ac90c1);
  [full preprint](https://arxiv.org/html/2202.05892v2).
- **[P2] Andrews et al. (2025).** *POSYDON Version 2: Population Synthesis with
  Detailed Binary-evolution Simulations across a Cosmological Range of
  Metallicities*, ApJS 281, 3.
  [DOI](https://doi.org/10.3847/1538-4365/adfb78);
  [full preprint](https://arxiv.org/html/2411.02376v2).

The section labels below were inspected in those preprints. Verify the version
when using a published equation to support a particular finding.

| Topic | P1 | P2 |
| --- | --- | --- |
| Architecture and evolution flow | II, IX, X | I, VI.1 |
| Stellar physics, mixing and winds | III | II.2–II.3 |
| Mass transfer and assumptions | IV.2 | II.2.2; IV.1 |
| Grids and processing | V–VI | II; IV |
| Classification/interpolation and constraints | VII, particularly VII.4 | III; Appendix A |
| Track matching and detached evolution | VIII.1 | V.1 |
| CE, collapse and compact remnants | VIII.2–VIII.3 | II.2.1; V.2–V.4 |
| Accretion luminosities | IV.2.3 | V.5 |
| Cosmic populations and detection rates | X for population context | VI.2–VI.4 |
| Limitations and model scope | VII.6; XI | VII |

### Selected original prescriptions

This is a starting map, not a claim that every original equation below was
independently revalidated in preparing this guide. Read the exact relevant
source/version when making a source-dependent PR finding. For engines and
prescriptions not listed here, follow the implementation's actual citation.

- **[R3] Eggleton (1983),** *Approximations to the radii of Roche lobes*,
  ApJ 268, 368–369.
  [Original paper](https://articles.adsabs.harvard.edu/pdf/1983ApJ...268..368E).
- **[R4] Peters (1964),** *Gravitational Radiation and the Motion of Two Point
  Masses*, Physical Review 136, B1224.
  [DOI](https://doi.org/10.1103/PhysRev.136.B1224).
- **[R5] Junker & Schäfer (1992),** *Binary systems: higher order gravitational
  radiation damping and wave emission*, MNRAS 254, 146–164.
  [DOI](https://doi.org/10.1093/mnras/254.1.146);
  [original paper](https://adsabs.harvard.edu/pdf/1992MNRAS.254..146J).
- **[R6] Hut (1981),** *Tidal evolution in close binary systems*, A&A 99, 126.
  [ADS record](https://ui.adsabs.harvard.edu/abs/1981A%26A....99..126H/abstract).
- **[R7] Webbink (1984),** *Double white dwarfs as progenitors of R Coronae
  Borealis stars and Type I supernovae*, ApJ 277, 355.
  [DOI](https://doi.org/10.1086/161701).
- **[R8] Kalogera (1996),** *Orbital Characteristics of Binary Systems after
  Asymmetric Supernova Explosions*, ApJ 471, 352.
  [DOI](https://doi.org/10.1086/177974);
  [preprint](https://arxiv.org/abs/astro-ph/9605186).
- **[R9] Wong, Valsecchi, Fragos & Kalogera (2012),** ApJ 747, 111,
  cited by POSYDON's eccentric pre-SN orbital treatment.
  [DOI](https://doi.org/10.1088/0004-637X/747/2/111).
- **[R10] Fryer et al. (2012),** *Compact Remnant Mass Function: Dependence on
  the Explosion Mechanism and Metallicity*, ApJ 749, 91.
  [DOI](https://doi.org/10.1088/0004-637X/749/1/91).
- **[R11] Hendriks et al. (2023),** *Pulsational pair-instability supernovae in
  gravitational-wave and electromagnetic transients*, MNRAS 526, 4130.
  [DOI](https://doi.org/10.1093/mnras/stad2857).
- **[R12] Xing et al. (2024),** *From ZAMS to merger: Detailed binary evolution
  models of coalescing neutron star-black hole systems at solar metallicity*,
  A&A 683, A144.
  [DOI](https://doi.org/10.1051/0004-6361/202347971).

### Repository and review configuration

All code-specific observations refer to the pinned source snapshot in Section 1.
In particular, inspect these primary files when revising this guide:

- [Interpolation implementation](https://github.com/POSYDON-code/POSYDON/blob/11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1/posydon/interpolation/IF_interpolation.py)
- [Physical constraints](https://github.com/POSYDON-code/POSYDON/blob/11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1/posydon/interpolation/constraints.py)
- [Orbital and stellar utilities](https://github.com/POSYDON-code/POSYDON/blob/11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1/posydon/utils/common_functions.py)
- [Population normalization](https://github.com/POSYDON-code/POSYDON/blob/11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1/posydon/popsyn/norm_pop.py)
- [Cosmic-rate utilities](https://github.com/POSYDON-code/POSYDON/blob/11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1/posydon/popsyn/rate_calculation.py)
- [Scientific comparison workflow](https://github.com/POSYDON-code/POSYDON/blob/11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1/.github/workflows/baseline-validation-on-comment.yml)

Official Codex documentation checked on 20 November 2026:

- [GitHub review and repository instructions](https://learn.chatgpt.com/docs/third-party/github)
- [AGENTS.md scope and discovery](https://learn.chatgpt.com/docs/agent-configuration/agents-md)
- [Companion review-guide pattern](https://learn.chatgpt.com/guides/best-practices)
- [Cloud internet-access configuration](https://learn.chatgpt.com/docs/cloud/internet-access)
