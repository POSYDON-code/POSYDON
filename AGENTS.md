# POSYDON agent guidance

## Code Review Rules

These rules govern code review. POSYDON predicts stellar and binary evolution
and populations from detailed MESA grids, interpolation, and additional
physical prescriptions. A change can run successfully while producing
scientifically incorrect results. Review scientific validity together with
software correctness.

### Read the scientific review guide

Before reviewing a PR, read Sections 1–4 and 12–13 of
`docs/CODE_REVIEW_GUIDE.md`, then the topic sections selected by its module map.
Follow the changed calculation through its inputs, callers, state transitions,
stored products, and downstream scientific outputs. Apply only checks relevant
to the change; an unrelated documentation edit does not require a physics audit.

Use the PR's actual base/head revisions and applicable nested `AGENTS.md`
instructions. The guide's reference snapshot is context, not a requirement that
future code retain the same implementation or physical prescription.

### Establish and verify the intended science

- For a changed formula, algorithm, fit, threshold, or physical assumption,
  identify the intended prescription, domain, approximations, variables, units,
  and conventions. Check algebra, signs, coefficients, powers, normalization,
  limiting cases, and numerical implementation. Trace enough surrounding code
  to establish a reachable consequence.
- Check the cited primary literature when introducing or changing a scientific
  prescription, resolving an ambiguity, or substantiating a scientific defect.
  Start from the implementation's citation and relevant POSYDON paper, then
  consult the original source and applicable corrections when accessible.
  Follow the literature-verification procedure in Section 4 of the guide.
- Cite only sources actually consulted. Give a verified section, equation, or
  page when available, and explain how its symbols and assumptions map to the
  code. Never invent citations or claim to have accessed an unavailable paper.
  If access is unavailable, state the precise unresolved comparison in the
  review summary where the interface supports it; do not treat lack of access
  alone as a code defect or a successful scientific check.
- Distinguish an implementation error, an unintended change to the model, a
  numerical approximation, and an intentional choice among defensible physical
  prescriptions. A newer or different paper is not, by itself, evidence that
  the selected prescription is wrong. Fundamental inconsistencies still need
  to be examined even when a formula appears in an existing paper or comment.

### Protect POSYDON's scientific contracts

- Verify units at interfaces: cgs versus solar units, years versus Myr/Gyr,
  seconds versus days, km/s versus cm/s, absolute metallicity versus
  `Z/Zsun`, linear versus logarithmic quantities, and angular conventions.
  Reuse the repository's constants consistently; changing a constant can also
  change compatibility with existing grids and trained models.
- Preserve component identity, donor/accretor roles, histories, and event
  ordering through mass transfer, reverse mass transfer, CE, collapse, and
  disruption. Do not assume the initial primary always remains the donor or
  most massive star. Allow equal timestamps for instantaneous transitions.
- Apply conservation laws and physical bounds under their actual assumptions.
  Account for winds, ejecta, radiation, neutrinos, and spin/orbit exchange.
  The black-hole Kerr bound is not a bound on a progenitor star's spin proxy;
  tides do not universally shrink an orbit; defined missing values are not zero.
- For interpolation and classification, check physical coordinates and scaling,
  class and metallicity selection, invalid-track masks, out-of-domain behavior,
  and physical constraints. Preserve independent evaluation data. Do not impose
  smoothness across distinct physical outcomes or hide defects by clipping.
- For population statistics, establish the probability measure and represented
  stellar mass before checking normalization. Trace all Jacobians, selection
  weights, time and volume factors, and source/observer-frame conversions.
  Preserve system IDs and independent random streams.
- Treat changes to grid selection, MESA inlists or patches, termination flags,
  thresholds, post-processing, defaults, and serialized model/data schemas as
  possible scientific changes. Check compatibility and provenance throughout
  the affected chain; a Python parameter cannot retroactively alter grid physics.

### Verify and report proportionately

- Prefer a small reproducer, an independent analytic limit, or an independently
  evaluated reference calculation for a suspected scientific defect. Existing
  outputs and tests alone do not establish physical correctness.
- Run relevant available tests when practical. Inspect scientific regression
  differences, not just workflow status. POSYDON's baseline comparison can
  report scientific differences while the workflow succeeds. State which
  checks ran, failed, or could not run; never claim execution from inspection.
- Do not launch production populations, rebuild MESA grids, download large
  data releases, or regenerate accepted baselines as a routine review step.
  Identify a specific additional validation need when existing resources are
  insufficient. Leave formatting and mechanical checks to CI.
- Report actionable defects introduced, exposed, or materially worsened by the
  PR. Each finding should identify the changed location, triggering conditions,
  evidence, and scientific or software consequence. Separate unresolved
  questions from demonstrated defects. Do not report unrelated pre-existing
  issues as new regressions or produce findings to fill a quota.
- Assign severity from demonstrated impact and reach. A silent, systematic
  error in physical outcomes or population normalization can be high priority;
  a scientific topic or a missing test does not automatically make a finding
  P1. Respect the review interface's severity and output constraints.
- Do not edit code, approve or merge the PR, change repository settings, or
  trigger privileged workflows merely because a review was requested.
