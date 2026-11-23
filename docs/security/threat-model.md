# POSYDON threat model for code review

## 1. Status and evidence scope

Prepared **23 November 2026** from read-only inspection of representative
source and workflow files in
[`POSYDON-code/POSYDON`, commit `11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1`](https://github.com/POSYDON-code/POSYDON/tree/11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1),
then the `main` head. This is review context, not a security audit, list of
confirmed vulnerabilities or statement that deployment controls are effective.
No exploit, suspicious artifact, simulation, cluster job or CI workflow was
executed to prepare it.

For each PR, use its actual revisions and relevant deployment evidence.
Update this model when entry points, privileges, artifact sources or execution
environments change. Read the repository's
[`SECURITY.md`](../../SECURITY.md) and scientific
[`CODE_REVIEW_GUIDE.md`](../CODE_REVIEW_GUIDE.md) alongside it.

## 2. System and deployment assumptions

The inspected code provides Python modules and command-line tools for stellar
and binary evolution, population synthesis, MESA grid processing, model
interpolation and analysis. It reads local scientific data, downloads selected
data releases, supports custom Python evolution steps, and generates or runs
local/HPC work. GitHub workflows test code, run selected scientific comparisons
and publish project artifacts.

| Context | Observed in the repository | Assumption requiring deployment confirmation |
| --- | --- | --- |
| Researcher workstation | Local paths, Python imports, model loading and command-line operations. | The operator controls the account; the intended data/configuration/model sources and directory ownership are known. |
| HPC work | Generated SLURM submission scripts and MESA execution/cleanup. | Scheduler accounts, filesystem sharing, quotas, job isolation and executable permissions are configured by the relevant site. |
| Scientific validation CI | A gated workflow executes PR-head code on a self-hosted runner. | Runner provisioning, reset/cleanup, mounted filesystems, network access, installed credentials and protection from cross-job persistence are not established by the YAML alone. |
| Publishing | Release/documentation workflows receive tokens or publishing credentials. | Maintainer roles, release authorization, secret management and any environment protection are configured outside the inspected source. |

This model does not establish an internet-facing, multi-tenant service or a
remote upload interface. If a PR introduces such a wrapper, identify its
authentication, input and privilege boundaries explicitly; local trusted-input
assumptions cannot simply be carried over to that deployment.

## 3. Assets and input origins

| Asset | Relevant failure or misuse |
| --- | --- |
| Source, packages, containers and documentation | Unauthorized modification, publication or replacement of trusted outputs. |
| Grids, interpolators, population files and provenance | Malicious replacement or silent corruption that changes scientific results or disguises their origin. |
| Researchers' and runners' files | Reads, overwrites, deletion or execution beyond the selected task and directories. |
| Credentials and publishing authority | Exposure in logs/artifacts, unintended inheritance by executed code, or use outside the authorized job. |
| HPC/runner resources | Work beyond authorized allocation or persistence that affects other jobs. |

Relevant input origins are repository/PR content, selected remote data releases,
local model and population files, user configuration/import hooks, command-line
arguments, environment paths and scheduler fields. They are not interchangeable
trust categories. For example, a trusted maintainer-selected pickle and a
pickle delivered by an arbitrary new upload channel have different implications
even when they reach the same loader.

A security finding should identify a plausible party or process controlling
the relevant input, the access it starts with, and the authority gained or
asset affected. Do not assume a compromised maintainer, repository, data host
or worker without explaining why that assumption is needed. Review improvements
to resilience against such compromise as hardening unless a concrete PR
regression establishes a narrower issue.

## 4. Trust boundaries and review checks

### A. Published data to local filesystem

**Observed.** `posydon/utils/datasets.py` supplies a curated registry of HTTPS
Zenodo URLs and expected MD5 values. In `data_download.py`,
`download_one_dataset` selects a registry entry, downloads it with
`urllib.request.urlretrieve`, optionally checks MD5, and extracts members with
`tarfile` into the configured data root. The downloaded archive is then removed.

**Review checks.** Trace changes to the registry, allowed sources, redirects,
destination selection, checking options and failure paths. Distinguish
accidental-corruption checks from an authenticated publisher/provenance
decision. An integrity-check exception or bypass must not be represented as
successful verification.

For extraction changes, assess whether members can write outside the intended
root through absolute names, `..`, symbolic/hard links or unusual member types;
whether unexpected existing files can be replaced; and whether size or member
counts can exceed the operation's intended resources. Inspect the Python
version and actual extraction controls rather than assuming library defaults.
These are desired review properties, not a claim that the current extractor
already enforces each one.

The observed downloader selects known registry entries. Do not report an
arbitrary-URL input problem without finding an actual path that allows it.

### B. Scientific artifacts to Python execution

**Observed.** `posydon/interpolation/IF_interpolation.py` loads interpolation
models and saved attributes with `pickle.load`.

**Review checks.** Require code-level trust in those model files before loading them.
Trace who creates or supplies them, where they are stored, and whether a PR
widens the accepted source, path or caller. Scientific usefulness, a familiar
filename or an HDF-like surrounding workflow does not establish trust in a
pickle. If a new service or automation accepts models from less-trusted sources,
review that boundary before calling the existing loader.

Keep source/data/model identity and the associated preprocessing metadata
consistent. A result-integrity concern should explain the unauthorized input
path; a mismatched training version or physics error alone may instead be a
scientific correctness defect.

### C. Configuration to custom code

**Observed.** `posydon/popsyn/io.py` parses literal values and can import
configured steps/hooks. An absolute-import path uses
`importlib.util.spec_from_file_location` and `exec_module`.
`posydon/utils/configfile.py` separately evaluates a restricted set of AST
nodes for arithmetic/configuration values; it is not an unrestricted `eval`.

**Review checks.** Configuration that selects custom Python implementations
requires code-level trust. Inspect changes that automatically consume another
user's configuration or move these operations into a privileged context.
Preserve parser restrictions, explicit import intent and path provenance.
Do not claim that parsing a literal makes the subsequent configured import
safe for arbitrary external input. Conversely, deliberate trusted local
extensions are part of the scientific workflow and are not automatically bugs.

### D. Paths and parameters to shell/HPC operations

**Observed.** `bin/posydon-run-grid:run_mesa` changes working directory,
constructs a shell invocation of the MESA executable and output path, and
cleans selected output patterns. `posydon/CLI/io.py` and
`bin/posydon-setup-pipeline` generate bash and SLURM submission scripts.

**Review checks.** Establish which fields are executable configuration and
which should remain literal paths or scheduler parameters. Follow quoting and
escaping through generated scripts as well as direct calls. Inspect shell
metacharacters/newlines, leading-option interpretation, spaces, missing paths,
symlink targets and working-directory failure behavior where relevant.

Cleanup must remain confined to intended work products; generated submissions
must retain intended dependencies, account and resource scope. Evaluate
overwrites, failed runs and retries for unintended effects on shared data or
jobs. A malicious-input finding needs a less-trusted input route, while an
accidental cleanup error may be independently reportable as data loss. A large
simulation requested by the operator is not itself an availability attack.

### E. PR content to self-hosted validation

**Observed.** `.github/workflows/baseline-validation-on-comment.yml` responds
to `/run-baseline-check`. It checks that the requester has `admin` or `write`
permission, rejects a PR whose head repository differs from its base repository,
records the head/base SHAs, and runs PR-head scripts on a self-hosted runner.
Checkout uses `persist-credentials: false`. Declared workflow permissions
include read access to contents and write access for issues and pull requests.

**Review checks.** Preserve the intended permission, origin and revision
boundary. Review any broadening of allowed callers, fork execution, checkout
refs or runner privileges together with its replacement controls. Permission
to trigger a particular run is not a reason to execute an unrelated revision.

Assess PR-controlled scripts and test definitions, cache/artifact provenance,
comment/report handling, credential availability and reusable host state.
Disabling persisted Git credentials is a useful specific control; it says
nothing by itself about other credentials, mounted resources or persistence
on that host. Ask for the exact missing deployment fact rather than assuming
the runner is either isolated or compromised.

Security review must not trigger this workflow merely to inspect it. Ordinary
unit CI and this privileged scientific-validation path have distinct inputs
and execution contexts.

### F. Builds to publication authority

**Observed.** `publish-to-anaconda.yml` runs on a published release, passes
Conda credentials to the publishing action, and includes a GHCR build/push job
with `packages: write`. Documentation workflows deploy with a GitHub token
from development-push or release-triggered work.

**Review checks.** Trace changed event filters, selected refs, code/actions
executed during builds, destinations and token permissions. Prevent a change
from giving an unintended source of input the ability to publish, overwrite
trusted artifacts, or expose credentials through build output. Check that
published artifacts correspond to the intended source and version.

Do not infer that all CI jobs have these credentials. Secret values,
organization settings and external release protections were not inspected.
Recommendations to alter those settings should identify the specific control
needed rather than claiming an unenforced policy already exists.

## 5. Evidence, validation and maintenance

For a candidate finding, record the changed location, input origin, violated
boundary, preconditions, existing controls and impact. A small inert archive
path fixture or a controlled script-generation test may establish a property
without touching real data or submitting a job. Use an isolated authorized
environment for any necessary execution; inspection alone must not be described
as a successful exploit or security test.

Separate proven behavior from assumptions and from proposed hardening. An
unknown runner configuration, inaccessible secret setting or missing test is
not proof of a vulnerability. Unrelated pre-existing issues do not become new
PR findings. Keep sensitive material out of public review examples.

Update this model when POSYDON adopts new schedulers, download mechanisms,
model formats, execution services or workflow privileges. A remote service,
multi-user artifact store or new external-model interface would require
additional boundaries and deployment evidence beyond this snapshot.

## 6. Primary source map

These sources establish the observed implementation, not the effectiveness of
all proposed review properties:

- [Dataset registry](https://github.com/POSYDON-code/POSYDON/blob/11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1/posydon/utils/datasets.py)
- [Download and extraction](https://github.com/POSYDON-code/POSYDON/blob/11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1/posydon/utils/data_download.py)
- [Interpolation serialization](https://github.com/POSYDON-code/POSYDON/blob/11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1/posydon/interpolation/IF_interpolation.py)
- [Configured population imports](https://github.com/POSYDON-code/POSYDON/blob/11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1/posydon/popsyn/io.py)
- [Configuration arithmetic parser](https://github.com/POSYDON-code/POSYDON/blob/11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1/posydon/utils/configfile.py)
- [MESA execution](https://github.com/POSYDON-code/POSYDON/blob/11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1/bin/posydon-run-grid)
- [Population job scripts](https://github.com/POSYDON-code/POSYDON/blob/11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1/posydon/CLI/io.py)
- [Pipeline job scripts](https://github.com/POSYDON-code/POSYDON/blob/11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1/bin/posydon-setup-pipeline)
- [Self-hosted baseline validation](https://github.com/POSYDON-code/POSYDON/blob/11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1/.github/workflows/baseline-validation-on-comment.yml)
- [Release publication](https://github.com/POSYDON-code/POSYDON/blob/11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1/.github/workflows/publish-to-anaconda.yml)
- [Development documentation publication](https://github.com/POSYDON-code/POSYDON/blob/11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1/.github/workflows/deploy-github-pages-development.yml)
- [Release documentation publication](https://github.com/POSYDON-code/POSYDON/blob/11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1/.github/workflows/deploy-github-pages-release.yml)
