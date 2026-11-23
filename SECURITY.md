# POSYDON security review guidance

## Purpose

Use this file for security review of POSYDON pull requests, together with
[`docs/security/threat-model.md`](docs/security/threat-model.md). Apply the
scientific correctness guidance in
[`docs/CODE_REVIEW_GUIDE.md`](docs/CODE_REVIEW_GUIDE.md) where the same change
affects physical results or population statistics.

This guidance describes review expectations for scientific Python, local
command-line and HPC workflows. It does not establish a vulnerability-reporting
address, response timetable, supported-version policy, or an automatically
enforced security gate. The threat model separates observed implementation
from deployment assumptions that maintainers must validate.

## Code Review Rules

### Establish the boundary and consequence

- Follow the changed input from its origin to file access, deserialization,
  imports, command execution, scheduler submission, CI execution or publication.
  Identify who can control it and the authority under which it is processed.
- Require a reachable path and a concrete consequence: unintended execution,
  unauthorized read/write, credential exposure, artifact or result corruption,
  privilege misuse, or exhaustion beyond the intended resource allocation.
  Identify the conditions needed for the concern to apply.
- Distinguish trusted researcher-supplied code/configuration from externally
  supplied data. A researcher deliberately choosing a local executable or
  custom step is not, by itself, an injection vulnerability. A change that
  feeds less-trusted input into the same mechanism can create a new boundary.
- Review defects introduced, exposed or materially worsened by the PR. Do not
  turn an unrelated existing design choice into a new regression. Separate a
  demonstrated vulnerability, a deployment-dependent concern and optional
  hardening. Do not infer a public server or remote upload endpoint that the
  code and deployment evidence do not establish.

### Review the relevant POSYDON surfaces

1. **Downloads and archives.** Inspect dataset provenance, download destination,
   integrity-check behavior and extraction paths. Consider absolute paths,
   parent traversal, links, unexpected file types, existing-file replacement
   and resource use. A checksum intended to detect accidental corruption is
   not sufficient evidence of publisher authenticity. Preserve the distinction
   between curated dataset entries and arbitrary external URLs.
2. **Models and configuration.** Treat pickle models and configurations that
   select executable imports/hooks as code-capable inputs. Check their origin,
   permissions and any newly introduced input channels. Preserve constrained
   literal/AST parsing; do not describe it as unrestricted `eval`. Do not assume
   all HDF files or all numeric arrays have pickle's execution semantics.
3. **Local execution and HPC.** Trace paths and arguments into shell commands,
   generated bash/SLURM files and cleanup operations. Check quoting, newlines,
   option interpretation, working-directory assumptions, symlinks, overwrite
   scope and job dependencies. Distinguish a control parameter from intentional
   user-supplied executable code, and a large authorized run from resource abuse.
4. **PR validation and self-hosted CI.** Inspect event and permission checks,
   repository origin, the exact executed revision, credentials, caches,
   artifacts and runner persistence. Preserve the intended authorization
   boundary when changing who can run PR-controlled code. A same-repository
   branch is not automatically equivalent to reviewed release code; disabling
   persisted checkout credentials does not establish complete host isolation.
5. **Release and documentation publishing.** Check triggers, refs, action and
   dependency provenance, token scope and publishing destinations. Keep
   unreviewed input from acquiring publishing authority or exposing credentials
   through build output, artifacts or logs. Establish the actual job's access;
   do not assume every CI job receives release credentials.
6. **Scientific data integrity.** Where an input can cross a real trust boundary,
   assess silent changes to models, grids, provenance, selected systems or
   results. Ordinary formula errors and model uncertainty remain scientific
   correctness concerns unless a relevant security mechanism is involved.

### Verify and communicate

Use the actual base/head code and current workflow configuration. Consult the
threat model's source map, then follow moved implementations and changed
deployment assumptions. Do not assume the documented snapshot is permanent.

Prefer a small inert test or precise data-flow explanation. Do not execute
untrusted pickles/imports, extract suspicious archives into a working data
tree, submit cluster jobs, trigger privileged workflows or publish artifacts
merely to substantiate a routine review. Any necessary execution should use
an appropriately isolated, authorized environment with bounded resources.

For each finding, state the changed location, controllable input, violated
boundary, necessary preconditions and concrete impact. Explain which existing
controls were checked. Assign severity from reach and consequence, not from
the presence of a risky-looking API alone. Follow the review surface's output
and severity constraints; do not inflate a finding to obtain a comment.

State what was inspected or tested and what deployment fact remains unknown.
Keep secrets, private datasets and sensitive infrastructure details out of
review evidence. This file does not authorize repository-setting changes,
workflow execution, external disclosure or code modifications during a review.
