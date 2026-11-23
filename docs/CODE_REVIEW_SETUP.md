# POSYDON Code Review setup

Updated 23 November 2026 for `POSYDON-code/POSYDON`.

These directions describe the repository files and settings needed for review.
They are not a deployment confirmation. Verify the committed revisions and
saved Codex settings when completing or updating the rollout.

## Repository files

| Repository file | Role |
| --- | --- |
| [`AGENTS.md`](../AGENTS.md) | Root instructions, including scientific `## Code Review Rules`. |
| [`docs/CODE_REVIEW_GUIDE.md`](CODE_REVIEW_GUIDE.md) | Scientific procedure, module checks, analytic examples, validation and references. |
| [`SECURITY.md`](../SECURITY.md) | Root security expectations, relevant risks, assumptions and reporting criteria. |
| [`docs/security/threat-model.md`](security/threat-model.md) | POSYDON assets, entry points, trust boundaries and deployment assumptions. |
| `docs/CODE_REVIEW_SETUP.md` | These installation and operational directions. |

The **Open file** rows for `AGENTS.md` and `SECURITY.md` in Codex settings
refer to files **in the repository**. They are not upload fields.
Commit the files through the normal contribution process. The threat-model
setting likewise takes a repository-relative path, not the document's text or
a local computer path.

The root `AGENTS.md` explicitly directs the reviewer to the scientific guide;
the companion filename alone does not cause automatic discovery. These Markdown
references are not automatically added to the existing Sphinx navigation.

The inspected `main` tree at commit
`11a6bb26e54b5a1736ba0a88bb2367f9b147c5f1` had no `AGENTS.md`. Recheck before
installation. Merge into existing instructions rather than overwriting them,
and preserve any existing vulnerability-reporting policy in `SECURITY.md`.
Keep the uppercase filenames. Root guidance covers the repository; nested
files can specialize guidance for the directories they cover.

## Repository settings

Open [Codex code review settings](https://chatgpt.com/settings/code-review) and
select `POSYDON-code/POSYDON`. Use the following configuration:

| Area | Setting | Value or action |
| --- | --- | --- |
| Code review | Automatic review | **On** |
| Code review | Review | **All PRs** |
| Code review | Trigger | **On every push** |
| Code review | Exhaustive review | Existing setting at preparation: **Follow personal preferences**. Inspect the live tooltip and options before changing it; see below. |
| Security review | Auto security review | **On** |
| Security review | Review | **All PRs** |
| Security review | Trigger | **Whenever code review runs**, to cover the same revisions as code review. |
| Security review | Threat model | Repository path: `docs/security/threat-model.md`, after that file is committed and available on the branches being reviewed. |
| Security review | Minimum automatic reporting severity | **High**, reporting High and Critical findings. |
| Security review | Minimum manual reporting severity | **Medium**, reporting Medium, High and Critical findings. |
| Security review | Path-based overrides | Leave empty initially. Add only for a documented, specific exception. |

The initial configuration supplied for this setup shows code review enabled
for all PRs on every push, and security review enabled for all PRs **On PR open**.
The recommended security trigger above changes that initial configuration: **On PR open** checks
the initial PR independently, while **Whenever code review runs** follows code
review, including its later push-triggered runs. Verify the saved value rather
than assuming the recommendation has been applied.

**Exhaustive review:** the official documentation and changelog checked for
this setup do not define this option's implementation, available choices,
severity behavior or pricing multiplier. The live UI establishes that the
control exists. Use its current help text before deciding whether to retain
personal preferences or select a repository-wide behavior; do not assume that
it enables all severities, a particular model or literature access.

The High/Medium thresholds and path overrides above apply to **security
reporting**. They do not change scientific Code Review priorities. The current
published GitHub Code Review guide describes P0/P1 findings; these instructions
do not guarantee P2 comments or override the review surface's output rules.
Use a separate broader review when moderate findings or a written scientific
assessment are required.

GitHub automatic reviews count as Code Review usage. Security reviews also
consume included Codex allowance or ChatGPT credits. Following every code
review can therefore add security-review usage compared with PR-open-only
review. Do not infer a fixed cost per PR or an exhaustive-review multiplier.

## Install, activate and verify

1. Add the four review-context files and these setup directions through a normal
   POSYDON contribution, preserving existing instructions and policies. Merge them into the normal PR target
   branch and carry them into other maintained target branches that need the
   same guidance.
2. Verify the exact threat-model path exists in the repository revisions being
   reviewed before saving that path in Codex settings. Check active PRs that
   predate the new files and update them through the repository's usual branch
   workflow where needed. If no threat-model source is configured, Codex can
   generate one for each security review; a local uncommitted file is not the
   configured repository source.
3. Save and recheck the settings above. The files customize a configured
   reviewer; they do not enable either integration themselves. Repository
   settings require the appropriate GitHub push or administrator permission.
4. On a representative PR, request `@codex review` if needed. For an explicit
   security check, use `@codex security review`; its findings use the manual
   threshold. If security is coupled to code review, avoid requesting an extra
   run unnecessarily. Check task status and the associated report for the
   reviewed revision, guidance used, evidence and actual verification limits.
5. Calibrate with a known consequential defect, a legitimate exception and an
   unrelated change, using existing examples or a controlled test branch.
   Do not merge deliberately incorrect science to test the reviewer. Refine
   rules that repeatedly create noise.

No inline findings does not prove that every check ran or that the review
found no lower-severity issues. Security findings posted on a PR have that
PR's visibility; the fuller Security Report remains available in the associated
Codex task. Review that report when assessing coverage.

Public documentation does not specify which exact base/head/merge checkout
supplies instructions for every hosted review. Committing files to maintained
branches is the practical deployment path; verify behavior on an actual review.
The separate `POSYDON-MESA-INLISTS` repository/submodule needs its own integration
and guidance if all PRs there should receive review.

## Literature checks and scientific validation

The scientific guide requires primary-source checks for material changed
claims when sources are accessible, and clear disclosure when they are not.
An `AGENTS.md` file does not enable browsing, ADS or journal access. General
cloud network documentation does not guarantee live retrieval in every hosted
review. Check the actual task record. Resolve a material inaccessible reference
through a maintainer or a separately configured review; do not describe an
unperformed search as verification.

For new prescriptions, commit concise source notes with the exact citation,
version, variable mapping, units, assumptions, validity range and an independently
checked example. This preserves useful evidence without requiring complete
papers to be stored in the repository.

POSYDON's existing scientific baseline workflow remains useful. An authorized
maintainer can use `/run-baseline-check` on supported same-repository PRs.
Inspect its comparison report: scientific differences can coexist with a
successful workflow status. Neither a green workflow nor an empty automated
review establishes complete scientific correctness. These instructions do not
authorize a reviewing agent to trigger privileged validation workflows.

## Basis and maintenance

The scientific guide identifies its source snapshot and evidence scope;
`docs/security/threat-model.md` separately records observed code and deployment
assumptions. Update path maps, test commands and source notes as the code and
physical prescriptions evolve. Keep root instructions concise. Preserve the
distinction between a physical law, an adopted model, a numerical approximation,
a security boundary and an unresolved question.

## Official configuration sources

- [GitHub Code Review and repository instructions](https://learn.chatgpt.com/docs/third-party/github)
- [AGENTS.md scope and discovery](https://learn.chatgpt.com/docs/agent-configuration/agents-md)
- [Companion review-guide recommendation](https://learn.chatgpt.com/guides/best-practices)
- [Security Review triggers, threat models and reporting thresholds](https://learn.chatgpt.com/docs/security/security-review)
- [AGENTS.md, SECURITY.md and threat-model roles](https://learn.chatgpt.com/use-cases/set-up-security-context)
- [Persistent SECURITY.md guidance](https://learn.chatgpt.com/docs/security/plugin/scans)
- [Usage and pricing](https://learn.chatgpt.com/docs/pricing)
- [Cloud internet-access controls](https://learn.chatgpt.com/docs/cloud/internet-access)
- [ChatGPT and Codex changelog](https://learn.chatgpt.com/docs/changelog)
