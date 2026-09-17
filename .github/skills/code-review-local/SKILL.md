---
name: code-review-local
description: >-
  Runs a local review of a change in OPM/opm-simulators, OPM/opm-common,
  OPM/opm-grid or OPM/opm-upscaling — an open PR or a local branch that has no
  PR yet. For an open PR already checked out and up to date, diffs locally
  instead of re-fetching; otherwise fetches with `gh`. Applies the
  code-review skill's domain checks and writes the findings up as a report,
  delivered in chat by default; posting to GitHub is only offered, never
  automatic, and only for an open PR. Invoked only after build and tests are
  already confirmed green; performs a static review only (no local build or
  ctest). Use when reviewing an OPM pull request or local branch, triaging
  merge readiness, or when the user mentions reviewing opm-simulators,
  opm-common, opm-grid, Flow, gpuistl, or an OPM PR number.
---

# OPM Pull Request Review (local)

This is the process for running a review locally: how to identify and fetch the change (an open PR, or a local branch with no PR yet), what domain checklist to apply, and how to write up the findings. It says nothing about *what* makes a good or bad change — that's in `code-review`, which this skill loads and applies.

Review as an automated reviewer, in your own voice: report what you checked, what you found, and what you could not determine. Do not imitate a maintainer's phrasing and do not present the review as if a person wrote it.

## Step 0 — Load the checklist skill

Before doing anything else, invoke the `code-review` skill. It defines severity/confidence, blocking-issue criteria, path-to-section routing, the domain passes, and the CI/merge conventions used below. This file assumes you have it loaded and does not repeat its content.

## Build and tests are already confirmed — do not re-run them

This skill is invoked only after a separate system, or the user, has already confirmed that the change **builds** and that the relevant **tests pass**. Treat compile success and green tests as given — for a local branch with no CI run yet, that means the user is telling you it already builds and passes locally; take that at face value rather than asking for proof.

**Do not** configure or run a local build, invoke `cmake` / `make` / `ctest`, or spend the review proving the code compiles. Focus on static review: diff, checkout greps, `Read`, and the checklist's domain passes. When a finding needs runtime evidence you cannot get from the diff (bit-identical output, deck-specific regression numbers), ask for it or escalate — do not try to produce it yourself.

You have a checkout locally (a build tree may exist alongside it): do a full static pass plus greps/`Read` against the checkout, but do not build or run `ctest`. Never claim bit-identical output, a specific regression PDF result, or a performance number unless that evidence is in the PR/commit description, CI artifacts, or the author's own comments — do not manufacture it by running the code.

## Step 1 — Identify and fetch the change

Figure out which of these applies — ask the user if it isn't obvious from how they asked:

**1. Local branch, no open PR.** Reviewing work before (or instead of) opening a PR.

```bash
git remote show origin | sed -n 's/.*HEAD branch: //p'   # the repo's base branch, e.g. master
git fetch origin <base>
git log origin/<base>..HEAD --oneline                     # commits on this branch
git diff origin/<base>...HEAD                              # merge-base diff — this is what you review
```

There is no PR title/body: build the Summary's context from the commit messages, and from whatever the user tells you about intent. If there are uncommitted or staged changes the user wants included, look at `git diff HEAD` (working tree) separately and say explicitly in the report that those are uncommitted — don't silently fold them into the committed diff.

**2. Open PR, already checked out and current — fast path.**

```bash
gh pr view <N> --repo OPM/<repo> --json title,body,baseRefName,headRefName,headRefOid,isDraft
git rev-parse HEAD
```

If local `HEAD` matches `headRefOid`, the checkout already *is* the PR head: skip `gh pr diff` and diff locally instead —

```bash
git diff origin/<baseRefName>...HEAD
```

This is the same diff `gh pr diff` would return, without the extra round-trip, and Step 2's greps/`Read` calls run directly against a checkout you already know matches the PR.

**3. Open PR, not checked out or stale locally.**

```bash
gh pr view <N> --repo OPM/<repo> --json title,body,files,additions,deletions,isDraft
gh pr diff <N> --repo OPM/<repo>
```

If a checkout would help — Step 2's checklist prefers verifying against real files over reasoning from the diff alone — and one is cheap to get, `gh pr checkout <N>` and switch to the fast path above. Otherwise review from the diff, and mark findings you could not verify against a checkout as lower confidence per the checklist's rubric.

Whichever mode applies, use the checklist's routing table to decide which domain passes apply, and match review depth to the change — a one-line fix does not get the whole checklist.

## Step 2 — Apply the checklist

Work through `code-review`: blocking issues first, then the routed domain passes, then the cross-cutting passes. Grade every finding with that skill's severity/confidence rubric.

## Step 3 — Write the review

**Tone.** Write factually, as an automated reviewer. Report findings; don't perform a persona, don't quote maintainers, and don't write anything that implies a human authored the review.

Tone follows confidence, not style. State a High-confidence finding as a fact and name the evidence ("`test_Serialization` has no entry for this type"). State a Medium-confidence finding conditionally ("if `X` can be null here, this dereference is unguarded"). Phrase a Low-confidence item as a question, because you genuinely do not know the answer. Never soften a verified finding to sound polite, and never harden a guess to sound authoritative. Always say what you did and did not check.

The **Summary** and **Readiness** lines are for the human deciding whether to merge or open the PR: state the verdict plainly and don't hedge it.

**Hand over the replacement code.** Paste a fenced C++ code block with the suggested rewrite rather than describing it in prose. GitHub `suggestion` blocks are effectively unused in these repos — don't introduce them.

Be specific and link to the line in master you're comparing against. Don't guess: when a finding is outside what you can verify, put it under Questions and say exactly what needs confirming — leave it to the human reviewer to decide who to ask, since no single person owns an area in this project.

Output template:

```markdown
## Summary
<2-3 sentences: what the change does, and the single thing that decides whether it's ready.>

**Readiness: Ready | Ready after fixes | Needs major rework**

## Blocking
- `path/file.cpp:123` — [confidence: High] <issue, why it matters, suggested fix>

## Should fix before merge
- `path/file.hpp:45` — [confidence: Medium] <issue, with a fenced C++ snippet of the suggested rewrite>

## Nits
<max 5 if anything is under Blocking, otherwise all of them, grouped by file>
- `path/file.cpp:200` — <...>

## Questions
- <what you need author knowledge for, e.g. mirror in MultisegmentWell? No severity here.>

## Test coverage assessment
- <what is covered, what is not, and whether the gap is acceptable>

## Verification
- Build/tests: assumed already confirmed (not re-run by this review)
- Checks verified by grep/read: <e.g. "new files present in CMakeLists_files.cmake">
- Recommended review state: `APPROVED` | `COMMENTED` | `CHANGES_REQUESTED` (open PR only — omit for a local branch with no PR to post to)
- Recommended CI: `jenkins build this <flags> please` (open PR only; see the checklist's CI and merge protocol section for which flags apply)
```

Write "None" under any empty section rather than dropping it — an empty Blocking section is information. For a local branch with no open PR, drop the two PR-only Verification lines entirely rather than writing "N/A" — everything else in the template applies unchanged.

## Step 4 — Offer to post (open PR only)

The review is always delivered in chat first — that's the outcome regardless of what happens next. Nothing is posted to GitHub automatically, ever.

If this was an open PR (Step 1 modes 2 or 3), after the write-up ask the user explicitly whether to post it to GitHub, and in what form (an issue comment with the report, inline comments, and/or a review state such as `COMMENTED` / `CHANGES_REQUESTED` / `APPROVED`). Only proceed on an explicit yes for *this* review — `gh pr review`, `gh pr comment`, or `gh api` write calls otherwise stay off the table, and so does any Jenkins trigger comment or `@`-mention. What you post must match what you already showed in chat, not a re-derived version.

If this was a local branch with no open PR (Step 1 mode 1), there is nothing to post — skip this step.

## Do not

- Build the project or run tests/`ctest` as part of the review — that work is done before this skill runs.
- Claim bit-identical output, regression numbers, or performance results you did not see in the PR/CI evidence.
- Write in a maintainer's voice, quote a maintainer, or imply that a person authored the review.
- Post to GitHub without an explicit yes for that specific review — approval on one run is not standing approval for the next.
