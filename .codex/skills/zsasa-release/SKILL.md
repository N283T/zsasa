---
name: zsasa-release
description: Use when preparing, merging, tagging, or troubleshooting a zsasa project release, including requests like release, version bump, changelog, create release PR, merge release PR, tag push, or publish vX.Y.Z.
---

# zsasa Release

Follow the **Release** section of `AGENTS.md` in the `N283T/zsasa` repository instead of a generic release flow. It is the single description of the procedure: version bump, checks, rehearsal of the publish workflow, merge and tag, and the packaging checksums afterwards.

Points that are easy to miss:

- The tag push is the publish trigger and cannot be undone. Merge the release PR and push the tag only after the user approves, and only with CI and the rehearsal green.
- Tag the merge commit on `main`, never the release branch.
- Do not edit `packaging/aur/` in the release PR; its checksums exist only after the release is published.
- A skipped `deploy` check on a pull request is normal.
- When the user asks to proceed after CI, do not stop at "the PR is open": carry on through merge, tag and the post-release checksum PR, asking for approval where the procedure says so.
