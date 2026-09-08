# Verified bubble publication

Published 2026-09-08 using the authenticated GitHub connector after the section 4.2 review.
Remote commit: `d9df8786553a24de220691e55ca06c4515ff092b`.
Local reviewed commit: `dce73d2e0534f1d9937ea505c0ad674e56ee6219`.
Both have the identical Git tree `8047d6bbb6b5084896b6619713a7640c28f1b722` (236 files).
All 176 unique blobs were checked against their local Git hashes. A short pipe read during binary transfer was detected by a hash mismatch and corrected with full-block reads and byte-length checks before updating main.

The existing remote MIT LICENSE and parent commit fe301db516317568a7984acedb08591626088666 were preserved. The branch was updated without force, then fetched again; the complete tree matched the local reviewed version.
Local intermediate commits remain on local-validation-work; the connector publication has a separate commit identity with the same files.
Local git push has no credential helper; the connector has verified push permission.

Only the reviewed bubble work was published. Porous development after this point is not a reproduction result.

## Porous development publication

Following the user’s additional publication instruction, the component-verified porous development is included. All 71 porous tests passed again before publication (logs/72-porous-prepublication-tests.log). Complete two-phase coupling, static contact-drop CFD and section 4.4 reproduction remain unfinished, as recorded in STATUS.toml and docs/POROUS_PROGRESS_JA.md.

The latest remote main, `cc8a27077398d6e082a9b82d999ffdfcd00a0083`, removes the opening local-handoff sentence from README. That remote edit is retained. Local source and evidence were reviewed at `5f85b1a`; this publication adds status documentation and the repeated test log. Earlier “local only” and pause notes are historical.
