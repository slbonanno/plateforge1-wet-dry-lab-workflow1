# 0001 — Modules communicate through registered artifacts

Date: 2026-09-12
Status: accepted

## Context

The project spans sequence handling, plate assays, reagent reasoning, and file
emission, and will be developed in bursts with long gaps. Direct function calls
between modules would produce a dependency tangle that is hard to re-enter.

## Decision

Modules never import each other. A module produces something, registers it in
the ledger (`core.artifacts.register`) with a typed `obj_id`, and returns that
id. A downstream module accepts the id and resolves it. Relationships between
artifacts are recorded as edges in `artifact_links`.

## Consequences

- Every module is independently runnable and testable.
- Provenance is free: `lineage()` walks either direction.
- The cost is indirection — you cannot follow a call graph from module to
  module. The ledger is the call graph.
- Adding a cross-module capability means adding an artifact type, not an import.
