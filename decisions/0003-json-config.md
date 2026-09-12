# 0003 — Experiment definitions are JSON, validated against a schema

Date: 2026-09-12
Status: accepted

## Context

Experiment setup (plate format, replicate scheme, control layout, reagents,
cloning notes) needs to be captured from the user today and written by an agent
later. Both paths must produce the same artifact.

## Decision

JSON, validated with `jsonschema`. The interactive prompt is a front end that
writes the JSON; an agent is a different front end writing the same JSON.

## Consequences

- No comments in the format; use a `"_comment"` key where needed.
- The schema is the single source of truth for what an experiment can express;
  extending the workflow means extending the schema first.
