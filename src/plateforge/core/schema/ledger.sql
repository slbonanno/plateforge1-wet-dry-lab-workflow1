CREATE TABLE IF NOT EXISTS artifacts (
    obj_id       TEXT PRIMARY KEY,
    obj_type     TEXT NOT NULL,
    label        TEXT,
    store        TEXT,
    path         TEXT,
    produced_by  TEXT NOT NULL,
    code_version TEXT,
    params_json  TEXT NOT NULL DEFAULT '{}',
    meta_json    TEXT NOT NULL DEFAULT '{}',
    notes        TEXT,
    created_at   TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS artifact_links (
    parent_id    TEXT NOT NULL,
    child_id     TEXT NOT NULL,
    relationship TEXT NOT NULL,
    PRIMARY KEY (parent_id, child_id, relationship)
);

CREATE INDEX IF NOT EXISTS ix_artifacts_type ON artifacts(obj_type);
CREATE INDEX IF NOT EXISTS ix_artifacts_producer ON artifacts(produced_by);
CREATE INDEX IF NOT EXISTS ix_links_parent ON artifact_links(parent_id);
CREATE INDEX IF NOT EXISTS ix_links_child ON artifact_links(child_id);
