CREATE TABLE IF NOT EXISTS reagents (
    rgt_id     TEXT PRIMARY KEY,
    name       TEXT NOT NULL,
    category   TEXT NOT NULL,
    attrs_json TEXT NOT NULL DEFAULT '{}',
    notes      TEXT,
    created_at TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS lots (
    lot_id     TEXT PRIMARY KEY,
    rgt_id     TEXT NOT NULL,
    lot_number TEXT,
    vendor     TEXT,
    catalog    TEXT,
    expires_on TEXT,
    meta_json  TEXT NOT NULL DEFAULT '{}',
    notes      TEXT,
    created_at TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS usage (
    run_id     TEXT NOT NULL,
    lot_id     TEXT NOT NULL,
    role       TEXT NOT NULL,
    meta_json  TEXT NOT NULL DEFAULT '{}',
    PRIMARY KEY (run_id, lot_id, role)
);

CREATE TABLE IF NOT EXISTS caveats (
    cav_id     TEXT PRIMARY KEY,
    run_id     TEXT NOT NULL,
    rule       TEXT NOT NULL,
    severity   TEXT NOT NULL,
    message    TEXT NOT NULL,
    meta_json  TEXT NOT NULL DEFAULT '{}',
    created_at TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS ix_lots_reagent ON lots(rgt_id);
CREATE INDEX IF NOT EXISTS ix_caveats_run ON caveats(run_id);
