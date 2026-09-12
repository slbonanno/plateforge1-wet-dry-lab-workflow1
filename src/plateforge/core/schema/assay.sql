CREATE TABLE IF NOT EXISTS plates (
    plate_id   TEXT PRIMARY KEY,
    barcode    TEXT UNIQUE,
    format     INTEGER NOT NULL,
    role       TEXT,
    is_virtual INTEGER NOT NULL DEFAULT 0,
    file_hash  TEXT,
    meta_json  TEXT NOT NULL DEFAULT '{}',
    notes      TEXT,
    created_at TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS wells (
    plate_id TEXT NOT NULL,
    well     TEXT NOT NULL,
    clone_id TEXT,
    content  TEXT,
    PRIMARY KEY (plate_id, well)
);

CREATE TABLE IF NOT EXISTS reads (
    plate_id   TEXT NOT NULL,
    well       TEXT NOT NULL,
    channel    TEXT NOT NULL,
    value      REAL,
    PRIMARY KEY (plate_id, well, channel)
);

CREATE INDEX IF NOT EXISTS ix_wells_clone ON wells(clone_id);
