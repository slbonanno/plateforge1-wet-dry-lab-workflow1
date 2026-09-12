CREATE TABLE IF NOT EXISTS sequences (
    seq_id     TEXT PRIMARY KEY,
    source     TEXT NOT NULL,
    source_ref TEXT,
    chain      TEXT,
    aa_seq     TEXT NOT NULL,
    meta_json  TEXT NOT NULL DEFAULT '{}',
    notes      TEXT,
    created_at TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS library_members (
    lib_id  TEXT NOT NULL,
    seq_id  TEXT NOT NULL,
    PRIMARY KEY (lib_id, seq_id)
);

CREATE INDEX IF NOT EXISTS ix_sequences_chain ON sequences(chain);
