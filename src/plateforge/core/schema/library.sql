-- Index side of the sequence library. Full OAS rows live in parquet
-- (core.bulk, table "oas_pool"); this table holds the columns we filter on.
CREATE TABLE IF NOT EXISTS sequences (
    seq_id        TEXT PRIMARY KEY,
    source        TEXT NOT NULL,
    source_ref    TEXT,
    chain         TEXT,
    aa_seq        TEXT NOT NULL,
    v_call        TEXT,
    v_gene        TEXT,
    v_family      TEXT,
    j_call        TEXT,
    j_gene        TEXT,
    cdr3_aa       TEXT,
    cdr3_len      INTEGER,
    redundancy    INTEGER,
    anarci_status TEXT,
    has_liability INTEGER DEFAULT 0,
    species       TEXT,
    btype         TEXT,
    disease       TEXT,
    isotype       TEXT,
    meta_json     TEXT NOT NULL DEFAULT '{}',
    notes         TEXT,
    created_at    TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS library_members (
    lib_id  TEXT NOT NULL,
    seq_id  TEXT NOT NULL,
    rank    INTEGER,
    PRIMARY KEY (lib_id, seq_id)
);

CREATE INDEX IF NOT EXISTS ix_sequences_chain ON sequences(chain);
CREATE INDEX IF NOT EXISTS ix_sequences_vgene ON sequences(v_gene);
CREATE INDEX IF NOT EXISTS ix_sequences_cdr3len ON sequences(cdr3_len);
CREATE INDEX IF NOT EXISTS ix_members_seq ON library_members(seq_id);
