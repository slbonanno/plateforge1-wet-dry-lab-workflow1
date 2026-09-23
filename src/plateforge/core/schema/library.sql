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
    n_ambiguous   INTEGER DEFAULT 0,
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

-- Clones: a sequence realised in a construct (decisions/0005). Minted at
-- construct assignment, never recomputed, and carried forward from here.
CREATE TABLE IF NOT EXISTS clones (
    clone_id            TEXT PRIMARY KEY,
    seq_id              TEXT NOT NULL,
    construct           TEXT NOT NULL,
    v_gene              TEXT,
    cdr3_aa             TEXT,
    linker              TEXT,
    orientation         TEXT,
    clone_aa_seq        TEXT NOT NULL,
    aa_length           INTEGER,
    clone_dna_seq       TEXT,
    dna_length          INTEGER,
    gc                  REAL,
    sites_removed       INTEGER,
    gc_swaps            INTEGER,
    synthesis_warnings  TEXT,
    germline            TEXT,
    germline_source     TEXT,
    well                TEXT,
    plate_barcode       TEXT,
    parent_clone_id     TEXT,
    meta_json           TEXT NOT NULL DEFAULT '{}',
    notes               TEXT,
    created_at          TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS ix_clones_seq ON clones(seq_id);
CREATE INDEX IF NOT EXISTS ix_clones_plate ON clones(plate_barcode);
