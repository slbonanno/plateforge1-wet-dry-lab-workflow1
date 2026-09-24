-- Plates, the wells in them, and the steps that made them.
--
-- The organising idea (decision 0019): the plate is NOT the sample. A clone
-- is the through-line; plates are containers it passes through. Every well
-- records which clone it descends from and which well it came from, so a
-- rearray is the same mechanism as a 1:1 transfer rather than a special case.
--
-- barcode is deliberately NOT UNIQUE. A barcode labels a physical container,
-- and several steps change what is in a container without moving it to a new
-- one -- adding beads, washing, eluting. Each of those states is its own PLT,
-- because artifacts are never modified (rule 10), and they all carry the same
-- barcode because it is the same piece of plastic. Making the column unique
-- meant INSERT OR REPLACE silently deleted the earlier state of the plate,
-- which broke lineage halfway along the chain and was found exactly that way.
-- Order states of one barcode by created_at.

CREATE TABLE IF NOT EXISTS plates (
    plate_id     TEXT PRIMARY KEY,
    barcode      TEXT,               -- NOT unique: see below
    format       INTEGER NOT NULL,
    role         TEXT,
    content_kind TEXT,              -- dna | cells | supernatant | beads |
                                    -- eluate | igg_prep | dilution | assay
    is_virtual   INTEGER NOT NULL DEFAULT 0,
    produced_by  TEXT,              -- the RUN that made it
    file_hash    TEXT,
    meta_json    TEXT NOT NULL DEFAULT '{}',
    notes        TEXT,
    created_at   TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS wells (
    plate_id        TEXT NOT NULL,
    well            TEXT NOT NULL,          -- always A01 (rule 3)
    clone_id        TEXT,                   -- NULL for controls and blanks
    content         TEXT,
    role            TEXT NOT NULL DEFAULT 'sample',  -- sample | control |
                                                     -- blank | empty
    volume_ul       REAL,
    concentration   REAL,
    conc_units      TEXT,
    source_plate_id TEXT,                   -- per-well, not per-plate: that is
    source_well     TEXT,                   -- what makes a rearray ordinary
    flags           TEXT,                   -- per-well deviations
    meta_json       TEXT NOT NULL DEFAULT '{}',
    PRIMARY KEY (plate_id, well)
);

-- One row per step that happened. Planning and recording are separate events
-- hours apart, so a row exists from the moment a step is planned and is
-- completed later.
CREATE TABLE IF NOT EXISTS steps (
    run_id       TEXT PRIMARY KEY,
    step_type    TEXT NOT NULL,
    status       TEXT NOT NULL DEFAULT 'planned',  -- planned | done | abandoned
    inputs_json  TEXT NOT NULL DEFAULT '[]',       -- plate ids in
    outputs_json TEXT NOT NULL DEFAULT '[]',       -- plate ids out
    params_json  TEXT NOT NULL DEFAULT '{}',
    lots_json    TEXT NOT NULL DEFAULT '[]',       -- reagent lots consumed
    operator     TEXT,
    instrument   TEXT,
    planned_at   TEXT NOT NULL,
    started_at   TEXT,
    finished_at  TEXT,
    deviations   TEXT,
    meta_json    TEXT NOT NULL DEFAULT '{}'
);

CREATE TABLE IF NOT EXISTS transfers (
    xfr_id          TEXT NOT NULL,
    run_id          TEXT,
    source_plate_id TEXT NOT NULL,
    source_well     TEXT NOT NULL,
    dest_plate_id   TEXT NOT NULL,
    dest_well       TEXT NOT NULL,
    volume_ul       REAL NOT NULL,
    ordinal         INTEGER NOT NULL DEFAULT 0,
    PRIMARY KEY (xfr_id, ordinal)
);

CREATE TABLE IF NOT EXISTS reads (
    plate_id   TEXT NOT NULL,
    well       TEXT NOT NULL,
    channel    TEXT NOT NULL,
    value      REAL,
    units      TEXT,
    run_id     TEXT,
    PRIMARY KEY (plate_id, well, channel)
);

CREATE INDEX IF NOT EXISTS ix_wells_clone ON wells(clone_id);
CREATE INDEX IF NOT EXISTS ix_wells_source ON wells(source_plate_id, source_well);
CREATE INDEX IF NOT EXISTS ix_plates_kind ON plates(content_kind);
CREATE INDEX IF NOT EXISTS ix_plates_barcode ON plates(barcode, created_at);
CREATE INDEX IF NOT EXISTS ix_steps_type ON steps(step_type, status);
CREATE INDEX IF NOT EXISTS ix_transfers_run ON transfers(run_id);
