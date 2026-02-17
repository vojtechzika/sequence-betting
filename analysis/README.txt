# Sequences 2026 — Analysis Scripts

This folder contains the reproducible data-prep pipeline for the Sequences 2026 project.

## Folder structure (expected)

Repository root:

- `data/raw/`  
  Raw oTree CSV exports. Do not edit these files.

- `data/clean/`  
  Generated, analysis-ready tables (created by scripts).

Analysis project:

- `analysis/scripts/`  
  R scripts implementing the ETL (data prep) steps.

## Conventions

- We keep **one raw folder** (`data/raw/`) containing all session exports.
- We write derived data to dataset-specific folders:
  - `data/clean/pilot/`
  - `data/clean/main/`
- Script numbering reflects pipeline order.

## Scripts

### `scripts/00_setup.R`
Loads packages, defines paths, and defines helper `msg()`.

Key path variables:
- `path_raw` = `../data/raw`
- `path_clean` = `../data/clean`
- `path_clean_ds(ds)` = `../data/clean/<ds>`

### `scripts/01_read_merge_filter.R`
Function: `run_read_merge_filter(dataset)`

- Finds raw CSV files in `data/raw/` by filename pattern:
  - `dataset == "pilot"` → filenames matching `pilot.*.csv`
  - `dataset == "main"`  → filenames matching `main.*.csv`
- Reads each file, filters rows where `participant._current_page_name == "Finished"` (case-insensitive).
- Merges rows across files.
- Saves to:
  - `data/clean/<dataset>/merged.csv`
- Prints per-file counts and total kept/dropped.

### `scripts/02_table_participants.R`
Function: `make_participants_table(dataset)`

- Reads: `data/clean/<dataset>/merged.csv`
- Keeps a specified participant-level column set (participant/session/demographics).
- Renames columns to shorter analysis-friendly names.
- Saves to:
  - `data/clean/<dataset>/participants.csv`

### `scripts/03_table_sequences.R`
Function: `make_sequences_table(dataset)`

- Reads: `data/clean/<dataset>/merged.csv`
- Extracts sequence-level variables from wide oTree format:
  - `sequences.<k>.player.<var>` for k = 1..64 (or whatever exists)
- Keeps **only** the variables listed in `vars` inside the script.
- Reshapes wide → long to produce one row per participant × sequence.
- Renames a subset of columns for readability (e.g., `pos_in_block -> pos`).
- Saves to:
  - `data/clean/<dataset>/sequences.csv`


### `scripts/run_etl.R`
Entry point to run the ETL steps.

- Sources the scripts
- Sets `dataset <- "pilot"` (switch to `"main"` when needed)
- Runs:
  - `run_read_merge_filter(dataset)`
  - `make_participants_table(dataset)`
  - `make_sequences_table(dataset)` (if included in the runner)

## How to run (RStudio)

Open the `analysis` R project (`.Rproj`), then:

```r
source("scripts/run_etl.R")
