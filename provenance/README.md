# Provenance

Git tracks the released files. The report from `acgut reproduce` compares recomputed numerical results with `data/published/`; `acgut validate` checks the released records and source data. The former repository-wide checksum file and its maintenance scripts have been removed.

`sampling_software.csv` records, for each of the 66 CREST runs, the CREST version and commit, backend, command, settings, normal termination and the SHA-256 of the CREST log (logs are not distributed).

Finer-grained checksums live next to the data: `data/calculations/manifest.csv` hashes each calculation record and names its source log by SHA-256; each record hashes its frequency and coordinate files; `data/source/manifest.csv` and `config/sampling.csv` hash the starting structures and CREST ensembles. Every `acgut reproduce` run writes its own `run_metadata.json` with the code, configuration and numerical-input hashes it used.
