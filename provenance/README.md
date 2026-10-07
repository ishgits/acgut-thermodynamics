# Provenance

`file_checksums.csv` records the SHA-256 and size of every released file except itself, ignored user runs and caches. `python workflows/check_integrity.py` verifies it; `python workflows/write_checksums.py` rewrites it after an intentional change.

`sampling_software.csv` records, for each of the 66 CREST runs, the CREST version and commit, backend, command, settings, normal termination and the SHA-256 of the CREST log (logs are not distributed).

Finer-grained checksums live next to the data: `data/calculations/manifest.csv` hashes each calculation record and names its source log by SHA-256; each record hashes its frequency and coordinate files; `data/source/manifest.csv` and `config/sampling.csv` hash the starting structures and CREST ensembles. Every `atcgu reproduce` run writes its own `run_metadata.json` with the code, configuration and numerical-input hashes it used.
