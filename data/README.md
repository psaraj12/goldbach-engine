# Data

Checkpoints, miss files, logs, provenance records and the run 6 binary are archived on
Zenodo ([doi:10.5281/zenodo.23082138](https://doi.org/10.5281/zenodo.23082138)). They
are not stored in this repository.

| Folder | Run |
| --- | --- |
| `1_campaign1_v4/` | campaign 1 over [4·10^18, 4.001·10^18]: 13 checkpoints (no logs) |
| `2_campaign_v5/` | v5 campaign over [4·10^18, 4.001·10^18], segments 1 and 2, with the campaign source and its `SHA256SUMS` |
| `3_rerun_v6_segment1/` | v6 recomputation of segment 1 |
| `rerun_seg2/` | v6.2 recomputation of segment 2 (run 3b) |
| `4_extension_v6/` | v6 extension over [4.001·10^18, 4.002·10^18] |
| `5_extension_v62/` | v6.2 extension over [4.002·10^18, 4.003·10^18], with the exact source and the build reproduction log |
| `06_rerun_ext/` | v5 recomputation of [4.001·10^18, 4.003·10^18], with source and binary |

The SHA-256 of every file is listed in
[`../docs/provenance/DATA_SHA256SUMS.txt`](../docs/provenance/DATA_SHA256SUMS.txt).
To check a downloaded copy, run this from the folder that contains the run folders:

```bash
sha256sum -c path/to/goldbach-engine/docs/provenance/DATA_SHA256SUMS.txt
```

To re-check the witnesses and compare two runs:

```bash
python3 scripts/verify_witnesses.py 5_extension_v62/goldbach_v62_checkpoint.csv --minimal
python3 scripts/compare_witnesses.py 06_rerun_ext/goldbach_v5_checkpoint.csv 5_extension_v62/goldbach_v62_checkpoint.csv
```
