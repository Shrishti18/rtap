# Priority 1 — JARVIS-WTB: Wannier tight-binding Hamiltonians

Garrity, K. F. & Choudhary, K. *Database of Wannier tight-binding Hamiltonians
using high-throughput density functional theory.* **Sci. Data 8, 106 (2021).**
<https://doi.org/10.1038/s41597-021-00885-z>

Data source: figshare collection
[10.6084/m9.figshare.c.5192276.v2](https://doi.org/10.6084/m9.figshare.c.5192276.v2)
(4 articles), which is the deposit referenced by the paper. Downloaded via the
figshare REST API; **every file was verified against the figshare-supplied MD5
and byte size.**

## What is here

| Path | Contents |
| --- | --- |
| `zips/` | **1,772 archives**, one per material, exactly as deposited. 12.6 GB |
| `extracted_small/<JVASP-id>/` | `POSCAR`, `wannier90.win`, `INCAR`, `KPOINTS`, `OSZICAR` unpacked for every material |
| `jarvis_wtb_manifest.csv` | 1,772 rows — full per-material manifest (see below) |
| `jvasp_id_formula_spacegroup.csv` | the compact **id → formula / space group** map |
| `figshare_filelist.json` | figshare file ids, sizes, MD5s, download URLs |
| `jarvis_dft_3d.json`, `jarvis_dft_2d.json` | JARVIS-DFT metadata dumps (93,902 + 1,103 entries) from `jarvis-tools` |
| `download_figshare.py` | the downloader (resumable, MD5-verified) |
| `build_manifest.py` | archive verification + metadata join |
| `extract_hr_dat.py` | unpack `wannier90_hr.dat` for all or selected materials |

## Coverage

* **1,772 archives — 1,407 3D + 365 2D.** The paper states 1,771 (1,406 3D +
  365 2D); the figshare deposit carries one extra 3D entry. All 1,772 are here.
* **1,772 / 1,772 contain `wannier90_hr.dat`** and `wannier90.win`;
  1,771 contain `wannier90.wout` (`JVASP-6145` has no `.wout`).
* **1,772 / 1,772 have formula + space group.** 1,755 came from the JARVIS-DFT
  `dft_3d`/`dft_2d` datasets; the remaining 17 were derived from the deposited
  `POSCAR` with spglib — the `metadata_source` column records which.

## Archive contents

Each `zips/JVASP-<n>.zip` holds the raw Wannier90 and VASP files:

```
INCAR  KPOINTS  KPOINTS.bands  OSZICAR  OUTCAR  POSCAR  JVASP-<n>.json
JVASP-<n>.png  wannier90.eig  wannier90.win  wannier90.wout  wannier90_hr.dat
```

`wannier90_hr.dat` is the standard Wannier90 real-space Hamiltonian: a comment
line, `num_wann`, `nrpts`, the R-point degeneracies (15 per line), then
`nrpts × num_wann²` rows of `R1 R2 R3 m n Re(H) Im(H)` in eV.

## Why the archives are not unpacked here

Unpacking just the `wannier90_hr.dat` files needs **69.5 GB**; the whole set is
larger still. This session's disk is 30 GB. The zips are the format figshare
ships and are byte-for-byte intact, so nothing is lost. To unpack:

```bash
python3 extract_hr_dat.py --out /big/disk/wtb                 # all 1,772
python3 extract_hr_dat.py --out ./sub --jid JVASP-5 JVASP-8   # selected ids
python3 extract_hr_dat.py --out ./sub --formula FeSe          # by formula
```

`jarvis-tools` reads these directly:

```python
from jarvis.io.wannier.outputs import WannierHam
WannierHam("JVASP-5_hr.dat").get_bandstructure_plot()
```

## `jarvis_wtb_manifest.csv` columns

`jid`, `formula`, `spg_number`, `spg_symbol`, `dimensionality` (3D/2D), `icsd`,
`mp_id` (Materials Project id, recovered from the calculation path),
`num_wann`, `num_bands`, `optb88vdw_bandgap`, `zip`, `zip_bytes`,
`n_files_in_zip`, `has_hr_dat`, `hr_dat_name`, `hr_dat_bytes`, `has_win`,
`has_wout`, `has_eig`, `has_poscar`, `figshare_article`, `error`,
`metadata_source`.
