# Priority 2 — flat-band candidate lists

Regnault, N. *et al.* *Catalogue of flat-band stoichiometric materials.*
**Nature 603, 824–828 (2022).** <https://doi.org/10.1038/s41586-022-04519-1>
(see also the [Author Correction](https://doi.org/10.1038/s41586-022-05065-6)).

Two independent sources are captured here, because they do not agree exactly:

1. **The paper's Supplementary Information** — the list of record, parsed out of
   Supplementary Tables XI and XII.
2. **The Materials Flatband Database**,
   <https://www.topologicalquantumchemistry.fr/flatbands/> — the live database,
   which is a slightly later revision of the same work.

## A. From the Supplementary Information (the paper of record)

`SI_41586_2022_4519_MOESM1_ESM.pdf` (25 MB, 166 pp.) is the Supplementary
Information; `SI_41586_2022_4519_MOESM2_ESM.pdf` is the peer-review file.
`SI1.txt` / `SI2.txt` are text extractions. `parse_si_tables.py` produces:

| File | Rows | Contents |
| --- | --- | --- |
| `tableXI_curated_flatband_materials.csv` | **2,379** | Supplementary Table XI — the curated flat-band materials |
| `tableXII_curated_flat_atomic_band_materials.csv` | **1,102** | Supplementary Table XII — curated flat *atomic* band materials |
| `top_candidates_best_flatbands.csv` | **339** | the subset of Table XI flagged as best candidates |

Columns: `num`, `formula`, `space_group_symbol`, `space_group_number`,
`topology_at_Ef` (TI / SM / Triv.), `is_top_candidate`, `best_figure`,
`best_figure_number`, `sublattices` (k = Kagome, p = pyrochlore, b = bipartite,
s = split, l = Lieb), `magnetic`, `superconductor`, `high_quality`, `n_icsd`,
`icsd_ids`, `icsd_extra_sublattice_labels`.

### Two numbers that do not match the paper's prose

* **Top candidates: 339, not 345.** Appendix H 3 states "the 345 best
  candidates (corresponding to 949 ICSD entries)", shown in Figs 26–69. The
  printed Table XI carries exactly **339** `Fig. NN` flags, spanning exactly
  those 44 figures. This is a discrepancy inside the published table, not a
  parsing artefact: a direct `pdftotext` scan of pages 48–94 also finds 339.
  The site's own `bestflatband` flag is carried in the site-derived files below
  as a second opinion.
* **ICSD entries: 6,262, not 6,338.** The prose says the 2,379 materials
  correspond to 6,338 ICSD entries. Each row's ICSD cell is a multi-line block
  vertically centred on its row, and where two tall blocks abut, the row
  boundary is genuinely ambiguous in the PDF. `finalize_flatbands.py` therefore
  re-derives the per-material ICSD sets from the live database and writes
  `tableXI_curated_flatband_materials_verified.csv`, which adds
  `material_id`, `icsd_ids_verified`, `n_icsd_verified`, `site_formula`,
  `site_space_group_number` and `verification`. Use the verified columns when
  the ICSD grouping matters; use `icsd_ids` if you want strictly what the PDF
  shows.

### How the tables were parsed

`pdftotext -layout` cannot recover these tables: the column x-positions are
re-typeset on every page, and each row's ICSD cell is a multi-line block
centred on the row rather than aligned to its baseline. The parser therefore
works from `-bbox-layout` word boxes — cells are found by intra-row x-gap
clustering, and ICSD text lines are clustered by vertical gap and assigned to
the nearest row anchor. Row numbering is asserted to be exactly 1..N, so a
mis-parse cannot pass silently.

## B. From the Materials Flatband Database (live site)

The site has no bulk export, but `showmaterial.php?ICSD=<id>` returns a JSON
record per ICSD entry. `fetch_flatband_records.py` walks every ICSD in the
curated lists — seeded from the SI parse and the site's own `Curated`/`Atomic`
search listings, then closed under each record's `otherICSDs` field so no entry
of a curated material is missed.

| File | Contents |
| --- | --- |
| `flatband_icsd_metadata.csv` | one row per ICSD entry: formula, space group, point group, topology with and without SOC, curated / best / atomic / high-quality flags, sublattice counts, cell parameters, Materials Project link, CIF path |
| `site_materials_curated_flatband.csv` | curated flat-band materials, grouped |
| `site_materials_best_flatband.csv` | the site's best flat-band materials |
| `site_materials_atomic_flatband.csv` | curated flat *atomic* band materials |
| `records/` | the trimmed JSON records (gitignored; band-structure plot payloads stripped) |
| `poscar/` | POSCAR for every ICSD entry |
| `site_curated.html`, `site_atomic.html` | the raw search listings |
| `seed_icsds.txt` | the seed ICSD list |

Materials are grouped by the site's own `otherICSDs` field, so the material
level does not depend on the PDF parse at all.

Counts differ from the paper: the site returns 2,055 curated and 1,081 atomic
materials against the SI's 2,379 and 1,102. The site is a later revision; both
sets are provided rather than silently reconciled.

## Not used

ICSD itself (<https://icsd.fiz-karlsruhe.de>) is subscription-only, so no ICSD
query was made. Everything here came from the open Nature supplementary files
and the public flat-band database.
