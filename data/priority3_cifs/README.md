# Priority 3 — CIFs for the Priority-2 candidates

**8,168 CIFs in `cif/<ICSD>.cif`, 37 MB. Every one contains an explicit
`_symmetry_equiv_pos_as_xyz` loop — not just a space-group number.**

This covers **every ICSD entry of every Priority-2 candidate**, with nothing
missing:

| Candidate list | ICSD entries | with CIF | with explicit symmetry ops |
| --- | ---: | ---: | ---: |
| Curated flat-band (2,379 materials) | 6,338 | **6,338** | **6,338** |
| — of which best candidates (345 materials) | 949 | **949** | **949** |
| Curated flat *atomic* band (1,102 materials) | 1,830 | **1,830** | **1,830** |
| **Total distinct ICSD entries** | **8,168** | **8,168** | **8,168** |

## Where these came from

Not from COD, Materials Project or OQMD. The Materials Flatband Database
(<https://www.topologicalquantumchemistry.fr/flatbands/>) serves the original
ICSD CIF for each entry through its `showmaterial.php?ICSD=<id>` endpoint, so
the structures are matched to the candidate lists **by ICSD id exactly** rather
than by a formula/space-group lookup that could land on the wrong polymorph.
They were harvested alongside the Priority-2 metadata by
`../priority2_flatbands/fetch_flatband_records.py`.

ICSD itself is subscription-only and was never queried.

### Licensing — read before redistributing

Each file begins with:

```
#(C) 2017 by FIZ Karlsruhe - Leibniz Institute for Information Infrastructure.  All rights reserved.
```

These are ICSD records made publicly readable by the flat-band database. That
copyright notice is FIZ Karlsruhe's, not a CC licence — fine to use for your
own research, but check with FIZ before redistributing the set. If you need
structures under an open licence instead, `../priority2_flatbands/flatband_icsd_metadata.csv`
carries a `materials_project` link for all 8,168 entries (a free MP API key
gets you MP's CIFs), and COD is an option where a matching entry exists.

## Format

Standard ICSD CIF: cell parameters, `_symmetry_space_group_name_H-M`,
`_symmetry_Int_Tables_number`, a `_symmetry_equiv_pos_as_xyz` loop, and the
asymmetric-unit `_atom_site_*` loop with Wyckoff symbols and occupancies.

The symmetry-operation counts are the true space-group orders — 2, 4, 8, 12,
16, 24, 48, 96, 192 and so on, up to 192 for the `Fd-3m` entries.
`../priority2_flatbands/flatband_icsd_metadata.csv` records the count per entry
in `cif_symmetry_ops`, alongside `cif_file`, formula, space group, the
curated / best / atomic flags and the cell parameters.

## Related

`../priority2_flatbands/poscar/<ICSD>.POSCAR` holds a VASP POSCAR for the same
8,168 entries (P1-expanded coordinates, no symmetry operations — use the CIFs
when you need the symmetry).
