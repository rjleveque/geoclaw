# Parametric TC Wind-Field Generator — Design & Verification Plan

> **Provenance — PARTIAL RECONSTRUCTION.** Rebuilt 2026-08-23. The original
> `.plans/met_wind_field_generator.md` was **never tracked in git and has no full
> vault mirror** — only a summary survives in the vault project note *"Implement
> parametric TC wind-field generator"* (last reviewed 2026-07-30). This file is a
> skeleton assembled from that project note, the M1 entry in
> `dev/design/met_forcing_roadmap.md`, and the current code stubs; it is **not the
> original document** and is known to be incomplete. Sections reconstructed from
> the summary are marked; anything not confirmed against code is `VERIFY:`.
> Code anchors below were confirmed present at reconstruction time. **Kernel
> territory:** the wind/pressure model numerics are implemented and verified by
> the maintainer — this doc scopes and plans, it does not supply kernel code.

## Status (from the project note, 2026-07-30)

Design + verification plan complete; the enabling surge→met refactor is merged to
`clawpack/master`. Implementation is **pinned (not started)**. Decision made:
build a **native** Python parametric wind-field generator in
`clawpack.geoclaw.met` (Option A), verified against the Fortran `aux` fields — use
`EmilioEchevarria/tcwindfields` (MIT) only as an independent cross-check, not
vendored. Immediate next decision: when to start phase G1 (Holland-1980 kernel)
and lock the first-model scope.

## Goal

Turn a `StormTrack` + a parametric model into gridded wind/pressure fields **in
Python**, and write netCDF (and/or OWI) that `GriddedMetForcing` already consumes
— no live-solver changes, no new Fortran forcing family. This makes the parametric
family available offline and, crucially, makes the Python and Fortran evaluators
**mutually verifying**.

## Code anchors (confirmed 2026-08-23)

Python stubs to fill (module-level functions in `met/storm.py`):
- `construct_fields(storm, r, t, model="holland_1980")` — `met/storm.py:447`
- `holland_1980(storm, r, t)` — `met/storm.py:457`
- `holland_2010(storm, r, t)` — `met/storm.py:463`
- `cle_2015(storm, r, t)` — `met/storm.py:469`

Proposed new API (does **not** exist yet — `VERIFY:` when built):
- `ParametricMetForcing.to_gridded(...)` on `met/parametric.py` (no such method at
  reconstruction time; `parametric.py` currently exposes `write_geoclaw` at :118
  and the 500 km `storm_radius` default fill at :150).

Fortran fidelity oracle — the per-cell field evaluators the Python generator must
reproduce, in `src/2d/shallow/surge/parametric_met_forcing_module.f90`:
- `set_holland_1980_fields` (:607), `set_holland_2008_fields` (:679),
  `set_holland_2010_fields` (:754), `set_CLE_fields` (:843), `set_SLOSH_fields`
  (:1208). *(Note: `met/storm.py` has no `holland_2008` stub, though the Fortran
  side implements `set_holland_2008_fields` — the Python phasing skips 2008 or
  reaches it via 2010; `VERIFY:` the intended first-pass scope.)*

## Phasing (G1–G5, from the roadmap M1 entry)

1. **G1 — Holland-1980.** First model; must match the existing
   `tests/regression/met_forcing/regression_data/aux_holland80.txt` golden.
2. **G2/G3 — Holland-2010 / Holland-2008.**
3. **G4/G5 — CLE15** (Chavas–Lin–Emanuel 2015): an ODE solve, the heaviest kernel;
   cross-check against the Fortran `set_CLE_fields` *and* the MIT `tcwindfields`
   CLE15 (reference only).

The **Chow** USACE PBL wind solver is the harder sub-case and is **deferred to
roadmap L2** (nested telescoping-grid iterative solver; verify vs gradient-wind
far-field, stationary-storm symmetry, published Chow output).

## Verification plan (the point of the exercise)

- **Python↔Fortran field equivalence** is the primary check: the Python generator
  must reproduce the Fortran forcing-aux goldens (e.g. `aux_holland80.txt`) on the
  same grid / frames / storm. Proposed tolerance (from the project note):
  `rtol ≈ 1e-9 / atol ≈ 1e-6` cross-language. `VERIFY:` tolerance against a real run.
- **Anti-circularity rule:** the equivalence golden always comes from the
  independent Fortran path, never regenerated from the Python generator — the
  generator never grades its own output.
- **Physical invariants:** gradient-wind far-field balance; stationary-storm
  symmetry; positivity / SI units.
- This harness is **shared with roadmap M3** (analytic/idealized forcing), which
  additionally supplies end-to-end manufactured ocean-response checks (Proudman
  resonance, Greenspan edge waves).

## Open questions (from the project note)

- First-pass model scope: Holland-1980 → CLE15 only, or also Holland-2010/2008 in
  between? (CLE15 is the heaviest kernel.)
- Is the Chow PBL solver part of this project or a separate follow-on (roadmap L2)?
- `StormTrack.to_xarray()` CF schema — match `create_era5_storm_file` /
  `create_nws13_storm_file`? (leaning yes)
- Gate the generator behind a soft `xarray` import? (xarray is already a
  netcdf/test dep.)
- Python↔Fortran equivalence tolerance target (proposed `rtol≈1e-9 / atol≈1e-6`).
- Timeline / appetite / who owns the kernel numerics.

## References

- Met-forcing roadmap (M1 / M2 / L2): `dev/design/met_forcing_roadmap.md`
- Fortran source of truth:
  `src/2d/shallow/surge/parametric_met_forcing_module.f90`,
  `met_forcing_module.f90`, `geoclaw_module.f90`
- Python home / stubs to fill: `src/python/geoclaw/met/storm.py`
  (`holland_1980` / `holland_2010` / `cle_2015` / `construct_fields`),
  `met/parametric.py`
- Verification oracle: `tests/regression/met_forcing/`
  (`regression_data/aux_holland80.txt`)
- External reference impl (MIT, cross-check only):
  https://github.com/EmilioEchevarria/tcwindfields
- Papers: Holland (1980); Holland, Belanger & Fritz (2010); Powell/Holland (2008);
  Chavas, Lin & Emanuel (2015); Emanuel & Rotunno (2011)
- Vault project note (tracks resumption): *Implement parametric TC wind-field
  generator*; sibling analytic thread: *Implement analytical meteorological forcing
  for GeoClaw*.
