# Meteorological Forcing — Post-Refactor Roadmap

> **Provenance.** Reconstructed 2026-08-23 from the Obsidian vault mirror
> *"Met Forcing Post-Refactor Roadmap"* (last synced 2026-07-30) + git history +
> merged code. The original `.plans/met_forcing_roadmap.md` was never tracked in
> git (unlike `met_forcing_refactor.md` / `met_forcing_docs_plan.md`, which were
> recovered verbatim from a blob), so the vault mirror is the only source; this is
> a reconstruction, not a recovered original. Every merged-status / PR claim below
> was verified against the geoclaw merge history: S1 = PR #715, S2 = PR #716,
> OWI/ASCII off-by-one = PR #717, Decision 1→B (met package) = PR #724 /
> commit `f39c0937` — all merged to `clawpack/master` by 2026-07-30. Obsidian
> `[[wikilinks]]` in the mirror have been rendered as plain text and `.plans/`
> paths repointed to `dev/design/`. Lines marked `VERIFY:` are not confirmed
> against code.

*Status: living document. Companion to the completed `met_forcing_refactor.md`
(object-model + Fortran rename) and `met_forcing_docs_plan.md` (documentation).*

**Progress (2026-07-30):** S1 (Fortran `met_*` de-prefix, PR #715), S2 (Python
OWI WIN/PRE reader/writer, PR #716) + the OWI/ASCII off-by-one fix (PR #717),
**Decision 1 → B** (the `clawpack.geoclaw.met` package + `surge` deprecation shim +
`MetData`/`rundata.met_data` aliases), and the **naming/docs conformance**
follow-up (geoclaw examples/tests + clawutil `met_data` + Sphinx docs) are all
merged to `clawpack/master`. **Next up: S3** (track-reader harvest), then the
**M1 parametric generator** — now set up as a tracked project, currently pinned.

**Progress (2026-09-01):** **S3a** (IBTrACS reader type/API hygiene) is done — see
below. It is the first of a short series driven by a downstream need to run a
45-year North Atlantic historical ensemble (IBTrACS v4 + HURDAT2, 1980–2025) through
GeoClaw, which surfaced that **the historical readers cannot currently produce a
runnable storm file at all**: `read_hurdat`/`read_ibtracs` mark missing radii `-1`
while `write_geoclaw` tests only `np.isnan`, so `fill_dict` never fires and `-1` is
written verbatim (a regression against v5.9.0, whose writer checked `== -1`). That
is the next item, as an S4 prerequisite; then multi-storm ingestion (**S5**, new —
both archives ship one file per basin and `read_hurdat` is single-storm-only), then
S3 proper and S4.

## Context

The met-forcing refactor is merged to `clawpack/master`. The object model is now
met-generic (`Track` / `StormTrack` / `ParametricMetForcing` /
`GriddedMetForcing`), the Fortran modules are renamed
(`met_forcing_module` / `parametric_met_forcing_module` /
`gridded_met_forcing_module`), and the Sphinx/tutorial documentation is done.

This roadmap triages the *next* wave of work: finishing the "surge → met"
naming transition, trimming verbose Fortran names, adding meteorological-format
support (OWI, richer track readers), and a Python wind-field generator
(the parametric family now; the USACE/Chow solver later).

**Guiding principle (from the 2026-07 audit):** GeoClaw's recurring failure mode
is *forking a parallel implementation and letting it drift*. Every item below
favors **extending the one merged object model** over introducing a second one.

**Kernel-territory note:** the wind/pressure models (M1's parametric family —
Holland, CLE, … — and the L2 Chow PBL PDE solver) are numerical kernel code. Any
wind model must be verified with a real check (Python↔Fortran field equivalence,
gradient-wind far-field balance, stationary-storm symmetry, or a match to
published Chow/SLOSH output), not a plausible-looking run.

---

## Decision 1 — Python "surge" → "met" naming  ✅ DONE (Option B, merged 2026-07-30)

**Option B shipped:** `met/` is now the real package (implementation moved from
`surge/`); `surge/` is a thin `DeprecationWarning` shim; `MetData = SurgeData` in
`clawpack.geoclaw.data`; `rundata.met_data` aliases `surge_data` (clawutil);
examples/tests and the Sphinx docs conform to the `met` name. The `surge.data`
wire filename and the `Storm` class name are unchanged. The full rename (Option C
— `SurgeData`→`MetData`, `surge.data`→`met` wire, dir/label renames) remains the
eventual major-version end state, tracked as **L1**. Option table kept below for
history.

Contract surface a rename touches: package path `clawpack.geoclaw.surge`
(imported by community `setrun.py`); the `Storm` class; `SurgeData` +
`rundata.surge_data` (the attr is registered in **clawutil**, a separate repo);
the `surge.data` filename (**cross-language contract** with `met_forcing_module`);
and Sphinx labels `surgedata` / `quick_surge` / `setrun_surge` (doc repo).

| Option | End state | Consequences | Effort |
|---|---|---|---|
| **A. Docstrings/comments only** | `surge/` stays; wording reframed | Zero API/wire/build/downstream churn; does **not** fix the core inconsistency (met objects inside a `surge` package). | S |
| **B. Add `met` aliases, deprecate `surge` (RECOMMENDED)** | Canonical `clawpack.geoclaw.met` re-exports the model; `surge` becomes a `DeprecationWarning` shim; optional `MetData = SurgeData` + accept `rundata.met_data` | New users get clean naming; **no forced downstream break**; `surge.data` wire + Fortran untouched; clawutil change only if `met_data` attr wanted. Two names coexist for a deprecation cycle. | M |
| **C. Full rename to `met`** | `surge`→`met`, `SurgeData`→`MetData`, `surge.data`→`met` wire, `surge_data`→`met_data` | Cleanest, fully coherent; **breaks every community `setrun.py`**; needs coordinated geoclaw + clawutil + doc release + Fortran wire change. | L (multi-repo) |

**Recommendation: B now, C as the eventual major-version end state.** B resolves
the packaging inconsistency without breaking the storm-surge user base, and sets
up C as the deprecation-cycle endpoint. Keep the class name `Storm` (it *is* a
storm). Keep the `surge.data` filename for now (cross-language + cross-repo cost
not yet justified). Internal docstring/comment reframing (A) is absorbed into B.

---

## Decision 2 — Fortran `met_*` de-prefix  ✅ DONE

The only redundant prefixing was four public names + one subroutine inside
`gridded_met_forcing_module.f90` (a module whose name already says `met_forcing`),
none referenced by any external `use` site: `has_met_crop`→`has_crop`,
`met_crop_extent`→`crop_extent`, `met_x_shift`/`met_y_shift`→`x_shift`/`y_shift`
(now an exact match to the topo/dtopo names they parallel), subroutine
`met_crop_indices`→`crop_indices`. Stale pre-rename `.o` files removed. The module
*names* themselves are unchanged.

*Verified behavior-neutral:* the gridded netCDF + Holland aux regression tests
pass against existing goldens (`rtol=1e-14/atol=1e-8`, no `GEOCLAW_REGEN`).
Delivered on branch `met-deprefix-gridded`.

---

## Work items

### SHORT TERM

- **S1. Fortran `met_*` de-prefix** — ✅ done (Decision 2).

- **S2. OWI (Oceanweather WIN/PRE) reader/writer — ✅ DONE (PR #716).**
  Rewrote the orphaned `surge/data_storms.py` into a clean field reader/writer
  (`read_owi` / `write_owi` / `read_owi_start_time` + the `OWIData` dataclass; SI
  in memory, mbar/f10.4 on disk), registered it in `meson.build`, and added thin
  `GriddedMetForcing.from_owi` / `to_owi` helpers. Tests: object round-trip +
  real-`isaac` sample reads in `test_met_objects.py`, plus a format-1 gridded
  regression case in `tests/regression/met_forcing/`. Documented the f10.4
  precision limit.
  - **OWI/ASCII coordinate off-by-one — ✅ FIXED (PR #717).** The Fortran OWI
    reader placed grid point *i* at `sw + i*d` instead of `sw + (i-1)*d`, shifting
    the whole grid one cell NE and disagreeing with the NetCDF path. Proven with a
    new OWI-vs-NetCDF equivalence regression test (failed pre-fix by ~1535 Pa),
    fixed, and the `aux_owi` + isaac `owi_ascii` gauge goldens regenerated. Note:
    `test_isaac.py`'s OWI coordinate reconstruction was updated to match.

- **S3. Harvest low-risk track-reader improvements** into `met/track.py`
  (do **not** adopt a parallel `CycloneDataset` object model): quadrant wind-radii
  (`r34/r50/r64` × NE/SE/SW/NW), and `Basin`/`StormStatus` enums + standardization
  maps. Port ideas, re-verify against real files, keep `StormTrack` as the single
  home.
  - *The original "…and an IBTrACS netCDF path if the current reader is CSV-only"
    clause was **stale** and has been struck: `read_ibtracs` has always been an
    xarray/netCDF reader. Confirmed against the real `IBTrACS.NA.v04r01.nc`
    (2,298 storms; dims `storm` / `date_time` / `quadrant`), which does carry the
    `quadrant` dimension this item needs — only the bundled
    `tests/data/storm/ibtracs.nc` fixture is a stripped single-storm slice without
    it, so offline coverage of quadrant radii will need a new fixture or a
    `remote` test.*

- **S3a. IBTrACS reader type/API hygiene — ✅ DONE (PR #NNN, verified against
  `7ef7a816`).** Prerequisite for S4: the fill path could not run at all, because
  `read_ibtracs` left `self.t` as an `xarray.DataArray`, so `storm.t[n]` was a 0-d
  DataArray and `fill_rad_w_other_source`'s `interp` raised
  `ValueError: Could not convert object to NumPy datetime`.
  - `self.t` is now a `datetime64[s]` ndarray (matching every other reader),
    `classification` is decoded from bytes (`b'TD'` → `'TD'`), `ds.dims.keys()` →
    `ds.sizes` (xarray `FutureWarning`; `Dataset.dims` becomes a set of names), and
    the storm-dimension squeeze is explicit so a single-observation storm no longer
    has `date_time` collapsed along with it.
  - Truncating to seconds also drops the sub-second roundoff IBTrACS decodes
    (`…T06:00:00.000039936`). **This moved two rows of the written storm file back
    into exact agreement with the committed `ibtracs_geoclaw.txt` baseline**
    (`-3.60000003e+03` → `-3.60000000e+03` and `7.19999997e+03` → `7.20000000e+03`
    at data rows 103 and 105); the v5.9.0-era golden holds the exact values, so the
    roundoff was a post-v5.9.0 regression rather than the baseline being wrong.
    Verification is that agreement plus four contract tests in
    `tests/test_met_objects.py` — types, whole-second times, a read under
    `FutureWarning`-as-error, and `fill_rad_w_other_source` accepting a scalar
    *and* a 0-d DataArray. All four fail on `7ef7a816`.
  - Only golden regenerated: `characterization/read_ibtracs.json`, which was
    storing the *xarray repr string* for `t` (including `Size: 952B`, coords and
    attrs) and `"b'TD'"` for `classification` — i.e. exactly the type leak being
    fixed, and a snapshot that would have broken on any xarray repr change.
  - `test_storm.py::test_storm_io[ibtracs]` stays `xfail`; the `-1`/NaN sentinel
    bug that causes it is untouched here and is the subject of the next item.

- **S4. Parametric geometry reconstruction (missing-data filling).** Restore and
  modernize the old code's "basic reconstruction from statistical and physical
  data" for the recurring pain of **missing storm geometry** in commonly-available
  best tracks (ATCF/HURDAT/IBTrACS): no radius of maximum winds (RMW), no
  outer/last-closed-isobar radius (ROCI), or only one of max-wind-speed /
  central-pressure. Today the readers mark these `-1`/`NaN`, warn, and refuse to
  fill physically — the only automatic default is a flat constant
  (`storm_radius = 500 km`, `met/parametric.py:150`). This item adds **opt-in,
  explicit** estimator functions with the **existing `fill_dict` signature**
  `fn(t, storm) -> value (SI)`, so they drop straight into
  `ParametricMetForcing.write_geoclaw(fill_dict=...)` and compose with the existing
  `fill_rad_w_other_source` time-interpolator (`met/track.py:1035`). **No new
  object model, no Fortran, no wire change** (extend, don't fork). The default
  write path stays warn/skip — nothing is silently reconstructed. New module
  `met/reconstruction.py`; optional `default_fill_dict(track, models={...})`
  convenience assembler. Estimator families (each cites its paper; **maintainer
  confirms coefficients & provenance**):
  - **Wind–pressure relationship** — fill `max_wind_speed` ↔ `central_pressure`:
    Knaff & Zehr (2007) (recommended; operational standard), alt Atkinson &
    Holliday (1977). *WPR is the essential companion not in the reference note —
    flagged for sign-off.*
  - **RMW from intensity + latitude** — fill `max_wind_radius`: Willoughby et al.
    (2006) (recommended), alts Kimball & Mulekar (2004) / Carrasco et al. (2014),
    optional Vickery & Wadhera (2008) for Holland-B-consistent RMW.
  - **Outer/storm radius (ROCI) from environment** — replace the 500 km constant:
    Chavas et al. (2016), with a **climatological fallback** so it degrades to a
    physically-motivated default when environment (SST/outflow) is unknown.

  *Not numerical-kernel code* — empirical regressions from cited papers, not
  solver/Riemann/limiter math, so it does not hit the kernel STOP rule; coefficient
  choices remain the maintainer's. **Verification = reproduce the published
  relationship**, not a manufactured PDE: match published reference values
  (Willoughby RMW at documented `Vmax`/lat; Knaff–Zehr WPR table) to tolerance;
  physical invariants (RMW↓ with `Vmax`, deeper `dp`⇒higher `Vmax`, positivity/SI,
  latitude direction); distributional check vs observed IBTrACS/EBTRK; and a
  sparse-ATCF (RMW/ROCI stripped) → `fill_dict` → `write_geoclaw` round-trip
  yielding a runnable 7-column storm file with all fields finite
  (`tests/test_met_objects.py`). Opt-in ⇒ existing aux + bowl-slosh goldens stay
  byte-identical (no `GEOCLAW_REGEN`). Sequence **after S3** (can also draw on the
  harvested quadrant radii / basin) and as a **companion feeding M1** (the
  generator consumes complete parameter sets). Reference note (TC-radius section):
  *Rework TC Geometry Estimates* (vault project note).

### MEDIUM TERM

- **M1. Parametric-family offline wind-field generator (Python).** Turn a
  `StormTrack` + parameterization into gridded wind/pressure fields in Python and
  write netCDF (and/or OWI) that `GriddedMetForcing` already consumes — no
  live-solver changes, no new Fortran family. Fills the `met/storm.py`
  `holland_1980`/`holland_2010`/`cle_2015`/`construct_fields` stubs; adds
  `ParametricMetForcing.to_gridded(...)`. Phasing (G1–G5): **Holland-1980**
  (matches the existing `aux_holland80.txt` golden) → Holland-2010/2008 →
  **CLE15** (ODE solve; heaviest kernel; cross-check vs the Fortran CLE *and*
  `EmilioEchevarria/tcwindfields` CLE15, MIT — reference only, not vendored).
  **Verification is the point** (see the discipline note below). *Kernel territory
  — the model numerics are implemented/verified by the maintainer.* The **Chow**
  USACE PBL wind solver is the harder sub-case, deferred to L2 (nested
  telescoping-grid iterative solver: advection + diffusion + surface drag +
  pressure-gradient + implicit Coriolis; verify vs gradient-wind far-field,
  stationary-storm symmetry, published Chow output).
  - **Status: designed, tracked, pinned (not started).** Durable design +
    Fortran-fidelity spec + verification plan: `dev/design/met_wind_field_generator.md`.
    Tracked as an Obsidian project (currently pinned):
    *Implement parametric TC wind-field generator* — this roadmap remains the
    technical index; the project note tracks resumption.
  - **Sibling: see M3** (analytic/idealized forcing) — both add new parametric
    subtypes via the same registration path and share the Python↔Fortran
    field-equivalence verification harness.

- **M2. `StormTrack.to_xarray()` / CF netCDF track I/O** — a single canonical
  serialization for tracks (idea harvested from the artifact's CF schema), giving
  multi-storm datasets a home *inside* `StormTrack` rather than a rival class.
  Delivered alongside M1's gridding layer (recommend the same CF schema as
  `create_era5_storm_file`/`create_nws13_storm_file`).

- **M3. Analytic / idealized meteorological forcing (new parametric subtypes).**
  Closed-form pressure + wind evaluated **per cell, live in Fortran** — new
  `ParametricMetForcing` subtypes on a **generic `Track`** (a moving front), not
  gridded and not a new family. Primary motivation: meteotsunami modeling
  (issue #694; the hydrocoast/Miyashita `meteotsunami` fork) generalized to a
  reusable idealized-forcing capability. Fields: closed-form pressure with
  along-track shape Φ (sine/hat/morlet/gaussian/jump) × cross-track apodization
  Ψ⊥ × temporal ramp R(t); **wind first-class** (`none`/`geostrophic`/`gradient`).
  Integration = the **same subtype path as M1's models**: register in
  `data.py` `forcing_subtype_registry` (`("parametric", <new int>)`) +
  `forcing_spec_type`/`set_storm` case → new `set_analytic_fields` in
  `parametric_met_forcing_module.f90`; mimics the `test_topography` analytic
  bathymetry pattern (`topo_module.f90` ~L375). Source coupling (src2
  pressure-gradient + wind-drag) is **untouched** — kernel-territory note applies.
  **Verification is a strength:** analytic ocean response — **Proudman resonance**
  (moving pressure jump over constant depth) and **Greenspan edge waves** — gives
  end-to-end *manufactured-solution* tests, and the Python↔Fortran
  field-equivalence harness (shared with M1) covers the field evaluator itself.
  Deliverable: a Proudman moving-jump canonical example + regression/CI. Open
  design decisions (see the spec): subtypes-vs-sibling-family, wind-coupling
  formulation, the moving-front AMR criterion (**seam 2**), and subtype naming +
  legacy signed-int alias on the `surge.data` wire. Full spec:
  *Analytic Met Forcing Plan* (vault note); tracked project (status **active**,
  GeoClaw-Users'-Workshop-2026 target):
  *Implement analytical meteorological forcing for GeoClaw*.
  **Sibling to M1** — both are new parametric subtypes sharing the
  subtype-registration path and the verification harness.

### LONG TERM

- **L1. Full `met` rename** (Decision 1, Option C) at a GeoClaw major version once
  B's deprecation cycle has run — coordinated with clawutil + doc repos.
- **L2. Chow PBL wind solver** — first as the Python generator's hardest sub-case
  (offline, per M1), then as a native Fortran forcing family *only if* verification
  shows the live-solver cost is warranted. Explicit kernel work.
- **L3. De-duplicate track readers** — retire scattered prototype readers once
  S3/M2 land.
- **L4. TC precipitation estimation (placeholder — separate scoping).** The second
  half of the *Rework TC Geometry Estimates* reference note: estimate rainfall
  from a parametric storm — from vertical velocity (Emanuel et al. 2006, Zhu et al.
  2013, Lu et al. 2018, Gori et al. 2022) or from the empirical TC-precip↔radii
  relationship (Tuleya et al. 2007, Mudd et al. 2017). A **distinct, larger
  capability** than S4's geometry filling: it needs a rainfall forcing field plus
  overland-flow coupling that does not exist yet. Listed here only so the reference
  note's precip thread is not lost; scope separately before committing.

---

## Suggested sequencing

1. **S1** — ✅ done (PR #715).
2. **S2** (OWI reader/writer) + off-by-one fix — ✅ done (PR #716 / #717).
3. **Decision 1 → B** (met package + `surge` shim + `MetData`/`met_data`) +
   naming/docs conformance (examples/tests + clawutil + Sphinx) — ✅ done.
4. **S3** (track-reader harvest) — **NEXT**.
5. **S4 — parametric geometry reconstruction** (opt-in missing-data filling).
   After S3 (can use harvested quadrant radii / basin) and a **companion to M1**
   (feeds it complete parameter sets); pure Python, opt-in, low-risk.
6. **M1 — parametric-family Python generator** (the substantive new capability;
   designed, tracked as a project, currently pinned). Recommended **after S3** so
   it is authored against the richest `StormTrack`, but it can start independently
   if appetite comes first. Chow is the harder sub-case → **L2**.
7. **M3 — analytic / idealized forcing** — the currently **active** analytic thread
   (Workshop-2026 target), **independent of S3 and M1** (it uses a generic `Track`,
   not the enriched `StormTrack`), so it can proceed now / in parallel. Sibling to
   M1; shares the subtype-registration path + verification harness.
8. **M2 / L1 / L2 / L3 / L4** as evidence and appetite allow.

## Verification discipline (all items)

Behavior-neutral changes must keep the met-forcing aux goldens and bowl-slosh
canary byte-identical (no `GEOCLAW_REGEN`). New capabilities add their own
regression cases under `tests/regression/met_forcing/` and/or object tests in
`tests/test_met_objects.py`. Cross-repo changes (clawutil / riemann / doc) are
called out explicitly and merged in lockstep.

**Parametric generator = verification engine (M1).** A faithful Python parametric
generator must reproduce the Fortran forcing-aux goldens (e.g. `aux_holland80.txt`)
on the same grid/frames/storm — a Python↔Fortran equivalence test that makes the
two implementations *mutually verifying* and is the path to retiring the hand-rolled
analytic-vortex fixtures (`_netcdf_forcing`/`_owi_forcing`). **Anti-circularity
rule:** the equivalence golden always comes from the independent Fortran path, never
regenerated from the generator, so it never grades its own output.

**Analytic forcing = end-to-end manufactured solutions (M3).** The analytic
subtypes add *closed-form ocean-response* checks — Proudman resonance (known
amplification as `c → √(gh)`) and Greenspan edge-wave amplification — that verify
the full forcing→SWE coupling, not just the field. Their per-cell field evaluator
reuses M1's Python↔Fortran field-equivalence harness (same anti-circularity rule).
Together M1 and M3 give the met-forcing subsystem both field-level and
system-level manufactured verification.
