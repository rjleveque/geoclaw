# File-Based Friction Field Input for GeoClaw — Implementation Plan

> **Provenance.** Relocated 2026-08-23 from the gitignored `.plans/friction-file-plan.md`
> (developer scratch) into tracked `dev/design/`, unchanged in content. It was
> never lost — this move just puts the durable design doc under version control
> alongside the met-forcing docs. Vault project notes previously referenced it as
> `.plans/friction_file_plan.md` (underscores); the canonical path is now
> `dev/design/friction-file-plan.md`. Implementation status below is unchanged and
> not re-verified against current `master` — heed its own "Verify first" section.

Handoff for Claude Code. Repo: `clawpack/geoclaw` (verified against master @ `1af3d3e`, 2026-08-01). **Double-check the "Verify first" section against current master before writing code**, since call sites and list formats may have moved.

## Goal

Implement the long-dormant file-based variable-friction path: read a gridded Manning's $n$ field from one or more files covering limited coastal regions, fill it into `aux(friction_index)`, and fall back to the existing depth/constant Manning elsewhere. Coverage is regional (like met forcing), not global (unlike topo), so uncovered cells must degrade gracefully to the depth-based coefficient. Also ship a Python preprocessing path that turns commonly available land-cover rasters into a friction file.

## Current state of the stub (precise)

**Fortran — `src/2d/shallow/friction_module.f90`**
- Documented precedence, top of file: `file > region > depth > constant`.
- `friction_file_type` derived type exists: `num_cells(2)`, `lower(2)`, `upper(2)`, `dx(2)`, `no_data_value`, `data(:,:)`.
- `setup_variable_friction` reads `friction.data`, parses `num_friction_regions` + regions (works), then reads `num_friction_files`, and for each prints `*** WARNING *** File based friction specification unimplemented.` and calls the stub reader.
  - **Bug:** inside the file loop it does `friction_files = read_friction_file(...)` (assigns scalar result to the whole array). Must be `friction_files(m) = ...`.
- `set_friction_field` applies the region path (depth-break lookup inside each rectangle), then `if (num_friction_files > 0) stop "*** ERROR *** File based friction specification not implemented."`. No interpolation of `file%data` onto aux exists.
  - **Bug/limitation:** depth-break comparison uses `aux(1,i,j) - sea_level` with an inline comment that `sea_level is not set properly here, just a test case`. Revisit while touching this routine.
- `read_friction_file` parses a header for `file_type` 2/3 (`num_cells(1)`, `num_cells(2)`, `lower(1)`, `lower(2)`, `dx(1)`, `no_data_value`) but every data-read branch is a commented-out port of the topo reader falling through to `stop`. `file_type` 1 (xyz) also stops. Header is GeoClaw topotype-3 minus `dy` (it sets `dx(2)=dx(1)`).

**Fortran — consumers / call sites**
- `set_friction_field` is called from `src/2d/shallow/setaux.f90:220` (and `src/2d/shallow/multilayer/setaux.f90:209`), *after* topo is cell-integrated into `aux(1)`. Runs at grid init and every regrid, so friction is static-filled there. **Both setaux paths need the fix.**
- `setup_variable_friction()` is called once at startup from `src/2d/shallow/amr2.f90` (and `src/2d/bouss/amr2.f90`).
- `friction_index` is consumed in `src/2d/shallow/src2.f90:101` and `src/2d/shallow/src1d.f90:63`.

**The precedence gap (important).** In `src2.f90` the friction source term does:
```fortran
if (.not. variable_friction) then
    do nman = num_manning, 1, -1                 ! depth-based
        if (aux(1,i,j) < manning_break(nman)) coeff = manning_coefficient(nman)
    enddo
else
    coeff = aux(friction_index,i,j)              ! variable: aux is the ONLY source
endif
```
When `variable_friction` is on there is **no per-cell fallback** — the coefficient comes solely from `aux(friction_index)`. But `set_friction_field` fills that slot only inside regions, so today the documented `depth > constant` levels are never realized for variable friction, and uncovered cells carry 0 (frictionless). This feature is the moment to fix the whole chain (see Task 2).

**Python — `src/python/geoclaw/data.py`, `FrictionData` (~line 1154)**
- `friction_files` attribute exists; `read()` hard-codes `self.friction_files = []  # Is not supported`.
- `write()` has a format-string bug: `self._out_file.write("'%s' %s\n " % fname)` — two specifiers, one value; raises the moment the list is non-empty. It also never writes a file_type, though Fortran reads one via `(a,i1)`.
- `friction_index` is written as `+1` for Fortran indexing (keep).

## Settled design decisions

1. **Do not fork a reader.** The stub's bespoke `read_friction_file` duplicates the topo reader, which is exactly the parallel-implementation failure mode we're avoiding elsewhere. Reuse the topo gridded-field machinery in `src/2d/shallow/topo_module.f90`: `read_topo_file` (line 746), the ASCII `case(2:3)` body, the netCDF `topo_type=4` path (`read_netcdf_descriptor`, `nf90_open` @ 959), and `cellgridintegrate` for cell-averaging. A friction file is a topo-like gridded field that lands in `aux(friction_index)` and is interpreted as Manning's $n$ instead of elevation. Prefer factoring a shared gridded-field reader over copy-paste; if that refactor is too large for one pass, at minimum make the friction file format *byte-identical* to topotype 2/3/4 so `clawpack.geoclaw.topotools` writers/readers produce and consume friction files for free.
2. **File format = topotype 2/3 (ASCII) and 4 (netCDF).** No new format. Extend the friction header to carry `dy` (drop the `dx(2)=dx(1)` assumption) so it matches topo exactly.
3. **Store resolved Manning's $n$ per cell, not land-cover class codes.** The class→$n$ lookup is a modeling choice and lives in Python preprocessing, versioned and auditable. Fortran stays a dumb reader.
4. **Regional coverage via no_data masking.** File overrides only where `data /= no_data_value`; elsewhere the cell keeps the region/depth/constant value. This is the "like met forcing" behavior.
5. **Static fill at regrid** via `set_friction_field` in `setaux`. No time dependence, simpler than met forcing (a degenerate single-time gridded field).
6. **Cell-average onto the target cell** (reuse `cellgridintegrate`), consistent with topo. Roughness is smooth relative to bathymetry, so nearest/bilinear is defensible, but cell-averaging composes with existing machinery and is the physically correct reduction of a fine $n$ field.

## Verify first (double-check before implementing)

- [ ] Confirm `src2.f90`/`src1d.f90` still have no per-cell fallback for `variable_friction` (drives Task 2 strategy).
- [ ] Confirm the current `topofiles` entry format in `setrun`/`TopographyData` (currently `[topo_type, path]`) and mirror it for `friction_files` as `[friction_type, path]`. Align the `friction.data` line format with what `read_friction_file`/`setup_variable_friction` parse.
- [ ] Confirm the `read_netcdf_descriptor` signature and the `topo.data` key=value descriptor block so the netCDF friction path can reuse it verbatim.
- [ ] Confirm both `setaux.f90` and `multilayer/setaux.f90` should receive identical treatment.
- [ ] Decide precedence-realization strategy (Task 2, A vs B).

## Tasks

### Task 1 — Fortran read path (reuse, don't fork)
- Replace the bespoke `read_friction_file` body with a call into the shared/topo gridded reader, directing output to a friction field object rather than `aux(1)`. Support ASCII types 2/3 and netCDF type 4.
- Fix the `friction_files(m) = ...` array-assignment bug in `setup_variable_friction`.
- Extend `friction_file_type` (or the reused topo field type) to carry `dy` and the registration convention (cell-centered vs corner) identically to topo.
- **Acceptance:** a topotype-3 file written by `topotools` loads into a `friction_file_type` with correct extents, `dx/dy`, and `no_data_value`; a netCDF (type 4) file loads equivalently.

### Task 2 — Precedence + interpolation in `set_friction_field`
Realize the full chain by seeding every cell, then overriding:
```text
for each cell (i,j) incl. ghosts:      ! aux(1)=topo already set by setaux
    coeff = depth_based_manning(aux(1,i,j))         ! (1) constant/depth from geoclaw_module
    if cell in friction_region r:                    ! (2) region override
        coeff = region_depth_break(r, aux(1,i,j))
    for each friction_file f covering (x,y), finest last:   ! (3) file override, highest priority
        v = cellgridintegrate/interp(f, cell)
        if v /= f%no_data_value: coeff = v
    aux(friction_index,i,j) = coeff
```
- Pull `manning_break`/`manning_coefficient`/`num_manning` from `geoclaw_module` to seed the depth/constant level (strategy **A**, recommended: keeps `src2.f90` untouched, fixes the current uncovered-cell gap as a side effect).
- Strategy **B** (fallback, lower effort): leave `set_friction_field` region-only, and in `src2.f90`/`src1d.f90` treat a sentinel `aux(friction_index)` (e.g. `<= 0`) as "use depth-based." More invasive to the hot loop; only if A proves awkward. Flag the choice for review.
- Apply the same fix in `multilayer/setaux.f90`'s path.
- Resolve the `sea_level` concern in the depth-break comparison while here.
- **Acceptance:** see Task 5 regression.

### Task 3 — Python `FrictionData`
- `write()`: emit each entry as `[friction_type, path]` mirroring `topofiles`; fix the format string; write the abspath resolved relative to `out_file`'s directory (existing intent).
- `read()`: parse `num_friction_files` and the entries instead of hard-coding empty.
- Add a short `setrun` usage example:
```python
fdata = rundata.friction_data
fdata.variable_friction = True
fdata.friction_index    = 3
fdata.friction_regions.append([lower, upper, depths, coeffs])   # existing
fdata.friction_files.append([3, 'friction_manning.tt3'])        # NEW: [friction_type, path]
```
- **Acceptance:** round-trip `write()`→`read()` reproduces `friction_files`; generated `friction.data` parses in `setup_variable_friction`.

### Task 4 — Python preprocessing: land-cover → friction file
New helper (e.g. `src/python/geoclaw/friction_tools.py`, or extend `topotools`). Optional deps `rasterio` (+ `rasterio.warp`, `rasterio.features`) or `rioxarray`/`geocube`; guard imports so core GeoClaw doesn't require them.

Pipeline (categorical-correct):
1. Read source raster/vector via GDAL/rasterio (COG/GeoTIFF/IMG/HDF; rasterize vector sources).
2. Reproject to lon/lat (EPSG:4326, GeoClaw grid) with **nearest** at native resolution (never bilinear on class codes).
3. Apply a two-stage, versioned lookup: `dataset_code → semantic_class → n`.
4. Area-average $n$ from native resolution to the target friction-tile resolution (mesh-scale averaging; the ADCIRC toolchain approach).
5. Write a topotype-3 or netCDF-4 friction file via the existing `Topography` writer.

Ship per-dataset legend crosswalks (NLCD-16, WorldCover-11, C-CAP, CORINE-44) plus one semantic-class→$n$ table.
- **Acceptance:** a small synthetic classified raster maps to a known $n$ field; a partial-coverage tile writes with `no_data` outside the classified footprint.

### Task 5 — Tests + example
- Fortran regression on a small domain with one friction tile partially covering it. Assert: (a) covered cells get the file $n$; (b) uncovered cells get the depth-based $n$, **not 0**; (c) `no_data` cells inside the tile fall to region/depth. Compare `aux(friction_index)` against an analytic expected field. Add to the geoclaw regression suite.
- Python unit test for Task 4 on a synthetic raster.
- Add an `examples/` case (storm-surge or tsunami) that ships a friction file end to end.

## Data sources the Python path should ingest (reference)

Format collapses to a few GDAL-readable families; the real work is CRS + class coding + raster/vector, not container parsing.

| Source | Container | CRS | Res | Classes |
|---|---|---|---|---|
| NLCD (Annual C1) | COG GeoTIFF | Albers/WGS84 | 30 m | 16-class Anderson II (+ RAT sidecar) |
| NOAA C-CAP (regional / 1 m) | GeoTIFF/IMG | Albers | 30 m, 1 m | ~25 coastal (+ RAT) |
| ESA WorldCover | COG GeoTIFF, 3°×3° tiles | EPSG:4326 | 10 m | 11-class (+ classes CSV) |
| Copernicus GLC (CGLS-LC100) | GeoTIFF | EPSG:4326 | 100 m | 23-class + fractional bands |
| CORINE (CLC) | GeoPackage/GDB **vector** or 100 m raster | ETRS89-LAEA 3035 | 100 m | 44-class |
| ESRI 10 m (Impact Observatory) | COG GeoTIFF, STAC | UTM/tile | 10 m | 9-class |
| MODIS MCD12Q1 | HDF-EOS (HDF4), or GeoTIFF via AppEEARS | Sinusoidal | 500 m | IGBP/UMD/PFT |
| NWI (wetlands) | GDB/shapefile/GPKG **vector** | various | vector | Cowardin |
| EMODnet Seabed Substrate / usSEABED (offshore) | **vector** / point CSV | 4326/3035 | vector | Folk grain-size |

Offshore friction is essentially unconstrained; use a constant open-water $n \approx 0.02\text{–}0.025$ and let the depth-based fallback carry the water column. That is why the file only needs to cover coastal tiles.

### Representative class → Manning's $n$ (surge/tsunami literature)
| Class | $n$ |
|---|---|
| Open water | 0.02 – 0.025 |
| Developed low→high | 0.05 – 0.15 |
| Forest (decid/evergreen) | 0.10 – 0.16 |
| Wetland / marsh | 0.045 – 0.07 |
| Cultivated / pasture | 0.035 – 0.05 |
| Barren | 0.025 |

Values are a starting table; keep them in the versioned lookup, not in code, so they can be swapped/cited (cf. Bunya et al. 2010, Mattocks & Forbes, Chow 1959).

## Out of scope
- Overland **wind reduction** (ADCIRC's `z0` / surface-canopy fields). GeoClaw's surge wind drag isn't keyed off land roughness; that's a separate aux field and a separate feature. Do not fold it into the bottom-friction file.
- Time-varying friction (morphodynamic roughness). Static only.

## How comparable codes do it (context)
- **ADCIRC:** `mannings_n_at_sea_floor` nodal attribute in `fort.13`; resolved $n$ computed offline by mesh-scale averaging a land-cover raster; negative-$n$ sentinel = fall back to depth-based Cd (our no_data mask is the gridded analog).
- **HEC-RAS 2D / TUFLOW / LISFLOOD-FP:** land-cover raster + lookup → gridded $n$ aligned to the DEM. Closest analog to this design.
- **Delft3D-FM:** trachytopes (class fractions → combined roughness); more than we need for a static $n$ field.
