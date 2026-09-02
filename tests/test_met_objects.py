#!/usr/bin/env python
# encoding: utf-8

"""Unit tests for the meteorological-forcing object model (Phase 1).

These tests exercise the new object model introduced by the met-forcing
refactor -- :class:`Track`, :class:`StormTrack`,
:class:`ParametricMetForcing`, and :class:`GriddedMetForcing` -- directly
(i.e. not only through the ``Storm`` compatibility wrapper).  They assert:

* the four objects construct;
* the readers (``read_geoclaw`` / ``read_atcf`` / ``read_data``) return the new
  objects;
* a new-object ``write_geoclaw`` / ``write_data`` produces bytes byte-identical
  to the ``Storm`` path *and* to the committed Phase-0 write goldens in
  ``tests/data/storm/characterization/``.
"""

from pathlib import Path
import sys
import warnings

import numpy as np
import pytest

import clawpack.geoclaw.met.storm as storm
from clawpack.geoclaw.met.track import Track, StormTrack, fill_rad_w_other_source
from clawpack.geoclaw.met.parametric import ParametricMetForcing
from clawpack.geoclaw.met.gridded import GriddedMetForcing

# ``tests/`` has no package __init__; make the sibling helpers importable.
sys.path.insert(0, str(Path(__file__).parent))
from test_storm import (_storm_input_path, _storm_check_path,  # noqa: E402
                        _make_storm_from_format)
from test_storm_characterization import _descriptor_head, golden_dir  # noqa: E402


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

@pytest.mark.python
@pytest.mark.storm
def test_objects_construct():
    """All four core objects construct with sensible empty defaults."""
    track = Track()
    assert track.t is None
    assert track.center is None
    # eye_location is an alias of center; ID an alias of id.
    track.eye_location = np.zeros((2, 2))
    assert track.center is track.eye_location

    storm_track = StormTrack()
    assert isinstance(storm_track, Track)
    for field in ("max_wind_speed", "max_wind_radius", "central_pressure",
                  "storm_radius", "classification", "basin", "wind_speeds"):
        assert getattr(storm_track, field) is None

    parametric = ParametricMetForcing()
    assert isinstance(parametric.track, StormTrack)
    assert parametric.file_paths == []

    gridded = GriddedMetForcing()
    assert gridded.scaling == [1.0, 1.0]
    assert gridded.crop_extent is None
    assert gridded.met_variable_map == {}


# ---------------------------------------------------------------------------
# Readers produce the new objects
# ---------------------------------------------------------------------------

@pytest.mark.python
@pytest.mark.storm
def test_read_geoclaw_returns_parametric():
    """ParametricMetForcing.read_geoclaw returns a ParametricMetForcing view."""
    forcing = ParametricMetForcing.read_geoclaw(_storm_check_path("tcvitals"))
    assert isinstance(forcing, ParametricMetForcing)
    assert isinstance(forcing.track, StormTrack)
    assert forcing.file_format == "geoclaw"
    assert forcing.eye_location.shape[1] == 2
    assert forcing.t.dtype.kind == "M"  # datetime64 track axis


@pytest.mark.python
@pytest.mark.storm
def test_read_atcf_returns_stormtrack():
    """StormTrack.read_atcf returns a populated StormTrack."""
    pytest.importorskip("pandas")
    track = StormTrack.read_atcf(_storm_input_path("atcf"))
    assert isinstance(track, StormTrack)
    assert track.file_format == "atcf"
    assert track.basin == "Atlantic"
    assert track.center.shape[1] == 2
    assert track.wind_speeds is not None


# ---------------------------------------------------------------------------
# IBTrACS reader output types.
#
# read_ibtracs is the only reader backed by xarray, and it used to let xarray
# types leak out: ``t`` stayed a DataArray and ``classification`` stayed bytes.
# These are contract tests rather than snapshots, so they keep holding as the
# fixture or the xarray version changes.
# ---------------------------------------------------------------------------

_IBTRACS_KWARGS = {"sid": "2008245N17323", "agency_pref": ["wmo", "usa"]}


@pytest.mark.python
@pytest.mark.storm
def test_ibtracs_reader_types():
    """read_ibtracs returns plain NumPy/Python types, like every other reader."""
    pytest.importorskip("xarray")
    track = StormTrack.read_ibtracs(_storm_input_path("ibtracs"),
                                    **_IBTRACS_KWARGS)

    # Times: a datetime64 ndarray, not an xarray DataArray.  Indexing a
    # DataArray yields 0-d DataArrays, which is what broke
    # fill_rad_w_other_source.
    assert isinstance(track.t, np.ndarray)
    assert track.t.dtype == np.dtype("datetime64[s]")
    assert isinstance(track.t[0], np.datetime64)

    # time_offset is drawn from t, so it follows.
    assert isinstance(track.time_offset, np.datetime64)

    # Classification decoded from the netCDF's bytes, so it compares against
    # ordinary string literals rather than b'TD'.
    assert track.classification.dtype.kind == "U"
    assert "TD" in track.classification
    assert not any(str(value).startswith("b'") for value in track.classification)

    assert track.event.dtype.kind == "U"


@pytest.mark.python
@pytest.mark.storm
def test_ibtracs_second_resolution_times():
    """IBTrACS times land on whole seconds, not a nanosecond roundoff tail.

    IBTrACS stores times that decode to values like
    ``2008-09-01T06:00:00.000039936``.  Carried through to ``write_geoclaw``
    those turn nominal hour offsets into ``-3.60000003e+03`` instead of
    ``-3.60000000e+03``; the committed ``ibtracs_geoclaw.txt`` baseline has the
    exact values, so truncating to seconds restores agreement with it.
    """
    pytest.importorskip("xarray")
    track = StormTrack.read_ibtracs(_storm_input_path("ibtracs"),
                                    **_IBTRACS_KWARGS)

    offsets = (track.t - track.time_offset) / np.timedelta64(1, "s")
    # Every IBTrACS observation is on a whole minute; no sub-second residue.
    assert np.allclose(offsets % 60.0, 0.0, atol=0.0)


@pytest.mark.python
@pytest.mark.storm
def test_ibtracs_reader_no_future_warnings():
    """read_ibtracs uses no deprecated xarray API.

    ``Dataset.dims`` returning a mapping is deprecated and becomes a set of
    dimension names; ``Dataset.sizes`` is the replacement.  Promoting the
    warning to an error keeps this from silently rotting on the next xarray
    release.
    """
    pytest.importorskip("xarray")
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        # The reader legitimately warns about missing RMW/ROCI; that is not the
        # class of warning under test here.
        warnings.filterwarnings("ignore", category=UserWarning)
        StormTrack.read_ibtracs(_storm_input_path("ibtracs"), **_IBTRACS_KWARGS)


@pytest.mark.python
@pytest.mark.storm
def test_fill_rad_accepts_scalar_or_zero_d_time():
    """fill_rad_w_other_source accepts a scalar and a 0-d DataArray time alike.

    Callers pass ``storm.t[n]``, so whichever of those a reader produces has to
    work.  A 0-d DataArray used to raise
    ``ValueError: Could not convert object to NumPy datetime``.
    """
    xr = pytest.importorskip("xarray")
    pytest.importorskip("pandas")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        target = StormTrack.read_ibtracs(_storm_input_path("ibtracs"),
                                         **_IBTRACS_KWARGS)
        source = StormTrack.read_atcf(_storm_input_path("atcf"))

    scalar_t = target.t[len(target.t) // 2]
    from_scalar = fill_rad_w_other_source(scalar_t, target, source,
                                          "max_wind_radius")

    zero_d = xr.DataArray(scalar_t)
    assert zero_d.ndim == 0
    from_zero_d = fill_rad_w_other_source(zero_d, target, source,
                                          "max_wind_radius")

    assert np.isclose(from_scalar, from_zero_d)


@pytest.mark.python
@pytest.mark.storm
def test_read_data_returns_gridded(tmp_path):
    """GriddedMetForcing.read_data returns a populated GriddedMetForcing."""
    # Produce an ASCII descriptor with the new object, then read it back.
    gridded = _make_ascii_gridded()
    descriptor = tmp_path / "ascii.storm"
    gridded.write_data(descriptor)

    read = GriddedMetForcing.read_data(descriptor)
    assert isinstance(read, GriddedMetForcing)
    assert read.file_format == 1
    assert read.crop_extent == [-100.0, -60.0, 10.0, 40.0]
    assert read.ramp_width == 3.0
    assert read.x_shift == 1.25
    assert read.y_shift == -0.5
    assert len(read.file_paths) == 4


# ---------------------------------------------------------------------------
# Missing-data contract
#
# np.nan is the in-memory marker in every reader; write_geoclaw additionally
# tolerates a negative for callers still on the v5.9.0 -1 convention; zero is
# never missing.  These lock all three halves of that contract.
# ---------------------------------------------------------------------------

def _synthetic_track(**overrides):
    """A minimal four-point ParametricMetForcing for writer tests."""
    n = 4
    forcing = ParametricMetForcing()
    forcing.t = np.array(["2020-01-01T00", "2020-01-01T06",
                          "2020-01-01T12", "2020-01-01T18"],
                         dtype="datetime64[s]")
    forcing.time_offset = forcing.t[0]
    forcing.eye_location = np.column_stack([np.linspace(-80.0, -83.0, n),
                                            np.linspace(25.0, 28.0, n)])
    forcing.max_wind_speed = np.full(n, 40.0)
    forcing.central_pressure = np.full(n, 98000.0)
    forcing.max_wind_radius = np.full(n, 40e3)
    forcing.storm_radius = np.full(n, 300e3)
    for name, value in overrides.items():
        setattr(forcing, name, value)
    return forcing


@pytest.mark.python
@pytest.mark.storm
@pytest.mark.parametrize("field", ["max_wind_speed", "central_pressure",
                                   "max_wind_radius", "storm_radius"])
def test_missing_marker_equivalence(tmp_path, field):
    """A NaN-marked and a -1-marked storm write byte-identical files.

    The readers now emit NaN, but objects built by hand and community fill
    functions may still use -1.  Both must reach the same fill/skip branch.
    """
    nan_values = np.full(4, np.nan)
    sentinel_values = np.full(4, -1.0)

    nan_path = tmp_path / "nan.storm"
    sentinel_path = tmp_path / "sentinel.storm"

    _synthetic_track(**{field: nan_values}).write_geoclaw(nan_path)
    with warnings.catch_warnings():
        # The -1 path warns that the convention is deprecated; that is the
        # documented behavior, not a failure.
        warnings.simplefilter("ignore", DeprecationWarning)
        _synthetic_track(**{field: sentinel_values}).write_geoclaw(
            sentinel_path)

    assert nan_path.read_bytes() == sentinel_path.read_bytes()


@pytest.mark.python
@pytest.mark.storm
def test_failed_fill_is_skipped_not_written(tmp_path):
    """A fill that cannot produce a value leaves the row skipped.

    Returning NaN or -1 from a fill used to be written into the file verbatim,
    producing a storm file GeoClaw cannot run.
    """
    for sentinel in (np.nan, -1.0):
        path = tmp_path / f"fill_{sentinel}.storm"
        forcing = _synthetic_track(max_wind_radius=np.full(4, np.nan))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            forcing.write_geoclaw(
                path, fill_dict={"max_wind_radius": lambda t, s: sentinel})
        # Header line is the cast count; every row should have been skipped.
        assert path.read_text().splitlines()[0].strip() == "0"


@pytest.mark.python
@pytest.mark.storm
def test_zero_is_not_missing(tmp_path):
    """Zero is real data, not a missing marker.

    HURDAT2 reports genuinely-zero quadrant wind radii for weak systems.  A
    'treat <= 0 as missing' rule would silently reinterpret those.
    """
    path = tmp_path / "zero.storm"
    forcing = _synthetic_track(max_wind_speed=np.zeros(4))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # No fill supplied: if zero were treated as missing every row would be
        # skipped and this would be a zero-cast file.
        forcing.write_geoclaw(path)

    lines = path.read_text().splitlines()
    assert lines[0].strip() == "4"
    assert float(lines[3].split()[3]) == 0.0


@pytest.mark.python
@pytest.mark.storm
def test_write_geoclaw_does_not_mutate_storm(tmp_path):
    """Writing twice with different fills honors the second fill.

    Fills used to be written back onto the storm object, so the second write saw
    a storm that was no longer missing and silently reused the first fill.
    """
    forcing = _synthetic_track(max_wind_radius=np.full(4, np.nan))

    first = tmp_path / "first.storm"
    second = tmp_path / "second.storm"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        forcing.write_geoclaw(
            first, fill_dict={"max_wind_radius": lambda t, s: 20e3})
        forcing.write_geoclaw(
            second, fill_dict={"max_wind_radius": lambda t, s: 60e3})

    assert float(first.read_text().splitlines()[3].split()[4]) == 20e3
    assert float(second.read_text().splitlines()[3].split()[4]) == 60e3
    # And the storm itself is untouched.
    assert np.all(np.isnan(forcing.max_wind_radius))


@pytest.mark.python
@pytest.mark.storm
def test_fill_dict_storm_radius_not_clobbered(tmp_path):
    """A caller-supplied storm_radius fill beats the built-in 500 km default.

    The default was applied with dict.update() *after* copying the caller's
    fill_dict, so a caller-supplied storm_radius fill was silently discarded.
    """
    path = tmp_path / "roci.storm"
    forcing = _synthetic_track(storm_radius=np.full(4, np.nan))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        forcing.write_geoclaw(
            path, fill_dict={"storm_radius": lambda t, s: 250e3})

    assert float(path.read_text().splitlines()[3].split()[6]) == 250e3


@pytest.mark.python
@pytest.mark.storm
def test_hurdat_sentinels_become_nan(tmp_path):
    """HURDAT2 -99 wind / -999 pressure become NaN, not scaled sentinels.

    Converting units before normalizing turned -999 mbar into -99900.0 Pa, a
    value that looks physical enough to be written out.
    """
    header = "AL011980,          UNNAMED,      2,\n"
    rows = ("19800101, 0000,  , TD, 25.0N,  80.0W, -99, -999,"
            + " -999," * 11 + " -999,\n"
            "19800101, 0600,  , TD, 25.5N,  80.5W,  35, 1005,"
            + " -999," * 11 + " -999,\n")
    path = tmp_path / "sentinel_hurdat.txt"
    path.write_text(header + rows)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        track = StormTrack.read_hurdat(path)

    assert np.isnan(track.max_wind_speed[0])
    assert np.isnan(track.central_pressure[0])
    # The valid second row is untouched and correctly converted.
    assert np.isclose(track.central_pressure[1], 100500.0)
    assert not np.isnan(track.max_wind_speed[1])


@pytest.mark.python
@pytest.mark.storm
@pytest.mark.parametrize("file_format", ["atcf", "tcvitals", "ibtracs"])
def test_written_radii_are_positive(tmp_path, file_format):
    """Every storm file written from a bundled input is runnable.

    A non-positive storm_radius zeros the forcing through the spatial ramp and a
    non-positive max_wind_radius divides by zero in the Holland profiles, so a
    file with either is not a valid GeoClaw input.  hurdat and jma are excluded
    until met.reconstruction supplies real RMW/ROCI estimators; their current
    baselines use placeholder zeros.
    """
    if file_format == "ibtracs":
        pytest.importorskip("xarray")
    pytest.importorskip("pandas")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        test_storm, fill_mwr, fill_rad = _make_storm_from_format(file_format)
        write_kwargs = {}
        if fill_mwr is not None:
            write_kwargs["max_wind_radius_fill"] = fill_mwr
        if fill_rad is not None:
            write_kwargs["storm_radius_fill"] = fill_rad
        path = tmp_path / f"{file_format}.storm"
        test_storm.write(path, file_format="geoclaw", **write_kwargs)

    values = np.loadtxt(path, skiprows=3)
    assert values.shape[0] > 0
    assert np.all(np.isfinite(values))
    assert np.all(values[:, 4] > 0.0), "max_wind_radius must be positive"
    assert np.all(values[:, 6] > 0.0), "storm_radius must be positive"


# ---------------------------------------------------------------------------
# Write byte-equality: new object == Storm path == committed golden
# ---------------------------------------------------------------------------

WRITE_GEOCLAW_FORMATS = ["atcf", "tcvitals"]


@pytest.mark.python
@pytest.mark.storm
@pytest.mark.parametrize("file_format", WRITE_GEOCLAW_FORMATS)
def test_write_geoclaw_matches_storm_and_golden(tmp_path, file_format):
    """New-object write_geoclaw == Storm-path write == committed golden."""
    if file_format == "atcf":
        pytest.importorskip("pandas")

    input_path = _storm_input_path(file_format)
    reader = getattr(StormTrack, "read_%s" % file_format)

    # New-object path.
    track = reader(input_path)
    forcing = ParametricMetForcing(track=track)
    new_out = tmp_path / f"{file_format}_new.storm"
    forcing.write_geoclaw(new_out)

    # Storm compatibility path.
    s = storm.Storm(input_path, file_format=file_format)
    storm_out = tmp_path / f"{file_format}_storm.storm"
    s.write(storm_out, file_format="geoclaw")

    new_bytes = new_out.read_text()
    golden = (golden_dir / f"write_geoclaw_{file_format}.txt").read_text()

    assert new_bytes == storm_out.read_text()
    assert new_bytes == golden


def _make_ascii_gridded():
    """A deterministic OWI/NWS12 ASCII GriddedMetForcing (fixed controls)."""
    gridded = GriddedMetForcing()
    gridded.time_offset = np.datetime64("2012-08-29")
    gridded.crop_extent = [-100.0, -60.0, 10.0, 40.0]
    gridded.ramp_width = 3
    gridded.x_shift = 1.25
    gridded.y_shift = -0.5
    gridded.file_format = "ascii"
    gridded.file_paths = [Path("storm_1.PRE"), Path("storm_1.WIN"),
                          Path("storm_2.PRE"), Path("storm_2.WIN")]
    return gridded


@pytest.mark.python
@pytest.mark.storm
def test_write_data_matches_storm_and_golden(tmp_path):
    """New-object write_data == Storm-path write == committed golden (ASCII)."""
    # New-object path.
    new_desc = tmp_path / "ascii_new.storm"
    _make_ascii_gridded().write_data(new_desc)

    # Storm compatibility path.
    s = storm.Storm()
    s.time_offset = np.datetime64("2012-08-29")
    s.crop_extent = [-100.0, -60.0, 10.0, 40.0]
    s.ramp_width = 3
    s.x_shift = 1.25
    s.y_shift = -0.5
    s.file_format = "ascii"
    s.file_paths = [Path("storm_1.PRE"), Path("storm_1.WIN"),
                    Path("storm_2.PRE"), Path("storm_2.WIN")]
    storm_desc = tmp_path / "ascii_storm.storm"
    s.write(storm_desc, file_format="data")

    new_head = _descriptor_head(new_desc.read_text())
    golden = (golden_dir / "write_data_ascii.txt").read_text()

    assert new_head == _descriptor_head(storm_desc.read_text())
    assert new_head == golden


@pytest.mark.python
@pytest.mark.storm
@pytest.mark.netcdf
def test_write_data_netcdf_matches_golden(tmp_path):
    """New-object write_data (NWS13 netCDF) == committed descriptor golden."""
    pytest.importorskip("xarray")
    pytest.importorskip("netCDF4")
    from test_storm import create_nws13_storm_file

    nc = tmp_path / "nws13.nc"
    create_nws13_storm_file(nc)

    gridded = GriddedMetForcing()
    gridded.time_offset = np.datetime64("2012-08-29")
    gridded.crop_extent = [-100.0, -60.0, 10.0, 40.0]
    gridded.ramp_width = 3
    gridded.x_shift = 1.25
    gridded.y_shift = -0.5
    gridded.file_format = "nws13"
    gridded.file_paths = [nc]
    descriptor = tmp_path / "nws13.storm"
    gridded.write_data(descriptor,
                       var_mapping={"wind_u": "uwnd", "wind_v": "vwnd",
                                    "pressure": "press"})

    new_head = _descriptor_head(descriptor.read_text())
    golden = (golden_dir / "write_data_nws13.txt").read_text()
    assert new_head == golden


# ---------------------------------------------------------------------------
# Regression coverage for previously-latent surge helper bugs (wrap-up fixes)
# ---------------------------------------------------------------------------

@pytest.mark.python
@pytest.mark.storm
def test_storm_str_datetime64():
    """Storm.__str__ renders a datetime64 track axis without error.

    Regression: the datetime64 branch used a typo (``np.datetiem64``) and
    ``.isoformat()``, which np.datetime64 does not provide.
    """
    s = storm.Storm()
    s.name = "TESTSTORM"
    s.t = np.array(["2020-08-01T00:00", "2020-08-01T06:00"],
                   dtype="datetime64[s]")
    s.file_paths = ["a.storm"]
    text = str(s)
    assert "TESTSTORM" in text
    assert "2020-08-01T00:00:00" in text
    assert "2020-08-01T06:00:00" in text


@pytest.mark.python
@pytest.mark.storm
def test_construct_fields_resolves_radius():
    """construct_fields resolves its radius argument (no NameError) and reaches
    the not-yet-implemented model stub.

    Regression: the call passed an undefined name ``x`` instead of ``r``.
    """
    s = storm.Storm()
    with pytest.raises(NotImplementedError):
        storm.construct_fields(s, 1.0e3, 0.0, model="holland_1980")


@pytest.mark.python
@pytest.mark.storm
def test_make_multi_structure_splits_by_storm(tmp_path):
    """make_multi_structure splits a multi-storm ATCF into per-storm Storms.

    Regression: the helper referenced ``os`` without importing it, and grouped
    by timestamp rather than storm identity (colliding on os.mkdir).
    """
    pytest.importorskip("pandas")
    from clawpack.geoclaw.met.tools import make_multi_structure

    # Build a 2-storm fixture from real ATCF records: relabel the cyclone
    # number (field 1) so the two synthetic storms are distinguishable.
    with open(_storm_input_path("atcf")) as data_file:
        records = [line for line in data_file
                   if len(line.split(",")) > 8][:4]
    assert records, "expected ATCF records in the test fixture"
    basin = records[0].split(",")[0].strip()

    def relabel(line, cyclone):
        fields = line.split(",")
        fields[1] = " %s" % cyclone
        return ",".join(fields)

    fixture = tmp_path / "multi.atcf"
    fixture.write_text("".join([relabel(r, "09") for r in records]
                               + [relabel(r, "11") for r in records]))

    storms = make_multi_structure(str(fixture),
                                  output_dir=str(tmp_path / "clipped"))
    assert list(storms.keys()) == [basin + "09", basin + "11"]
    for split in storms.values():
        assert isinstance(split, storm.Storm)
        assert len(split.t) == len(records)


@pytest.mark.python
@pytest.mark.storm
def test_surgedata_lives_in_geoclaw_data():
    """SurgeData/FrictionData live in clawpack.geoclaw.data, not surge.data.

    Regression: isaac/setplot_kml.py imported the nonexistent
    ``clawpack.geoclaw.surge.data`` module.
    """
    import clawpack.geoclaw.data as geodata
    assert hasattr(geodata, "SurgeData")
    assert hasattr(geodata, "FrictionData")
    with pytest.raises(ImportError):
        import clawpack.geoclaw.surge.data  # noqa: F401


@pytest.mark.python
@pytest.mark.storm
def test_isaac_setplot_kml_imports():
    """isaac/setplot_kml.py imports cleanly after the surge.data import fix."""
    pytest.importorskip("matplotlib")
    pytest.importorskip("clawpack.visclaw")
    import importlib.util

    kml_path = (Path(__file__).parents[1] / "examples" / "storm-surge"
                / "isaac" / "setplot_kml.py")
    if not kml_path.exists():
        pytest.skip("isaac/setplot_kml.py not present")
    spec = importlib.util.spec_from_file_location("isaac_setplot_kml", kml_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert hasattr(module, "setplot")


# ---------------------------------------------------------------------------
# met / surge package aliasing (Decision 1 -> B)
# ---------------------------------------------------------------------------

@pytest.mark.python
@pytest.mark.storm
def test_met_package_reexports_model():
    """``clawpack.geoclaw.met`` re-exports the storm/met object model."""
    import clawpack.geoclaw.met as met
    for name in ("Storm", "Track", "StormTrack", "ParametricMetForcing",
                 "GriddedMetForcing", "OWIData", "SurgeData", "MetData",
                 "construct_fields", "available_formats", "available_models"):
        assert hasattr(met, name), f"clawpack.geoclaw.met is missing {name}"


@pytest.mark.python
@pytest.mark.storm
def test_metdata_is_surgedata_alias():
    """``clawpack.geoclaw.data`` exposes ``MetData`` as an alias of ``SurgeData``."""
    import clawpack.geoclaw.data as geodata
    assert geodata.MetData is geodata.SurgeData


@pytest.mark.python
@pytest.mark.storm
def test_surge_import_emits_deprecation_warning():
    """Importing the deprecated ``surge`` package warns and points at ``met``."""
    import importlib
    import warnings
    import clawpack.geoclaw.surge as surge
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        importlib.reload(surge)   # re-run the package body to observe the warning
    dep = [w for w in caught if issubclass(w.category, DeprecationWarning)]
    assert len(dep) >= 1, "surge import did not emit a DeprecationWarning"
    assert "clawpack.geoclaw.met" in str(dep[0].message)


@pytest.mark.python
@pytest.mark.storm
def test_surge_shim_delegates_to_met():
    """``surge.*`` names resolve to the identical objects in ``met.*``."""
    import clawpack.geoclaw.met.storm as met_storm
    import clawpack.geoclaw.met.track as met_track
    import clawpack.geoclaw.surge.storm as surge_storm
    import clawpack.geoclaw.surge.track as surge_track
    assert surge_storm.Storm is met_storm.Storm
    assert surge_track.StormTrack is met_track.StormTrack
    # Private names delegate too (gridded.py imports _Meta from track).
    assert surge_track._Meta is met_track._Meta


# ---------------------------------------------------------------------------
# OWI (Oceanweather WIN/PRE) field I/O
# ---------------------------------------------------------------------------

def _make_owi_data(nt=3, ny=4, nx=5):
    """Build a small deterministic OWIData object."""
    from clawpack.geoclaw.met.data_storms import OWIData

    lon = -99.0 + np.arange(nx) * 0.25
    lat = 8.0 + np.arange(ny) * 0.25
    t = (np.datetime64("2012-08-20T12:00:00")
         + (np.arange(nt) * 3600).astype("timedelta64[s]"))
    ramp = np.arange(nt * ny * nx, dtype=float).reshape(nt, ny, nx)
    return OWIData(time=t, longitude=lon, latitude=lat,
                   wind_u=ramp, wind_v=-ramp, pressure=101300.0 - ramp)


@pytest.mark.python
@pytest.mark.storm
def test_owi_roundtrip(tmp_path):
    """write_owi -> read_owi reproduces every field exactly."""
    from clawpack.geoclaw.met import data_storms

    d = _make_owi_data()
    pre, win = tmp_path / "x.PRE", tmp_path / "x.WIN"
    data_storms.write_owi(d, pre, win)
    r = data_storms.read_owi(pre, win)

    assert np.array_equal(d.time, r.time)
    np.testing.assert_allclose(d.longitude, r.longitude)
    np.testing.assert_allclose(d.latitude, r.latitude)
    np.testing.assert_allclose(d.wind_u, r.wind_u)
    np.testing.assert_allclose(d.wind_v, r.wind_v)
    # Pressure survives the Pa<->mbar conversion to the f10.4 precision.
    np.testing.assert_allclose(d.pressure, r.pressure, atol=1e-2)
    # The cheap start-time peek agrees with the first record.
    assert data_storms.read_owi_start_time(pre) == d.time[0]


@pytest.mark.python
@pytest.mark.storm
def test_owi_written_format_is_80_col(tmp_path):
    """Written OWI lines match the fixed 80-column contract the Fortran reads."""
    from clawpack.geoclaw.met import data_storms

    pre, win = tmp_path / "x.PRE", tmp_path / "x.WIN"
    data_storms.write_owi(_make_owi_data(), pre, win)
    lines = win.read_text().splitlines()
    assert len(lines[0]) == 80                       # file header
    assert lines[0].startswith("Oceanweather WIN/PRE Format")
    assert len(lines[1]) == 80                       # grid-spec header
    assert lines[1].startswith("iLat=") and "iLong=" in lines[1]
    assert "DT=" in lines[1]


@pytest.mark.python
@pytest.mark.storm
def test_owi_read_real_isaac():
    """The reader parses the committed real isaac.WIN/PRE sample."""
    from clawpack.geoclaw.met import data_storms

    isaac = (Path(__file__).parents[1] / "examples" / "storm-surge" / "isaac")
    if not (isaac / "isaac.PRE").exists():
        pytest.skip("isaac OWI sample not present")

    d = data_storms.read_owi(isaac / "isaac.PRE", isaac / "isaac.WIN")
    assert d.shape == (96, 116)                       # (ny, nx)
    assert d.wind_u.shape == (d.num_times, 96, 116)
    np.testing.assert_allclose(d.longitude[0], -99.0)
    np.testing.assert_allclose(d.latitude[0], 8.0)
    # Ambient far-field pressure: 1013 mbar -> 101300 Pa.
    np.testing.assert_allclose(d.pressure[0].max(), 101300.0)


@pytest.mark.python
@pytest.mark.storm
def test_gridded_from_owi_descriptor(tmp_path):
    """from_owi builds a format-1 descriptor pointing at the WIN/PRE pair."""
    from clawpack.geoclaw.met import data_storms

    pre, win = tmp_path / "x.PRE", tmp_path / "x.WIN"
    data_storms.write_owi(_make_owi_data(), pre, win)

    forcing = GriddedMetForcing.from_owi(pre, win)
    assert forcing.file_format == "owi"
    assert forcing.file_paths == [pre, win]           # [PRE, WIN] order
    assert forcing.time_offset == np.datetime64("2012-08-20T12:00:00")

    # A format-1 (ascii/owi) descriptor round-trips through read_data.
    descriptor = tmp_path / "met.storm"
    forcing.write_data(descriptor)
    back = GriddedMetForcing.read_data(descriptor)
    assert back.file_format == 1
    assert [p.name for p in back.file_paths] == ["x.PRE", "x.WIN"]


@pytest.mark.python
@pytest.mark.storm
def test_gridded_to_owi(tmp_path):
    """to_owi writes the pair and points the forcing at it."""
    from clawpack.geoclaw.met import data_storms

    d = _make_owi_data()
    pre, win = tmp_path / "x.PRE", tmp_path / "x.WIN"
    forcing = GriddedMetForcing().to_owi(d, pre, win)

    assert pre.exists() and win.exists()
    assert forcing.file_format == "owi"
    assert forcing.file_paths == [pre, win]
    # The files it wrote read back to the same fields.
    r = data_storms.read_owi(pre, win)
    np.testing.assert_allclose(d.wind_u, r.wind_u)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
