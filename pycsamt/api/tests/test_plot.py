from __future__ import annotations

import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

from pycsamt.api import plot as plot_mod
from pycsamt.api.plot import (
    PLOT_CONFIG,
    PlotConfig,
    add_colorbar,
    add_polar_colorbar,
    load_plot_config,
    reset_plot_config,
    save_fig,
    set_dpi,
    set_fmt,
    set_savedir,
    write_default_config,
)
from pycsamt.api.plot import (
    _load_config_file,
    _load_env,
    _normalise_token,
    _parse_fmt_tokens,
    _resolve_formats,
    _safe_print,
    _validate_fmt,
    _build_singleton,
)


@pytest.fixture(autouse=True)
def _restore_global_plot_config():
    """PLOT_CONFIG is a process-wide singleton; snapshot/restore it so
    tests that mutate it (set_fmt, set_dpi, ...) never leak into other
    tests."""
    snapshot = PLOT_CONFIG.to_dict()
    yield
    for k, v in snapshot.items():
        setattr(PLOT_CONFIG, k, v)


def _figure_with_mappable():
    fig, ax = plt.subplots()
    mappable = ax.imshow([[0, 1], [1, 0]])
    return fig, ax, mappable


# ─────────────────────────────────────────────────────────────────────────
# _safe_print
# ─────────────────────────────────────────────────────────────────────────


def test_safe_print_normal_message(capsys):
    _safe_print("hello")
    assert capsys.readouterr().out.strip() == "hello"


def test_safe_print_falls_back_on_unicode_encode_error(monkeypatch):
    class _BoomOnceStream:
        encoding = "cp1252"

        def __init__(self):
            self.calls = 0
            self.written = []

        def write(self, s):
            self.calls += 1
            if self.calls == 1:
                raise UnicodeEncodeError("cp1252", s, 0, 1, "boom")
            self.written.append(s)

    stream = _BoomOnceStream()
    monkeypatch.setattr(plot_mod.sys, "stdout", stream)
    monkeypatch.setattr(
        "builtins.print",
        lambda msg: stream.write(msg),
    )
    _safe_print("bad ✔ char")
    assert stream.written == ["bad ? char\n"]


# ─────────────────────────────────────────────────────────────────────────
# token helpers
# ─────────────────────────────────────────────────────────────────────────


def test_normalise_token_strips_whitespace_and_dot():
    assert _normalise_token("  .svg  ") == "svg"


def test_parse_fmt_tokens_from_comma_string():
    assert _parse_fmt_tokens("png, +svg ,pdf") == ["png", "+svg", "pdf"]


def test_parse_fmt_tokens_from_iterable_with_embedded_commas():
    assert _parse_fmt_tokens(["png,svg", "+pdf"]) == ["png", "svg", "+pdf"]


def test_resolve_formats_plain_and_additive():
    assert _resolve_formats(["png"], "png") == ["png"]
    assert _resolve_formats(["+svg"], "png") == ["png", "svg"]
    assert _resolve_formats(["svg", "pdf"], "png") == ["svg", "pdf"]
    assert _resolve_formats([], "png") == ["png"]


def test_resolve_formats_skips_empty_and_duplicate_tokens():
    # "+" alone -> fmt becomes "" and is skipped; "png" repeated is deduped
    assert _resolve_formats(["+", "png", "png"], "png") == ["png"]


def test_validate_fmt_warns_on_unknown_format():
    with pytest.warns(UserWarning, match="not in the known list"):
        _validate_fmt(["bogus"])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _validate_fmt(["png", "svg"])  # no warning raised


# ─────────────────────────────────────────────────────────────────────────
# _load_config_file / _load_env
# ─────────────────────────────────────────────────────────────────────────


def test_load_config_file_returns_empty_when_no_file_present(monkeypatch, tmp_path):
    monkeypatch.setattr(plot_mod, "_CFG_LOCAL", tmp_path / "missing_local.cfg")
    monkeypatch.setattr(plot_mod, "_CFG_GLOBAL", tmp_path / "missing_global.cfg")
    assert _load_config_file() == {}


def test_load_config_file_reads_local_cfg(monkeypatch, tmp_path):
    cfg_path = tmp_path / "pycsamt_plot.cfg"
    cfg_path.write_text(
        "[plot]\n"
        "fmt = png, +svg\n"
        "base_fmt = png\n"
        "dpi = 222\n"
        "bbox_inches = standard\n"
        "transparent = true\n"
        "facecolor = black\n"
        "savedir = ~/figs\n"
        "close_after_save = true\n"
        "verbose = false\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(plot_mod, "_CFG_LOCAL", cfg_path)
    monkeypatch.setattr(plot_mod, "_CFG_GLOBAL", tmp_path / "unused.cfg")
    out = _load_config_file()
    assert out["fmt"] == ["png", "+svg"]
    assert out["dpi"] == 222
    assert out["transparent"] is True
    assert out["close_after_save"] is True
    assert out["verbose"] is False
    assert out["savedir"].endswith("figs")


def test_load_config_file_partial_fields_only(monkeypatch, tmp_path):
    cfg_path = tmp_path / "pycsamt_plot.cfg"
    cfg_path.write_text("[plot]\ndpi = 200\n", encoding="utf-8")
    monkeypatch.setattr(plot_mod, "_CFG_LOCAL", cfg_path)
    monkeypatch.setattr(plot_mod, "_CFG_GLOBAL", tmp_path / "unused.cfg")
    out = _load_config_file()
    assert out == {"dpi": 200}


def test_load_config_file_single_non_dpi_field(monkeypatch, tmp_path):
    cfg_path = tmp_path / "pycsamt_plot.cfg"
    cfg_path.write_text("[plot]\nbase_fmt = svg\n", encoding="utf-8")
    monkeypatch.setattr(plot_mod, "_CFG_LOCAL", cfg_path)
    monkeypatch.setattr(plot_mod, "_CFG_GLOBAL", tmp_path / "unused.cfg")
    out = _load_config_file()
    assert out == {"base_fmt": "svg"}


def test_load_config_file_without_plot_section_returns_empty(monkeypatch, tmp_path):
    cfg_path = tmp_path / "pycsamt_plot.cfg"
    cfg_path.write_text("[other]\nkey = value\n", encoding="utf-8")
    monkeypatch.setattr(plot_mod, "_CFG_LOCAL", cfg_path)
    monkeypatch.setattr(plot_mod, "_CFG_GLOBAL", tmp_path / "unused.cfg")
    assert _load_config_file() == {}


def test_load_env_reads_every_pycsamt_variable(monkeypatch):
    monkeypatch.setenv("PYCSAMT_FMT", "svg,+pdf")
    monkeypatch.setenv("PYCSAMT_BASE_FMT", "tiff")
    monkeypatch.setenv("PYCSAMT_DPI", "333")
    monkeypatch.setenv("PYCSAMT_SAVEDIR", "~/envdir")
    monkeypatch.setenv("PYCSAMT_TRANSPARENT", "true")
    monkeypatch.setenv("PYCSAMT_BBOX", "standard")
    out = _load_env()
    assert out["fmt"] == ["svg", "+pdf"]
    assert out["base_fmt"] == "tiff"
    assert out["dpi"] == 333
    assert out["savedir"].endswith("envdir")
    assert out["transparent"] is True
    assert out["bbox_inches"] == "standard"


def test_load_env_empty_when_unset(monkeypatch):
    for var in (
        "PYCSAMT_FMT",
        "PYCSAMT_BASE_FMT",
        "PYCSAMT_DPI",
        "PYCSAMT_SAVEDIR",
        "PYCSAMT_TRANSPARENT",
        "PYCSAMT_BBOX",
    ):
        monkeypatch.delenv(var, raising=False)
    assert _load_env() == {}


def test_build_singleton_layers_file_then_env(monkeypatch):
    monkeypatch.setattr(plot_mod, "_load_config_file", lambda: {"dpi": 111})
    monkeypatch.setattr(
        plot_mod, "_load_env", lambda: {"dpi": 222, "base_fmt": "svg"}
    )
    cfg = _build_singleton()
    assert cfg.dpi == 222  # env wins over file
    assert cfg.base_fmt == "svg"


# ─────────────────────────────────────────────────────────────────────────
# PlotConfig
# ─────────────────────────────────────────────────────────────────────────


def test_plot_config_resolve_formats_uses_own_fmt_by_default():
    cfg = PlotConfig(fmt="+svg")
    assert cfg.resolve_formats() == ["png", "svg"]
    assert cfg.resolve_formats(fmt="pdf") == ["pdf"]


def test_plot_config_save_accepts_axes_and_figure(tmp_path):
    fig, ax, _ = _figure_with_mappable()
    cfg = PlotConfig(dpi=80, verbose=False)
    paths = cfg.save(ax, tmp_path / "out_axes")
    assert paths[0].exists()
    assert paths[0].suffix == ".png"

    paths2 = cfg.save(fig, tmp_path / "out_fig", fmt="svg")
    assert paths2[0].suffix == ".svg"
    assert paths2[0].exists()
    plt.close(fig)


def test_plot_config_save_rejects_non_figure_object():
    cfg = PlotConfig()
    with pytest.raises(TypeError, match="Expected a matplotlib Figure or Axes"):
        cfg.save(object(), "whatever")


def test_plot_config_save_resolves_relative_path_under_savedir(tmp_path):
    fig, ax, _ = _figure_with_mappable()
    cfg = PlotConfig(savedir=tmp_path / "figs", verbose=False)
    paths = cfg.save(ax, "nested/out")
    assert paths[0] == (tmp_path / "figs" / "nested" / "out").with_suffix(".png")
    assert paths[0].exists()
    plt.close(fig)


def test_plot_config_save_absolute_path_ignores_savedir(tmp_path):
    fig, ax, _ = _figure_with_mappable()
    cfg = PlotConfig(savedir=tmp_path / "figs", verbose=False)
    target = tmp_path / "elsewhere" / "out"
    paths = cfg.save(ax, target)
    assert paths[0] == target.with_suffix(".png")
    plt.close(fig)


def test_plot_config_save_verbose_prints_each_file(tmp_path, capsys):
    fig, ax, _ = _figure_with_mappable()
    cfg = PlotConfig(verbose=True, fmt=["png", "svg"])
    cfg.save(ax, tmp_path / "out")
    out = capsys.readouterr().out
    assert "out.png" in out
    assert "out.svg" in out
    plt.close(fig)


def test_plot_config_save_close_after_save_closes_figure(tmp_path):
    fig, ax, _ = _figure_with_mappable()
    cfg = PlotConfig(close_after_save=True, verbose=False)
    cfg.save(ax, tmp_path / "out")
    assert not plt.fignum_exists(fig.number)


def test_plot_config_configure_sets_multiple_attrs():
    cfg = PlotConfig()
    cfg.configure(dpi=500, verbose=False)
    assert cfg.dpi == 500
    assert cfg.verbose is False


def test_plot_config_configure_rejects_unknown_attr():
    cfg = PlotConfig()
    with pytest.raises(AttributeError, match="has no attribute"):
        cfg.configure(bogus=1)


def test_plot_config_reset_restores_defaults():
    cfg = PlotConfig(dpi=999, verbose=False, savedir="/x")
    cfg.reset()
    assert cfg.dpi == 150
    assert cfg.verbose is True
    assert cfg.savedir is None


def test_plot_config_context_reverts_on_success():
    cfg = PlotConfig(dpi=150)
    with cfg.context(dpi=999) as ctx:
        assert ctx is cfg
        assert cfg.dpi == 999
    assert cfg.dpi == 150


def test_plot_config_context_reverts_on_exception():
    cfg = PlotConfig(dpi=150)
    with pytest.raises(ValueError):
        with cfg.context(dpi=999):
            raise ValueError("boom")
    assert cfg.dpi == 150


def test_plot_config_to_dict_and_summary_and_repr():
    cfg = PlotConfig(dpi=222)
    d = cfg.to_dict()
    assert d["dpi"] == 222
    assert set(d) == set(cfg._fields())
    summary = cfg.summary()
    assert "dpi" in summary
    assert "222" in summary
    assert repr(cfg) == summary


# ─────────────────────────────────────────────────────────────────────────
# module-level convenience functions
# ─────────────────────────────────────────────────────────────────────────


def test_save_fig_uses_global_config(tmp_path):
    fig, ax, _ = _figure_with_mappable()
    PLOT_CONFIG.verbose = False
    paths = save_fig(ax, tmp_path / "global_out")
    assert paths[0].exists()
    plt.close(fig)


def test_set_fmt_single_and_multiple():
    set_fmt("svg")
    assert PLOT_CONFIG.fmt == "svg"
    set_fmt("svg", "pdf")
    assert PLOT_CONFIG.fmt == ["svg", "pdf"]


def test_set_dpi_and_set_savedir():
    set_dpi("321")
    assert PLOT_CONFIG.dpi == 321
    set_savedir("~/somewhere")
    assert PLOT_CONFIG.savedir.endswith("somewhere")


def test_reset_plot_config_restores_defaults():
    PLOT_CONFIG.dpi = 999
    reset_plot_config()
    assert PLOT_CONFIG.dpi == 150


def test_load_plot_config_missing_explicit_path_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_plot_config(tmp_path / "missing.cfg")


def test_load_plot_config_explicit_path_applies_settings(tmp_path):
    cfg_path = tmp_path / "custom.cfg"
    cfg_path.write_text(
        "[plot]\n"
        "fmt = png\n"
        "base_fmt = png\n"
        "dpi = 444\n"
        "bbox_inches = tight\n"
        "transparent = true\n"
        "facecolor = white\n"
        "savedir = ~/x\n"
        "close_after_save = true\n",
        encoding="utf-8",
    )
    load_plot_config(cfg_path)
    assert PLOT_CONFIG.dpi == 444
    assert PLOT_CONFIG.transparent is True
    assert PLOT_CONFIG.close_after_save is True


def test_load_plot_config_explicit_path_partial_fields(tmp_path):
    cfg_path = tmp_path / "partial.cfg"
    cfg_path.write_text("[plot]\nbase_fmt = svg\n", encoding="utf-8")
    load_plot_config(cfg_path)
    assert PLOT_CONFIG.base_fmt == "svg"


def test_load_plot_config_explicit_path_without_plot_section(tmp_path):
    cfg_path = tmp_path / "no_plot_section.cfg"
    cfg_path.write_text("[other]\nkey = value\n", encoding="utf-8")
    before = PLOT_CONFIG.dpi
    load_plot_config(cfg_path)
    assert PLOT_CONFIG.dpi == before  # nothing applied, no error


def test_load_plot_config_default_search_uses_env_override(monkeypatch, tmp_path):
    monkeypatch.setattr(plot_mod, "_CFG_LOCAL", tmp_path / "missing.cfg")
    monkeypatch.setattr(plot_mod, "_CFG_GLOBAL", tmp_path / "missing2.cfg")
    monkeypatch.setenv("PYCSAMT_DPI", "777")
    load_plot_config()
    assert PLOT_CONFIG.dpi == 777
    monkeypatch.delenv("PYCSAMT_DPI", raising=False)


def test_write_default_config_writes_readable_ini(tmp_path, capsys):
    PLOT_CONFIG.savedir = None
    out = write_default_config(tmp_path / "written.cfg")
    assert out.exists()
    text = out.read_text(encoding="utf-8")
    assert "[plot]" in text
    assert "dpi" in text
    assert "written.cfg" in capsys.readouterr().out


# ─────────────────────────────────────────────────────────────────────────
# add_colorbar / add_polar_colorbar
# ─────────────────────────────────────────────────────────────────────────


def test_add_colorbar_default_side_and_label():
    fig, ax, mappable = _figure_with_mappable()
    cbar = add_colorbar(mappable, ax, label="value")
    assert cbar.ax.get_ylabel() == "value"
    plt.close(fig)


@pytest.mark.parametrize("side", ["left", "top", "bottom"])
def test_add_colorbar_other_sides(side):
    fig, ax, mappable = _figure_with_mappable()
    cbar = add_colorbar(mappable, ax, side=side)
    assert cbar is not None
    plt.close(fig)


def test_add_colorbar_rejects_bad_side():
    fig, ax, mappable = _figure_with_mappable()
    with pytest.raises(ValueError, match="colorbar side must be"):
        add_colorbar(mappable, ax, side="diagonal")
    plt.close(fig)


def test_add_colorbar_tick_format_and_max_ticks():
    fig, ax, mappable = _figure_with_mappable()
    cbar = add_colorbar(mappable, ax, max_ticks=3, tick_format="%.1f")
    assert cbar.formatter is not None
    plt.close(fig)


def test_add_colorbar_max_ticks_none_skips_locator_override():
    fig, ax, mappable = _figure_with_mappable()
    cbar = add_colorbar(mappable, ax, max_ticks=None)
    assert cbar is not None
    plt.close(fig)


def test_add_polar_colorbar_label_and_ticks():
    fig = plt.figure()
    ax = fig.add_subplot(projection="polar")
    mappable = ax.scatter([0, 1], [1, 2], c=[0.1, 0.9])
    cbar = add_polar_colorbar(
        mappable, ax, label="phase", max_ticks=4, tick_format="%.2f"
    )
    assert cbar.ax.get_ylabel() == "phase"
    assert cbar.formatter is not None
    plt.close(fig)


def test_add_polar_colorbar_defaults_without_label_or_format():
    fig = plt.figure()
    ax = fig.add_subplot(projection="polar")
    mappable = ax.scatter([0, 1], [1, 2], c=[0.1, 0.9])
    cbar = add_polar_colorbar(mappable, ax, max_ticks=None)
    assert cbar is not None
    plt.close(fig)
