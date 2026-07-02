# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for the ``pycsamt map`` command group."""

from __future__ import annotations

import csv
import json
from io import StringIO

from click.testing import CliRunner

from pycsamt.cli import main
from pycsamt.cli.tests.conftest import make_fake_sites


def test_map_group_help(runner: CliRunner) -> None:
    result = runner.invoke(main, ["map", "--help"])

    assert result.exit_code == 0
    assert "stations" in result.output
    assert "plot" in result.output


def test_map_registered_on_root_help(runner: CliRunner) -> None:
    result = runner.invoke(main, ["--help"])

    assert result.exit_code == 0
    assert "map" in result.output


def test_map_stations_help(runner: CliRunner) -> None:
    result = runner.invoke(main, ["map", "stations", "--help"])

    assert result.exit_code == 0
    assert "--drop-missing" in result.output
    assert "--format" in result.output
    assert "--output" in result.output


def test_map_plot_help(runner: CliRunner) -> None:
    result = runner.invoke(main, ["map", "plot", "--help"])

    assert result.exit_code == 0
    assert "--output" in result.output
    assert "--no-label" in result.output
    assert "--dpi" in result.output


def test_map_stations_json(runner: CliRunner, monkeypatch) -> None:
    fake = make_fake_sites(3)
    monkeypatch.setattr(
        "pycsamt.cli.commands.map.stations._get_sites",
        lambda *a, **kw: fake,
    )

    result = runner.invoke(
        main,
        ["map", "stations", "--survey", ".", "--format", "json"],
    )

    assert result.exit_code == 0
    data = json.loads(result.output)
    assert len(data) == 3
    assert data[0]["station"] == "S01"
    assert data[0]["lat"] == 0.0
    assert data[0]["lon"] == 100.0


def test_map_stations_csv(runner: CliRunner, monkeypatch) -> None:
    fake = make_fake_sites(2)
    monkeypatch.setattr(
        "pycsamt.cli.commands.map.stations._get_sites",
        lambda *a, **kw: fake,
    )

    result = runner.invoke(
        main,
        ["map", "stations", "--survey", ".", "--format", "csv"],
    )

    assert result.exit_code == 0
    rows = list(csv.DictReader(StringIO(result.output)))
    assert len(rows) == 2
    assert rows[1]["station"] == "S02"
    assert rows[1]["lat"] == "1.0"
    assert rows[1]["lon"] == "101.0"

