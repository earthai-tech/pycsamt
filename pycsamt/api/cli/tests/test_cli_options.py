from __future__ import annotations

import click
from click.testing import CliRunner

from pycsamt.api.cli.options import (
    common_options,
    format_option,
    fresh_option,
    n_jobs_option,
    no_cache_option,
    no_color_option,
    output_dir_option,
    overwrite_option,
    survey_option,
    verbose_option,
)


def test_verbose_option_counts_occurrences():
    @click.command()
    @verbose_option
    def cmd(verbose):
        click.echo(f"verbose={verbose}")

    result = CliRunner().invoke(cmd, ["-vv"])
    assert result.exit_code == 0
    assert "verbose=2" in result.output


def test_no_color_option_flag():
    @click.command()
    @no_color_option
    def cmd(no_color):
        click.echo(f"no_color={no_color}")

    result = CliRunner().invoke(cmd, ["--no-color"])
    assert "no_color=True" in result.output
    result_default = CliRunner().invoke(cmd, [])
    assert "no_color=False" in result_default.output


def test_output_dir_option_default_and_override(tmp_path):
    @click.command()
    @output_dir_option
    def cmd(output_dir):
        click.echo(f"output_dir={output_dir}")

    result = CliRunner().invoke(cmd, ["-o", str(tmp_path)])
    assert str(tmp_path) in result.output


def test_format_option_choice_and_default():
    @click.command()
    @format_option
    def cmd(output_format):
        click.echo(f"format={output_format}")

    assert "format=text" in CliRunner().invoke(cmd, []).output
    assert "format=json" in CliRunner().invoke(cmd, ["-f", "json"]).output


def test_format_option_rejects_invalid_choice():
    @click.command()
    @format_option
    def cmd(output_format):
        click.echo("ok")

    result = CliRunner().invoke(cmd, ["-f", "xml"])
    assert result.exit_code != 0


def test_overwrite_option_flag():
    @click.command()
    @overwrite_option
    def cmd(overwrite):
        click.echo(f"overwrite={overwrite}")

    assert "overwrite=True" in CliRunner().invoke(cmd, ["--overwrite"]).output


def test_survey_option_accepts_existing_path(tmp_path):
    @click.command()
    @survey_option
    def cmd(survey_path):
        click.echo(f"survey_path={survey_path}")

    result = CliRunner().invoke(cmd, ["-S", str(tmp_path)])
    assert str(tmp_path) in result.output


def test_survey_option_defaults_to_none():
    @click.command()
    @survey_option
    def cmd(survey_path):
        click.echo(f"survey_path={survey_path}")

    result = CliRunner().invoke(cmd, [])
    assert "survey_path=None" in result.output


def test_fresh_option_flag():
    @click.command()
    @fresh_option
    def cmd(fresh):
        click.echo(f"fresh={fresh}")

    assert "fresh=True" in CliRunner().invoke(cmd, ["--fresh"]).output


def test_n_jobs_option_default_and_range():
    @click.command()
    @n_jobs_option
    def cmd(n_jobs):
        click.echo(f"n_jobs={n_jobs}")

    assert "n_jobs=1" in CliRunner().invoke(cmd, []).output
    assert "n_jobs=4" in CliRunner().invoke(cmd, ["-j", "4"]).output


def test_n_jobs_option_rejects_below_minimum():
    @click.command()
    @n_jobs_option
    def cmd(n_jobs):
        click.echo("ok")

    result = CliRunner().invoke(cmd, ["-j", "0"])
    assert result.exit_code != 0


def test_no_cache_option_flag():
    @click.command()
    @no_cache_option
    def cmd(no_cache):
        click.echo(f"no_cache={no_cache}")

    assert "no_cache=True" in CliRunner().invoke(cmd, ["--no-cache"]).output


def test_common_options_attaches_all_four():
    @click.command()
    @common_options
    def cmd(verbose, no_color, output_format, output_dir):
        click.echo(
            f"verbose={verbose} no_color={no_color} "
            f"format={output_format} output_dir={output_dir}"
        )

    result = CliRunner().invoke(cmd, ["-v", "--no-color", "-f", "csv"])
    assert result.exit_code == 0
    assert "verbose=1" in result.output
    assert "no_color=True" in result.output
    assert "format=csv" in result.output


def test_option_docstrings_are_assigned():
    assert "verbose" in verbose_option.__doc__.lower()
    assert "no-color" in no_color_option.__doc__
    assert "output-dir" in output_dir_option.__doc__
    assert "format" in format_option.__doc__.lower()
    assert "overwrite" in overwrite_option.__doc__.lower()
    assert "jobs" in n_jobs_option.__doc__.lower()
    assert "no-cache" in no_cache_option.__doc__
