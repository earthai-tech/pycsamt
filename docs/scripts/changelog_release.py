#!/usr/bin/env python
"""Assemble docs/changelog.d/ fragments into a versioned changelog entry.

Usage
-----
    python docs/scripts/changelog_release.py --version 2.2.0 --date 2026-08-15

Groups every fragment in ``docs/changelog.d/`` by section, inserts the
assembled block at the top of the matching ``docs/source/changelog/vN.rst``
series file, and deletes the consumed fragments (git history keeps them).
See ``docs/changelog.d/README.rst`` for the fragment format.
"""
from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from datetime import date
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
FRAGMENTS_DIR = REPO_ROOT / "docs" / "changelog.d"
CHANGELOG_DIR = REPO_ROOT / "docs" / "source" / "changelog"

# Badge order matches the `.. |Badge| replace::` definitions in
# docs/source/conf.py's rst_epilog.
BADGE_ORDER = [
    "Feature", "New", "Fix", "Enhancement", "Perf", "Breaking",
    "Security", "API Change", "Deprecated", "Docs", "Build", "Tests",
]
# fragment type -> (changelog section, badge)
TYPE_INFO = {
    "feature": ("Added", "Feature"),
    "new": ("Added", "New"),
    "enhancement": ("Changed", "Enhancement"),
    "perf": ("Changed", "Perf"),
    "fix": ("Fixed", "Fix"),
    "api": ("Changed", "API Change"),
    "deprecated": ("Changed", "Deprecated"),
    "breaking": ("Changed", "Breaking"),
    "removed": ("Removed", "Breaking"),
    "security": ("Security", "Security"),
    "docs": ("Docs & tooling", "Docs"),
    "build": ("Docs & tooling", "Build"),
    "tests": ("Docs & tooling", "Tests"),
}
SECTION_ORDER = ["Added", "Fixed", "Changed", "Removed", "Security", "Docs & tooling"]

FRAGMENT_RE = re.compile(r"^(?P<ref>\+?[^.]+)\.(?P<type>[a-z]+)\.rst$")
VERSION_RE = re.compile(r"^\d+\.\d+\.\d+([a-z]+\d+)?$")
VERSION_ANCHOR_RE = re.compile(r"^\.\. _changelog-\d.*$", re.MULTILINE)


@dataclass
class Fragment:
    path: Path
    ref: str
    type: str
    body: str


def _sort_key(ref: str):
    if ref.startswith("+"):
        return (1, ref)
    try:
        return (0, int(ref))
    except ValueError:
        return (1, ref)


def load_fragments() -> list[Fragment]:
    if not FRAGMENTS_DIR.is_dir():
        raise SystemExit(f"fragments directory not found: {FRAGMENTS_DIR}")

    fragments = []
    for path in sorted(FRAGMENTS_DIR.iterdir()):
        if not path.is_file() or path.name == "README.rst":
            continue
        match = FRAGMENT_RE.match(path.name)
        if not match:
            raise SystemExit(
                f"changelog fragment {path.name!r} doesn't match "
                "'<ref>.<type>.rst' -- see docs/changelog.d/README.rst"
            )
        type_ = match.group("type")
        if type_ not in TYPE_INFO:
            raise SystemExit(
                f"changelog fragment {path.name!r} has unknown type "
                f"{type_!r}; valid types: {', '.join(sorted(TYPE_INFO))}"
            )
        body = path.read_text(encoding="utf-8").strip()
        if not body:
            raise SystemExit(f"changelog fragment {path.name!r} is empty")
        fragments.append(
            Fragment(path=path, ref=match.group("ref"), type=type_, body=body)
        )

    fragments.sort(key=lambda f: _sort_key(f.ref))
    return fragments


def _badges_for(fragments: list[Fragment]) -> str:
    present = {TYPE_INFO[f.type][1] for f in fragments}
    ordered = [b for b in BADGE_ORDER if b in present]
    return " ".join(f"|{b}|" for b in ordered)


def build_block(version: str, date_str: str, summary: str, fragments: list[Fragment]) -> str:
    anchor = version.replace(".", "-")
    heading = f"{version} {_badges_for(fragments)}"

    lines = [
        f".. _changelog-{anchor}:",
        "",
        heading,
        "-" * len(heading),
        "",
        f"*Released {date_str}.*",
        "",
        summary,
        "",
    ]

    by_section: dict[str, list[Fragment]] = {}
    for fragment in fragments:
        section = TYPE_INFO[fragment.type][0]
        by_section.setdefault(section, []).append(fragment)

    for section in SECTION_ORDER:
        entries = by_section.get(section)
        if not entries:
            continue
        lines.append(section)
        lines.append("~" * len(section))
        lines.append("")
        for fragment in entries:
            badge = TYPE_INFO[fragment.type][1]
            body_lines = fragment.body.splitlines()
            lines.append(f"* |{badge}| {body_lines[0]}")
            for extra in body_lines[1:]:
                lines.append(f"  {extra}" if extra else "")
        lines.append("")

    lines.append("----")
    lines.append("")
    # A trailing "\n" beyond what "\n".join gives turns the final "----", ""
    # pair into "----\n\n" -- a blank line is needed before whatever follows
    # (the next existing entry, or end of file for a brand-new series).
    return "\n".join(lines) + "\n"


def _series_header(series: str) -> str:
    title = f"{series}.x changelog"
    return (
        f".. _changelog-{series}:\n\n"
        f"{title}\n"
        f"{'=' * len(title)}\n\n"
        "Every notable change in this series, newest first. See the\n"
        ":ref:`Changelog index <changelog>` for the badge legend and links\n"
        "to other series.\n\n"
        "----\n\n"
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", required=True, help="e.g. 2.2.0")
    parser.add_argument(
        "--date", default=date.today().isoformat(), help="YYYY-MM-DD, default today"
    )
    parser.add_argument(
        "--series", help="e.g. v2; default derived from --version's major number"
    )
    parser.add_argument(
        "--summary", help="release summary paragraph; default a TODO placeholder"
    )
    args = parser.parse_args(argv)

    if not VERSION_RE.match(args.version):
        raise SystemExit(f"--version {args.version!r} doesn't look like X.Y.Z or X.Y.ZrcN")
    try:
        date.fromisoformat(args.date)
    except ValueError:
        raise SystemExit(f"--date {args.date!r} must be YYYY-MM-DD")

    series = args.series or f"v{args.version.split('.')[0]}"
    series_path = CHANGELOG_DIR / f"{series}.rst"

    fragments = load_fragments()
    if not fragments:
        raise SystemExit(f"no fragments found in {FRAGMENTS_DIR} -- nothing to release")

    summary = args.summary or "*TODO -- write the one-paragraph release summary here.*"
    block = build_block(args.version, args.date, summary, fragments)
    anchor_label = f".. _changelog-{args.version.replace('.', '-')}:"

    created = not series_path.exists()
    if created:
        new_text = _series_header(series) + block
    else:
        text = series_path.read_text(encoding="utf-8")
        if anchor_label in text:
            raise SystemExit(f"{series_path} already has an entry for {args.version}")
        match = VERSION_ANCHOR_RE.search(text)
        if match:
            new_text = text[: match.start()] + block + text[match.start() :]
        else:
            new_text = text.rstrip("\n") + "\n\n----\n\n" + block

    series_path.write_text(new_text, encoding="utf-8")
    for fragment in fragments:
        fragment.path.unlink()

    rel_series = series_path.relative_to(REPO_ROOT)
    print(f"{'Created' if created else 'Updated'} {rel_series} with {args.version}")
    print(f"Consumed {len(fragments)} fragment(s) from {FRAGMENTS_DIR.relative_to(REPO_ROOT)}")
    print("Next steps:")
    print(f"  1. Hand-edit the release summary in {rel_series}.")
    print(f"  2. Write docs/source/release_notes/v{args.version}.rst.")
    if created:
        print(f"  3. Add 'changelog/{series}' to the toctree in docs/source/changelog.rst.")


if __name__ == "__main__":
    sys.exit(main())
