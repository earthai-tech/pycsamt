"""Validate that a pytest JUnit report represents a complete passing run."""

from __future__ import annotations

import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def passing_report(path: Path) -> bool:
    """Return True only for a non-empty report with no failures or errors."""

    try:
        root = ET.parse(path).getroot()
    except (OSError, ET.ParseError):
        return False

    suites = (
        [root] if root.tag == "testsuite" else list(root.findall("testsuite"))
    )
    if not suites:
        return False
    tests = sum(int(suite.get("tests", "0")) for suite in suites)
    failures = sum(int(suite.get("failures", "0")) for suite in suites)
    errors = sum(int(suite.get("errors", "0")) for suite in suites)
    return tests > 0 and failures == 0 and errors == 0


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 1:
        print("usage: verify_pytest_junit.py REPORT.xml", file=sys.stderr)
        return 2
    report = Path(args[0])
    if not passing_report(report):
        print(
            f"incomplete or failing pytest report: {report}", file=sys.stderr
        )
        return 1
    print(f"verified complete passing pytest report: {report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
