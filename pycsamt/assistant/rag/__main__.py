# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Entry point: ``python -m pycsamt.assistant.rag <command>``."""

from __future__ import annotations

import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
