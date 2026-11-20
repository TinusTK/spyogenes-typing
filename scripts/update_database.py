#!/usr/bin/env python3
"""
update_database.py - rebuild database/ for spyogenes-typing from PubMLST
(and, with --references or when missing, the NCBI reference parts), then
write a change summary.

    python scripts/update_database.py --db database --summary CHANGES_DB.md

The code lives in tool/spyogenes_typing_update.py so that a pip install
carries it (`spyogenes-typing --setup` uses it to build a user's own
database). This file only calls it; options and exit codes are unchanged
(0 = updated, 3 = no changes, other = error).
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tool"))

from spyogenes_typing_update import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
