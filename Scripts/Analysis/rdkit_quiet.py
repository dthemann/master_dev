"""Suppress RDKit's harmless 2D/3D coordinate warning, keeping all other logs.

RDKit emits one warning per molecule for SDFs tagged 2D but carrying real
(non-zero Z) coordinates:

    Warning: molecule is tagged as 2D, but at least one Z coordinate is not
    zero. Marking the mol as 3D.

It's harmless (RDKit just promotes the mol to 3D, which is what we want) but
floods the console when reading thousands of docked poses. ``warnings.filterwarnings``
does NOT catch it — it's emitted by RDKit's C++ logger, not Python's ``warnings``.

This routes RDKit's C++ logs through Python's ``sys.stderr`` and drops ONLY this
one message; every other RDKit warning/error (unparseable elements, sanitization
failures, ...) still gets through. Import this module and call
``silence_rdkit_2d3d_warning()`` once per process before reading molecules.
"""

from __future__ import annotations

import sys

_SUPPRESS_SUBSTR = "is tagged as 2D, but at least one Z coordinate is not zero"


class _RDKitStderrFilter:
    """stderr wrapper that swallows only the harmless 2D/3D RDKit warning."""

    def __init__(self, stream):
        self._stream = stream

    def write(self, msg):
        if _SUPPRESS_SUBSTR not in msg:
            self._stream.write(msg)
        return len(msg)

    def flush(self):
        self._stream.flush()

    def __getattr__(self, name):
        return getattr(self._stream, name)


def silence_rdkit_2d3d_warning() -> None:
    """Install the 2D/3D-warning filter on this process's stderr (idempotent).

    Call once per process that loads molecules. For ``multiprocessing`` /
    ``ProcessPoolExecutor`` workers, calling this at module import is enough:
    fork inherits the parent's installed filter and spawn re-imports the module.
    """
    from rdkit import rdBase

    rdBase.LogToPythonStderr()
    if not isinstance(sys.stderr, _RDKitStderrFilter):
        sys.stderr = _RDKitStderrFilter(sys.stderr)
