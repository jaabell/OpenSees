"""Reader of the SWMR HDF5 sink written by `recorder Monitor`.

Layout (ladruno::MonitorSink, SRC/recorder/LadrunoMonitorSink.h):

    /            attrs: FORMAT="ladruno-monitor", FORMAT_VERSION, GENERATOR, units?
    COLUMNS      [nCols]      variable-length strings (self-describing channels)
    STEP         [T]          int32   (commitTag per frame)
    TIME         [T]          float64 (pseudo-time per frame)
    FRAMES       [T, nCols]   float64 (one row per emitted frame)

The file is opened in SWMR read mode and refreshed on each read, so it can be
read while another process is still writing it. A file that never armed SWMR
is read normally. Depends on h5py only.
"""

from __future__ import annotations

import os

# Read a file another process may still hold open: disable libhdf5 locking.
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import h5py  # noqa: E402


def _s(v):
    return v.decode() if isinstance(v, (bytes, bytearray)) else v


class MonitorReader:
    """Stateless reader: opens the file per call, so a growing SWMR sink is
    always seen at its current length (no stale handle, no long-lived lock)."""

    def __init__(self, path: str):
        self.path = path

    def _open(self):
        # a missing sink raises FileNotFoundError whatever the h5py version
        if not os.path.exists(self.path):
            raise FileNotFoundError(self.path)
        try:
            return h5py.File(self.path, "r", swmr=True)
        except Exception:
            # not (yet) an SWMR file, or the writer has not armed it: plain read
            return h5py.File(self.path, "r")

    @staticmethod
    def _refresh(*dsets):
        for d in dsets:
            try:
                d.id.refresh()
            except Exception:
                pass

    def meta(self) -> dict:
        with self._open() as f:
            if _s(f.attrs.get("FORMAT")) != "ladruno-monitor":
                raise ValueError(f"{self.path}: not a ladruno-monitor sink "
                                 f"(FORMAT={_s(f.attrs.get('FORMAT'))!r})")
            cols = [_s(c) for c in f["COLUMNS"][:]]
            n = int(f["STEP"].shape[0])
            return {
                "columns": cols,
                "nColumns": len(cols),
                "nframes": n,
                "format": _s(f.attrs.get("FORMAT")),
                "format_version": int(f.attrs.get("FORMAT_VERSION", 0)),
                "generator": _s(f.attrs.get("GENERATOR", "")),
                "units": _s(f.attrs.get("units", "")),
                "path": os.path.abspath(self.path),
            }

    def frames_since(self, since: int = 0) -> dict:
        """Return frames [since:n). `n` is the current total, so a caller that
        polls advances its cursor to `n` for the next read."""
        with self._open() as f:
            step_ds, time_ds, frames_ds = f["STEP"], f["TIME"], f["FRAMES"]
            self._refresh(step_ds, time_ds, frames_ds)
            n = int(step_ds.shape[0])
            since = max(0, int(since))
            if since >= n:
                return {"since": since, "n": n, "step": [], "t": [], "rows": []}
            step = step_ds[since:n].astype("int64").tolist()
            t = time_ds[since:n].astype("float64").tolist()
            rows = frames_ds[since:n, :].astype("float64").tolist()
            return {"since": since, "n": n, "step": step, "t": t, "rows": rows}

    def all(self) -> dict:
        """meta + every frame."""
        m = self.meta()
        fr = self.frames_since(0)
        return {**m, **fr}
