#!/usr/bin/env python3
"""
Edit Ladder Rungs (Ink/Stitch satin column helper)

Select one ladder object (2 rail subpaths + rungs), or two separate rail
paths, and pick a mode:

 - Add: with no rungs yet, builds the full set (same as the old Generate
   Ladder Rungs). With rungs already present, inserts one new rung midway
   between each adjacent pair of unprotected rungs (halves the spacing).
 - Decimate: removes every other unprotected rung (keeps the 1st, 3rd, ...).

The first and last PROTECTED_END_COUNT rungs (in order along the rails) are
never touched by either mode, so tapers stay as generated. Both modes can be
run repeatedly.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import inkex
from inkex import EffectExtension, ShapeElement

from rung_rail_satin import (
    Rail, make_csp, unit_normal, apply_rungs, resolve_path, stroke_width,
    STROKE_DEFAULT_MULT,
)

PROTECTED_END_COUNT = 4


def s_of_nearest(rail, q):
    best_s, best_d = 0.0, float("inf")
    for i, (ts, arclen, pts) in enumerate(rail._tables):
        d = np.linalg.norm(pts - q, axis=1)
        j = int(np.argmin(d))
        if d[j] < best_d:
            best_d, best_s = d[j], rail.seg_bounds[i] + arclen[j]
    return best_s


def a_end_s(rail, rung):
    """Arc position on `rail` of whichever rung endpoint lies nearer to it."""
    cands = [s_of_nearest(rail, np.array(rung[k][1])) for k in (0, -1)]
    def gap(k):
        q = np.array(rung[k][1])
        return min(np.linalg.norm(t[2] - q, axis=1).min() for t in rail._tables)
    return cands[0] if gap(0) <= gap(-1) else cands[1]


class EditLadderRungs(EffectExtension):
    def add_arguments(self, pars):
        pars.add_argument("--mode", default="add")

    def effect(self):
        sel = list(self.svg.selection.get(ShapeElement))
        el = resolve_path(sel)
        csp = [list(sp) for sp in el.path.to_absolute().to_superpath()]

        if len(csp) < 2:
            raise inkex.AbortExtension("Need a path with two rail subpaths.")
        if len(csp) == 2:
            if self.options.mode == "decimate":
                inkex.errormsg("No rungs to decimate.")
                return
            apply_rungs(el)
            return

        railA, railB = Rail(csp[0]), Rail(csp[1])
        rungs = csp[2:]
        rungs.sort(key=lambda r: s_of_nearest(
            railA, (np.array(r[0][1]) + np.array(r[-1][1])) / 2.0))
        n = len(rungs)
        P = PROTECTED_END_COUNT

        if self.options.mode == "decimate":
            keep = [r for i, r in enumerate(rungs)
                    if i < P or i >= n - P or (i - P) % 2 == 0]
            if len(keep) == n:
                inkex.errormsg("Nothing to decimate: no unprotected rungs to remove.")
                return
            new_rungs = keep
        else:
            if n < 2 * P:
                inkex.errormsg("Nothing to add: no gap between the protected end rungs.")
                return
            B = stroke_width(el)
            new_rungs = []
            for i, r in enumerate(rungs):
                new_rungs.append(r)
                if P - 1 <= i < n - P:
                    nxt = rungs[i + 1]
                    s_a = a_end_s(railA, r)
                    s_b = a_end_s(railA, nxt)
                    pA, tanA = railA.point_tangent((s_a + s_b) / 2.0)
                    res = railB.line_cross(np.array(pA), unit_normal(tanA))
                    if res is None:
                        continue
                    pB = res[0]
                    dvec = np.array(pB) - np.array(pA)
                    length = np.linalg.norm(dvec)
                    if length < 1e-9:
                        continue
                    u = dvec / length
                    endA = np.array(pA) - u * (STROKE_DEFAULT_MULT * B)
                    endB = np.array(pB) + u * (STROKE_DEFAULT_MULT * B)
                    new_rungs.append([[list(endA)] * 3, [list(endB)] * 3])

        el.path = make_csp(csp[:2] + new_rungs).to_path()


if __name__ == "__main__":
    EditLadderRungs().run()
