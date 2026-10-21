#!/usr/bin/env python3
"""
Edit Ladder Rungs (Ink/Stitch satin column helper)

Select one or more ladder objects (2 rail subpaths + rungs), and/or separate
rail paths in pairs (Shift+click order: rails 1+2 are one ladder, 3+4 the
next, ...), and pick a mode:

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
from inkex import EffectExtension, ShapeElement, PathElement

from rung_rail_satin import (
    Rail, make_csp, unit_normal, apply_rungs, resolve_path, stroke_width,
    STROKE_DEFAULT_MULT,
)

PROTECTED_END_COUNT = 4
PAIR_GAP_RATIO = 0.5


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


def name_of(el):
    return el.label or el.get_id()


def rail_points(el):
    csp = el.path.to_absolute().to_superpath()
    rail = Rail(csp[0])
    return rail, np.vstack([t[2] for t in rail._tables])


def plausible_pair(a, b):
    """False when two rails picked as a pair sit far apart for their length:
    the average gap from each rail's points to the other rail is more than
    PAIR_GAP_RATIO times the longer rail's length."""
    ra, pa = rail_points(a)
    rb, pb = rail_points(b)
    def mean_gap(p, q):
        return float(np.mean([np.min(np.linalg.norm(q - x, axis=1)) for x in p]))
    gap = (mean_gap(pa, pb) + mean_gap(pb, pa)) / 2.0
    return gap <= PAIR_GAP_RATIO * max(ra.length, rb.length)


class EditLadderRungs(EffectExtension):
    def add_arguments(self, pars):
        pars.add_argument("--mode", default="add")

    def effect(self):
        # Selection order is click order (Shift+click). A path with 2+
        # subpaths is a ladder on its own; single-subpath rails pair up in
        # the order they were selected: 1+2, 3+4, ...
        sel = list(self.svg.selection.get(ShapeElement))
        if not sel:
            raise inkex.AbortExtension(
                "Select one or more ladders, or rail paths in pairs (Shift+click rail 1, rail 2, rail 3, rail 4, ...).")
        ladders, rails = [], []
        for e in sel:
            if not isinstance(e, PathElement):
                if not hasattr(e, "to_path_element"):
                    raise inkex.AbortExtension(f"{name_of(e)} is not a path or a shape.")
                path_el = e.to_path_element()
                e.replace_with(path_el)
                e = path_el
            n = len(e.path.to_absolute().to_superpath())
            (ladders if n >= 2 else rails).append(e)
        if len(rails) % 2:
            raise inkex.AbortExtension(
                f"{len(rails)} separate rails are selected; rails pair up in selection order, so the count must be even. "
                f"Unpaired: {name_of(rails[-1])}. Nothing was changed.")
        pairs = [rails[i:i + 2] for i in range(0, len(rails), 2)]
        bad = [p for p in pairs if not plausible_pair(*p)]
        if bad:
            raise inkex.AbortExtension(
                "These rails were paired by selection order but are far apart for their length, "
                "so the click order is probably off: "
                + "; ".join(f"{name_of(a)} + {name_of(b)}" for a, b in bad) + ". Nothing was changed.")
        targets = ladders + [resolve_path(p) for p in pairs]
        self.multi = len(targets) > 1
        for el in targets:
            self.edit_one(el)

    def note(self, el, msg):
        inkex.errormsg(f"{name_of(el)}: {msg}" if getattr(self, "multi", False) else msg)

    def edit_one(self, el):
        csp = [list(sp) for sp in el.path.to_absolute().to_superpath()]

        if len(csp) < 2:
            raise inkex.AbortExtension("Need a path with two rail subpaths.")
        if len(csp) == 2:
            if self.options.mode == "decimate":
                self.note(el, "No rungs to decimate.")
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
                self.note(el, "Nothing to decimate: no unprotected rungs to remove.")
                return
            new_rungs = keep
        else:
            if n < 2 * P:
                self.note(el, "Nothing to add: no gap between the protected end rungs.")
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
