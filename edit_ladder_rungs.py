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
 - Remove all: removes every rung, keeping the two rails.
 - Fit to rails: every rung keeps its position and angle, and its ends are
   moved to just past each rail (RUNG_OVERSHOOT_MM, 0.125 mm, on each
   side, the same overshoot new rungs get). Rungs whose line
   does not cross both rails are left as they are. All rungs are fitted,
   including the protected ends.

The first and last PROTECTED_END_COUNT rungs (in order along the rails) are
never touched by Add or Decimate, so tapers stay as generated. Both modes can be
run repeatedly.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import inkex
from inkex import EffectExtension, ShapeElement, PathElement

import subpath_numbers
from angle_rungs import AngleRungs
from rung_rail_satin import (
    Rail, make_csp, unit_normal, apply_rungs, resolve_path, stroke_width, rails_first,
    STROKE_DEFAULT_MULT, overshoot_uu,
    normalize_rail_directions, measure_A, pair_rails, rung_between,
    pair_progress, pair_of_rung, longer_first, clear_rung_ends, column_outline, rung_span,
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


def fitted(rung, railA, railB, over):
    """The rung re-ended at its own line's rail crossings plus overshoot, or None."""
    e0, e1 = np.array(rung[0][1], dtype=float), np.array(rung[-1][1], dtype=float)
    length = np.linalg.norm(e1 - e0)
    if length < 1e-9:
        return None
    u = (e1 - e0) / length
    mid = (e0 + e1) / 2.0
    reach = max(3.0 * length, 1.0)
    hits = []
    for rail in (railA, railB):
        res = rail.line_cross(mid, u)
        if res is None or np.linalg.norm(np.array(res[0]) - mid) > reach:
            return None
        hits.append(np.array(res[0], dtype=float))
    # order the crossings along the rung's own direction, so its ends keep their order
    pa, pb = sorted(hits, key=lambda q: float(np.dot(q - mid, u)))
    width = float(np.dot(pb - pa, u))
    if width < 1e-9:
        return None
    return [[list(pa - u * over)] * 3, [list(pb + u * over)] * 3]


class EditLadderRungs(EffectExtension):
    def add_arguments(self, pars):
        pars.add_argument("--mode", default="add")
        pars.add_argument("--angle", type=float, default=30.0)

    def effect(self):
        mode = self.options.mode
        if mode == "angle":
            tool = AngleRungs()
            tool.svg = self.svg
            tool.options = self.options
            tool.effect()
            return
        if mode in ("number", "unnumber"):
            paths = [e for e in self.svg.selection.get(ShapeElement) if isinstance(e, PathElement)]
            if mode == "unnumber":
                subpath_numbers.unnumber(self.svg, paths)
                return
            if not paths:
                raise inkex.AbortExtension("Select one or more paths to number their subpaths.")
            for note in subpath_numbers.number(self.svg, paths):
                inkex.errormsg(note)
            return
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
        csp = rails_first([list(sp) for sp in el.path.to_absolute().to_superpath()])

        if len(csp) < 2:
            raise inkex.AbortExtension("Need a path with two rail subpaths.")
        if len(csp) == 2:
            if self.options.mode in ("decimate", "fit", "remove"):
                self.note(el, "No rungs to " + {"decimate": "decimate.", "fit": "fit.", "remove": "remove."}[self.options.mode])
                return
            apply_rungs(el)
            return

        if self.options.mode == "remove":
            # every rung goes, the two rails stay
            el.path = make_csp(csp[:2]).to_path()
            return

        railA, railB = Rail(csp[0]), Rail(csp[1])
        over = overshoot_uu(el)
        rungs = csp[2:]

        if self.options.mode == "fit":
            new = [fitted(r, railA, railB, over) for r in rungs]
            missed = sum(1 for r in new if r is None)
            if missed:
                self.note(el, f"{missed} rung(s) left as they were: their line does not cross both rails.")
            new = clear_rung_ends([o if r is None else r for o, r in zip(rungs, new)], railA, railB)
            el.path = make_csp(csp[:2] + new).to_path()
            return
        # Order rungs along the column by where each sits in the both-rails
        # pairing (pair_rails), not by its nearest point on one rail: on a
        # curl a rung's end can sit right next to a different turn of a rail.
        c0, c1 = normalize_rail_directions(*longer_first(csp[0], csp[1]))
        oA, oB = Rail(c0), Rail(c1)
        pairs = pair_rails(oA, oB, measure_A(oA, oB))
        outline = column_outline(oA, oB)
        progress = pair_progress(pairs)
        place = {id(r): progress[pair_of_rung(pairs, r)] for r in rungs}
        rungs.sort(key=lambda r: place[id(r)])
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
            # A new rung goes halfway along the column between two
            # neighbors, at the same both-rails pairing the builder uses, so
            # it is as near square to both rails as the column allows.
            new_rungs, added = [], []
            fits = None
            for i, r in enumerate(rungs):
                new_rungs.append(r)
                if P - 1 <= i < n - P:
                    mid = (place[id(r)] + place[id(rungs[i + 1])]) / 2.0
                    if fits is None:
                        fits = {}
                        for j in range(len(pairs)):
                            sp = rung_span(outline, pairs[j][2], pairs[j][3])
                            if sp is not None:
                                fits[j] = sp
                    if not fits:
                        continue
                    k = min(fits, key=lambda j: abs(progress[j] - mid))
                    rung = rung_between(*fits[k], over)
                    if rung is not None:
                        added.append(rung)
                        new_rungs.append(rung)

        if self.options.mode != "decimate" and added:
            # only the new rungs: existing ones are the user's to place
            swap = {id(a): c for a, c in zip(added, clear_rung_ends(added, oA, oB))}
            new_rungs = [swap.get(id(r), r) for r in new_rungs]
        el.path = make_csp(csp[:2] + new_rungs).to_path()


if __name__ == "__main__":
    EditLadderRungs().run()
