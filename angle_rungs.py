#!/usr/bin/env python3
"""
Angle Satin Rungs

Select one or more already-built ladder objects (2 rail subpaths followed by
rung subpaths). Every rung is set to the given angle across the column,
measured from square across: 0 = perpendicular to the column, positive leans
the rung clockwise on screen, negative counterclockwise. The angle is
absolute, not added to what is there, so running it again with another value
just re-aims the rungs.

Each rung turns about the middle of the column where it sits now (halfway
between its two rail crossings), and keeps the overshoot it had past each
rail. The column direction there is the average of the two rails' directions
at the points square across from that middle.

Rungs that are left as they were:
 - rungs inside a taper, where the two rails are not running alongside each
   other (more than TAPER_ANGLE_DEG apart), so there is no column direction
   to lean against;
 - rungs whose new line would miss a rail (leaning off the end of the
   column) or would cross another rung. A message says how many of these.
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import inkex
from inkex import EffectExtension, PathElement

from rung_rail_satin import Rail, make_csp, unit_normal, rails_first

TAPER_ANGLE_DEG = 20.0

TAPER, BLOCKED = "taper", "blocked"


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else None


def nearest_on(rail, q):
    """(point, s) of the rail sample nearest to q."""
    best = (None, 0.0, float("inf"))
    for i, (ts, arclen, pts) in enumerate(rail._tables):
        d = np.linalg.norm(pts - q, axis=1)
        j = int(np.argmin(d))
        if d[j] < best[2]:
            best = (pts[j], rail.seg_bounds[i] + arclen[j], d[j])
    return np.array(best[0]), best[1]


def crossing(rail, mid, u, reach):
    """Where the rung's own line crosses the rail, nearest to the rung."""
    res = rail.line_cross(mid, u)
    if res is not None and np.linalg.norm(res[0] - mid) <= reach:
        return np.array(res[0])
    return nearest_on(rail, mid)[0]


def ends(rung):
    return np.array(rung[0][1], dtype=float), np.array(rung[-1][1], dtype=float)


def angled(rung, railA, railB, cos_a, sin_a, taper_cos):
    """The rung re-aimed, or TAPER / BLOCKED when it has to stay as it is."""
    e0, e1 = ends(rung)
    u = unit(e1 - e0)
    if u is None:
        return BLOCKED
    mid = (e0 + e1) / 2.0
    reach = np.linalg.norm(e1 - e0)
    pA = crossing(railA, mid, u, reach)
    pB = crossing(railB, mid, u, reach)
    across = unit(pB - pA)
    if across is None:
        return TAPER

    # Rung ends, and how far each runs past its own rail.
    flipped = np.dot(e1 - e0, across) < 0
    endA, endB = (e1, e0) if flipped else (e0, e1)
    overA = float(np.dot(pA - endA, across))
    overB = float(np.dot(endB - pB, across))

    pivot = (pA + pB) / 2.0
    _, sA = nearest_on(railA, pivot)
    _, sB = nearest_on(railB, pivot)
    tanA = unit(np.array(railA.point_tangent(sA)[1], dtype=float))
    tanB = unit(np.array(railB.point_tangent(sB)[1], dtype=float))
    if tanA is None or tanB is None:
        return TAPER
    if np.dot(tanA, tanB) < 0:
        tanB = -tanB
    if np.dot(tanA, tanB) < taper_cos:
        return TAPER
    along = unit(tanA + tanB)

    square = unit_normal(along)
    aim = np.array([cos_a * square[0] - sin_a * square[1],
                    sin_a * square[0] + cos_a * square[1]])
    if np.dot(aim, across) < 0:
        aim = -aim

    resA = railA.line_cross(pivot, aim, s_guess=sA, _extrapolate=False)
    resB = railB.line_cross(pivot, aim, s_guess=sB, _extrapolate=False)
    if resA is None or resB is None:
        return BLOCKED
    newA = np.array(resA[0]) - aim * overA
    newB = np.array(resB[0]) + aim * overB
    first, last = (newB, newA) if flipped else (newA, newB)
    return [[list(first)] * 3, [list(last)] * 3]


def rungs_cross(r1, r2):
    a0, a1 = ends(r1)
    b0, b1 = ends(r2)

    def side(p, q, r):
        return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])

    return (side(a0, a1, b0) * side(a0, a1, b1) < 0
            and side(b0, b1, a0) * side(b0, b1, a1) < 0)


class AngleRungs(EffectExtension):
    def add_arguments(self, pars):
        pars.add_argument("--angle", type=float, default=30.0)

    def effect(self):
        sel = [e for e in self.svg.selection if isinstance(e, PathElement)]
        if not sel:
            raise inkex.AbortExtension("Select one or more satin column ladder objects.")

        a = math.radians(self.options.angle)
        cos_a, sin_a = math.cos(a), math.sin(a)
        taper_cos = math.cos(math.radians(TAPER_ANGLE_DEG))
        kept = []

        for el in sel:
            csp = rails_first([list(sp) for sp in el.path.to_absolute().to_superpath()])
            if len(csp) < 3:
                raise inkex.AbortExtension(
                    f"{el.get_id()} needs two rail subpaths plus at least one rung.")
            railA, railB = Rail(csp[0]), Rail(csp[1])

            old = csp[2:]
            new = [angled(r, railA, railB, cos_a, sin_a, taper_cos) for r in old]
            stay = sum(1 for r in new if r == BLOCKED)
            rungs = [o if isinstance(r, str) else r for o, r in zip(old, new)]
            changed = [not isinstance(r, str) for r in new]

            # A re-aimed rung that now crosses another rung goes back to how
            # it was. Repeat, since putting one back can free its neighbor.
            while True:
                bad = [i for i in range(len(rungs)) if changed[i] and any(
                    j != i and rungs_cross(rungs[i], rungs[j]) for j in range(len(rungs)))]
                if not bad:
                    break
                for i in bad:
                    rungs[i] = old[i]
                    changed[i] = False
                stay += len(bad)

            if stay:
                kept.append(f"{el.get_id()}: {stay}")
            el.path = make_csp(csp[:2] + rungs).to_path()

        if kept:
            inkex.errormsg(
                "Some rungs were left as they were, because at this angle they "
                "would run off the end of the column or cross another rung. "
                "Count per object: " + "; ".join(kept))


if __name__ == "__main__":
    AngleRungs().run()
