#!/usr/bin/env python3
"""
Touch Rails at Crossing (Ink/Stitch satin column helper)

Follow-on to Weave Paths: after a woven strand has been converted to Ink/
Stitch satin columns, the "under" piece's cut end just stops at a gap
sized to the covering stroke's width - it doesn't actually reach the
"over" column's rails. This tool closes that gap by extending (or
trimming) the under piece's two rail tips, independently, until each one
touches the near rail of the over column, matching a hand-edit done and
confirmed as the target look (rails touching, not floating in a gap).

Select exactly 2 satin-column paths (each already exactly 2 rails, plus
any rungs, which are identified as the two longest subpaths and left
untouched) that meet at one crossing. Whichever of the two has an end
sitting nearest the other's rails is the "under" piece needing fixed;
the other is the "over" piece and is never modified. Run once per
crossing that needs closing (same one-crossing-at-a-time selection model
as Cross-Splice Rails).

Extension method (confirmed with the user against a real hand-edited
example before building this): each of the under end's two rail tips is
extended by a straight line continuing its own local tangent, out to
wherever that line crosses the over column's nearest rail - the two tips
are NOT forced to a shared point, each reaches its own crossing
independently. If that crossing point lies further out than the current
tip (a real gap), a new straight segment is appended to reach it. If it
lies back toward the rail's existing body instead (the cut already
overshot), the tip node is simply moved back to that point instead
(flattening its handles) rather than fully re-deriving the trimmed curve
- a straight-line approximation in both directions, consistent with the
extension method chosen. The little rung capping the cut end is left
exactly as it was - not regenerated - per the user's explicit call.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import inkex
from inkex import EffectExtension, ShapeElement, PathElement

import rung_rail_satin as rrs

NEAR_ROUND = 1e-6


def _rails_and_rungs(csp):
    """Split a satin column's subpaths into its 2 rails (the longest 2, by
    arc length) and everything else (rungs), preserved as-is."""
    rails_info = sorted(
        ((rrs.Rail(sp), sp) for sp in csp),
        key=lambda pair: -pair[0].length,
    )
    if len(rails_info) < 2:
        raise inkex.AbortExtension("Selected path needs at least 2 subpaths (2 rails).")
    rails = [rails_info[0], rails_info[1]]
    other_subpaths = [sp for _, sp in rails_info[2:]]
    return rails, other_subpaths


def _min_dist_point_to_rails(point, rails):
    point = np.array(point, dtype=float)
    best = None
    for rail, _ in rails:
        for (ts, arclen, pts) in rail._tables:
            d = float(np.min(np.linalg.norm(pts - point, axis=1)))
            if best is None or d < best:
                best = d
    return best


def _end_tangent(rail, at_end):
    """Outward-pointing unit tangent at the given end of the rail (the
    direction continuing past the cut, away from the rail body)."""
    seg = rail.segs[-1] if at_end else rail.segs[0]
    d = rrs.Rail._boundary_tangent(seg, at_end=at_end)
    if not at_end:
        # _boundary_tangent(at_end=False) is the curve's forward derivative
        # at t=0, which points INTO the rail body (toward its interior) -
        # the outward direction continuing past a start-end cut is the
        # reverse of that.
        d = -d
    n = np.linalg.norm(d)
    if n < 1e-9:
        return np.array([1.0, 0.0])
    return d / n


def _extend_or_trim(subpath, rail, at_end, hit_point):
    """Return a new subpath with the given end's tip node moved/extended
    to hit_point, per the straight-line approximation described above."""
    tip = np.array(rail.end() if at_end else rail.start())
    tangent = _end_tangent(rail, at_end)
    hit = np.array(hit_point, dtype=float)
    t = float(np.dot(hit - tip, tangent))

    sp = [list(map(list, node)) for node in subpath]
    if t >= 0:
        # Extend: append one new straight segment out to the hit point.
        if at_end:
            sp[-1][2] = list(sp[-1][1])  # flatten this node's out-handle
            sp.append([list(hit), list(hit), list(hit)])
        else:
            sp[0][0] = list(sp[0][1])  # flatten this node's in-handle
            sp.insert(0, [list(hit), list(hit), list(hit)])
    else:
        # Trim: move the tip node itself back to the hit point, flattening
        # its handles (approximates the true curve-trim with a straight
        # local edit, consistent with the extension method).
        if at_end:
            sp[-1] = [list(hit), list(hit), list(hit)]
        else:
            sp[0] = [list(hit), list(hit), list(hit)]
    return sp


def _line_hit(anchor, direction, rails):
    """Nearest valid crossing of the infinite line (anchor, direction)
    against either of `rails`, or None."""
    best = None
    for rail, _ in rails:
        hit = rail.line_cross(np.array(anchor, dtype=float), np.array(direction, dtype=float))
        if hit is None:
            continue
        pt, s = hit
        d = float(np.linalg.norm(pt - np.array(anchor)))
        if best is None or d < best[0]:
            best = (d, pt)
    return None if best is None else best[1]


class TouchRails(EffectExtension):

    def effect(self):
        sel = list(self.svg.selection.get(ShapeElement))
        if len(sel) != 2:
            raise inkex.AbortExtension(
                "Select exactly 2 satin-column paths that meet at one crossing: "
                "the under piece with a cut end, and the over piece it should touch."
            )

        els = []
        for e in sel:
            if not isinstance(e, PathElement):
                if hasattr(e, "to_path_element"):
                    e = e.to_path_element()
                else:
                    raise inkex.AbortExtension(
                        "Both selected objects must be paths (or shapes convertible to paths)."
                    )
            els.append(e)

        data = []
        for e in els:
            csp = e.path.to_absolute().transform(e.composed_transform()).to_superpath()
            rails, rungs = _rails_and_rungs(csp)
            data.append({"el": e, "rails": rails, "rungs": rungs})

        # Find the (element, end) pair whose tip sits nearest the OTHER
        # element's rails - that identifies the under piece and which end
        # of it needs fixing.
        best = None  # (dist, under_idx, at_end)
        for i in (0, 1):
            other = data[1 - i]["rails"]
            for rail, _ in data[i]["rails"]:
                for at_end in (True, False):
                    tip = rail.end() if at_end else rail.start()
                    d = _min_dist_point_to_rails(tip, other)
                    if best is None or d < best[0]:
                        best = (d, i, at_end)

        _, under_idx, at_end = best
        over_idx = 1 - under_idx
        under = data[under_idx]
        over = data[over_idx]

        new_rail_subpaths = []
        for rail, subpath in under["rails"]:
            tip = rail.end() if at_end else rail.start()
            tangent = _end_tangent(rail, at_end)
            hit = _line_hit(tip, tangent, over["rails"])
            if hit is None:
                new_rail_subpaths.append(subpath)
                continue
            new_rail_subpaths.append(_extend_or_trim(subpath, rail, at_end, hit))

        under_el = under["el"]
        under_el.path = rrs.make_csp(new_rail_subpaths + under["rungs"]).to_path()
        under_el.transform = -under_el.getparent().composed_transform()


if __name__ == "__main__":
    TouchRails().run()
