#!/usr/bin/env python3
"""
Weave Paths (Ink/Stitch satin column helper)

Select 2+ open paths that cross each other (e.g. strands of a decorative
knot/braid) and run. Produces the visual effect of the strands physically
weaving over and under each other at every crossing, by cutting a gap in
whichever strand is "under" at each one - the same effect as manually
double-clicking to add nodes and deleting the hidden segment, done for
every crossing at once.

Algorithm (confirmed with the user against a real 3-strand test file
before building this):

- Paths are processed one at a time, in SELECTION order (the order
  `self.svg.selection` iterates them - same unverified click-order
  assumption already relied on elsewhere in this extension set, e.g.
  Cross-Splice Rails' `--id` handling).
- For the path currently being processed, its crossings with every other
  selected path are visited in the order they occur along ITS OWN length
  (start to end), and alternate over -> under -> over -> under, starting
  at "over".
- A crossing is one shared fact between exactly two paths, not two
  independent decisions: if a crossing was already resolved while an
  earlier path was being processed, the current path inherits the
  (guaranteed opposite) outcome for free and this does NOT consume a
  step in its own toggle - only genuinely fresh crossings do. Verified
  against knot-test.svg's 3-strand case: treating every crossing as
  fresh (no inheritance) produces direct contradictions after the very
  first crossing; skipping already-resolved ones does not, for any
  selection order.
- "over": the OTHER path gets the gap, sized to 1.5x the CURRENT path's
  own stroke width (the current path is the one staying visually on
  top, so the gap in the other path must be wide enough that its own
  stroke doesn't peek out from under it).
- "under": the CURRENT path gets the gap instead, sized to 1.5x the
  OTHER path's stroke width, for the same reason in reverse.

Each cut path is split into separate open PathElements (one per surviving
piece), replacing the original.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import inkex
from inkex import EffectExtension, ShapeElement, PathElement
from inkex.paths import CubicSuperPath

import bezier_geom as bg
import rung_rail_satin as rrs

GAP_WIDTH_MULT = 1.5
POINT_DEDUP_TOL = 1e-2  # doc units - merges an intersection reported twice
                          # (e.g. right at a segment boundary) into one


def stroke_width(el):
    return rrs.stroke_width(el)


def find_crossings(railA, railB):
    """All intersections between two Rail objects, as (sA, sB, point) in
    each rail's own arc-length parametrization."""
    hits = []
    for i, segA in enumerate(railA.segs):
        tsA, arclenA, _ = railA._tables[i]
        for j, segB in enumerate(railB.segs):
            tsB, arclenB, _ = railB._tables[j]
            raw = bg.cubic_intersect(
                tuple(map(tuple, segA)), tuple(map(tuple, segB))
            )
            for t1, t2, pt in raw:
                sA = railA.seg_bounds[i] + float(np.interp(t1, tsA, arclenA))
                sB = railB.seg_bounds[j] + float(np.interp(t2, tsB, arclenB))
                hits.append((sA, sB, pt))
    return _dedup_by_point(hits)


def _dedup_by_point(hits):
    out = []
    for sA, sB, pt in hits:
        if any(abs(pt[0] - o[2][0]) < POINT_DEDUP_TOL and
               abs(pt[1] - o[2][1]) < POINT_DEDUP_TOL for o in out):
            continue
        out.append((sA, sB, pt))
    return out


def _crossing_key(i, j, pt):
    a, b = (i, j) if i < j else (j, i)
    return (a, b, round(pt[0], 1), round(pt[1], 1))


def resolve_weave(items):
    """items: list of dicts each with a 'rail' (Rail). Returns cuts: list
    (parallel to items) of [(s, half_width), ...] to remove from each
    item's own arc length."""
    n = len(items)
    cuts = [[] for _ in range(n)]
    resolved = {}

    for i in range(n):
        entries = []
        for j in range(n):
            if i == j:
                continue
            for (s_self, s_other, pt) in find_crossings(items[i]["rail"], items[j]["rail"]):
                entries.append({"other": j, "s_self": s_self, "point": pt})
        entries.sort(key=lambda e: e["s_self"])

        toggle = "over"
        for e in entries:
            j = e["other"]
            key = _crossing_key(i, j, e["point"])
            if key in resolved:
                outcome = resolved[key][i]
            else:
                outcome = toggle
                resolved[key] = {i: outcome, j: ("under" if outcome == "over" else "over")}
                toggle = "under" if toggle == "over" else "over"
            if outcome == "under":
                half_width = GAP_WIDTH_MULT * items[j]["sw"] / 2.0
                cuts[i].append((e["s_self"], half_width))

    return cuts


def _s_to_seg_t(rail, s):
    i, local_s = rail.seg_index_at(s)
    ts, arclen, _ = rail._tables[i]
    t = float(np.interp(local_s, arclen, ts))
    return i, t


def _node(pt_in, pt, pt_out):
    return [list(pt_in), list(pt), list(pt_out)]


def crop_subpath(rail, s0, s1):
    """Return a new CSP subpath spanning arc-length [s0, s1] of `rail`,
    with proper De Casteljau splitting at both ends (not a polyline
    approximation) so curved rails stay exactly curved."""
    s0 = max(0.0, min(s0, rail.length))
    s1 = max(0.0, min(s1, rail.length))
    i0, t0 = _s_to_seg_t(rail, s0)
    i1, t1 = _s_to_seg_t(rail, s1)

    if i0 == i1:
        _, right0 = bg.split_cubic(*rail.segs[i0], t0)
        t1_local = (t1 - t0) / (1 - t0) if t0 < 1.0 else 0.0
        left_final, _ = bg.split_cubic(*right0, t1_local)
        p0, p1, p2, p3 = left_final
        return [_node(p0, p0, p1), _node(p2, p3, p3)]

    _, right0 = bg.split_cubic(*rail.segs[i0], t0)
    left1, _ = bg.split_cubic(*rail.segs[i1], t1)

    nodes = [_node(right0[0], right0[0], right0[1])]
    prev_out = right0[2]
    prev_pt = right0[3]
    for k in range(i0 + 1, i1):
        p0, p1, p2, p3 = rail.segs[k]
        nodes.append(_node(prev_out, prev_pt, p1))
        prev_out = p2
        prev_pt = p3
    p0, p1, p2, p3 = left1
    nodes.append(_node(prev_out, prev_pt, p1))
    nodes.append(_node(p2, p3, p3))
    return nodes


def apply_weave_cuts(rail, cuts):
    """cuts: list of (s, half_width). Returns list of CSP subpaths for the
    surviving pieces after removing each cut interval (merging overlaps)."""
    intervals = []
    for s, hw in cuts:
        intervals.append((max(0.0, s - hw), min(rail.length, s + hw)))
    intervals.sort()
    merged = []
    for s0, s1 in intervals:
        if merged and s0 <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], s1))
        else:
            merged.append((s0, s1))

    pieces = []
    cursor = 0.0
    for s0, s1 in merged:
        if s0 > cursor + 1e-6:
            pieces.append(crop_subpath(rail, cursor, s0))
        cursor = s1
    if cursor < rail.length - 1e-6:
        pieces.append(crop_subpath(rail, cursor, rail.length))
    return pieces


class WeavePaths(EffectExtension):

    def effect(self):
        sel = list(self.svg.selection.get(ShapeElement))
        if len(sel) < 2:
            raise inkex.AbortExtension(
                "Select 2 or more open paths that cross each other."
            )

        items = []
        for e in sel:
            if not isinstance(e, PathElement):
                if hasattr(e, "to_path_element"):
                    e = e.to_path_element()
                else:
                    raise inkex.AbortExtension(
                        "All selected objects must be paths (or shapes convertible to paths)."
                    )
            abs_path = e.path.to_absolute().transform(e.composed_transform())
            csp = abs_path.to_superpath()
            if len(csp) != 1:
                raise inkex.AbortExtension(
                    f"Each selected path must be a single open subpath (found {len(csp)} on one selection)."
                )
            items.append({
                "el": e,
                "sw": stroke_width(e),
                "rail": rrs.Rail(csp[0]),
            })

        cuts = resolve_weave(items)

        for item, item_cuts in zip(items, cuts):
            pieces = apply_weave_cuts(item["rail"], item_cuts)
            el = item["el"]
            parent = el.getparent()
            idx = list(parent).index(el)
            if not pieces:
                el.delete()
                continue
            for k, subpath in enumerate(pieces):
                new_el = el.copy() if k > 0 else el
                new_el.path = rrs.make_csp([subpath]).to_path()
                new_el.transform = -parent.composed_transform()
                if k > 0:
                    new_el.set_id(f"{el.get_id()}_w{k}")
                    parent.insert(idx + k, new_el)


if __name__ == "__main__":
    WeavePaths().run()
