#!/usr/bin/env python3
"""
Clip Under to Over (Ink/Stitch satin column helper)

Replaces Touch Rails at Crossing and Taper Join at Crossing. Confirmed
against a 4-figure reference from the user:

  fig.1: a satin column is rails + rungs; Ink/Stitch computes the actual
  zigzag stitch fill BETWEEN them (fig.2) - fill can only exist where
  BOTH rails exist at a given point.

  fig.3: independently clipping each of the under piece's two rails at
  its own entry point into the over piece looks reasonable on the rails
  alone, but at any non-perpendicular crossing angle the two rails enter
  at different points along the crossing - so past wherever the SHORTER
  one stops, there's no second rail to pair with, and Ink/Stitch can't
  fill that stretch even though the longer rail still has material
  there. A visible unstitched gap right at the tip.

  fig.4: the fix is NOT to converge both rails to one shared point.
  Instead, the shorter rail is extended past its own entry point,
  tracking ALONG the over piece's near rail (following its curve, not a
  straight line to a point) until it reaches the position that lines up
  with - is incident with - the longer rail's own endpoint. Both
  under-rails then end at corresponding positions, so Ink/Stitch can
  fill the satin all the way to the longer rail's tip.

Select exactly 2 satin-column paths at one crossing. Under/over and
which end of the under piece is involved are auto-detected the same way
as the earlier tools (whichever end sits nearest the other piece's
rails).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import inkex
from inkex import EffectExtension, ShapeElement, PathElement

import bezier_geom as bg
import rung_rail_satin as rrs
import touch_rails as tr


# bg.point_in_polygon is O(polygon size) per query - the rails' full dense
# tables (N_PER_SEG=400 points/segment in rung_rail_satin) are far denser
# than needed for a gently-curving satin rail, and using them directly for
# both the polygon AND the per-point scan made this O(n*m) with n,m in the
# thousands (minutes, not seconds). Decimating both to a coarse stride
# keeps the shape accurate enough for this purpose and is fast.
POLY_STRIDE = 10
SCAN_STRIDE = 10


def _strided_points(rail, stride):
    """[(s, point), ...] across the whole rail, decimated to every
    `stride`-th sample (always including each segment's first/last point
    so segment joins aren't skipped)."""
    out = []
    for i, (ts, arclen, pts) in enumerate(rail._tables):
        base = rail.seg_bounds[i]
        idxs = list(range(0, len(pts), stride))
        if idxs[-1] != len(pts) - 1:
            idxs.append(len(pts) - 1)
        for j in idxs:
            out.append((base + float(arclen[j]), pts[j]))
    return out


def _ribbon_polygon(rails):
    """Closed polygon tracing the ribbon bounded by the 2 rails: forward
    along the first, backward along the second."""
    (railA, _), (railB, _) = rails
    ptsA = [p for (_, p) in _strided_points(railA, POLY_STRIDE)]
    ptsB = [p for (_, p) in _strided_points(railB, POLY_STRIDE)]
    return [tuple(p) for p in ptsA] + [tuple(p) for p in ptsB[::-1]]


def _entry_s(rail, at_end, polygon, known_outside_s=None):
    """Walking from the OUTSIDE end inward, the arc-length s where this
    rail first enters `polygon` - the over piece's ribbon. None if it
    never does within the rail's current (already weave-cut) length.

    The coarse decimated scan only brackets the crossing (between the
    last outside sample and the first inside one) - returning that inside
    sample directly let the clip land a couple mm PAST the real boundary,
    into the over piece's interior, instead of exactly on it. Bisect
    within the bracket to land on the boundary itself.

    `known_outside_s`: the original (pre-extension) tip's position in
    this rail's own s parametrization - a point ALWAYS known to be
    outside the polygon (that's the premise of extending at all). Seeds
    the bisection bracket for the edge case where the very first scanned
    sample is already inside (no earlier outside sample to bracket
    against) - without this, that case fell back to returning the raw
    coarse sample verbatim, which could be many mm past the true
    boundary."""
    pts = _strided_points(rail, SCAN_STRIDE)
    if at_end:
        pts = list(reversed(pts))
    prev_s = known_outside_s
    for s, p in pts:
        if bg.point_in_polygon(tuple(p), polygon):
            if prev_s is None:
                return s
            lo_s, hi_s = prev_s, s
            for _ in range(20):
                mid_s = (lo_s + hi_s) / 2.0
                mid_p, _ = rail.point_tangent(mid_s)
                if bg.point_in_polygon(tuple(mid_p), polygon):
                    hi_s = mid_s
                else:
                    lo_s = mid_s
            return lo_s  # last point confirmed still outside/on the boundary
        prev_s = s
    return None


def _exit_s(rail, at_end, polygon):
    """Mirror of _entry_s for the OPPOSITE situation: the rail's raw tip
    is already INSIDE the over piece's ribbon (Weave Paths' cut is
    computed for one centerline, but Ink/Stitch then offsets that into
    two separate rails - one can clear the real boundary while the
    other still overlaps it). Walking from the tip INTO the rail's own
    existing body, the arc-length s where it first exits back outside -
    where to trim the overlapping material back to. None if it never
    exits within the rail's own length."""
    pts = _strided_points(rail, SCAN_STRIDE)
    if at_end:
        pts = list(reversed(pts))
    prev_s = None
    for s, p in pts:
        inside = bg.point_in_polygon(tuple(p), polygon)
        if prev_s is not None and not inside:
            lo_s, hi_s = prev_s, s  # lo=inside, hi=outside
            for _ in range(20):
                mid_s = (lo_s + hi_s) / 2.0
                mid_p, _ = rail.point_tangent(mid_s)
                if bg.point_in_polygon(tuple(mid_p), polygon):
                    lo_s = mid_s
                else:
                    hi_s = mid_s
            return hi_s  # first point confirmed outside/on the boundary
        prev_s = s
    return None


def _resolve_rail_boundary(rail0, subpath0, at_end, polygon):
    """Returns (rail, subpath, s_cut) in a consistent form ready for
    _clip_subpath/_project_s/_extend_along, handling both real cases:
    the rail's raw tip is outside the over piece (the normal case -
    extend forward via a synthetic continuation to find where it
    enters), or already inside it (trim backward into the rail's own
    existing material to find where it exits) - see _exit_s."""
    tip0 = rail0.end() if at_end else rail0.start()
    if bg.point_in_polygon(tuple(tip0), polygon):
        s_cut = _exit_s(rail0, at_end, polygon)
        if s_cut is None:
            raise inkex.AbortExtension(
                "A rail's tip is already inside the over piece's shape and "
                "never exits back out within its own length - check this "
                "is really the right pair."
            )
        return rail0, subpath0, s_cut
    else:
        rail_ext, sp_ext = _extended_rail(rail0, subpath0, at_end)
        known_outside = rail0.length if at_end else EXTRAPOLATE_LEN
        s_cut = _entry_s(rail_ext, at_end, polygon, known_outside)
        return rail_ext, sp_ext, s_cut


# Weave Paths' own cut leaves the under piece's raw tip short of the over
# piece by design (a real gap, typically a few mm) - _entry_s alone would
# never find a crossing there since it only looks at the EXISTING rail.
# Extend a synthetic straight continuation past the tip first, generous
# enough to clear any such gap, then search that too.
EXTRAPOLATE_LEN = 20.0


def _extended_rail(rail, subpath, at_end):
    tangent = tr._end_tangent(rail, at_end)
    tip = rail.end() if at_end else rail.start()
    new_pt = np.array(tip, dtype=float) + tangent * EXTRAPOLATE_LEN
    sp = [list(map(list, node)) for node in subpath]
    if at_end:
        sp[-1][2] = list(sp[-1][1])
        sp.append([list(new_pt), list(new_pt), list(new_pt)])
    else:
        sp[0][0] = list(sp[0][1])
        sp.insert(0, [list(new_pt), list(new_pt), list(new_pt)])
    return rrs.Rail(sp), sp


def _project_s(point, rail):
    """Nearest arc-length position of `point` on `rail`."""
    point = np.array(point, dtype=float)
    best_s, best_d = None, None
    for i, (ts, arclen, pts) in enumerate(rail._tables):
        d = np.linalg.norm(pts - point, axis=1)
        j = int(np.argmin(d))
        if best_d is None or d[j] < best_d:
            best_d = d[j]
            best_s = rail.seg_bounds[i] + float(arclen[j])
    return best_s


def _nearest_over_rail(point, over_rails):
    best = None
    for rail, sp in over_rails:
        s = _project_s(point, rail)
        p, _ = rail.point_tangent(s)
        d = float(np.linalg.norm(np.array(p) - np.array(point)))
        if best is None or d < best[0]:
            best = (d, rail, sp, s)
    return best[1], best[2]


def _clip_subpath(subpath, rail, at_end, s_cut):
    """Trim `subpath` so the given end now stops exactly at arc-length
    s_cut, via De Casteljau split at that point."""
    i, local_s = rail.seg_index_at(s_cut)
    ts, arclen, _ = rail._tables[i]
    t = float(np.interp(local_s, arclen, ts))
    left, right = bg.split_cubic(*rail.segs[i], t)

    sp = [list(map(list, node)) for node in subpath]
    if at_end:
        keep = sp[: i + 1]
        keep[-1][2] = list(left[1])
        keep.append([list(left[2]), list(left[3]), list(left[3])])
        return keep
    else:
        keep = sp[i + 1:]
        keep[0][0] = list(right[2])
        keep.insert(0, [list(right[0]), list(right[0]), list(right[1])])
        return keep


EXTEND_TRACE_POINTS = 8


def _extend_along(subpath, at_end, near_rail, s_from, s_to):
    """Append a light polyline from the given end of `subpath` out to the
    point on `near_rail` at arc-length s_to, sampling near_rail's own
    curve at a handful of evenly-spaced points along the way - the
    "shorter rail extended along the over rail" step from fig.4.

    A single straight chord (tried first, matching a short g10 example)
    is wrong whenever s_from/s_to straddle a real bend in the over rail -
    over any real separation it cuts across the bend toward the ribbon's
    other edge instead of following it, which is the "extends through
    the interior" symptom. A dense per-sample trace (tried before that)
    swung too far the other way, hundreds of nearly-collinear nodes for
    a few mm. A handful of evenly-spaced points on the real curve is
    the middle ground: follows the bend, stays small."""
    sp = [list(map(list, node)) for node in subpath]
    steps = np.linspace(s_from, s_to, EXTEND_TRACE_POINTS)[1:]
    trace_pts = [near_rail.point_tangent(s)[0] for s in steps]
    if at_end:
        sp[-1][2] = list(sp[-1][1])
        for pt in trace_pts:
            sp.append([list(pt), list(pt), list(pt)])
    else:
        sp[0][0] = list(sp[0][1])
        for pt in reversed(trace_pts):
            sp.insert(0, [list(pt), list(pt), list(pt)])
    return sp


def _drop_orphaned_rungs(rungs, polygon):
    return [sp for sp in rungs if not any(bg.point_in_polygon(tuple(n[1]), polygon) for n in sp)]


class ClipToOver(EffectExtension):

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
            rails, rungs = tr._rails_and_rungs(csp)
            data.append({"el": e, "rails": rails, "rungs": rungs})

        best = None
        for i in (0, 1):
            other = data[1 - i]["rails"]
            for rail, _ in data[i]["rails"]:
                for at_end in (True, False):
                    tip = rail.end() if at_end else rail.start()
                    d = tr._min_dist_point_to_rails(tip, other)
                    if best is None or d < best[0]:
                        best = (d, i, at_end)

        _, under_idx, at_end = best
        over_idx = 1 - under_idx
        under = data[under_idx]
        over = data[over_idx]

        polygon = _ribbon_polygon(over["rails"])

        (railA0, spA0), (railB0, spB0) = under["rails"]
        # Each rail is resolved independently - one can legitimately start
        # outside the over piece (the normal weave-cut gap) while the
        # other already overlaps it (Ink/Stitch's rail offset can eat the
        # gap on one side even though it was fine on the centerline) - see
        # _resolve_rail_boundary.
        railA, spA, entry_A = _resolve_rail_boundary(railA0, spA0, at_end, polygon)
        railB, spB, entry_B = _resolve_rail_boundary(railB0, spB0, at_end, polygon)

        tipA = railA.point_tangent(entry_A)[0]
        tipB = railB.point_tangent(entry_B)[0]

        near_rail, _ = _nearest_over_rail(((np.array(tipA) + np.array(tipB)) / 2.0), over["rails"])
        proj_A = _project_s(tipA, near_rail)
        proj_B = _project_s(tipB, near_rail)

        tanA = tr._end_tangent(railA, at_end)
        tanB = tr._end_tangent(railB, at_end)
        avg_tan = tanA + tanB
        n = np.linalg.norm(avg_tan)
        avg_tan = avg_tan / n if n > 1e-9 else tanA
        mid_s = (proj_A + proj_B) / 2.0
        _, over_tan = near_rail.point_tangent(mid_s)
        forward = np.dot(avg_tan, over_tan) > 0

        a_ahead = (proj_A > proj_B) if forward else (proj_A < proj_B)

        new_subpaths = {}
        if a_ahead:
            new_subpaths["A"] = _clip_subpath(spA, railA, at_end, entry_A) if entry_A is not None else spA
            base = _clip_subpath(spB, railB, at_end, entry_B) if entry_B is not None else spB
            new_subpaths["B"] = _extend_along(base, at_end, near_rail, proj_B, proj_A)
        else:
            new_subpaths["B"] = _clip_subpath(spB, railB, at_end, entry_B) if entry_B is not None else spB
            base = _clip_subpath(spA, railA, at_end, entry_A) if entry_A is not None else spA
            new_subpaths["A"] = _extend_along(base, at_end, near_rail, proj_A, proj_B)

        kept_rungs = _drop_orphaned_rungs(under["rungs"], polygon)

        under_el = under["el"]
        under_el.path = rrs.make_csp([new_subpaths["A"], new_subpaths["B"]] + kept_rungs).to_path()
        under_el.transform = -under_el.getparent().composed_transform()


if __name__ == "__main__":
    ClipToOver().run()
