#!/usr/bin/env python3
"""
Generate Ladder Rungs (Ink/Stitch satin column helper)

Select either:
 - two open "rail" paths (they will be combined automatically), or
 - one already-combined path with exactly two rail subpaths

and run. Builds the cross-tie "rung" segments (zig-zag stitches) between the
two rails, following a rule set worked out by hand on real satin designs:

 - A = measured rail-to-rail distance, B = stroke width.
 - Rungs are perpendicular to the local rail tangent, EXCEPT inside a taper,
   where they run parallel to the last non-taper rung instead.
 - A taper is detected structurally: at a given end of the ladder, if one
   rail has more than one segment and its LAST segment lands at that end
   (regardless of whether that segment is straight or curved - a node
   marks the boundary either way), that segment is the taper.
 - Taper zone: 3 rungs (before/after the taper's start node, and as close
   to the convergence point as possible while staying >= 2B from it on
   BOTH rails). A 4th "middle of taper" rung only if the taper's arc
   length is clearly longer than A (> 1.15*A).
 - An end with no taper on either rail (gap there is close to A, not a
   near-convergence) gets one final rung parallel to the straight line
   between the two raw endpoints, pulled back >=4*B from them on both
   rails.
 - Interior spacing: 3*A on gentle curves/straights, 1*A on tighter curves.
 - Every rung ends RUNG_OVERSHOOT_MM (0.125 mm) past each rail, unless
   clear_rung_ends shortens it where the column curls back on itself.
 - Interior rungs join PAIRED points on the two rails, not a line square to
   one rail: pair_rails() matches the rails end to end so each rung is as
   near square to both rails as it can be, never going back along either
   rail. Where both rails curve around a common center (a curl, a spiral)
   the rungs fan like spokes, and none can cross another or jump to a
   different turn of a curling rail. On parallel rails this gives the same
   square rungs as before.
 - An end where the rails meet at a point without a corner (a smooth tip,
   e.g. the tail of a curl) is a "point" end: the paired rungs run right
   down into the tip, the last one where the column is still at least
   POINT_MIN_WIDTH_MULT * B wide and MIN_FROM_POINT_MULT * B from the point.

This is a from-scratch implementation of that rule set, checked against a
small number of real objects and not independently verified beyond that. Treat its output as a draft to
check, especially the taper-vs-open classification and the 3-vs-4 rung
threshold, both of which are heuristics.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import math
import numpy as np
import inkex
from inkex import EffectExtension, ShapeElement, PathElement
from inkex.paths import CubicSuperPath

import bezier_geom as bg

STROKE_TAPER_MULT = 4.0
# Overshoot beyond each rail's CENTERLINE crossing point. A rail's outer
# edge sits 0.5*B from its centerline, so 0.5 (to reach the edge) + 1.0
# (requested overshoot past the edge) = 1.5 gives 1*B of clearance past
# the actual rail stroke's outer edge, matching the node markers to the
# rail edges as intended.
STROKE_DEFAULT_MULT = 1.5
STROKE_OPEN_MULT = 3.0
MIN_FROM_POINT_MULT = 2.0
TAPER_EXTRA_RUNG_MARGIN = 1.15
TAPER_MAX_A_MULT = 3.0
SPACING_GENTLE_MULT = 3.0
SPACING_TIGHT_MULT = 1.0
STRAIGHTNESS_THRESHOLD = 0.9
N_PER_SEG = 400

# Rule (user, 2026-09-24): no rung ever sits exactly at a rail end. The
# last (closest-to-end) rung on a non-taper open end is inset this many
# column-widths (local rail separation, not stroke width) from the end -
# EXCEPT a genuine taper's near-the-convergence-point rung, which is
# still required close to the tip to keep the zigzag perpendicular to the
# rails there (governed by MIN_FROM_POINT_MULT instead, unchanged).
END_INSET_MULT = 0.5

# Every rung ends this many mm past each rail crossing, the same on both
# sides whatever the column width - replacing the old asymmetric/B-based
# overshoot multipliers above for actual rung placement (those stay only
# as historical defaults/unused otherwise). clear_rung_ends can still
# shorten an end where a column curls back on itself.
RUNG_OVERSHOOT_MM = 0.125
# A rung ending on its own rail at a curl-back must cross that rail at
# 30 degrees or more at both ends (sine of the angle), or it is a graze.
CURL_BACK_MIN_SIN = 0.5

# A point end's closest rung sits where the column is at least this many
# stroke widths wide, so it still clearly crosses both rails.
POINT_MIN_WIDTH_MULT = 4.0

# pair_rails: samples per rail are about 1 per this fraction of A (column
# width), between PAIR_MIN_SAMPLES and PAIR_MAX_SAMPLES. PAIR_DIST_WEIGHT is
# a small pull toward nearer partners, so a choice between equally square
# pairings goes to the shorter rung.
PAIR_SAMPLE_A_FRACTION = 0.1
PAIR_MIN_SAMPLES = 60
PAIR_MAX_SAMPLES = 400
PAIR_DIST_WEIGHT = 0.02


def reverse_subpath(subpath):
    """Reverse node order and each node's in/out handles, so the subpath
    traces the same geometry in the opposite direction."""
    rev = []
    for node in reversed(subpath):
        rev.append([list(node[2]), list(node[1]), list(node[0])])
    return rev


def normalize_rail_directions(csp0, csp1):
    """Rails are often hand-drawn anti-parallel (one rail's start sits near
    the other's end, not its start) - as parallel as not, going by the
    examples seen so far. Detect whichever pairing (start-start/end-end vs
    start-end/end-start) is shorter and reverse csp1 if needed, so all
    downstream code can safely assume start pairs with start."""
    s0, e0 = np.array(csp0[0][1]), np.array(csp0[-1][1])
    s1, e1 = np.array(csp1[0][1]), np.array(csp1[-1][1])
    same_dir = np.linalg.norm(s0 - s1) + np.linalg.norm(e0 - e1)
    opp_dir = np.linalg.norm(s0 - e1) + np.linalg.norm(e0 - s1)
    if opp_dir < same_dir:
        return csp0, reverse_subpath(csp1)
    return csp0, csp1


def find_rail_indices(csp):
    """Rails aren't reliably subpaths 0 and 1 - a hand-edited ladder, or one
    another tool rewrote, can have rungs ahead of or between the rails. A
    rung is never longer than a rail, so the rails are the two LONGEST
    subpaths by arc length. Returned in their order in the path."""
    lengths = [Rail(sp).length for sp in csp]
    order = sorted(range(len(csp)), key=lambda i: lengths[i], reverse=True)
    return sorted(order[:2])


def rails_first(csp):
    """The same subpaths with the two rails moved to the front (in their own
    order) and every rung after them, so csp[0], csp[1] are the rails and
    csp[2:] the rungs. Written back, Ink/Stitch and every tool see that order."""
    csp = list(csp)
    if len(csp) <= 2:
        return csp
    ia, ib = find_rail_indices(csp)
    return [csp[ia], csp[ib]] + [sp for i, sp in enumerate(csp) if i not in (ia, ib)]


def make_csp(subpaths):
    obj = CubicSuperPath.__new__(CubicSuperPath)
    list.__init__(obj, subpaths)
    return obj


def unit_normal(d):
    d = np.array(d, dtype=float)
    n = np.linalg.norm(d)
    if n < 1e-12:
        return np.array([0.0, 0.0])
    d = d / n
    return np.array([-d[1], d[0]])


class Rail:
    """A rail = one CSP subpath, generically N segments. Provides arc-length
    parametrized point/tangent lookup and per-segment crossing search."""

    def __init__(self, subpath):
        self.segs = []
        for i in range(len(subpath) - 1):
            p0 = np.array(subpath[i][1], dtype=float)
            p1 = np.array(subpath[i][2], dtype=float)
            p2 = np.array(subpath[i + 1][0], dtype=float)
            p3 = np.array(subpath[i + 1][1], dtype=float)
            self.segs.append((p0, p1, p2, p3))
        self.n = len(self.segs)
        # per-segment dense sample tables for arc length <-> t
        self._tables = []
        bounds = [0.0]
        for (p0, p1, p2, p3) in self.segs:
            ts = np.linspace(0, 1, N_PER_SEG)
            pts = np.array([bg.split_point(p0, p1, p2, p3, t) for t in ts])
            seglen = np.linalg.norm(np.diff(pts, axis=0), axis=1)
            arclen = np.concatenate([[0], np.cumsum(seglen)])
            self._tables.append((ts, arclen, pts))
            bounds.append(bounds[-1] + arclen[-1])
        self.seg_bounds = bounds  # len n+1, seg i spans [bounds[i], bounds[i+1]]
        self.length = bounds[-1]

    def seg_index_at(self, s):
        s = min(max(s, 0.0), self.length)
        for i in range(self.n):
            if s <= self.seg_bounds[i + 1] or i == self.n - 1:
                return i, s - self.seg_bounds[i]
        return self.n - 1, self.seg_bounds[-1] - self.seg_bounds[-2]

    def point_tangent(self, s):
        i, local_s = self.seg_index_at(s)
        ts, arclen, _ = self._tables[i]
        t = float(np.interp(local_s, arclen, ts))
        p0, p1, p2, p3 = self.segs[i]
        tangent = bg.bez_deriv(p0, p1, p2, p3, t) \
            if hasattr(bg, "bez_deriv") else self._deriv(p0, p1, p2, p3, t)
        # Same degenerate-boundary case _boundary_tangent guards against:
        # a straight segment's collapsed control points make the cubic
        # derivative evaluate to exactly zero at t=0 or t=1 (found live
        # 2026-09-24 - a short single-segment reference rail's t=0 end
        # produced a zero tangent here, silently killing a taper rung's
        # line_cross search with no error). Fall back to a local secant.
        if np.linalg.norm(tangent) < 1e-9 and (t < 1e-9 or t > 1 - 1e-9):
            near_t = 0.9 if t > 0.5 else 0.1
            near_pt = np.array(bg.split_point(p0, p1, p2, p3, near_t))
            end_pt = np.array(p3) if t > 0.5 else np.array(p0)
            tangent = end_pt - near_pt if t > 0.5 else near_pt - end_pt
        return bg.split_point(p0, p1, p2, p3, t), tangent

    @staticmethod
    def _deriv(p0, p1, p2, p3, t):
        p0, p1, p2, p3 = (np.array(p) for p in (p0, p1, p2, p3))
        mt = 1 - t
        return 3 * (mt ** 2) * (p1 - p0) + 6 * mt * t * (p2 - p1) + 3 * (t ** 2) * (p3 - p2)

    def start(self):
        return self.segs[0][0]

    def end(self):
        return self.segs[-1][3]

    @staticmethod
    def _boundary_tangent(seg, at_end):
        """Tangent of a segment at t=1 (at_end) or t=0, falling back to a
        LOCAL secant (the last/first 10% of the curve, not the whole-segment
        chord) when the cubic derivative degenerates to (near) zero - which
        happens whenever a control point collapses onto its neighboring
        endpoint (the same pattern seen earlier this session for plain `L`
        commands, but also for a curve deliberately easing to zero exit
        velocity at a node, which is a real, meaningful shape - not a
        straight line - so the whole-segment chord is too coarse an
        approximation of the true approach direction there. A short local
        secant tracks what the curve is actually doing near the boundary
        instead."""
        p0, p1, p2, p3 = seg
        t = 1.0 if at_end else 0.0
        d = Rail._deriv(p0, p1, p2, p3, t)
        if np.linalg.norm(d) < 1e-9:
            near_t = 0.9 if at_end else 0.1
            near_pt = np.array(bg.split_point(p0, p1, p2, p3, near_t))
            end_pt = np.array(p3) if at_end else np.array(p0)
            d = end_pt - near_pt if at_end else near_pt - end_pt
        return d

    def corner_cos_at_start(self):
        """cos(angle) between the tangent leaving segment[0] and the tangent
        entering segment[1], both evaluated at their shared node. 1.0 =
        perfectly smooth continuation (not a corner); near -1 to ~0.9 =
        a real direction change (a genuine taper corner)."""
        if self.n < 2:
            return 1.0
        out = self._boundary_tangent(self.segs[0], at_end=True)
        inn = self._boundary_tangent(self.segs[1], at_end=False)
        no, ni = np.linalg.norm(out), np.linalg.norm(inn)
        if no < 1e-9 or ni < 1e-9:
            return 1.0
        return float(np.dot(out, inn) / (no * ni))

    def corner_cos_at_end(self):
        if self.n < 2:
            return 1.0
        out = self._boundary_tangent(self.segs[-2], at_end=True)
        inn = self._boundary_tangent(self.segs[-1], at_end=False)
        no, ni = np.linalg.norm(out), np.linalg.norm(inn)
        if no < 1e-9 or ni < 1e-9:
            return 1.0
        return float(np.dot(out, inn) / (no * ni))

    def last_seg_length(self):
        return self.seg_bounds[-1] - self.seg_bounds[-2]

    def line_cross(self, anchor, direction, s_guess=None, _extrapolate=True, _local_window=None):
        """Where a fixed-direction line through anchor crosses this rail.
        Returns (point, s) or None."""
        L = 1e4
        line = (anchor - direction * L, anchor - direction * L,
                anchor + direction * L, anchor + direction * L)

        def segs_in_window(lo, hi):
            for i in range(self.n):
                if self.seg_bounds[i + 1] >= lo and self.seg_bounds[i] <= hi:
                    yield i

        def search(seg_indices):
            found = []
            for i in seg_indices:
                p0, p1, p2, p3 = self.segs[i]
                hits = bg.dedup_intersections(
                    bg.cubic_intersect(tuple(map(tuple, line)), (tuple(p0), tuple(p1), tuple(p2), tuple(p3)))
                )
                ts, arclen, _ = self._tables[i]
                for t1, t2, pt in hits:
                    s = self.seg_bounds[i] + float(np.interp(t2, ts, arclen))
                    found.append((s, np.array(pt)))
            return found

        candidates = []
        # Found live 2026-09-24 (g6): a rail that loops back near itself can
        # have a spatially-close-but-arclength-far crossing that wins purely
        # because it's the only one *found*, even with s_guess set (s_guess
        # only ranks candidates post-hoc, it can't rescue a search that
        # never looked locally). When s_guess is given, search a bounded
        # local window FIRST and use it if anything turns up there at all -
        # only fall through to the full rail (and extrapolation) when the
        # local window is genuinely empty.
        if s_guess is not None:
            window = _local_window if _local_window is not None else max(10.0 * (self.length / max(self.n, 1)), 0.15 * self.length)
            candidates = search(segs_in_window(s_guess - window, s_guess + window))
            if not candidates and _extrapolate:
                # Found live 2026-09-24: a short reference rail (e.g. the
                # unmodified side of a tight double taper) can have its true
                # geometric crossing fall just outside its own bounded end -
                # by construction (not a real gap in the object), not
                # because no reasonable crossing exists. cubic_intersect's
                # r1/r2 only relabel output parameters, they don't
                # extrapolate the actual curve - so build a real (straight)
                # extension segment beyond whichever boundary s_guess is
                # nearest, using the boundary tangent, and test against
                # that instead. Only that ONE boundary (not both) - s_guess
                # means we already know which end we're near.
                at_end = s_guess > self.length / 2
                i = self.n - 1 if at_end else 0
                seg = self.segs[i]
                d = self._boundary_tangent(seg, at_end)
                n = np.linalg.norm(d)
                if n > 1e-9:
                    u = d / n
                    MARGIN = 0.3 * self.length
                    boundary_pt = np.array(seg[3]) if at_end else np.array(seg[0])
                    far_pt = boundary_pt + u * MARGIN if at_end else boundary_pt - u * MARGIN
                    ext = (tuple(boundary_pt), tuple(boundary_pt), tuple(far_pt), tuple(far_pt)) if at_end \
                        else (tuple(far_pt), tuple(far_pt), tuple(boundary_pt), tuple(boundary_pt))
                    hits = bg.dedup_intersections(bg.cubic_intersect(tuple(map(tuple, line)), ext))
                    for t1, t2, pt in hits:
                        dist = np.linalg.norm(np.array(pt) - boundary_pt)
                        s = (self.seg_bounds[i + 1] + dist) if at_end else (self.seg_bounds[i] - dist)
                        candidates.append((s, np.array(pt)))
            # Deliberately NOT falling through to a blind full-rail search
            # here: s_guess means the caller has a specific local
            # expectation (a rung's local partner point), and a rail that
            # loops back near itself can otherwise "successfully" find a
            # technically-valid but arclength-nonsensical far crossing -
            # worse than just skipping this one optional rung.
            if not candidates:
                return None
            best = min(candidates, key=lambda c: abs(c[0] - s_guess))
            return best[1], best[0]

        candidates = search(range(self.n))
        if not candidates:
            return None
        best = min(candidates, key=lambda c: np.linalg.norm(c[1] - anchor))
        return best[1], best[0]


def pair_rails(railA, railB, A):
    """Match points along railA to points along railB, start to end, so the
    line joining each matched pair is as near square to BOTH rails as
    possible. Both rails only ever move forward (dynamic time warping), so
    two rungs taken from these pairs can never cross, and a rung can't jump
    to another turn of a rail that curls around. Expects the rails to run in
    the same direction (normalize_rail_directions). Returns a list of
    (sA, sB, pA, pB) from the start end to the finish end."""
    def samples(rail):
        m = int(rail.length / max(A * PAIR_SAMPLE_A_FRACTION, 1e-9)) + 1
        m = max(PAIR_MIN_SAMPLES, min(PAIR_MAX_SAMPLES, m))
        ss = np.linspace(0.0, rail.length, m)
        pts, tans = [], []
        for v in ss:
            p, t = rail.point_tangent(v)
            t = np.array(t, dtype=float)
            n = np.linalg.norm(t)
            pts.append(np.array(p, dtype=float))
            tans.append(t / n if n > 1e-12 else np.array([1.0, 0.0]))
        return ss, np.array(pts), np.array(tans)

    sA, PA, TA = samples(railA)
    sB, PB, TB = samples(railB)
    V = PB[None, :, :] - PA[:, None, :]
    D = np.linalg.norm(V, axis=2)
    U = V / np.maximum(D, 1e-9)[:, :, None]
    # (cos of the angle to each rail)^2: 0 when the rung is square to it.
    cost = (np.einsum('ijk,ik->ij', U, TA) ** 2 + np.einsum('ijk,jk->ij', U, TB) ** 2
            + PAIR_DIST_WEIGHT * (D / max(A, 1e-9)) ** 2)
    cost[D < 1e-6] = 0.0

    nA, nB = len(sA), len(sB)
    acc = np.full((nA, nB), np.inf)
    back = np.zeros((nA, nB), dtype=np.int8)  # 0 both moved, 1 only A, 2 only B
    acc[0, 0] = cost[0, 0]
    for i in range(nA):
        row, prev = acc[i], acc[i - 1] if i else None
        for j in range(nB):
            if i == 0 and j == 0:
                continue
            best, b = np.inf, 0
            if i and j and prev[j - 1] < best:
                best, b = prev[j - 1], 0
            if i and prev[j] < best:
                best, b = prev[j], 1
            if j and row[j - 1] < best:
                best, b = row[j - 1], 2
            row[j] = best + cost[i, j]
            back[i, j] = b
    pairs = []
    i, j = nA - 1, nB - 1
    while True:
        pairs.append((float(sA[i]), float(sB[j]), PA[i], PB[j]))
        if i == 0 and j == 0:
            break
        b = back[i, j]
        if b == 0:
            i, j = i - 1, j - 1
        elif b == 1:
            i -= 1
        else:
            j -= 1
    pairs.reverse()
    return pairs


def pair_progress(pairs):
    """Distance along the column at each pair from pair_rails, counted along
    whichever rail moved more at each step (the outside of a curve)."""
    progress = [0.0]
    for k in range(1, len(pairs)):
        progress.append(progress[-1] + max(abs(pairs[k][0] - pairs[k - 1][0]),
                                           abs(pairs[k][1] - pairs[k - 1][1])))
    return progress


def pair_of_rung(pairs, rung):
    """Index of the pair whose two points both lie closest to this rung's
    line segment: the rung's place along the column. Works on a curl, where a
    rung's end can sit next to a different turn of a rail."""
    a = np.array(rung[0][1], dtype=float)
    b = np.array(rung[-1][1], dtype=float)
    ab = b - a
    denom = float(np.dot(ab, ab)) or 1e-12

    def dist(P):
        t = np.clip((P - a) @ ab / denom, 0.0, 1.0)
        return np.linalg.norm(P - (a + np.outer(t, ab)), axis=1)

    PA = np.array([q[2] for q in pairs])
    PB = np.array([q[3] for q in pairs])
    return int(np.argmin(dist(PA) + dist(PB)))


def overshoot_uu(el):
    """RUNG_OVERSHOOT_MM in el's own path units: document mm to user units,
    then divided by the scale of every transform above the path.
    viewport_to_unit, not unittouu: unittouu falls back to px when the
    width/viewBox ratio is not a whole unit (e.g. 1 user unit = 0.5 mm)."""
    over = el.root.viewport_to_unit(f"{RUNG_OVERSHOOT_MM}mm")
    (a, c, _), (b, d, _) = el.composed_transform().matrix
    scale = math.sqrt(abs(a * d - b * c))
    return over / scale if scale > 1e-12 else over


def rung_between(pA, pB, over):
    """A rung subpath joining pA and pB, running `over` (path units) past
    each. None if the points coincide."""
    dvec = np.array(pB, dtype=float) - np.array(pA, dtype=float)
    n = np.linalg.norm(dvec)
    if n < 1e-9:
        return None
    u = dvec / n
    endA = np.array(pA, dtype=float) - u * over
    endB = np.array(pB, dtype=float) + u * over
    return [[list(endA)] * 3, [list(endB)] * 3]


def longer_first(csp0, csp1):
    """The two rails with the longer one first. The width (measure_A) is
    measured along the first rail, so a short inner rail hugging one turn of
    a long curl would otherwise read as a hairline width and pack the column
    with hundreds of rungs; the result must not depend on the path's order."""
    if Rail(csp1).length > Rail(csp0).length:
        return csp1, csp0
    return csp0, csp1


def measure_A(railA, railB):
    ts = np.linspace(0.15, 0.85, 200)
    pts_all_B = []
    for (p0, p1, p2, p3) in railB.segs:
        pts_all_B.append(np.array([bg.split_point(p0, p1, p2, p3, t) for t in np.linspace(0, 1, 200)]))
    pts_all_B = np.concatenate(pts_all_B, axis=0)
    dists = []
    for s in ts * railA.length:
        p, _ = railA.point_tangent(s)
        dists.append(np.min(np.linalg.norm(pts_all_B - np.array(p), axis=1)))
    return float(np.median(dists))


def stroke_width(el):
    style = el.style
    w = style.get("stroke-width")
    try:
        return float(w)
    except (TypeError, ValueError):
        return 0.254


def rail_polyline(rail, step):
    """Points along a rail, about `step` apart."""
    n = max(2, int(math.ceil(rail.length / step)) + 1)
    return np.array([rail.point_tangent(s)[0] for s in np.linspace(0.0, rail.length, n)], dtype=float)


def _ray_hits(origin, u, segs):
    """Signed distances along the line origin + t*u where it crosses segs
    (an (N, 2, 2) array), sorted."""
    p, q = segs[:, 0, :], segs[:, 1, :]
    e = q - p
    den = u[0] * e[:, 1] - u[1] * e[:, 0]
    ok = np.abs(den) > 1e-12
    w = p - origin
    t = np.where(ok, (w[:, 0] * e[:, 1] - w[:, 1] * e[:, 0]) / np.where(ok, den, 1.0), 0.0)
    v = np.where(ok, (w[:, 0] * u[1] - w[:, 1] * u[0]) / np.where(ok, den, 1.0), -1.0)
    hit = ok & (v >= 0.0) & (v <= 1.0)
    return np.unique(np.round(t[hit], 9))


def column_outline(railA, railB):
    """The column as one closed outline: railA out, railB back. Both rails
    must run the same way (normalize_rail_directions). Returns the outline
    points and, for each outline edge, which rail it belongs to: 0 for
    railA, 1 for railB, -1 for the two closing edges across the ends."""
    step = max(min(railA.length, railB.length) / 2000.0, 1e-6)
    pa = rail_polyline(railA, step)
    pb = rail_polyline(railB, step)[::-1]
    pts = np.concatenate([pa, pb], axis=0)
    owner = np.concatenate([np.zeros(len(pa) - 1, int), [-1], np.ones(len(pb) - 1, int), [-1]])
    along = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
    return pts, owner, along


def inside_outline(outline, pts):
    """Which of pts lie inside the closed outline (even-odd rule)."""
    pts = np.asarray(pts, dtype=float)
    a = outline[0]
    b = np.roll(a, -1, axis=0)
    x, y = pts[:, 0][:, None], pts[:, 1][:, None]
    ay, by = a[:, 1][None, :], b[:, 1][None, :]
    ax, bx = a[:, 0][None, :], b[:, 0][None, :]
    straddle = (ay > y) != (by > y)
    dy = np.where(by == ay, 1.0, by - ay)
    xcross = ax + (y - ay) * (bx - ax) / dy
    return (np.sum(straddle & (x < xcross), axis=1) % 2) == 1


def _hits_with_owner(outline, p, u):
    """(distance along p + t*u, rail) for every outline edge the line crosses."""
    pts, owner = outline[0], outline[1]
    q = np.roll(pts, -1, axis=0)
    e = q - pts
    den = u[0] * e[:, 1] - u[1] * e[:, 0]
    ok = np.abs(den) > 1e-12
    w = pts - p
    safe = np.where(ok, den, 1.0)
    t = (w[:, 0] * e[:, 1] - w[:, 1] * e[:, 0]) / safe
    v = (w[:, 0] * u[1] - w[:, 1] * u[0]) / safe
    hit = ok & (v >= 0.0) & (v <= 1.0)
    order = np.argsort(t[hit])
    return t[hit][order], owner[hit][order]


def rung_span(outline, pA, pB):
    """The rung for rail points pA (railA) and pB (railB), as its two rail
    points, or None if it does not cross the column cleanly.

    Normally nothing lies between pA and pB. At a curl-back a rail bends
    across the rung's way, so the rung meets that same rail a second time
    before it reaches the other one: the rung then ends at that second
    crossing. Anything else in the way (the other rail first, or the gap
    across the rail ends) fails, as does a stretch that is not inside the
    column."""
    pA = np.array(pA, dtype=float)
    pB = np.array(pB, dtype=float)
    n = np.linalg.norm(pB - pA)
    if n < 1e-9:
        return None
    u = (pB - pA) / n
    t, who = _hits_with_owner(outline, pA, u)
    eps = 0.02 * n
    between = (t > eps) & (t < n - eps)
    t, who = t[between], who[between]
    if len(t):
        if who[0] == 0:
            pB = pA + t[0] * u           # railA curls back: stop at it
            rail = 0
        elif who[-1] == 1:
            pA = pA + t[-1] * u          # railB curls back: start at it
            rail = 1
        else:
            return None
        m = np.linalg.norm(pB - pA)
        if m < 1e-9:
            return None
        # A real curl-back crosses the rail squarely enough at both ends; a
        # line that only grazes the outside of a bend meets it at a slant.
        pts, owner = outline[0], outline[1]
        idx = np.nonzero(owner == rail)[0]
        w = (pB - pA) / m
        for q in (pA, pB):
            i = idx[np.argmin(np.linalg.norm(pts[idx] - q, axis=1))]
            tan = pts[(i + 1) % len(pts)] - pts[i]
            tn = np.linalg.norm(tan)
            if tn > 1e-12 and abs(w[0] * tan[1] - w[1] * tan[0]) / tn < CURL_BACK_MIN_SIN:
                return None
        t2, _ = _hits_with_owner(outline, pA, (pB - pA) / m)
        if np.any((t2 > 0.02 * m) & (t2 < 0.98 * m)):
            return None
    ts = np.linspace(0.15, 0.85, 8)
    if not np.all(inside_outline(outline, pA[None, :] + ts[:, None] * (pB - pA)[None, :])):
        return None
    return pA, pB


def rung_in_column(outline, pA, pB):
    """True when rung_span finds a clean rung for these rail points."""
    return rung_span(outline, pA, pB) is not None


def clear_rung_ends(rungs, railA, railB):
    """Keep each rung's ends off the neighboring turn of a curl. A rung keeps
    the stretch between its own two rails; past each rail its end may run
    at most halfway across the gap to the next rail its line meets, so it
    crosses only its own two rails. Leaving enough gap for the ends is the
    designer's job (the help text suggests 0.6 mm)."""
    step = max(min(railA.length, railB.length) / 4000.0, 1e-6)
    pts = [rail_polyline(railA, step), rail_polyline(railB, step)]
    segs = np.concatenate([np.stack([p[:-1], p[1:]], axis=1) for p in pts], axis=0)
    out = []
    for r in rungs:
        e0 = np.array(r[0][1], dtype=float)
        e1 = np.array(r[-1][1], dtype=float)
        n = np.linalg.norm(e1 - e0)
        if len(r) != 2 or n < 1e-9:
            out.append(r)
            continue
        u = (e1 - e0) / n
        mid = (e0 + e1) / 2.0
        ts = _ray_hits(mid, u, segs)
        neg, pos = ts[ts < 0], ts[ts > 0]
        if not len(neg) or not len(pos):
            out.append(r)
            continue
        ends = []
        for side, hits in ((-1, -neg[::-1]), (1, pos)):
            own = hits[0]
            over = max(0.0, n / 2.0 - own)
            beyond = hits[1:][hits[1:] > own + 1e-9]
            if len(beyond):
                over = min(over, (beyond[0] - own) / 2.0)
            ends.append(mid + side * (own + over) * u)
        out.append([[list(ends[0])] * 3, [list(ends[1])] * 3])
    return out


def apply_rungs(path_el):
    """Given a path element with exactly 2 rail subpaths, build and append
    the rung geometry and mark it as an Ink/Stitch satin column. Shared by
    the Generate Ladder Rungs extension and Cross-Splice Rails (which feeds
    it a freshly-spliced 2-rail path)."""
    csp = path_el.path.to_absolute().to_superpath()
    if len(csp) != 2:
        raise inkex.AbortExtension(
            f"Expected exactly 2 rails (subpaths) on the combined path, found {len(csp)}."
        )

    csp0, csp1 = normalize_rail_directions(*longer_first(csp[0], csp[1]))
    railA = Rail(csp0)
    railB = Rail(csp1)
    B = stroke_width(path_el)
    A = measure_A(railA, railB)
    over = overshoot_uu(path_el)

    new_subpaths = []

    # ---- classify both ends ----
    end0 = _classify_end(railA, railB, at_start=True, A=A)
    end1 = _classify_end(railA, railB, at_start=False, A=A)

    for end in (end0, end1):
        if end["type"] == "taper" and end["dual"]:
            new_subpaths.extend(_build_dual_taper(end, A, B, over))
        elif end["type"] == "taper":
            new_subpaths.extend(_build_taper(end, A, B, over))
        elif end["type"] == "open":
            new_subpaths.extend(_build_open(end, A, B, over))
        # a point end's rungs come from _build_interior's pairing

    new_subpaths.extend(_build_interior(railA, railB, end0, end1, A, B, over))
    new_subpaths = clear_rung_ends(new_subpaths, railA, railB)

    path_el.path = make_csp(list(csp) + new_subpaths).to_path()

    # Mark as a satin column, same convention Ink/Stitch itself uses
    # (confirmed against its own bundled example files and observed
    # live behavior - inkex doesn't know the "inkstitch" prefix, so
    # this must use the fully-qualified Clark-notation name rather
    # than a bare "inkstitch:satin_column", or it writes the wrong
    # namespace silently).
    path_el.set("{http://inkstitch.org/namespace}satin_column", "True")
    path_el.set("{http://inkstitch.org/namespace}start_at_nearest_point", "False")
    path_el.set("{http://inkstitch.org/namespace}end_at_nearest_point", "False")


def resolve_path(sel):
    if len(sel) == 1 and isinstance(sel[0], PathElement):
        return sel[0]
    if len(sel) == 2:
        els = []
        for e in sel:
            if not isinstance(e, PathElement):
                if hasattr(e, "to_path_element"):
                    e = e.to_path_element()
                else:
                    raise inkex.AbortExtension("Both selected objects must be paths (or shapes convertible to paths).")
            els.append(e)
        els.sort(key=lambda e: list(e.getparent()).index(e))
        keep = els[-1]
        other = els[0]
        csp_keep = keep.path.to_absolute().to_superpath()
        csp_other = other.path.to_absolute().to_superpath()
        combined = make_csp(list(csp_other) + list(csp_keep))
        keep.path = combined.to_path()
        other.delete()
        return keep
    raise inkex.AbortExtension(
        "Select either one combined 2-rail path, or exactly 2 separate rail paths to combine."
    )


def _classify_end(railA, railB, at_start, A):
        pA = railA.start() if at_start else railA.end()
        pB = railB.start() if at_start else railB.end()
        gap = np.linalg.norm(pA - pB)

        # Candidate taper segment per rail = the segment adjacent to THIS
        # end (first segment for a start, last for an end). The real signal
        # for "is this actually a taper" is a genuine tangent discontinuity
        # (a corner) where it meets the rest of the rail - NOT segment
        # length. A multi-segment rail is often just one smooth curve an
        # artist split into several bezier pieces; those joins are
        # tangent-continuous and must NOT be mistaken for a taper, even if
        # one piece happens to be short. Length is kept as a secondary
        # sanity check only, on top of the corner requirement. If both
        # rails have a genuine corner here (a taper can converge from both
        # sides at once), prefer the sharper corner as the reference for
        # direction-setting - the crossing search itself still walks every
        # segment of the other rail regardless, so this choice only affects
        # which rail sets the "parallel" reference direction.
        # Confirmed real tapers seen so far range cos 0.42-0.923 (sharp to
        # fairly gentle); confirmed smooth non-taper joints have all been
        # 0.9998+ (essentially parallel). 0.97 sits safely in that gap.
        CORNER_COS_THRESHOLD = 0.97
        candidates = []
        for name, rail, other in (("A", railA, railB), ("B", railB, railA)):
            if rail.n < 2:
                continue
            corner_cos = rail.corner_cos_at_start() if at_start else rail.corner_cos_at_end()
            if corner_cos > CORNER_COS_THRESHOLD:
                continue  # smooth continuation, not a real taper corner
            seg_len = (rail.seg_bounds[1] - rail.seg_bounds[0]) if at_start else rail.last_seg_length()
            # A single sharp corner is visible from BOTH ends of a 2-segment
            # rail (each end's classification checks the same joint from its
            # own side), so a real taper must ALSO be disambiguated from the
            # long main body sharing that same corner. Comparing length to
            # the rail's own total length is too strict (a near-50/50 split
            # taper is still a real taper) - compare to A instead, since a
            # taper converges over a distance on the order of the ladder's
            # width, regardless of how long the rest of the rail is.
            if seg_len > TAPER_MAX_A_MULT * A:
                continue
            candidates.append((corner_cos, name, seg_len))
        taper_rail, taper_len = None, None
        if candidates:
            candidates.sort(key=lambda c: c[0])  # sharpest corner (lowest cos) first
            _, taper_rail, taper_len = candidates[0]

        is_taper = gap < 0.35 * A and taper_rail is not None
        # Both rails can genuinely taper into the same point at once. Rather
        # than pick one as "the" reference and treat the other as if it
        # were the plain rail (which it isn't), draw rungs across a dual
        # taper by connecting proportional points along each taper's own
        # chord directly - i.e. treating both as straight regardless of
        # their actual curvature, since with two independently-converging
        # tapers there's no single well-defined perpendicular reference
        # direction to fall back on anyway.
        is_dual = is_taper and len(candidates) == 2
        # The rails meet (or nearly) without a corner: a smooth point, like
        # the tail of a curl. The paired interior rungs run into it.
        is_point = not is_taper and gap < 0.35 * A
        return {
            "type": "taper" if is_taper else ("point" if is_point else "open"),
            "dual": is_dual,
            "at_start": at_start,
            "railA": railA, "railB": railB,
            "pA": pA, "pB": pB, "gap": gap,
            "taper_rail": taper_rail, "taper_len": taper_len,
            "candidates": candidates,
        }

def _rail_s_for_end(rail, at_start):
    return 0.0 if at_start else rail.length

def _step_from_end(rail, at_start, delta):
    s0 = _rail_s_for_end(rail, at_start)
    return s0 + delta if at_start else s0 - delta

def _build_taper(end, A, B, over):
    railA, railB, at_start = end["railA"], end["railB"], end["at_start"]
    which_taper = end["taper_rail"]
    taper_rail = railA if which_taper == "A" else railB
    plain_rail = railB if which_taper == "A" else railA
    taper_len = end["taper_len"]
    subpaths = []

    transition_s = taper_rail.seg_bounds[1] if at_start else taper_rail.seg_bounds[-2]
    # "before" = still on the main-body/non-taper side of the transition;
    # for at_start the taper occupies s in [0, transition_s] with the
    # point at s=0, so main-body is s>transition_s. For at_end the taper
    # occupies s in [transition_s, length] with the point at s=length,
    # so main-body is s<transition_s. "after" is the mirror of this.
    s_before = transition_s + (1.0 if at_start else -1.0)
    pT, _ = taper_rail.point_tangent(s_before)
    # Nearest point on plain_rail to set the perpendicular reference.
    # Window the search to the end plain_rail corresponds to (at_start ->
    # its own low-s end, at_end -> its own high-s end) instead of scanning
    # its FULL length - found live 2026-09-24 (g6): a rail that loops back
    # near itself can have a point near one end sit spatially close to a
    # point near the OTHER end, so a blind whole-rail nearest-Euclidean
    # search picks the wrong (far-in-arclength) region, seeding every
    # downstream line_cross call in this function with a bad locality no
    # s_guess can recover from. A generous window (relative to A, the rail
    # spacing) stays local while still comfortably covering any real taper.
    WINDOW = max(10.0 * A, 0.25 * plain_rail.length)
    if at_start:
        s_lo, s_hi = 0.0, min(plain_rail.length, WINDOW)
    else:
        s_lo, s_hi = max(0.0, plain_rail.length - WINDOW), plain_rail.length
    best_s, best_d = None, None
    for i in range(60):
        s = s_lo + (s_hi - s_lo) * i / 59
        p, _ = plain_rail.point_tangent(s)
        d = np.linalg.norm(np.array(p) - np.array(pT))
        if best_d is None or d < best_d:
            best_d, best_s = d, s
    pP, tanP = plain_rail.point_tangent(best_s)
    ref_dir = unit_normal(tanP)

    def emit(pA_taper, pA_plain, offset_dir):
        # Symmetric overshoot: each side runs `over` (RUNG_OVERSHOOT_MM in
        # path units) past the crossing point it connects - centered over
        # the column, not the old asymmetric taper/plain split.
        d = np.array(pA_taper) - np.array(pA_plain)
        n = np.linalg.norm(d)
        u = d / n if n > 1e-9 else np.array([1.0, 0.0])
        overshoot = over
        end_plain = np.array(pA_plain) - u * overshoot
        end_taper = np.array(pA_taper) + u * overshoot
        pt1, pt2 = (end_taper, end_plain) if which_taper == "A" else (end_plain, end_taper)
        subpaths.append([[list(pt1), list(pt1), list(pt1)], [list(pt2), list(pt2), list(pt2)]])

    res = plain_rail.line_cross(pT, ref_dir, s_guess=best_s)
    if res is not None:
        pP_before, _ = res
        emit(pT, pP_before, ref_dir)

    s_after = transition_s + (-1.0 if at_start else 1.0)
    p_after, _ = taper_rail.point_tangent(s_after)
    # s_guess=best_s throughout (found live 2026-09-24, g6): without it,
    # line_cross falls back to Euclidean-nearest, which for a rail that
    # loops back near itself can find a far-away-in-arclength crossing
    # that just happens to sit physically close - a long spurious rung
    # across the whole loop instead of the local one.
    res = plain_rail.line_cross(p_after, ref_dir, s_guess=best_s)
    if res is not None:
        pP_after, _ = res
        emit(p_after, pP_after, ref_dir)

    if taper_len > TAPER_EXTRA_RUNG_MARGIN * A:
        s_mid = (transition_s + _rail_s_for_end(taper_rail, at_start)) / 2.0
        p_mid, _ = taper_rail.point_tangent(s_mid)
        res = plain_rail.line_cross(p_mid, ref_dir, s_guess=best_s)
        if res is not None:
            pP_mid, _ = res
            emit(p_mid, pP_mid, ref_dir)

    # Near-the-convergence-point rung: the sanctioned exception to "never
    # at the end" / the 0.5*A end-inset rule - an extreme taper still
    # needs a rung close to its tip to keep the zigzag perpendicular to
    # the rails there, so this one is governed by its own minimum
    # clearance (MIN_FROM_POINT_MULT * B) instead of END_INSET_MULT * A.
    min_from_point = MIN_FROM_POINT_MULT * B
    point = (railA.start() if at_start else railA.end())  # approx convergence
    offset = min_from_point
    for _ in range(400):
        s_try = _step_from_end(taper_rail, at_start, offset)
        p_try, _ = taper_rail.point_tangent(s_try)
        res = plain_rail.line_cross(p_try, ref_dir, s_guess=best_s)
        if res is not None:
            pP_try, _ = res
            if (np.linalg.norm(np.array(p_try) - end["pA" if which_taper == "A" else "pB"])
                    >= min_from_point and
                    np.linalg.norm(np.array(pP_try) - end["pB" if which_taper == "A" else "pA"])
                    >= min_from_point):
                emit(p_try, pP_try, ref_dir)
                break
        offset += 0.02

    return subpaths

def _build_dual_taper(end, A, B, over):
    """Both rails genuinely taper into the same point. There's no
    single well-defined plain-rail tangent to use as a perpendicular
    reference here, so instead: treat each taper as a straight chord
    (point end to its own transition end, regardless of its actual
    curvature) and connect the two chords at matching fractions along
    their length - "as if both are straight". Fraction 0 = at the
    convergence point, fraction 1 = at each rail's own transition into
    its main body."""
    railA, railB, at_start = end["railA"], end["railB"], end["at_start"]
    cand = {name: seg_len for _, name, seg_len in end["candidates"]}
    len_A, len_B = cand["A"], cand["B"]

    def chord_for(rail, taper_len):
        if at_start:
            p_point, _ = rail.point_tangent(0.0)
            p_trans, _ = rail.point_tangent(taper_len)
        else:
            p_point, _ = rail.point_tangent(rail.length)
            p_trans, _ = rail.point_tangent(rail.length - taper_len)
        return np.array(p_point), np.array(p_trans)

    point_a, trans_a = chord_for(railA, len_A)
    point_b, trans_b = chord_for(railB, len_B)

    taper_len = max(len_A, len_B)
    n_rungs = 4 if taper_len > TAPER_EXTRA_RUNG_MARGIN * A else 3

    min_from_point = MIN_FROM_POINT_MULT * B
    f_near = min(max(min_from_point / len_A, min_from_point / len_B), 0.9)
    fractions = np.linspace(f_near, 1.0, n_rungs)

    subpaths = []
    for f in fractions:
        pa = point_a + (trans_a - point_a) * f
        pb = point_b + (trans_b - point_b) * f
        d = pb - pa
        n = np.linalg.norm(d)
        if n < 1e-9:
            continue
        u = d / n
        overshoot = over  # symmetric, RUNG_OVERSHOOT_MM past each rail
        end_a = pa - u * overshoot
        end_b = pb + u * overshoot
        subpaths.append([[list(end_a)] * 3, [list(end_b)] * 3])
    return subpaths

def _build_open(end, A, B, over):
    railA, railB, at_start = end["railA"], end["railB"], end["at_start"]
    # Never-at-the-end rule: inset the last rung 0.5 * column width (A,
    # rail separation) from the rail end, not a stroke-width-based margin.
    pullback = END_INSET_MULT * A
    n0 = np.linalg.norm(np.array(end["pB"]) - np.array(end["pA"]))
    if n0 < 1e-9:
        return []
    # Search is bounded by railA's own remaining length at this end - a
    # short open segment (e.g. a curved stub) has no length to walk past,
    # and stepping beyond it previously ran the whole 600-step search
    # against an out-of-range s with no chance of success.
    max_offset = railA.length

    def try_offset(offset, require_B_pullback):
        if offset > max_offset:
            return None
        s_try = _step_from_end(railA, at_start, offset)
        pA, tanA = railA.point_tangent(s_try)
        # Local perpendicular at THIS point, not a single direction fixed
        # at the raw endpoints - same reasoning as _build_taper's ref_dir
        # (line ~526): on a curved open end the true perpendicular rotates
        # as you step inward, so a fixed endpoint-to-endpoint direction
        # can miss railB entirely even though the rails are genuinely
        # running parallel to each other locally.
        local_dir = unit_normal(tanA)
        # same s_guess reasoning as _build_interior/_build_taper: bias
        # toward railB's locally-corresponding position, not whatever
        # happens to be Euclidean-nearest on a self-proximate rail.
        res = railB.line_cross(pA, local_dir, s_guess=_step_from_end(railB, at_start, offset))
        if res is None:
            return None
        pB, sB = res
        s_on_B_from_end = sB if at_start else (railB.length - sB)
        if require_B_pullback and s_on_B_from_end < pullback:
            return None
        d = np.array(pB) - np.array(pA)
        n = np.linalg.norm(d)
        if n < 1e-9:
            return None
        u = d / n
        # Symmetric overshoot: RUNG_OVERSHOOT_MM past each rail.
        overshoot = over
        endA = np.array(pA) - u * overshoot
        endB = np.array(pB) + u * overshoot
        return [[[list(endA)] * 3, [list(endB)] * 3]]

    # Primary search: require both rails to clear the pullback.
    offset = pullback
    for _ in range(600):
        hit = try_offset(offset, require_B_pullback=True)
        if hit is not None:
            return hit
        offset += 0.02

    # Fallback: this end must always get a rung - relax the requirement
    # that railB also clears the pullback (its crossing may just be
    # geometrically closer to its own endpoint here) rather than
    # silently placing nothing.
    offset = pullback
    for _ in range(600):
        hit = try_offset(offset, require_B_pullback=False)
        if hit is not None:
            return hit
        offset += 0.02

    raise inkex.AbortExtension(
        "Could not place the required rung on a non-tapered end "
        f"(near {tuple(round(x, 3) for x in end['pA'])}) - the rails may not "
        "actually run roughly parallel there. Check this end by hand."
    )
    return []

def _build_interior(railA, railB, end0, end1, A, B, over):
    """Rungs along the body of the column, from the both-rails pairing
    (pair_rails), spaced along whichever rail is moving faster at each point
    (the outside of a curve). Taper and open ends keep the zone their own
    builders use; a point end's rungs run right down into the tip."""
    lo = 0.0
    hi = railA.length
    for end in (end0, end1):
        if end["type"] == "point":
            continue
        if end["type"] == "taper":
            if end["taper_rail"] == "A":
                # reserved zone measured on railA directly (it's the taper rail)
                reserved = end["taper_len"] + 2.0 * B
            else:
                # taper is on railB; approximate railA's reserve via the
                # "before" rung's crossing distance from this end (rough)
                reserved = A * 0.5 + 2.0 * B
        else:
            # Keep interior rungs clear of the open-end rung, which now
            # sits at END_INSET_MULT * A from the rail end.
            reserved = END_INSET_MULT * A + 1.0
        if end["at_start"]:
            lo = max(lo, reserved)
        else:
            hi = min(hi, railA.length - reserved)

    usable = hi - lo
    if usable <= 0:
        return []
    chord = np.linalg.norm(np.array(railA.point_tangent(hi)[0]) - np.array(railA.point_tangent(lo)[0]))
    straightness = chord / usable if usable > 0 else 1.0
    spacing = (SPACING_GENTLE_MULT if straightness > STRAIGHTNESS_THRESHOLD else SPACING_TIGHT_MULT) * A

    pairs = pair_rails(railA, railB, A)
    progress = pair_progress(pairs)

    near = MIN_FROM_POINT_MULT * B
    min_w = POINT_MIN_WIDTH_MULT * B
    start_point = end0["type"] == "point"
    end_point = end1["type"] == "point"

    def usable_pair(k):
        sA, sB, pA, pB = pairs[k]
        if not (lo <= sA <= hi):
            return False
        if (start_point or end_point) and np.linalg.norm(np.array(pB) - np.array(pA)) < min_w:
            return False
        return True

    outline = column_outline(railA, railB)

    def rung_points(k):
        """Pair k's rail points, pulled back off a point end's tip on each rail."""
        sA, sB, pA, pB = pairs[k]
        if start_point:
            if sA < near:
                pA = railA.point_tangent(near)[0]
            if sB < near:
                pB = railB.point_tangent(near)[0]
        if end_point:
            if railA.length - sA < near:
                pA = railA.point_tangent(railA.length - near)[0]
            if railB.length - sB < near:
                pB = railB.point_tangent(railB.length - near)[0]
        return pA, pB

    # A rung must cross the column itself: near a curl's center the pairing
    # can reach across the gap to the next turn, and those pairs are skipped;
    # at a curl-back the rung ends where its own rail crosses its way again.
    spans = {}
    for k in range(len(pairs)):
        if usable_pair(k):
            sp = rung_span(outline, *rung_points(k))
            if sp is not None:
                spans[k] = sp

    def rung_at(k):
        return rung_between(*spans[k], over)

    ks = sorted(spans)
    if not ks:
        return []
    k0, k1 = ks[0], ks[-1]

    chosen = []
    if start_point:
        chosen.append(k0)  # the rung nearest the tail tip
    target = progress[k0] + spacing
    stop = progress[k1] - spacing * 0.5
    for k in ks:
        if progress[k] >= target and progress[k] < stop:
            chosen.append(k)
            target = progress[k] + spacing
    if end_point:
        # the rung nearest this tip; drop a regular rung crowding it
        if chosen and progress[k1] - progress[chosen[-1]] < spacing * 0.5 and chosen[-1] != k0:
            chosen.pop()
        chosen.append(k1)

    if not chosen:
        # Short object: normal spacing didn't fit even one rung between the
        # two end zones. Guarantee one, halfway along the usable span, so no
        # rail stretch is left with no cross-ties at all.
        mid = (progress[k0] + progress[k1]) / 2.0
        chosen = [min(ks, key=lambda k: abs(progress[k] - mid))]

    subpaths = []
    for k in chosen:
        r = rung_at(k)
        if r is not None:
            subpaths.append(r)
    return subpaths


class RungRailSatin(EffectExtension):

    def effect(self):
        sel = list(self.svg.selection.get(ShapeElement))
        path_el = resolve_path(sel)
        apply_rungs(path_el)


if __name__ == "__main__":
    RungRailSatin().run()
