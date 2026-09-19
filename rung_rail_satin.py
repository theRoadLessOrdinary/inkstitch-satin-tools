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
 - Taper rungs overshoot 1*B on the non-taper rail's crossing, and 8*B on
   the taper rail's crossing (so the node marker there doesn't hide it).
 - An end with no taper on either rail (gap there is close to A, not a
   near-convergence) gets one final rung parallel to the straight line
   between the two raw endpoints, pulled back >=4*B from them on both
   rails, overshooting 4*B on both sides.
 - Interior spacing: 2*A on gentle curves/straights, 1*A on tighter curves,
   default 1*B overshoot both sides.

This is a from-scratch implementation of that rule set, checked against a
small number of real objects and not independently verified beyond that. Treat its output as a draft to
check, especially the taper-vs-open classification and the 3-vs-4 rung
threshold, both of which are heuristics.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import inkex
from inkex import EffectExtension, ShapeElement, PathElement
from inkex.paths import CubicSuperPath

import bezier_geom as bg

STROKE_TAPER_MULT = 4.0
STROKE_DEFAULT_MULT = 1.0
STROKE_OPEN_MULT = 3.0
MIN_FROM_POINT_MULT = 2.0
TAPER_EXTRA_RUNG_MARGIN = 1.15
TAPER_MAX_A_MULT = 3.0
SPACING_GENTLE_MULT = 2.0
SPACING_TIGHT_MULT = 1.0
STRAIGHTNESS_THRESHOLD = 0.9
N_PER_SEG = 400


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
        return bg.split_point(p0, p1, p2, p3, t), bg.bez_deriv(p0, p1, p2, p3, t) \
            if hasattr(bg, "bez_deriv") else self._deriv(p0, p1, p2, p3, t)

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

    def line_cross(self, anchor, direction, s_guess=None):
        """Where a fixed-direction line through anchor crosses this rail.
        Returns (point, s) or None."""
        L = 1e4
        line = (anchor - direction * L, anchor - direction * L,
                anchor + direction * L, anchor + direction * L)
        candidates = []
        for i, (p0, p1, p2, p3) in enumerate(self.segs):
            hits = bg.dedup_intersections(
                bg.cubic_intersect(tuple(map(tuple, line)), (tuple(p0), tuple(p1), tuple(p2), tuple(p3)))
            )
            ts, arclen, _ = self._tables[i]
            for t1, t2, pt in hits:
                s = self.seg_bounds[i] + float(np.interp(t2, ts, arclen))
                candidates.append((s, np.array(pt)))
        if not candidates:
            return None
        if s_guess is None:
            best = min(candidates, key=lambda c: np.linalg.norm(c[1] - anchor))
        else:
            best = min(candidates, key=lambda c: abs(c[0] - s_guess))
        return best[1], best[0]


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

    csp0, csp1 = normalize_rail_directions(csp[0], csp[1])
    railA = Rail(csp0)
    railB = Rail(csp1)
    B = stroke_width(path_el)
    A = measure_A(railA, railB)

    new_subpaths = []

    # ---- classify both ends ----
    end0 = _classify_end(railA, railB, at_start=True, A=A)
    end1 = _classify_end(railA, railB, at_start=False, A=A)

    for end in (end0, end1):
        if end["type"] == "taper" and end["dual"]:
            new_subpaths.extend(_build_dual_taper(end, A, B))
        elif end["type"] == "taper":
            new_subpaths.extend(_build_taper(end, A, B))
        else:
            new_subpaths.extend(_build_open(end, A, B))

    new_subpaths.extend(_build_interior(railA, railB, end0, end1, A, B))

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
        return {
            "type": "taper" if is_taper else "open",
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

def _build_taper(end, A, B):
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
    # nearest point on plain_rail to set the perpendicular reference
    best_s, best_d = None, None
    for i in range(60):
        s = plain_rail.length * i / 59
        p, _ = plain_rail.point_tangent(s)
        d = np.linalg.norm(np.array(p) - np.array(pT))
        if best_d is None or d < best_d:
            best_d, best_s = d, s
    pP, tanP = plain_rail.point_tangent(best_s)
    ref_dir = unit_normal(tanP)

    def emit(pA_taper, pA_plain, overshoot_taper, overshoot_plain, offset_dir):
        d = np.array(pA_taper) - np.array(pA_plain)
        n = np.linalg.norm(d)
        u = d / n if n > 1e-9 else np.array([1.0, 0.0])
        end_plain = np.array(pA_plain) - u * overshoot_plain
        end_taper = np.array(pA_taper) + u * overshoot_taper
        pt1, pt2 = (end_taper, end_plain) if which_taper == "A" else (end_plain, end_taper)
        subpaths.append([[list(pt1), list(pt1), list(pt1)], [list(pt2), list(pt2), list(pt2)]])

    res = plain_rail.line_cross(pT, ref_dir, s_guess=best_s)
    if res is not None:
        pP_before, _ = res
        emit(pT, pP_before, B * STROKE_TAPER_MULT, B * STROKE_DEFAULT_MULT, ref_dir)

    s_after = transition_s + (-1.0 if at_start else 1.0)
    p_after, _ = taper_rail.point_tangent(s_after)
    res = plain_rail.line_cross(p_after, ref_dir)
    if res is not None:
        pP_after, _ = res
        emit(p_after, pP_after, B * STROKE_TAPER_MULT, B * STROKE_DEFAULT_MULT, ref_dir)

    if taper_len > TAPER_EXTRA_RUNG_MARGIN * A:
        s_mid = (transition_s + _rail_s_for_end(taper_rail, at_start)) / 2.0
        p_mid, _ = taper_rail.point_tangent(s_mid)
        res = plain_rail.line_cross(p_mid, ref_dir)
        if res is not None:
            pP_mid, _ = res
            emit(p_mid, pP_mid, B * STROKE_TAPER_MULT, B * STROKE_DEFAULT_MULT, ref_dir)

    min_from_point = MIN_FROM_POINT_MULT * B
    point = (railA.start() if at_start else railA.end())  # approx convergence
    offset = min_from_point
    for _ in range(400):
        s_try = _step_from_end(taper_rail, at_start, offset)
        p_try, _ = taper_rail.point_tangent(s_try)
        res = plain_rail.line_cross(p_try, ref_dir)
        if res is not None:
            pP_try, _ = res
            if (np.linalg.norm(np.array(p_try) - end["pA" if which_taper == "A" else "pB"])
                    >= min_from_point and
                    np.linalg.norm(np.array(pP_try) - end["pB" if which_taper == "A" else "pA"])
                    >= min_from_point):
                emit(p_try, pP_try, B * STROKE_TAPER_MULT, B * STROKE_DEFAULT_MULT, ref_dir)
                break
        offset += 0.02

    return subpaths

def _build_dual_taper(end, A, B):
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

    overshoot = STROKE_TAPER_MULT * B
    subpaths = []
    for f in fractions:
        pa = point_a + (trans_a - point_a) * f
        pb = point_b + (trans_b - point_b) * f
        d = pb - pa
        n = np.linalg.norm(d)
        if n < 1e-9:
            continue
        u = d / n
        end_a = pa - u * overshoot
        end_b = pb + u * overshoot
        subpaths.append([[list(end_a)] * 3, [list(end_b)] * 3])
    return subpaths

def _build_open(end, A, B):
    railA, railB, at_start = end["railA"], end["railB"], end["at_start"]
    pullback = STROKE_OPEN_MULT * B
    overshoot = STROKE_OPEN_MULT * B
    open_dir = np.array(end["pB"]) - np.array(end["pA"])
    n = np.linalg.norm(open_dir)
    if n < 1e-9:
        return []
    open_dir = open_dir / n

    def try_offset(offset, require_B_pullback):
        s_try = _step_from_end(railA, at_start, offset)
        pA, _ = railA.point_tangent(s_try)
        res = railB.line_cross(pA, open_dir)
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

def _build_interior(railA, railB, end0, end1, A, B):
    lo = 0.0
    hi = railA.length
    for end in (end0, end1):
        reserved = None
        if end["type"] == "taper":
            if end["taper_rail"] == "A":
                # reserved zone measured on railA directly (it's the taper rail)
                reserved = end["taper_len"] + 2.0 * B
            else:
                # taper is on railB; approximate railA's reserve via the
                # "before" rung's crossing distance from this end (rough)
                reserved = A * 0.5 + 2.0 * B
        else:
            reserved = STROKE_OPEN_MULT * B + 1.0
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

    subpaths = []
    s = lo + spacing
    while s < hi - spacing * 0.5:
        pA, tanA = railA.point_tangent(s)
        d = unit_normal(tanA)
        res = railB.line_cross(pA, d)
        if res is not None:
            pB, _ = res
            dvec = np.array(pB) - np.array(pA)
            u = dvec / np.linalg.norm(dvec)
            endA = np.array(pA) - u * (STROKE_DEFAULT_MULT * B)
            endB = np.array(pB) + u * (STROKE_DEFAULT_MULT * B)
            subpaths.append([[list(endA)] * 3, [list(endB)] * 3])
        s += spacing

    if not subpaths:
        # Short object: normal spacing didn't fit even one rung between
        # the two taper/open zones (e.g. an object that ends up as
        # taper-rungs-only). Guarantee at least one, at the midpoint of
        # the usable span, so no rail segment is left with zero
        # cross-ties at all.
        s_mid = (lo + hi) / 2.0
        pA, tanA = railA.point_tangent(s_mid)
        d = unit_normal(tanA)
        res = railB.line_cross(pA, d)
        if res is not None:
            pB, _ = res
            dvec = np.array(pB) - np.array(pA)
            n = np.linalg.norm(dvec)
            if n > 1e-9:
                u = dvec / n
                endA = np.array(pA) - u * (STROKE_DEFAULT_MULT * B)
                endB = np.array(pB) + u * (STROKE_DEFAULT_MULT * B)
                subpaths.append([[list(endA)] * 3, [list(endB)] * 3])

    return subpaths


class RungRailSatin(EffectExtension):

    def effect(self):
        sel = list(self.svg.selection.get(ShapeElement))
        path_el = resolve_path(sel)
        apply_rungs(path_el)


if __name__ == "__main__":
    RungRailSatin().run()
