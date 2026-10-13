#!/usr/bin/env python3
"""
Taper Join at Crossing (Ink/Stitch satin column helper)

Follow-on to Touch Rails at Crossing. Touching each of the under piece's
two rail tips to the over column INDEPENDENTLY (Touch Rails' behavior)
looks fine at a near-tangent graze, but at any real crossing angle the
two independent touch points land noticeably far apart along the over
column's rail - not the clean pointed tip a real satin crossing should
have. Where the crossing angle is more than ~5 degrees off a graze, the
under piece's two rails should instead taper together to ONE shared
point on the over column, the same way a satin column's own natural end
tapers to a point (see rung_rail_satin._build_taper /
_build_dual_taper, which this reuses the rung-placement idea from) -
except here the two rails end at that shared point as separate,
merely-incident subpaths (their endpoints coincide) rather than being
fused into one continuous/closed path, plus one or two new rungs near
the tip.

Select exactly 2 satin-column paths at one crossing, same as Touch
Rails - under/over and which end needs fixing are auto-detected the
same way (proximity). If the crossing angle is <= 5 degrees, this just
falls back to Touch Rails' independent-touch behavior (a near-graze
doesn't need a converging taper).

Refinement from the user: the new rungs must run PARALLEL to the first
of the two selected paths (selection order - the same click-order
assumption already relied on elsewhere in this extension set, e.g.
Cross-Splice Rails' --id handling, flagged there as unverified) "if it
can be determined" - i.e. this is a best-effort reference direction, not
a guaranteed one. Concretely: take the first selected path's rail
nearest the convergence point, use ITS local tangent there as a single
fixed direction, and build each rung as a line in that exact direction,
clipped between the under piece's two taper chords - mirroring how
rung_rail_satin._build_taper holds ONE fixed reference direction for an
entire taper zone rather than recomputing per rung. If that direction
turns out to be parallel to the taper chords themselves (degenerate, no
intersection), each rung falls back to a direct connection between its
two chord points instead.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import inkex
from inkex import EffectExtension, ShapeElement, PathElement

import rung_rail_satin as rrs
import touch_rails as tr

ANGLE_THRESHOLD_DEG = 5.0
TAPER_LEN_MULT = 1.5          # how far back each rail tapers, x the rail separation (A) -
                               # real tapers elsewhere in this codebase run roughly 1-2x A
                               # (rung_rail_satin's own TAPER_MAX_A_MULT=3.0 is the upper
                               # bound for something to still count as a taper, not a target)
MAX_TAPER_LEN_FRAC = 0.4       # never consume more than this fraction of a rail's own length
MIN_FROM_POINT_MULT = rrs.MIN_FROM_POINT_MULT
STROKE_TAPER_MULT = rrs.STROKE_TAPER_MULT
N_RUNGS = 2


def _crossing_angle_deg(under_tangent, over_tangent):
    u = under_tangent / np.linalg.norm(under_tangent)
    v = over_tangent / np.linalg.norm(over_tangent)
    cos_a = abs(float(np.dot(u, v)))
    cos_a = min(1.0, max(-1.0, cos_a))
    return np.degrees(np.arccos(cos_a))


def _nearest_tangent_on_rails(point, rails):
    """Tangent of whichever of `rails` passes closest to `point`."""
    point = np.array(point, dtype=float)
    best = None
    for rail, _ in rails:
        for i, (ts, arclen, pts) in enumerate(rail._tables):
            d = np.linalg.norm(pts - point, axis=1)
            j = int(np.argmin(d))
            if best is None or d[j] < best[0]:
                s = rail.seg_bounds[i] + float(arclen[j])
                _, tan = rail.point_tangent(s)
                best = (d[j], tan)
    return best[1]


def _crop_to_s(subpath, rail, at_end, s_cut):
    """Trim `subpath` back to arc-length s_cut (keeping the body-side
    portion), via De Casteljau splitting - reuses the same approach as
    weave_paths.crop_subpath for a single-sided crop."""
    s_cut = max(0.0, min(s_cut, rail.length))
    i, local_s = rail.seg_index_at(s_cut)
    ts, arclen, _ = rail._tables[i]
    t = float(np.interp(local_s, arclen, ts))

    import bezier_geom as bg
    left, right = bg.split_cubic(*rail.segs[i], t)

    sp = [list(map(list, node)) for node in subpath]
    if at_end:
        # keep segments [0, i-1] intact, replace segment i with `left`
        keep = sp[: i + 1]
        keep[-1][2] = list(left[1])
        keep.append([list(left[2]), list(left[3]), list(left[3])])
        return keep, left[3]
    else:
        keep = sp[i + 1:]
        keep[0][0] = list(right[2])
        keep.insert(0, [list(right[0]), list(right[0]), list(right[1])])
        return keep, right[0]


def _append_point(subpath, at_end, point):
    sp = [list(map(list, node)) for node in subpath]
    if at_end:
        sp[-1][2] = list(sp[-1][1])
        sp.append([list(point), list(point), list(point)])
    else:
        sp[0][0] = list(sp[0][1])
        sp.insert(0, [list(point), list(point), list(point)])
    return sp


def _line_line_intersect(p, d, a, b):
    """Point where the infinite line (p, direction d) crosses the finite
    segment a-b, or None if parallel/no crossing."""
    p, d, a, b = (np.array(x, dtype=float) for x in (p, d, a, b))
    seg = b - a
    denom = d[0] * seg[1] - d[1] * seg[0]
    if abs(denom) < 1e-9:
        return None
    diff = a - p
    t = (diff[0] * seg[1] - diff[1] * seg[0]) / denom
    u = (diff[0] * d[1] - diff[1] * d[0]) / denom
    if u < -0.05 or u > 1.05:
        return None
    return p + d * t


def _build_taper_join_rungs(conv_pt, trans_a, trans_b, ref_dir, overshoot):
    """1-2 rungs between the two straight taper chords (conv_pt->trans_a,
    conv_pt->trans_b), each drawn parallel to ref_dir where that's
    possible, falling back to a direct connection otherwise."""
    min_from_point = MIN_FROM_POINT_MULT
    len_a = np.linalg.norm(np.array(trans_a) - np.array(conv_pt))
    len_b = np.linalg.norm(np.array(trans_b) - np.array(conv_pt))
    f_near = min(max(min_from_point / max(len_a, 1e-6), min_from_point / max(len_b, 1e-6)), 0.9)
    fractions = np.linspace(f_near, 1.0, N_RUNGS)

    subpaths = []
    for f in fractions:
        pa = np.array(conv_pt) + (np.array(trans_a) - np.array(conv_pt)) * f
        pb_direct = np.array(conv_pt) + (np.array(trans_b) - np.array(conv_pt)) * f

        hit = None
        if ref_dir is not None:
            hit = _line_line_intersect(pa, ref_dir, conv_pt, trans_b)
        pb = hit if hit is not None else pb_direct

        d = pb - pa
        n = np.linalg.norm(d)
        if n < 1e-9:
            continue
        u = d / n
        end_a = pa - u * overshoot
        end_b = pb + u * overshoot
        subpaths.append([[list(end_a)] * 3, [list(end_b)] * 3])
    return subpaths


class TaperJoin(EffectExtension):

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

        data = []  # kept in ORIGINAL selection order - data[0] is "the first selected"
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
        first_selected = data[0]

        (railA, spA), (railB, spB) = under["rails"]
        tipA = railA.end() if at_end else railA.start()
        tipB = railB.end() if at_end else railB.start()
        mid_tip = (np.array(tipA) + np.array(tipB)) / 2.0
        tanA = tr._end_tangent(railA, at_end)
        tanB = tr._end_tangent(railB, at_end)
        avg_tan = tanA + tanB
        n = np.linalg.norm(avg_tan)
        avg_tan = avg_tan / n if n > 1e-9 else tanA

        conv_pt = tr._line_hit(mid_tip, avg_tan, over["rails"])

        angle = _crossing_angle_deg(avg_tan, _nearest_tangent_on_rails(mid_tip, over["rails"])) \
            if conv_pt is not None else 0.0

        if conv_pt is None or angle <= ANGLE_THRESHOLD_DEG:
            # Near-graze (or no clean crossing found) - fall back to Touch
            # Rails' independent per-rail touch, unchanged.
            new_rail_subpaths = []
            for rail, subpath in under["rails"]:
                tip = rail.end() if at_end else rail.start()
                tangent = tr._end_tangent(rail, at_end)
                hit = tr._line_hit(tip, tangent, over["rails"])
                new_rail_subpaths.append(
                    subpath if hit is None else tr._extend_or_trim(subpath, rail, at_end, hit)
                )
            under_el = under["el"]
            under_el.path = rrs.make_csp(new_rail_subpaths + under["rungs"]).to_path()
            under_el.transform = -under_el.getparent().composed_transform()
            return

        # A = real rail-to-rail separation (measure_A) - the scale a taper's
        # LENGTH (how far back it trims) should be based on. B = the
        # element's own cosmetic stroke-width style - these satin columns
        # draw a thin preview line (often "1px") unrelated to A, but it's
        # still the right scale for a small rung OVERSHOOT clearance, same
        # as rung_rail_satin's own B usage elsewhere.
        A_under = rrs.measure_A(railA, railB)
        B_under = rrs.stroke_width(under["el"])
        taper_len = TAPER_LEN_MULT * A_under
        taper_len = min(taper_len, MAX_TAPER_LEN_FRAC * min(railA.length, railB.length))

        new_rail_subpaths = []
        trans_pts = []
        for rail, subpath in under["rails"]:
            s_cut = (rail.length - taper_len) if at_end else taper_len
            cropped, trans_pt = _crop_to_s(subpath, rail, at_end, s_cut)
            new_rail_subpaths.append(_append_point(cropped, at_end, conv_pt))
            trans_pts.append(trans_pt)

        # Reference direction for the new rungs: the first SELECTED path's
        # rail nearest the convergence point, if one can be found.
        ref_dir = None
        try:
            ref_dir = _nearest_tangent_on_rails(conv_pt, first_selected["rails"])
            rn = np.linalg.norm(ref_dir)
            ref_dir = ref_dir / rn if rn > 1e-9 else None
        except Exception:
            ref_dir = None

        overshoot = STROKE_TAPER_MULT * B_under
        new_rungs = _build_taper_join_rungs(conv_pt, trans_pts[0], trans_pts[1], ref_dir, overshoot)

        under_el = under["el"]
        under_el.path = rrs.make_csp(new_rail_subpaths + under["rungs"] + new_rungs).to_path()
        under_el.transform = -under_el.getparent().composed_transform()


if __name__ == "__main__":
    TaperJoin().run()
