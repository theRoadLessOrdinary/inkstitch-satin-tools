#!/usr/bin/env python3
"""
Cross-Splice Rails (Ink/Stitch satin column helper)

Select 2 or 4 open strokes at a satin-column crossing whose direction
needs changing - the point of the tool is to swap which path continues
into which past an intersection, not to guess at long-range topology.
Because the selection is always scoped to just the strokes around that
one crossing (the user selects them that way), the closest pairing among
the open ends is always the correct one - no directional/tangent
analysis is needed or wanted. (An earlier version of this tool scored
candidate pairs by tangent alignment instead of distance, on the theory
that the nearest endpoint could be the wrong one - that version mismatched
strokes across an entire Celtic-knot design, latching onto a distant loop
that happened to have a similar crossing angle by the design's own
rotational symmetry. Plain nearest-distance matching doesn't have that
failure mode, and is correct precisely because the selection is already
scoped to one local crossing.)

With exactly 4 strokes, of the 3 ways to split them into two pairs,
whichever split has the smallest TOTAL distance (summed over both pairs'
closest end combo) wins.

Once paired, each connection is spliced into one continuous open path.
A genuine taper at the joining end (same corner test
`rung_rail_satin._classify_end` itself uses) is discarded - the whole
taper segment is dropped, not just bridged past - and replaced by a
single straight line to the new joint point; a non-tapered end is simply
extended with a new straight segment. No attempt is made to match the
old curve's arc at the joint - per the source process, a plain straight
connector is enough, and the rung generator's own tighter interior
spacing on low-straightness runs already produces sane rungs there. (A
human doing this by hand would go on to add nodes a couple mm either
side of the joint and delete the joint nodes so Inkscape's own curve
fit smooths it cosmetically - that's a manual polish step on top of
this, not something the rung math needs.)

With 4 strokes selected, splicing always produces exactly 2 rails, so
`rung_rail_satin.apply_rungs` is run on the result directly - one
operation, start to finish. With 2 strokes selected, this just produces
the single spliced open path (no rung generation) - useful as a
standalone splice to feed into Generate Ladder Rungs by hand.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import inkex
from inkex import EffectExtension, ShapeElement, PathElement

import rung_rail_satin as rrs

# Same threshold rung_rail_satin._classify_end uses to separate a genuine
# taper corner from a smooth multi-segment join - kept in sync by hand.
CORNER_COS_THRESHOLD = 0.97


def _is_taper_end(rail, at_start):
    """A genuine taper at this end: a real corner (per the same tangent
    test rung_rail_satin uses) AND the segment adjacent to this end is
    the SHORTER side of that corner. With exactly 2 segments, the one
    internal joint is structurally "adjacent" to both ends at once -
    without this length check, a long main-body curve ending in one
    short corner-turn tail (e.g. path501c) would get its long curve
    dropped when joining at the FAR end, not just its actual short
    taper - same disambiguation problem rung_rail_satin's own
    _classify_end solves by comparing segment length to the rail
    spacing A; there's no such external reference here, so the segment
    is compared to the rest of its own rail instead."""
    if rail.n < 2:
        return False
    cos = rail.corner_cos_at_start() if at_start else rail.corner_cos_at_end()
    if cos >= CORNER_COS_THRESHOLD:
        return False
    adjacent_len = (rail.seg_bounds[1] - rail.seg_bounds[0]) if at_start else rail.last_seg_length()
    return adjacent_len < (rail.length - adjacent_len)


def _end_point(rail, at_start):
    return rail.start() if at_start else rail.end()


def _gap(railA, a_start, railB, b_start):
    pA = np.array(_end_point(railA, a_start))
    pB = np.array(_end_point(railB, b_start))
    return float(np.linalg.norm(pB - pA))


def _best_end_combo(railA, railB):
    """Of the 4 ways to pick one end from each rail, the closest pair:
    (gap, a_start, b_start)."""
    best = None
    for a_start in (True, False):
        for b_start in (True, False):
            gap = _gap(railA, a_start, railB, b_start)
            if best is None or gap < best[0]:
                best = (gap, a_start, b_start)
    return best


def _splice(railA, a_start, csp_A, railB, b_start, csp_B):
    """Merge into one CSP subpath connecting railA's chosen end to
    railB's chosen end - discarding a genuine taper segment at either
    joining end, and bridging the two remaining ends with ONE new
    straight segment (not an extra intermediate point - the join is
    exactly two adjacent nodes with their facing handles flattened, so
    the one segment between them renders straight)."""
    left = rrs.reverse_subpath(csp_A) if a_start else list(csp_A)
    right = list(csp_B) if b_start else rrs.reverse_subpath(csp_B)

    if _is_taper_end(railA, a_start):
        left = left[:-1]
    if _is_taper_end(railB, b_start):
        right = right[1:]

    left = list(left)
    left[-1] = [list(left[-1][0]), list(left[-1][1]), list(left[-1][1])]
    right = list(right)
    right[0] = [list(right[0][1]), list(right[0][1]), list(right[0][2])]

    return left + right


class CrossSpliceSatin(EffectExtension):

    def effect(self):
        sel = list(self.svg.selection.get(ShapeElement))
        els = []
        for e in sel:
            if not isinstance(e, PathElement):
                if hasattr(e, "to_path_element"):
                    e = e.to_path_element()
                else:
                    raise inkex.AbortExtension(
                        "All selected objects must be paths (or shapes convertible to paths)."
                    )
            els.append(e)

        if len(els) not in (2, 4):
            raise inkex.AbortExtension(
                "Select exactly 2 open strokes (to splice into one rail) or "
                "4 (two rail pairs - spliced and fed straight into rung generation)."
            )

        # Selected strokes commonly carry DIFFERENT transforms (seen live:
        # path500a/d have none, path500b/c carry a translate()) - comparing
        # or splicing their raw local path data directly would silently
        # corrupt both the alignment scoring and the final geometry. Bake
        # each one's full composed (root-SVG) transform into its path data
        # up front, so everything downstream operates in one shared space.
        csps, rails = [], []
        for e in els:
            abs_path = e.path.to_absolute().transform(e.composed_transform())
            csp = abs_path.to_superpath()
            if len(csp) != 1:
                raise inkex.AbortExtension(
                    f"Each selected stroke must be a single open subpath (found {len(csp)} on one selection)."
                )
            csps.append(csp[0])
            rails.append(rrs.Rail(csp[0]))

        if len(els) == 2:
            gap, a_start, b_start = _best_end_combo(rails[0], rails[1])
            pairs = [(0, 1, a_start, b_start)]
        else:
            pairs, best_total = None, None
            for partition in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
                combos, total = [], 0.0
                for (i, j) in partition:
                    gap, a_start, b_start = _best_end_combo(rails[i], rails[j])
                    total += gap
                    combos.append((i, j, a_start, b_start))
                if best_total is None or total < best_total:
                    best_total, pairs = total, combos

        merged_subpaths = [
            _splice(rails[i], a_start, csps[i], rails[j], b_start, csps[j])
            for (i, j, a_start, b_start) in pairs
        ]

        # The merged path data is now in root-SVG coordinates (see above),
        # so `keep` must lose its own transform, and gain the exact
        # inverse of whatever its parent contributes, so nothing is
        # transformed a second time on render.
        keep = sorted(els, key=lambda e: list(e.getparent()).index(e))[-1]
        keep.path = rrs.make_csp(merged_subpaths).to_path()
        keep.transform = -keep.getparent().composed_transform()
        for e in els:
            if e is not keep:
                e.delete()

        if len(els) == 4:
            rrs.apply_rungs(keep)


if __name__ == "__main__":
    CrossSpliceSatin().run()
